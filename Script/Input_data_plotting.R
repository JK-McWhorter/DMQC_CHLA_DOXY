### =========================================================================
### BGC-Argo Raw Input Data Visualization Pipeline
### Direct Section Contour & QC Plots (Enhanced Missing-Data Visualization)
### =========================================================================

# =========================================================================
# 1. ENVIRONMENT CLEANUP & PACKAGE INITIALIZATION
# =========================================================================

cat("\014")
rm(list = ls())

required_pkgs <- c("dplyr", "ggplot2", "ncdf4", "viridis", "patchwork", "interp")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# =========================================================================
# 2. HELPER FUNCTIONS FOR MISSING VALUE DETECTION
# =========================================================================

is_argo_missing_val <- function(x) {
  is.na(x) | is.nan(x) | is.infinite(x) | (x == -99999) | 
    (x == 99999) | (x == -99999.0) | (x == 99999.0)
}

standardize_to_na <- function(x) {
  x[is_argo_missing_val(x)] <- NA_real_
  return(x)
}

# Python-identical QC Color Palette (Flag 9 shown as Black)
qc_color_scale <- c(
  "0" = "#b3b3b3",
  "1" = "#33A02C",
  "2" = "#B2DF8A",
  "3" = "#FDBF6F",
  "4" = "#E31A1C",
  "5" = "#1F78B4",
  "8" = "#CAB2D6",
  "9" = "#000000"
)

# =========================================================================
# 3. DIRECTORY CONFIGURATION & FLOAT IDS
# =========================================================================

dir_main_root  <- "C:/Users/Jennifer.McWhorter/Documents/GitHub/DMQC_CHLA_DOXY"
dir_dmqc_input <- "C:/Users/Jennifer.McWhorter/Documents/Data/2026/DMQC_INPUT"
dir_plot_base  <- file.path(dir_main_root, "Plots")

today_str <- format(Sys.Date(), "%Y-%m-%d")

float_ids <- c(
  4903622, 2904010, 2904011, 4903624, 4903625, 4903904, 
  6999992, 7902327, 3902693, 1902800, 7901009, 4903911, 
  4903910, 5907206
)

# =========================================================================
# 4. PLOTTING & GRID INTERPOLATION ROUTINES
# =========================================================================

# Linear interpolation with distance-bounded nearest fallback to prevent unrealistic fill across large gaps
create_interpolated_grid_r <- function(clean_df, val_col, grid_nx = 120, grid_ny = 120, max_gap_norm = 0.08) {
  if (is.null(clean_df) || nrow(clean_df) < 5) return(NULL)
  if (!val_col %in% names(clean_df)) return(NULL)
  
  sub_df <- clean_df %>%
    filter(
      !is_argo_missing_val(CYCLE_NUMBER),
      !is_argo_missing_val(PRES),
      !is_argo_missing_val(.data[[val_col]])
    )
  
  grouped <- sub_df %>%
    group_by(CYCLE_NUMBER, PRES) %>%
    summarise(z_val = mean(.data[[val_col]], na.rm = TRUE), .groups = "drop")
  
  if (nrow(grouped) < 5 || n_distinct(grouped$CYCLE_NUMBER) < 2 || n_distinct(grouped$PRES) < 2) {
    return(NULL)
  }
  
  val_range <- range(grouped$z_val, na.rm = TRUE)
  if (!all(is.finite(val_range)) || diff(val_range) == 0) return(NULL)
  
  xo_seq <- seq(min(grouped$CYCLE_NUMBER), max(grouped$CYCLE_NUMBER), length.out = grid_nx)
  yo_seq <- seq(min(grouped$PRES), max(grouped$PRES), length.out = grid_ny)
  
  interp_lin <- tryCatch({
    suppressWarnings(
      interp::interp(
        x = grouped$CYCLE_NUMBER,
        y = grouped$PRES,
        z = grouped$z_val,
        xo = xo_seq,
        yo = yo_seq,
        linear = TRUE,
        extrap = FALSE,
        duplicate = "mean"
      )
    )
  }, error = function(e) NULL)
  
  if (is.null(interp_lin) || is.null(interp_lin$z) || !is.matrix(interp_lin$z)) {
    return(NULL)
  }
  
  z_mat <- interp_lin$z
  
  # Distance-bounded nearest-neighbor fill for internal contour continuity
  if (any(is.na(z_mat))) {
    na_coords <- which(is.na(z_mat), arr.ind = TRUE)
    c_scale <- diff(range(xo_seq))
    p_scale <- diff(range(yo_seq))
    if (c_scale == 0) c_scale <- 1
    if (p_scale == 0) p_scale <- 1
    
    norm_obs_x <- grouped$CYCLE_NUMBER / c_scale
    norm_obs_y <- grouped$PRES / p_scale
    
    for (idx in seq_len(nrow(na_coords))) {
      r <- na_coords[idx, 1]
      c <- na_coords[idx, 2]
      
      pt_x <- xo_seq[r] / c_scale
      pt_y <- yo_seq[c] / p_scale
      
      dist_sq <- (norm_obs_x - pt_x)^2 + (norm_obs_y - pt_y)^2
      min_dist <- min(dist_sq)
      
      # Only fill if the point is within proximity to an observation
      if (sqrt(min_dist) <= max_gap_norm) {
        nearest_idx <- which.min(dist_sq)
        z_mat[r, c] <- grouped$z_val[nearest_idx]
      }
    }
  }
  
  grid_df <- expand.grid(x = xo_seq, y = yo_seq)
  grid_df$z <- as.vector(z_mat)
  grid_df <- grid_df %>% filter(!is.na(z) & is.finite(z))
  
  if (nrow(grid_df) == 0) return(NULL)
  return(grid_df)
}

# 2-panel figure generation with enhanced missing value overlays
generate_section_qc_plot <- function(df_plot, var_name, qc_var_name, title_prefix, file_suffix, 
                                     output_dir, float_id, col_palette = "viridis", 
                                     unit_expr = "", is_chla = FALSE) {
  if (!var_name %in% names(df_plot)) return(NULL)
  
  if (qc_var_name %in% names(df_plot)) {
    raw_qc <- as.character(df_plot[[qc_var_name]])
    match_res <- regexpr("[0-9]", raw_qc)
    matched_qc <- ifelse(match_res != -1, regmatches(raw_qc, match_res), "9")
    df_plot$QC_VAL <- ifelse(is.na(matched_qc) | matched_qc == "", "9", matched_qc)
  } else {
    df_plot$QC_VAL <- "9"
  }
  
  v_raw <- df_plot[[var_name]]
  is_missing_mask <- is_argo_missing_val(v_raw) | (df_plot$QC_VAL == "9")
  
  df_valid <- df_plot[!is_missing_mask, ] %>%
    filter(!is_argo_missing_val(PRES) & !is_argo_missing_val(CYCLE_NUMBER))
  
  df_missing <- df_plot[is_missing_mask, ] %>%
    filter(!is_argo_missing_val(PRES) & !is_argo_missing_val(CYCLE_NUMBER))
  
  if (nrow(df_valid) == 0 && nrow(df_missing) == 0) return(NULL)
  
  pres_all <- df_plot$PRES[!is_argo_missing_val(df_plot$PRES)]
  max_pres <- if (length(pres_all) > 0) max(pres_all, na.rm = TRUE) else 2000.0
  if (is.na(max_pres) || max_pres <= 0) max_pres <- 2000.0
  
  grid_df <- if (nrow(df_valid) > 0) create_interpolated_grid_r(df_valid, var_name) else NULL
  
  # Cycle completeness summary for top rug indicator (flags entirely dead cycles)
  cycle_status <- df_plot %>%
    group_by(CYCLE_NUMBER) %>%
    summarise(
      has_valid_data = any(!is_argo_missing_val(.data[[var_name]]) & QC_VAL != "9"),
      .groups = "drop"
    )
  dead_cycles <- cycle_status %>% filter(!has_valid_data)
  
  # --- Upper Plot: Section Contour + Micro-Stippled Missing Cloud + Dead Cycle Rug ---
  p_data <- ggplot()
  
  if (!is.null(grid_df)) {
    if (is_chla) {
      chla_breaks <- c(-0.1, 0.01, 0.03, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.2, 2.0, 5.0)
      p_data <- p_data + 
        geom_contour_filled(data = grid_df, aes(x = x, y = y, z = z), breaks = chla_breaks) +
        scale_fill_viridis_d(option = "turbo", name = unit_expr, drop = FALSE)
    } else {
      p_data <- p_data + 
        geom_contour_filled(data = grid_df, aes(x = x, y = y, z = z), bins = 15) +
        scale_fill_viridis_d(option = col_palette, name = unit_expr)
    }
  } else if (nrow(df_valid) > 0) {
    p_data <- p_data + 
      geom_point(data = df_valid, aes(x = CYCLE_NUMBER, y = PRES, color = !!sym(var_name)), 
                 size = 1.0, alpha = 0.85) +
      scale_color_viridis_c(option = col_palette, name = unit_expr)
  }
  
  # High-density, low-opacity micro-stippling for missing/QC 9 measurements
  if (nrow(df_missing) > 0) {
    p_data <- p_data + 
      geom_point(
        data = df_missing, 
        aes(x = CYCLE_NUMBER, y = PRES), 
        shape = 16, 
        size = 0.5, 
        alpha = 0.20, 
        color = "#2E2E2E"
      )
  }
  
  # Red rug tick on the surface axis indicating a cycle completely devoid of valid data
  if (nrow(dead_cycles) > 0) {
    p_data <- p_data + 
      geom_rug(
        data = dead_cycles, 
        aes(x = CYCLE_NUMBER), 
        sides = "t", 
        color = "#D73027", 
        linewidth = 1.0
      )
  }
  
  p_data <- p_data + 
    scale_y_reverse(limits = c(max_pres, 0), expand = c(0, 0)) +
    labs(
      title = paste0("Float ", float_id, " - ", title_prefix, " Raw Input Section Profile"),
      subtitle = if (nrow(dead_cycles) > 0) "Top red ticks denote cycles with 100% missing data" else NULL,
      x = NULL,
      y = "Pressure (dbar)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "#7F0000"),
      panel.grid.major = element_line(linetype = "dashed", color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  # --- Lower Plot: Raw QC Scatter Map ---
  df_qc_plot <- df_plot %>% 
    filter(!is_argo_missing_val(PRES) & !is_argo_missing_val(CYCLE_NUMBER))
  
  u_qcs <- unique(df_qc_plot$QC_VAL)
  unique_qcs <- u_qcs[order(ifelse(grepl("^[0-9]+$", u_qcs), as.numeric(u_qcs), 99))]
  
  p_qc <- ggplot()
  
  for (q_flag in unique_qcs) {
    sub_df <- df_qc_plot %>% filter(QC_VAL == q_flag)
    p_qc <- p_qc + 
      geom_point(
        data = sub_df, 
        aes(x = CYCLE_NUMBER, y = PRES, color = QC_VAL), 
        size = 1.1, 
        alpha = if (q_flag == "9") 0.45 else 0.8, 
        shape = 16
      )
  }
  
  p_qc <- p_qc + 
    scale_y_reverse(limits = c(max_pres, 0), expand = c(0, 0)) +
    scale_color_manual(values = qc_color_scale, limits = unique_qcs, name = "QC Flag") +
    labs(
      title = paste0(title_prefix, " Input QC Flags Section Map"),
      x = "Cycle Number",
      y = "Pressure (dbar)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(linetype = "dashed", color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1)))
  
  # Combined two-panel section graphic
  p_combined <- (p_data / p_qc) + 
    plot_layout(heights = c(1, 1), guides = "keep") +
    plot_annotation(
      title = paste0("BGC-Argo Float ", float_id, " - ", title_prefix, " Input Data"),
      subtitle = paste0("Plotted on: ", today_str),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  out_file <- file.path(output_dir, sprintf("Float_%d_%s_Raw_Input_Summary.png", float_id, file_suffix))
  ggsave(out_file, plot = p_combined, width = 12, height = 8, dpi = 150)
  message(sprintf("[%s] Saved Enhanced Raw QC Plot -> %s", float_id, out_file))
}

# =========================================================================
# 5. INGESTION & PLOTTING LOOP
# =========================================================================

for (woko in float_ids) {
  float_id   <- as.numeric(woko)
  active_wmo <- as.character(woko)
  
  dir_float_plots <- file.path(dir_plot_base, today_str, active_wmo, "RAW_INPUT_PLOTS")
  if (!dir.exists(dir_float_plots)) dir.create(dir_float_plots, recursive = TRUE)
  
  message(paste0("\n=========================================="))
  message(paste0("  PLOTTING INPUT DATA FOR FLOAT: ", active_wmo))
  message(paste0("==========================================\n"))
  
  br_pattern <- paste0("^BR", active_wmo, "_.*\\.nc$")
  local_br_files <- list.files(dir_dmqc_input, pattern = br_pattern, full.names = FALSE)
  local_br_files <- local_br_files[!grepl("D\\.nc$", local_br_files)]
  
  if (length(local_br_files) == 0) {
    warning(paste("No matching BR files found for WMO:", active_wmo))
    next
  }
  
  float_data <- NULL
  
  for (prof_name in sort(local_br_files)) {
    www <- sub("^BR", "", sub("\\.nc$", "", prof_name))
    br_file <- file.path(dir_dmqc_input, paste0("BR", www, ".nc"))
    d_file  <- file.path(dir_dmqc_input, paste0("D",  www, ".nc"))
    r_file  <- file.path(dir_dmqc_input, paste0("R",  www, ".nc"))
    
    physic_file <- if (file.exists(d_file)) d_file else if (file.exists(r_file)) r_file else NULL
    if (!file.exists(br_file) || is.null(physic_file)) next
    
    profile_bio    <- nc_open(br_file, readunlim = FALSE, write = FALSE)
    profile_physic <- nc_open(physic_file, readunlim = FALSE, write = FALSE)
    
    pres_raw <- as.vector(ncvar_get(profile_physic, "PRES"))
    pres <- standardize_to_na(pres_raw)
    
    pres_qc_raw <- ncvar_get(profile_physic, "PRES_QC")
    pres_qc_str <- if (is.matrix(pres_qc_raw)) apply(pres_qc_raw, 2, function(x) paste(x, collapse = "")) else paste(pres_qc_raw, collapse = "")
    pres_qc <- as.numeric(unlist(strsplit(pres_qc_str[1], "")))[1:length(pres)]
    pres_qc[is.na(pres) | is.na(pres_qc)] <- 9
    
    chl <- rep(NA_real_, length(pres));  chl_qc  <- rep("9", length(pres))
    fluo <- rep(NA_real_, length(pres)); fluo_qc <- rep("9", length(pres))
    bbp <- rep(NA_real_, length(pres));  bbp_qc  <- rep("9", length(pres))
    doxy <- rep(NA_real_, length(pres)); doxy_qc <- rep("9", length(pres))
    
    if ("CHLA" %in% names(profile_bio$var)) {
      chl <- standardize_to_na(as.vector(ncvar_get(profile_bio, "CHLA")))
      chl_qc_v <- paste(ncvar_get(profile_bio, "CHLA_QC"), collapse = "")
      chl_qc <- unlist(strsplit(chl_qc_v, ""))[1:length(pres)]
    }
    
    if ("FLUORESCENCE_CHLA" %in% names(profile_bio$var)) {
      fluo <- standardize_to_na(as.vector(ncvar_get(profile_bio, "FLUORESCENCE_CHLA")))
      fluo_qc_v <- paste(ncvar_get(profile_bio, "FLUORESCENCE_CHLA_QC"), collapse = "")
      fluo_qc <- unlist(strsplit(fluo_qc_v, ""))[1:length(pres)]
    }
    
    bbp_var_name <- intersect(c("BBP700", "BBP_700"), names(profile_bio$var))
    if (length(bbp_var_name) > 0) {
      bbp_raw <- standardize_to_na(as.vector(ncvar_get(profile_bio, bbp_var_name[1])))
      if (length(bbp_raw) == length(pres)) bbp <- bbp_raw
      qc_name <- paste0(bbp_var_name[1], "_QC")
      if (qc_name %in% names(profile_bio$var)) {
        bbp_qc_v <- paste(ncvar_get(profile_bio, qc_name), collapse = "")
        bbp_qc <- unlist(strsplit(bbp_qc_v, ""))[1:length(pres)]
      }
    }
    
    if ("DOXY" %in% names(profile_bio$var)) {
      doxy <- standardize_to_na(as.vector(ncvar_get(profile_bio, "DOXY")))
      if ("DOXY_QC" %in% names(profile_bio$var)) {
        doxy_qc_raw <- ncvar_get(profile_bio, "DOXY_QC")
        doxy_qc_str <- if (is.matrix(doxy_qc_raw)) apply(doxy_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(doxy_qc_raw, collapse = "")
        doxy_qc <- unlist(strsplit(doxy_qc_str, ""))[1:length(pres)]
      }
    }
    
    data_sub <- data.frame(
      CYCLE_NUMBER = as.numeric(sub(".*_([0-9]+)$", "\\1", www)),
      PRES = as.numeric(pres),
      PRES_QC = pres_qc,
      CHLA = chl,
      CHLA_QC = chl_qc,
      FLUO_CHLA = fluo,
      FLUO_CHLA_QC = fluo_qc,
      BBP700 = bbp,
      BBP700_QC = bbp_qc,
      DOXY = doxy,
      DOXY_QC = doxy_qc,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(PRES) & !is_argo_missing_val(PRES)) %>%
      arrange(PRES)
    
    if (nrow(data_sub) > 0) {
      float_data <- rbind(float_data, data_sub)
    }
    
    nc_close(profile_bio)
    nc_close(profile_physic)
  }
  
  if (is.null(float_data) || nrow(float_data) == 0) {
    warning(paste("No valid profiles found for float:", active_wmo))
    next
  }
  
  # --- Generate Section Plots for Raw Inputs ---
  generate_section_qc_plot(
    df_plot = float_data, 
    var_name = "DOXY", 
    qc_var_name = "DOXY_QC", 
    title_prefix = "DOXY Raw", 
    file_suffix = "DOXY_Raw", 
    output_dir = dir_float_plots, 
    float_id = float_id, 
    col_palette = "plasma", 
    unit_expr = expression(Oxygen ~ (mu * "mol/kg"))
  )
  
  generate_section_qc_plot(
    df_plot = float_data, 
    var_name = "BBP700", 
    qc_var_name = "BBP700_QC", 
    title_prefix = "BBP700 Raw", 
    file_suffix = "BBP700_Raw", 
    output_dir = dir_float_plots, 
    float_id = float_id, 
    col_palette = "turbo", 
    unit_expr = expression(b[bp](700) ~ (m^{-1}))
  )
  
  generate_section_qc_plot(
    df_plot = float_data, 
    var_name = "CHLA", 
    qc_var_name = "CHLA_QC", 
    title_prefix = "CHLA Raw", 
    file_suffix = "CHLA_Raw", 
    output_dir = dir_float_plots, 
    float_id = float_id, 
    col_palette = "turbo", 
    unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), 
    is_chla = TRUE
  )
  
  generate_section_qc_plot(
    df_plot = float_data, 
    var_name = "FLUO_CHLA", 
    qc_var_name = "FLUO_CHLA_QC", 
    title_prefix = "FLUO_CHLA Raw", 
    file_suffix = "FLUO_CHLA_Raw", 
    output_dir = dir_float_plots, 
    float_id = float_id, 
    col_palette = "mako", 
    unit_expr = "Fluorescence (RU)"
  )
}