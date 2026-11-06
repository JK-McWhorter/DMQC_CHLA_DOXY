### =========================================================================
### BGC-Argo CHLA, BBP700 & DOXY Corrections DMQC Pipeline
### Integrated Multi-Mode Processing: R-mode, A-mode & D-mode Preservation
### QC Summary Layout & Missing Data Visualization Exactly Matching Python
### Filtering Depth Levels Consistently Missing Across All BGC Parameters
### =========================================================================

# =========================================================================
# 1. ENVIRONMENT CLEANUP, PACKAGE INITIALIZATION & BUG HOTPATCH
# =========================================================================

cat("\014")
rm(list = ls())

required_pkgs <- c("dplyr", "ggplot2", "lubridate", "gridExtra", "tidyverse", 
                   "ggeffects", "mapdata", "suncalc", "ncdf4", "oce", "terra", 
                   "sf", "stars", "gstat", "rnaturalearth", "rnaturalearthdata", 
                   "viridis", "knitr", "patchwork", "stringr", "interp", "scales")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# HOTPATCH: Fixes ggplot2 3.5.0 bug with scales >= 1.3.0 in stat_contour_filled
# -------------------------------------------------------------------------
safe_zero_range <- function(x, tol = 1000 * .Machine$double.eps) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  if (length(x) > 2) x <- range(x, na.rm = TRUE)
  if (length(x) == 1) return(TRUE)
  if (anyNA(x)) return(NA)
  abs(x[1] - x[2]) < tol
}

try(assignInNamespace("zero_range", safe_zero_range, ns = "scales"), silent = TRUE)
try(assignInNamespace("zero_range", safe_zero_range, ns = "ggplot2"), silent = TRUE)

# =========================================================================
# 2. HELPER FUNCTIONS FOR MISSING VALUES, QC & CONSISTENCY AUDITING
# =========================================================================

is_argo_missing_val <- function(x) {
  is.na(x) | is.nan(x) | is.infinite(x) | (x == -99999) | 
    (x == 99999) | (x == -99999.0) | (x == 99999.0)
}

standardize_to_na <- function(x) {
  x[is_argo_missing_val(x)] <- NA_real_
  return(x)
}

standardize_missing_val <- function(x) {
  x[is_argo_missing_val(x)] <- -99999.0
  return(x)
}

standardize_missing_qc <- function(qc) {
  qc_chr <- as.character(qc)
  qc_chr[is.na(qc_chr) | is.nan(qc) | qc_chr == "" | qc_chr == "NA" | qc_chr %in% c("-99999", "99999")] <- "9"
  return(qc_chr)
}

safe_max <- function(x) {
  x_clean <- x[!is.na(x) & is.finite(x) & !is_argo_missing_val(x)]
  if (length(x_clean) > 0) max(x_clean) else NA_real_
}

safe_min <- function(x) {
  x_clean <- x[!is.na(x) & is.finite(x) & !is_argo_missing_val(x)]
  if (length(x_clean) > 0) min(x_clean) else NA_real_
}

audit_missing_mismatches <- function(df_export, float_wmo, output_dir) {
  is_qc_missing <- function(qc_val) {
    if (is.na(qc_val) || qc_val == "" || qc_val == " ") return(TRUE)
    qc_str <- as.character(qc_val)
    match_res <- regmatches(qc_str, regexpr("[0-9]", qc_str))
    if (length(match_res) == 0) return(TRUE)
    return(match_res == "9")
  }
  
  mismatches <- list()
  
  var_stems <- c("DOXY", "DOXY_ADJUSTED", "DOXY_FINAL",
                 "CHLA", "CHLA_ADJUSTED", "CHLA_FINAL",
                 "BBP700", "BBP700_ADJUSTED", "BBP700_FINAL",
                 "CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_FINAL")
  
  for (stem in var_stems) {
    qc_name <- paste0(stem, "_QC")
    if (all(c(stem, qc_name) %in% names(df_export))) {
      vals <- df_export[[stem]]
      qcs  <- df_export[[qc_name]]
      val_is_miss <- is_argo_missing_val(vals)
      qc_is_miss  <- sapply(qcs, is_qc_missing)
      
      bad_idx <- which(val_is_miss != qc_is_miss)
      if (length(bad_idx) > 0) {
        for (idx in bad_idx) {
          mismatches[[length(mismatches) + 1]] <- data.frame(
            WMO             = as.character(float_wmo),
            CYCLE           = df_export$CYCLE_NUMBER[idx],
            PRES            = df_export$PRES[idx],
            LEVEL_IDX       = idx,
            AUDIT_CATEGORY  = "INTERNAL_VALUE_VS_QC",
            PAIR_OR_VAR     = paste(stem, "vs", qc_name),
            DESCRIPTION     = if (val_is_miss[idx]) "Value is Missing (-99999/NA) but QC is NOT 9" else "QC is 9 but Value is NOT Missing",
            VAL_1           = as.character(vals[idx]),
            VAL_2           = as.character(qcs[idx]),
            IS_MISSING_1    = val_is_miss[idx],
            IS_MISSING_2    = qc_is_miss[idx],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  cross_pairs <- list(
    list(v1 = "DOXY", v2 = "DOXY_ADJUSTED", q1 = "DOXY_QC", q2 = "DOXY_ADJUSTED_QC", label = "DOXY: R vs A"),
    list(v1 = "DOXY", v2 = "DOXY_FINAL",    q1 = "DOXY_QC", q2 = "DOXY_FINAL_QC",    label = "DOXY: R vs D"),
    list(v1 = "DOXY_ADJUSTED", v2 = "DOXY_FINAL", q1 = "DOXY_ADJUSTED_QC", q2 = "DOXY_FINAL_QC", label = "DOXY: A vs D"),
    list(v1 = "CHLA", v2 = "CHLA_ADJUSTED", q1 = "CHLA_QC", q2 = "CHLA_ADJUSTED_QC", label = "CHLA: R vs A"),
    list(v1 = "CHLA", v2 = "CHLA_FINAL",    q1 = "CHLA_QC", q2 = "CHLA_FINAL_QC",    label = "CHLA: R vs D"),
    list(v1 = "CHLA_ADJUSTED", v2 = "CHLA_FINAL", q1 = "CHLA_ADJUSTED_QC", q2 = "CHLA_FINAL_QC", label = "CHLA: A vs D"),
    list(v1 = "BBP700", v2 = "BBP700_ADJUSTED", q1 = "BBP700_QC", q2 = "BBP700_ADJUSTED_QC", label = "BBP700: R vs A"),
    list(v1 = "BBP700", v2 = "BBP700_FINAL",    q1 = "BBP700_QC", q2 = "BBP700_FINAL_QC",    label = "BBP700: R vs D"),
    list(v1 = "BBP700_ADJUSTED", v2 = "BBP700_FINAL", q1 = "BBP700_ADJUSTED_QC", q2 = "BBP700_FINAL_QC", label = "BBP700: A vs D"),
    list(v1 = "CHLA_FLUORESCENCE", v2 = "CHLA_FLUORESCENCE_FINAL", 
         q1 = "CHLA_FLUORESCENCE_QC", q2 = "CHLA_FLUORESCENCE_FINAL_QC", label = "CHLA_FLUO: R vs D")
  )
  
  for (pair in cross_pairs) {
    if (all(c(pair$v1, pair$v2, pair$q1, pair$q2) %in% names(df_export))) {
      v1_vals <- df_export[[pair$v1]]
      v2_vals <- df_export[[pair$v2]]
      q1_vals <- df_export[[pair$q1]]
      q2_vals <- df_export[[pair$q2]]
      
      v1_miss <- is_argo_missing_val(v1_vals)
      v2_miss <- is_argo_missing_val(v2_vals)
      q1_miss <- sapply(q1_vals, is_qc_missing)
      q2_miss <- sapply(q2_vals, is_qc_missing)
      
      val_diff_idx <- which(v1_miss != v2_miss)
      if (length(val_diff_idx) > 0) {
        for (idx in val_diff_idx) {
          mismatches[[length(mismatches) + 1]] <- data.frame(
            WMO             = as.character(float_wmo),
            CYCLE           = df_export$CYCLE_NUMBER[idx],
            PRES            = df_export$PRES[idx],
            LEVEL_IDX       = idx,
            AUDIT_CATEGORY  = "CROSS_MODE_VALUE_DIVERGENCE",
            PAIR_OR_VAR     = pair$label,
            DESCRIPTION     = sprintf("Value missing in %s but present in %s", 
                                      ifelse(v1_miss[idx], pair$v1, pair$v2),
                                      ifelse(v1_miss[idx], pair$v2, pair$v1)),
            VAL_1           = as.character(v1_vals[idx]),
            VAL_2           = as.character(v2_vals[idx]),
            IS_MISSING_1    = v1_miss[idx],
            IS_MISSING_2    = v2_miss[idx],
            stringsAsFactors = FALSE
          )
        }
      }
      
      qc_diff_idx <- which(q1_miss != q2_miss)
      if (length(qc_diff_idx) > 0) {
        for (idx in qc_diff_idx) {
          mismatches[[length(mismatches) + 1]] <- data.frame(
            WMO             = as.character(float_wmo),
            CYCLE           = df_export$CYCLE_NUMBER[idx],
            PRES            = df_export$PRES[idx],
            LEVEL_IDX       = idx,
            AUDIT_CATEGORY  = "CROSS_MODE_QC_DIVERGENCE",
            PAIR_OR_VAR     = paste(pair$q1, "vs", pair$q2),
            DESCRIPTION     = sprintf("QC 9 in %s but valid in %s",
                                      ifelse(q1_miss[idx], pair$q1, pair$q2),
                                      ifelse(q1_miss[idx], pair$q2, pair$q1)),
            VAL_1           = as.character(q1_vals[idx]),
            VAL_2           = as.character(q2_vals[idx]),
            IS_MISSING_1    = q1_miss[idx],
            IS_MISSING_2    = q2_miss[idx],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  df_mismatches <- if (length(mismatches) > 0) bind_rows(mismatches) else data.frame()
  
  cat("\n", rep("=", 95), "\n", sep = "")
  cat(sprintf("MISSING DATA & INCONSISTENCY AUDIT SUMMARY: FLOAT %s\n", float_wmo))
  cat(rep("=", 95), "\n", sep = "")
  
  if (nrow(df_mismatches) == 0) {
    cat(sprintf("[%s] PASSED: All missing value masks and QC flags align consistently internally & across all modes.\n", float_wmo))
  } else {
    cat(sprintf("[%s] Found %d total level-by-level inconsistencies.\n\n", float_wmo, nrow(df_mismatches)))
    summary_table <- df_mismatches %>%
      group_by(AUDIT_CATEGORY, PAIR_OR_VAR) %>%
      summarise(MISMATCH_COUNT = n(), AFFECTED_CYCLES = n_distinct(CYCLE), .groups = "drop") %>%
      arrange(desc(MISMATCH_COUNT))
    print(as.data.frame(summary_table))
    
    audit_csv_path <- file.path(output_dir, sprintf("%s_missing_data_mismatches.csv", float_wmo))
    write.csv(df_mismatches, audit_csv_path, row.names = FALSE)
    cat(sprintf("\n>>> SAVED DETAILED AUDIT REPORT -> %s\n", audit_csv_path))
  }
  cat(rep("=", 95), "\n\n", sep = "")
  return(df_mismatches)
}

# =========================================================================
# 3. CONFIGURATION & DIRECTORY MANAGEMENT
# =========================================================================

dir_main_root    <- "C:/Users/Jennifer.McWhorter/Documents/GitHub/DMQC_CHLA_DOXY"
dir_dmqc_input   <- "C:/Users/Jennifer.McWhorter/Documents/Data/2026/DMQC_INPUT"

dir_functions    <- file.path(dir_main_root, "Functions")
dir_output       <- file.path(dir_main_root, "Output")
dir_plot_base    <- file.path(dir_main_root, "Plots")
dir_script       <- file.path(dir_main_root, "Script")
dir_tables       <- file.path(dir_main_root, "Tables")
dir_profs        <- file.path(dir_main_root, "Tables/profs")
dir_data_root    <- file.path(dir_main_root, "Data")
combo_output_dir <- file.path(dir_output, "COMBO_CHLA_DOXY_DMODE")

for (d in c(dir_functions, dir_output, dir_plot_base, dir_script, dir_tables, dir_profs, dir_data_root, combo_output_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

func.sources <- list.files(dir_functions, pattern = "\\.R$", full.names = TRUE)
for (f in func.sources) {
  source(f)
}

launch_date_dict <- tibble::tribble(
  ~float_ids, ~launch_date_str,
  4903622,    "2021-10-04 23:36:00",
  2904010,    "2025-05-22 04:32:00",
  2904011,    "2025-05-19 12:15:00",
  4903624,    "2021-09-25 09:51:00",
  4903625,    "2021-09-22 10:53:00",
  4903904,    "2025-06-24 07:00:00",
  6999992,    "2025-06-24 11:41:00",
  7902327,    "2026-05-21 04:14:00",
  3902693,    "2026-05-20 18:58:00",
  1902800,    "2026-05-29 09:33:00",
  7901009,    "2023-06-05 10:30:00",
  4903911,    "2026-05-15 21:52:00", 
  4903910,    "2026-05-20 18:00:00", 
  5907206,    "2026-05-29 09:33:00"
)

float_type_dict <- tibble::tribble(
  ~float_ids, ~float_type,
  4903622,    "aoml_apex",
  2904010,    "aoml_apex",
  2904011,    "aoml_apex",
  4903624,    "aoml_apex",
  4903625,    "aoml_apex",
  4903904,    "aoml_navis",
  6999992,    "aoml_navis",
  7902327,    "aoml_navis",
  3902693,    "aoml_apex",
  1902800,    "aoml_apex",
  7901009,    "aoml_navis",
  4903911,    "aoml_navis",
  4903910,    "aoml_navis",
  5907206,    "aoml_apex"
)

write_DOXY_adjusted_error <- function(nc, profile_idx, err_mbar, psal, temp, pres, dens, DOXY_adj) {
  valid_idx <- !is.na(psal) & !is.na(DOXY_adj) & !is.na(dens) & (dens > 0) & 
    !(psal %in% c(-99999, 99999)) & !(DOXY_adj %in% c(-99999, 99999, 0))
  DOXY_adj_error <- rep(NA_real_, length(psal))
  err_umol_L <- err_mbar * 1.00
  DOXY_adj_error[valid_idx] <- err_umol_L
  DOXY_adj_error_umol_kg <- (DOXY_adj_error * 1000) / dens
  
  if (!is.null(nc)) {
    DOXY_adj_full <- ncvar_get(nc, "DOXY_ADJUSTED")
    if (is.matrix(DOXY_adj_full)) {
      ncvar_put(nc, "DOXY_ADJUSTED_ERROR", DOXY_adj_error_umol_kg, 
                start = c(1, profile_idx), count = c(length(DOXY_adj_error_umol_kg), 1))
    } else {
      ncvar_put(nc, "DOXY_ADJUSTED_ERROR", DOXY_adj_error_umol_kg)
    }
  }
  return(DOXY_adj_error_umol_kg)
}

# Exact Python QC_COLOR_MAP
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

mode_levels      <- c("R-mode", "A-mode", "D-mode", "Bottle")
color_mapping    <- c("R-mode" = "#E31A1C", "A-mode" = "#8C510A", "D-mode" = "#33A02C", "Bottle" = "#000000")
linetype_mapping <- c("R-mode" = "dotted",  "A-mode" = "dashed",  "D-mode" = "solid",   "Bottle" = "blank")
linewidth_map    <- c("R-mode" = 0.8,       "A-mode" = 1.0,       "D-mode" = 1.2,       "Bottle" = 0.0)
shape_mapping    <- c("R-mode" = NA,        "A-mode" = NA,        "D-mode" = NA,        "Bottle" = 16)

calib_csv_path <- file.path(dir_tables, "Jul2026_inair_output_binflags5.csv")

calib_df <- read.csv(calib_csv_path) %>%
  mutate(float_ids = as.numeric(float_ids)) %>%
  select(-any_of("LAUNCH_DATE")) %>%
  left_join(launch_date_dict, by = "float_ids") %>%
  left_join(float_type_dict, by = "float_ids") %>%
  mutate(
    launch_date = lubridate::ymd_hms(launch_date_str, tz = "UTC")
  ) %>%
  rename(
    doxy_slope = SLOPE,
    doxy_drift = DRIFT
  ) %>%
  select(float_ids, doxy_slope, doxy_drift, launch_date, float_type)

float_ids     <- calib_df$float_ids 
target_floats <- c(6999992, 4903904, 1902800, 5907206, 4903624)
today_str     <- format(Sys.Date(), format = "%Y-%m-%d")

path_lut_nc     <- file.path(dir_data_root, "LUT/127630.nc")
path_bottle_val <- file.path(dir_data_root, "Bottle_samples/Master_Vali_24_July_2025.csv")

# =========================================================================
# 4. GLOBAL DATA INGESTION & PREPARATION
# =========================================================================

LUT <- rast(path_lut_nc)
df_LUT <- as.data.frame(LUT, xy = TRUE, na.rm = TRUE)
lut_sf <- st_as_sf(df_LUT, coords = c("x", "y"), crs = 4326)

if (file.exists(path_bottle_val)) {
  Bottle_Vali <- read.csv(path_bottle_val)
} else {
  warning("Bottle validation file not found at path: ", path_bottle_val)
  Bottle_Vali <- data.frame(Float_num = numeric(0), DOXY_bottle = numeric(0), PRES = numeric(0))
}

# =========================================================================
# 5. INDIVIDUAL PER-CYCLE PROFILE PLOTTING ROUTINE
# =========================================================================

generate_cycle_profile_plot <- function(df_cyc, var_cols, var_modes, var_title, unit_expr, 
                                        out_file, float_id, cyc_num, bottle_df = NULL) {
  plot_rows <- list()
  for (i in seq_along(var_cols)) {
    col_name <- var_cols[i]
    mode_name <- var_modes[i]
    if (col_name %in% names(df_cyc)) {
      vals <- standardize_to_na(as.numeric(df_cyc[[col_name]]))
      pres <- standardize_to_na(as.numeric(df_cyc$PRES))
      valid <- !is.na(vals) & !is_argo_missing_val(vals) & !is.na(pres) & !is_argo_missing_val(pres)
      if (any(valid)) {
        plot_rows[[length(plot_rows) + 1]] <- data.frame(
          PRES  = pres[valid],
          Value = vals[valid],
          Mode  = mode_name,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(plot_rows) == 0 && (is.null(bottle_df) || nrow(bottle_df) == 0)) {
    return(invisible(NULL))
  }
  
  df_plot_lines <- if (length(plot_rows) > 0) bind_rows(plot_rows) else data.frame()
  if (nrow(df_plot_lines) > 0) {
    df_plot_lines$Mode <- factor(df_plot_lines$Mode, levels = mode_levels)
    df_plot_lines <- df_plot_lines %>% arrange(Mode, PRES)
  }
  
  df_paths   <- df_plot_lines %>% group_by(Mode) %>% filter(n() >= 2) %>% ungroup()
  df_singles <- df_plot_lines %>% group_by(Mode) %>% filter(n() == 1) %>% ungroup()
  
  p <- ggplot()
  
  if (nrow(df_paths) > 0) {
    p <- p + 
      geom_path(
        data = df_paths,
        aes(x = Value, y = PRES, color = Mode, linetype = Mode, linewidth = Mode, group = Mode),
        na.rm = TRUE
      )
  }
  
  if (nrow(df_singles) > 0) {
    p <- p + 
      geom_point(
        data = df_singles,
        aes(x = Value, y = PRES, color = Mode),
        size = 2.4, shape = 18, na.rm = TRUE
      )
  }
  
  if (!is.null(bottle_df) && nrow(bottle_df) > 0) {
    bottle_df$Mode <- factor("Bottle", levels = mode_levels)
    p <- p + 
      geom_point(
        data = bottle_df,
        aes(x = Value, y = PRES, color = Mode, shape = Mode),
        size = 2.4, na.rm = TRUE
      )
  }
  
  cyc_time <- ""
  if ("DATE_TIME" %in% names(df_cyc) && !is.na(df_cyc$DATE_TIME[1]) && df_cyc$DATE_TIME[1] != "-99999") {
    dt_val <- df_cyc$DATE_TIME[1]
    if (inherits(dt_val, "POSIXt")) {
      cyc_time <- format(dt_val, format = "%Y-%m-%d %H:%M")
    } else {
      parsed_dt <- suppressWarnings(lubridate::as_datetime(as.character(dt_val)))
      if (!is.na(parsed_dt)) {
        cyc_time <- format(parsed_dt, format = "%Y-%m-%d %H:%M")
      } else {
        cyc_time <- as.character(dt_val)
      }
    }
  }
  
  cyc_lat <- if ("LATITUDE" %in% names(df_cyc)) round(df_cyc$LATITUDE[1], 3) else NA
  cyc_lon <- if ("LONGITUDE" %in% names(df_cyc)) round(df_cyc$LONGITUDE[1], 3) else NA
  
  sub_text <- paste0(
    if (!is.na(cyc_time) && cyc_time != "") paste0("Date: ", cyc_time, " UTC | ") else "",
    if (!is.na(cyc_lat) && !is.na(cyc_lon)) sprintf("Lat: %.3f, Lon: %.3f", cyc_lat, cyc_lon) else ""
  )
  
  p <- p + 
    scale_y_reverse(expand = expansion(mult = c(0.02, 0.05))) +
    scale_color_manual(name = "Data Mode", values = color_mapping, breaks = mode_levels, drop = FALSE) +
    scale_linetype_manual(name = "Data Mode", values = linetype_mapping, breaks = mode_levels, drop = FALSE) +
    scale_linewidth_manual(name = "Data Mode", values = linewidth_map, breaks = mode_levels, drop = FALSE) +
    scale_shape_manual(name = "Data Mode", values = shape_mapping, breaks = mode_levels, drop = FALSE, na.translate = FALSE) +
    labs(
      title = sprintf("Float %s - Cycle %03d: %s Profile", float_id, cyc_num, var_title),
      subtitle = sub_text,
      x = unit_expr,
      y = "Pressure (dbar)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      panel.grid.major = element_line(linetype = "dashed", color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80")
    ) +
    guides(
      color = guide_legend(override.aes = list(
        linetype = linetype_mapping,
        linewidth = c(0.8, 1.0, 1.2, 0.0),
        shape = c(NA, NA, NA, 16)
      ))
    )
  
  ggsave(out_file, plot = p, width = 7, height = 8, dpi = 130)
  return(invisible(TRUE))
}

# =========================================================================
# 6. BATCH EXECUTION LOOP (PROCESSING, EXPORT, AUDIT & DUAL PLOTTING)
# =========================================================================

for (woko in float_ids) {
  
  float_id   <- as.numeric(woko)
  active_wmo <- as.character(woko)
  
  dir_date_base   <- file.path(dir_plot_base, today_str)
  dir_float_plots <- file.path(dir_date_base, active_wmo)
  
  dir_chla_plots  <- file.path(dir_float_plots, "CHLA")
  dir_fluo_plots  <- file.path(dir_float_plots, "CHLA_FLUORESCENCE")
  dir_bbp_plots   <- file.path(dir_float_plots, "BBP")
  dir_lut_plots   <- file.path(dir_float_plots, "CHLA", "LUT")
  dir_doxy_plots  <- file.path(dir_float_plots, "DOXY")
  dir_qc_summary  <- file.path(dir_float_plots, "QC_SUMMARIES")
  
  for (p_dir in c(dir_chla_plots, dir_fluo_plots, dir_bbp_plots, dir_lut_plots, dir_doxy_plots, dir_qc_summary)) {
    if (!dir.exists(p_dir)) dir.create(p_dir, recursive = TRUE)
  }
  
  message(paste0("\n=========================================="))
  message(paste0("  STARTING PROCESSING FOR FLOAT: ", active_wmo))
  message(paste0("==========================================\n"))
  
  float_data <- NULL
  
  br_pattern <- paste0("^BR", active_wmo, "_.*\\.nc$")
  local_br_files <- list.files(dir_dmqc_input, pattern = br_pattern, full.names = FALSE)
  local_br_files <- local_br_files[!grepl("D\\.nc$", local_br_files)]
  
  if (length(local_br_files) == 0) {
    warning(paste("No matching BR files found in local directory for WMO:", active_wmo))
    next
  }
  
  meta_path <- file.path(dir_dmqc_input, paste0(active_wmo, "_meta.nc"))
  if (!file.exists(meta_path)) {
    meta_path <- file.path(dir_profs, paste0(active_wmo, "_meta.nc"))
  }
  
  # --- STEP 1: EXTRACT LOCAL PROFILES ---
  for (prof_name in sort(local_br_files)) {
    www <- sub("^BR", "", sub("\\.nc$", "", prof_name))
    
    br_file <- file.path(dir_dmqc_input, paste0("BR", www, ".nc"))
    d_file  <- file.path(dir_dmqc_input, paste0("D",  www, ".nc"))
    r_file  <- file.path(dir_dmqc_input, paste0("R",  www, ".nc"))
    
    physic_file <- if (file.exists(d_file)) d_file else if (file.exists(r_file)) r_file else NULL
    
    if (!file.exists(br_file) || is.null(physic_file)) {
      message(sprintf("Skipping %s: Missing corresponding physics file (D or R).", prof_name))
      next
    }
    
    profile_bio    <- nc_open(br_file, readunlim = FALSE, write = FALSE)
    profile_physic <- nc_open(physic_file, readunlim = FALSE, write = FALSE)
    
    lat_val <- ncvar_get(profile_physic, "LATITUDE")[1]
    lon_val <- ncvar_get(profile_physic, "LONGITUDE")[1]
    if (is.na(lat_val) || is.na(lon_val)) {
      nc_close(profile_bio); nc_close(profile_physic); next
    }
    
    jd <- ncvar_get(profile_physic, "JULD")[1]
    origin <- as.POSIXct("1950-01-01 00:00:00", tz = "UTC") 
    time_val <- origin + jd * 3600 * 24
    if (is.na(time_val)) {
      nc_close(profile_bio); nc_close(profile_physic); next
    }
    
    pres_raw <- as.vector(ncvar_get(profile_physic, "PRES"))
    pres <- standardize_to_na(pres_raw)
    
    pres_qc_raw <- ncvar_get(profile_physic, "PRES_QC")
    pres_qc_str <- if (is.matrix(pres_qc_raw)) apply(pres_qc_raw, 2, function(x) paste(x, collapse = "")) else paste(pres_qc_raw, collapse = "")
    pres_qc <- as.numeric(unlist(strsplit(pres_qc_str[1], "")))[1:length(pres)]
    pres_qc[is.na(pres) | is.na(pres_qc)] <- 9
    
    temp <- rep(NA_real_, length(pres)); temp_qc <- rep(9, length(pres))
    sal  <- rep(NA_real_, length(pres)); sal_qc  <- rep(9, length(pres))
    
    if ("TEMP" %in% names(profile_physic$var)) {
      temp <- standardize_to_na(as.vector(ncvar_get(profile_physic, "TEMP")))
      temp_qc_raw <- ncvar_get(profile_physic, "TEMP_QC")
      temp_qc_str <- if (is.matrix(temp_qc_raw)) apply(temp_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(temp_qc_raw, collapse = "")
      temp_qc <- as.numeric(unlist(strsplit(temp_qc_str, "")))[1:length(pres)]
      temp_qc[is.na(temp) | is.na(temp_qc)] <- 9
    }
    
    if ("PSAL" %in% names(profile_physic$var)) {
      sal <- standardize_to_na(as.vector(ncvar_get(profile_physic, "PSAL")))
      sal_qc_raw <- ncvar_get(profile_physic, "PSAL_QC")
      sal_qc_str <- if (is.matrix(sal_qc_raw)) apply(sal_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(sal_qc_raw, collapse = "")
      sal_qc <- as.numeric(unlist(strsplit(sal_qc_str, "")))[1:length(pres)]
      sal_qc[is.na(sal) | is.na(sal_qc)] <- 9
    }
    
    MLD <- NA_real_
    if (length(which(!is.na(pres))) > 15 && safe_max(pres) > 250) {
      sigma_all <- swSigmaTheta(sal, temp, pres)
      valid_mld <- !is.na(pres) & !is.na(sigma_all)
      if (sum(valid_mld) > 5) {
        dep_sigma <- swDepth(pres[valid_mld], lat_val)
        MLD <- MLD_calc(sigma_all[valid_mld], dep_sigma)
      }
    }
    
    chl <- rep(NA_real_, length(pres)); chl_qc <- rep(9, length(pres))
    chl_adj <- rep(NA_real_, length(pres)); chl_adj_qc <- rep(9, length(pres))
    bbp <- rep(NA_real_, length(pres)); bbp_qc <- rep(9, length(pres))
    fluo <- rep(NA_real_, length(pres)); fluo_qc <- rep(9, length(pres))
    doxy <- rep(NA_real_, length(pres)); doxy_qc <- rep(9, length(pres))
    doxy_adj <- rep(NA_real_, length(pres)); doxy_adj_qc <- rep(9, length(pres))
    scale_chla <- NA_real_; dark_chla <- NA_real_
    
    if ("CHLA" %in% names(profile_bio$var)) {
      chl <- standardize_to_na(as.vector(ncvar_get(profile_bio, "CHLA")))
      chl_qc_v <- paste(ncvar_get(profile_bio, "CHLA_QC"), collapse = "")
      chl_qc <- unlist(strsplit(chl_qc_v, ""))[1:length(pres)]
      chl_qc[is.na(chl) | is.na(chl_qc)] <- 9
    }
    
    if ("CHLA_ADJUSTED" %in% names(profile_bio$var)) {
      chl_adj <- standardize_to_na(as.vector(ncvar_get(profile_bio, "CHLA_ADJUSTED")))
      if ("CHLA_ADJUSTED_QC" %in% names(profile_bio$var)) {
        chl_adj_qc_v <- paste(ncvar_get(profile_bio, "CHLA_ADJUSTED_QC"), collapse = "")
        chl_adj_qc <- unlist(strsplit(chl_adj_qc_v, ""))[1:length(pres)]
        chl_adj_qc[is.na(chl_adj) | is.na(chl_adj_qc)] <- 9
      }
    }
    
    if ("FLUORESCENCE_CHLA" %in% names(profile_bio$var)) {
      fluo <- standardize_to_na(as.vector(ncvar_get(profile_bio, "FLUORESCENCE_CHLA")))
      fluo_qc_v <- paste(ncvar_get(profile_bio, "FLUORESCENCE_CHLA_QC"), collapse = "")
      fluo_qc <- unlist(strsplit(fluo_qc_v, ""))[1:length(pres)]
      fluo_qc[is.na(fluo) | is.na(fluo_qc)] <- 9
      
      if (file.exists(meta_path)) {
        try({
          meta_file <- nc_open(meta_path, readunlim = FALSE, write = FALSE)
          params <- ncvar_get(meta_file, "PARAMETER") 
          index_meta <- grep("^CHLA\\s*$", params)
          if (length(index_meta) > 0) {
            coeff_str <- ncvar_get(meta_file, "PREDEPLOYMENT_CALIB_COEFFICIENT")[index_meta]
            scale_chla <- as.numeric(sub(".*SCALE_CHLA=([-+]?[0-9]*\\.?[0-9]+);.*", "\\1", coeff_str))
            dark_chla  <- as.numeric(sub(".*DARK_CHLA=([-+]?[0-9]*\\.?[0-9]+);.*", "\\1", coeff_str))
          }
          nc_close(meta_file)
        }, silent = TRUE)
      }
    }
    
    bbp_var_name <- intersect(c("BBP700", "BBP_700"), names(profile_bio$var))
    if (length(bbp_var_name) > 0) {
      bbp_raw <- standardize_to_na(as.vector(ncvar_get(profile_bio, bbp_var_name[1])))
      if (length(bbp_raw) == length(pres)) bbp <- bbp_raw
      qc_name <- paste0(bbp_var_name[1], "_QC")
      if (qc_name %in% names(profile_bio$var)) {
        bbp_qc_v <- paste(ncvar_get(profile_bio, qc_name), collapse = "")
        parsed_qc <- unlist(strsplit(bbp_qc_v, ""))[1:length(pres)]
        bbp_qc <- ifelse(is.na(parsed_qc) | is.na(bbp), 9, as.numeric(parsed_qc))
      }
    }
    
    if ("DOXY" %in% names(profile_bio$var)) {
      doxy <- standardize_to_na(as.vector(ncvar_get(profile_bio, "DOXY")))
      if ("DOXY_QC" %in% names(profile_bio$var)) {
        doxy_qc_raw <- ncvar_get(profile_bio, "DOXY_QC")
        doxy_qc_str <- if (is.matrix(doxy_qc_raw)) apply(doxy_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(doxy_qc_raw, collapse = "")
        doxy_qc <- as.numeric(unlist(strsplit(doxy_qc_str, "")))[1:length(pres)]
        doxy_qc[is.na(doxy) | is.na(doxy_qc)] <- 9
      }
    }
    
    if ("DOXY_ADJUSTED" %in% names(profile_bio$var)) {
      doxy_adj <- standardize_to_na(as.vector(ncvar_get(profile_bio, "DOXY_ADJUSTED")))
      if ("DOXY_ADJUSTED_QC" %in% names(profile_bio$var)) {
        doxy_adj_qc_raw <- ncvar_get(profile_bio, "DOXY_ADJUSTED_QC")
        doxy_adj_qc_str <- if (is.matrix(doxy_adj_qc_raw)) apply(doxy_adj_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(doxy_adj_qc_raw, collapse = "")
        doxy_adj_qc <- as.numeric(unlist(strsplit(doxy_adj_qc_str, "")))[1:length(pres)]
        doxy_adj_qc[is.na(doxy_adj) | is.na(doxy_adj_qc)] <- 9
      }
    }
    
    data_sub <- data.frame(
      PRES = as.numeric(pres), PRES_QC = as.numeric(pres_qc),
      TEMP = as.numeric(temp), PSAL = as.numeric(sal),
      CHLA = chl, CHLA_QC = as.numeric(chl_qc), 
      CHLA_ADJUSTED = chl_adj, CHLA_ADJUSTED_QC = as.numeric(chl_adj_qc),
      FLUO_CHLA = fluo, FLUO_CHLA_QC = as.numeric(fluo_qc),
      BBP700 = bbp, BBP700_QC = as.numeric(bbp_qc),
      DOXY = doxy, DOXY_QC = as.numeric(doxy_qc),
      DOXY_ADJUSTED = doxy_adj, DOXY_ADJUSTED_QC = as.numeric(doxy_adj_qc)
    ) %>%
      filter(!is.na(PRES) & !is_argo_missing_val(PRES)) %>%
      arrange(PRES)
    
    if (nrow(data_sub) > 0) {
      data_sub$TIME <- time_val
      data_sub$LONGITUDE <- lon_val
      data_sub$LATITUDE  <- lat_val
      data_sub$CYCLE_NUMBER <- as.numeric(sub(".*_([0-9]+)$", "\\1", www))
      data_sub$float_num <- substr(www, 1, 7)
      data_sub$MLD <- MLD
      data_sub$DARK_FLUO  <- dark_chla
      data_sub$SCALE_FLUO <- scale_chla
      
      float_data <- rbind(float_data, data_sub)
    }
    
    nc_close(profile_bio)
    nc_close(profile_physic)
  }
  
  if (is.null(float_data) || nrow(float_data) == 0) {
    warning(paste("Skipping WMO", active_wmo, "- No valid bio, oxy, or bbp profiles collected."))
    next
  }
  
  # =========================================================================
  # STEP 1.5: BBP700 MANUAL QUALITY CONTROL TESTS
  # =========================================================================
  message("\nProcessing BBP700 Manual Quality Control Tests...")
  float_data$BBP700_QC_MANUAL <- ifelse(is.na(float_data$BBP700), 9, 1)
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    data_sub <- float_data[which(float_data$CYCLE_NUMBER == i), ]
    bbp <- data_sub$BBP700
    if (all(is.na(bbp))) next
    
    try({
      bbp_qc_test <- bbp_offset(bbp)
      for (j in seq_along(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    try({
      bbp_qc_test <- Range_test_bbp(bbp)$bbp_qc
      for (j in seq_along(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
  }
  
  float_data <- float_data %>%
    mutate(
      BBP700_QC_NUM      = suppressWarnings(as.numeric(as.character(BBP700_QC))),
      BBP700_FINAL_QC    = pmax(coalesce(BBP700_QC_NUM, 9), coalesce(BBP700_QC_MANUAL, 9), na.rm = TRUE),
      BBP700_FINAL_QC    = ifelse(is.na(BBP700) | PRES_QC %in% c("9", 9), 9, BBP700_FINAL_QC),
      BBP700_ADJUSTED_QC = BBP700_FINAL_QC,
      BBP700_ADJUSTED    = ifelse(BBP700_FINAL_QC == 9, NA_real_, BBP700),
      BBP700_FINAL       = ifelse(is.na(BBP700) | BBP700_FINAL_QC %in% c(4, 9), NA_real_, BBP700)
    ) %>%
    select(-BBP700_QC_NUM)
  
  # =========================================================================
  # STEP 3: ARGO CHLA COOKBOOK IN SITU DARK CORRECTION
  # =========================================================================
  float_data$FLUO_CHLA_RANGE <- ifelse(
    !is.na(float_data$FLUO_CHLA) & 
      !is_argo_missing_val(float_data$FLUO_CHLA) & 
      float_data$FLUO_CHLA > 0 & 
      float_data$FLUO_CHLA <= 4095,
    float_data$FLUO_CHLA,
    NA_real_
  )
  
  all_avail_cycs <- sort(unique(float_data$CYCLE_NUMBER[!is.na(float_data$CYCLE_NUMBER)]))
  cycle_deep_mins   <- c()
  qualifying_cycles <- c()
  
  for (cyc_i in all_avail_cycs) {
    if (length(cycle_deep_mins) >= 5) break
    
    data_sub <- float_data %>% filter(CYCLE_NUMBER == cyc_i)
    if (nrow(data_sub) == 0) next
    
    max_p <- safe_max(data_sub$PRES)
    if (is.na(max_p) || max_p < 900) next 
    
    deep_layer <- data_sub %>% 
      filter(!is.na(PRES) & 
               PRES >= 900 & 
               PRES <= 1000 & 
               !is.na(FLUO_CHLA_RANGE) & 
               !is_argo_missing_val(FLUO_CHLA_RANGE))
    
    if (nrow(deep_layer) < 3) next
    
    min_val <- safe_min(deep_layer$FLUO_CHLA_RANGE)
    if (!is.na(min_val) && is.finite(min_val)) {
      cycle_deep_mins   <- c(cycle_deep_mins, min_val)
      qualifying_cycles <- c(qualifying_cycles, cyc_i)
    }
  }
  
  MED <- if (length(cycle_deep_mins) > 0) median(cycle_deep_mins, na.rm = TRUE) else NA_real_
  if (is.na(MED) || is.infinite(MED) || MED < 0) {
    MED <- coalesce(float_data$DARK_FLUO[1], 0)
  }
  
  scale_factor <- coalesce(float_data$SCALE_FLUO[1], 1)
  
  float_data$CHLA_DARK              <- (float_data$FLUO_CHLA_RANGE - MED) * scale_factor
  float_data$CHLA_NoLUT             <- float_data$CHLA_DARK
  float_data$CHLA_ADJUSTED          <- coalesce(float_data$CHLA_ADJUSTED, float_data$CHLA_DARK)
  float_data$CHLA_FLUORESCENCE       <- (float_data$FLUO_CHLA - coalesce(float_data$DARK_FLUO, 0)) * scale_factor
  float_data$CHLA_FLUORESCENCE_FINAL <- (float_data$FLUO_CHLA - MED) * scale_factor
  
  # =========================================================================
  # STEP 4: PRE-NPQ AUTOMATED QC & DESPIKING
  # =========================================================================
  float_data$CHLA_QC_AUTO <- 1
  for (i in unique(float_data$CYCLE_NUMBER)) {
    idx <- which(float_data$CYCLE_NUMBER == i)
    chla_v <- float_data$CHLA_NoLUT[idx]
    n_obs  <- length(chla_v)
    qc_vec <- rep(1, n_obs)
    
    range_fail <- !is.na(chla_v) & (chla_v < -0.1 | chla_v > 50.0)
    qc_vec[range_fail] <- 4
    
    if (sum(!is.na(chla_v)) >= 5) {
      med_5 <- stats::runmed(ifelse(is.na(chla_v), 0, chla_v), k = 5, endrule = "median")
      res   <- chla_v - med_5
      valid_res <- res[!is.na(res)]
      if (length(valid_res) > 0) {
        q10 <- quantile(valid_res, probs = 0.10, na.rm = TRUE)
        negative_spike <- !is.na(res) & (res < 2 * q10)
        qc_vec[negative_spike] <- 4
      }
    }
    
    valid_v <- chla_v[!is.na(chla_v)]
    if (length(valid_v) > 5 && length(unique(valid_v)) == 1) {
      qc_vec[!is.na(chla_v)] <- 4
    }
    
    qc_vec[is.na(chla_v) | is_argo_missing_val(chla_v)] <- 9
    qc_vec[float_data$PRES_QC[idx] %in% c("9", 9)] <- 9
    float_data$CHLA_QC_AUTO[idx] <- qc_vec
  }
  
  # =========================================================================
  # STEP 4.5: SOLAR ANGLE & COOKBOOK NPQ
  # =========================================================================
  float_data <- tibble::rowid_to_column(float_data, "ID")
  float_data$SOLAR_ANGLE <- NA_real_
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    idx <- which(float_data$CYCLE_NUMBER == i)
    float_data$SOLAR_ANGLE[idx] <- getSunlightPosition(
      unique(float_data$TIME[idx])[1],
      unique(float_data$LATITUDE[idx])[1],
      unique(float_data$LONGITUDE[idx])[1]
    )$altitude * (180 / pi)
  }
  
  float_data$sun_test <- ifelse(float_data$SOLAR_ANGLE >= 0, "day", "night")
  
  float_data$NPQ_APPLIED <- FALSE
  for (i in unique(float_data$CYCLE_NUMBER)) {
    idx <- which(float_data$CYCLE_NUMBER == i)
    cyc_sun <- float_data$sun_test[idx][1]
    cyc_mld <- float_data$MLD[idx][1]
    
    if (!is.na(cyc_sun) && cyc_sun == "day" && !is.na(cyc_mld) && cyc_mld > 5) {
      npq_layer_depth <- 0.9 * cyc_mld
      mld_sub <- idx[which(float_data$PRES[idx] <= npq_layer_depth & 
                             !is.na(float_data$CHLA_NoLUT[idx]) & 
                             !is_argo_missing_val(float_data$CHLA_NoLUT[idx]) & 
                             float_data$CHLA_QC_AUTO[idx] < 4)]
      
      if (length(mld_sub) > 1) {
        max_chla_idx <- mld_sub[which.max(float_data$CHLA_NoLUT[mld_sub])]
        max_chla_val <- float_data$CHLA_NoLUT[max_chla_idx]
        pres_at_max  <- float_data$PRES[max_chla_idx]
        shallowest_p <- safe_min(float_data$PRES[mld_sub])
        
        if (!is.na(shallowest_p) && pres_at_max > shallowest_p && is.finite(max_chla_val)) {
          upward_idx <- idx[which(float_data$PRES[idx] < pres_at_max & float_data$CHLA_QC_AUTO[idx] < 4)]
          if (length(upward_idx) > 0) {
            float_data$CHLA_NoLUT[upward_idx] <- max_chla_val
            float_data$NPQ_APPLIED[upward_idx] <- TRUE
          }
        }
      }
    }
  }
  
  # --- STEP 5: SPATIAL LUT ALIGNMENT & RATIOS ---
  float_sf <- st_as_sf(float_data, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  nearest_indices <- st_nearest_feature(float_sf, lut_sf)
  
  matched_result <- float_data %>%
    mutate(
      nearest_indices = nearest_indices,
      physio_ratio = df_LUT$fluorescence_chlorophyll_ratio[nearest_indices],
      CHLA_LUT = CHLA_NoLUT / physio_ratio,
      matched_lat = df_LUT$y[nearest_indices],
      matched_lon = df_LUT$x[nearest_indices]
    )
  
  if (float_id %in% target_floats) {
    message(paste0("--> Bypassing CHLA calibration adjustments for Float: ", active_wmo))
    matched_result <- matched_result %>%
      mutate(CHLA_LUT = CHLA_NoLUT, physio_ratio = 1)
  }
  
  matched_result_v2 <- matched_result %>%
    mutate(
      SCALE_CHLA   = scale_chla,
      DARK_CHLA    = dark_chla,
      physio_ratio = ifelse(float_num %in% target_floats, 1, physio_ratio),
      CHLA_FINAL   = ifelse(float_num %in% target_floats, CHLA_NoLUT, CHLA_LUT)
    ) %>%
    group_by(CYCLE_NUMBER) %>%
    group_modify(~ {
      df_cyc <- .x
      chla_v <- df_cyc$CHLA_FINAL
      fluo_v <- df_cyc$CHLA_FLUORESCENCE
      n_obs  <- length(chla_v)
      
      qc_flags_chla <- df_cyc$CHLA_QC_AUTO
      qc_flags_chla[df_cyc$NPQ_APPLIED] <- 5
      qc_flags_chla[df_cyc$CHLA_QC %in% c("4", 4)] <- 4
      qc_flags_chla[is_argo_missing_val(chla_v)]   <- 9
      pres_qc_9_mask <- df_cyc$PRES_QC %in% c("9", 9)
      qc_flags_chla[pres_qc_9_mask] <- 9
      
      df_cyc$CHLA_ADJUSTED_QC <- ifelse(is.na(df_cyc$CHLA_ADJUSTED) | pres_qc_9_mask, 9, qc_flags_chla)
      
      qc_flags_fluo <- rep(1, n_obs)
      range_fail_fluo <- !is.na(fluo_v) & (fluo_v < -0.2 | fluo_v > 100)
      qc_flags_fluo[range_fail_fluo] <- 4
      qc_flags_fluo[is_argo_missing_val(fluo_v)] <- 9
      qc_flags_fluo[pres_qc_9_mask]              <- 9
      
      qc_flags_fluo_final <- ifelse(length(cycle_deep_mins) >= 5, 1, 2)
      qc_flags_fluo_final <- rep(qc_flags_fluo_final, n_obs)
      qc_flags_fluo_final[qc_flags_fluo >= 4] <- qc_flags_fluo[qc_flags_fluo >= 4]
      qc_flags_fluo_final[is_argo_missing_val(df_cyc$CHLA_FLUORESCENCE_FINAL)] <- 9
      qc_flags_fluo_final[pres_qc_9_mask] <- 9
      
      chla_v[qc_flags_chla == 9] <- NA_real_
      fluo_v[qc_flags_fluo == 9] <- NA_real_
      df_cyc$CHLA_FLUORESCENCE_FINAL[qc_flags_fluo_final == 9] <- NA_real_
      
      df_cyc$CHLA_FINAL                 <- chla_v
      df_cyc$CHLA_FINAL_QC              <- qc_flags_chla
      df_cyc$CHLA_FLUORESCENCE_QC       <- qc_flags_fluo
      df_cyc$CHLA_FLUORESCENCE_FINAL_QC <- qc_flags_fluo_final
      return(df_cyc)
    }) %>%
    ungroup() %>%
    select(-any_of(c(
      "SOLAR_ANGLE", "FLUO_CHLA_QC", 
      "match_dist_m", "DARK_FLUO", "SCALE_FLUO", 
      "CHLA_LUT", "base_qc", "NPQ_APPLIED", "CHLA_QC_AUTO"
    )))
  
  # ==============================================================================
  # MASTER BGC-ARGO DMQC OXYGEN ADJUSTMENT PIPELINE
  # ==============================================================================
  current_float_type <- float_type_dict$float_type[float_type_dict$float_ids == float_id]
  if (length(current_float_type) == 0) current_float_type <- "aoml_apex"
  DOXY_adj_err <- if (current_float_type == "aoml_navis") 5 else 2
  iprof         <- if (current_float_type == "aoml_navis") 1 else 2
  
  df_all <- matched_result_v2 %>%
    mutate(
      float_ids     = as.numeric(float_num),
      WMOID         = as.numeric(float_num),
      CYCLE_NUMBER  = as.numeric(CYCLE_NUMBER),
      DATE_TIME     = TIME,
      LATITUDE      = LATITUDE,
      LONGITUDE     = LONGITUDE,
      PRES_ADJUSTED = PRES,
      DOXY          = if ("DOXY" %in% names(.)) DOXY else NA_real_,
      DOXY_QC       = if ("DOXY_QC" %in% names(.)) DOXY_QC else 9,
      DOXY_ADJUSTED = if ("DOXY_ADJUSTED" %in% names(.)) DOXY_ADJUSTED else NA_real_,
      DOXY_ADJUSTED_QC = if ("DOXY_ADJUSTED_QC" %in% names(.)) DOXY_ADJUSTED_QC else 9,
      CHLA_ADJUSTED = coalesce(CHLA_ADJUSTED, CHLA_NoLUT)
    ) %>%
    left_join(calib_df, by = "float_ids")
  
  delta_days <- as.numeric(difftime(df_all$DATE_TIME, df_all$launch_date, units = "days"))
  gain       <- coalesce(df_all$doxy_slope, 1) * (1 + (coalesce(df_all$doxy_drift, 0) / 100) * (delta_days / 365.25))
  df_all$DOXY_CALCULATED_ADJ <- df_all$DOXY * gain
  
  df_all$DOXY_ADJUSTED    <- coalesce(df_all$DOXY_ADJUSTED, df_all$DOXY_CALCULATED_ADJ)
  df_all$DOXY_ADJUSTED_QC <- ifelse(is.na(df_all$DOXY_ADJUSTED) | df_all$PRES_QC %in% c("9", 9), 9, df_all$DOXY_QC)
  
  df_all <- df_all %>%
    mutate(
      DENSITY = oce::swRho(salinity = PSAL, temperature = TEMP, pressure = PRES_ADJUSTED),
      DOXY_ADJUSTED_ERROR = write_DOXY_adjusted_error(
        nc = NULL, profile_idx = iprof, err_mbar = DOXY_adj_err, 
        psal = PSAL, temp = TEMP, pres = PRES_ADJUSTED, 
        dens = DENSITY, DOXY_adj = DOXY_CALCULATED_ADJ
      )
    )
  
  # Vector-safe hook detection
  df_all <- df_all %>%
    group_by(CYCLE_NUMBER) %>%
    mutate(
      is_valid_pres = !is.na(PRES_ADJUSTED) & !is.na(DOXY_CALCULATED_ADJ) & !is_argo_missing_val(DOXY_CALCULATED_ADJ),
      max_pres_val  = if (any(is_valid_pres)) max(PRES_ADJUSTED[is_valid_pres], na.rm = TRUE) else NA_real_,
      sec_pres_val  = if (sum(is_valid_pres) >= 2) safe_max(PRES_ADJUSTED[is_valid_pres & PRES_ADJUSTED < max_pres_val]) else NA_real_,
      doxy_deepest  = ifelse(!is.na(max_pres_val) & PRES_ADJUSTED == max_pres_val, DOXY_CALCULATED_ADJ, NA_real_),
      doxy_second   = ifelse(!is.na(sec_pres_val) & PRES_ADJUSTED == sec_pres_val, DOXY_CALCULATED_ADJ, NA_real_),
      deepest_diff  = abs(safe_max(doxy_deepest) - safe_max(doxy_second)),
      HOOK_FLAG_4   = !is.na(PRES_ADJUSTED) & !is.na(max_pres_val) & (PRES_ADJUSTED == max_pres_val) & 
        (!is.na(deepest_diff) & deepest_diff > 2)
    ) %>%
    ungroup()
  
  df_all <- df_all %>%
    mutate(DOXY_FINAL = coalesce(DOXY_CALCULATED_ADJ, DOXY_ADJUSTED, DOXY)) %>%
    group_by(CYCLE_NUMBER) %>%
    group_modify(~ {
      df_cyc <- .x
      v <- df_cyc$DOXY_FINAL
      n_obs <- length(v)
      qc_flags <- rep(1, n_obs)
      
      range_fail <- !is.na(v) & (v < -5 | v > 600)
      qc_flags[range_fail] <- 4
      
      if (n_obs >= 3) {
        v_prev <- c(NA, v[1:(n_obs - 1)])
        v_next <- c(v[2:n_obs], NA)
        spike_val <- abs(v - 0.5 * (v_prev + v_next)) - 0.5 * abs(v_next - v_prev)
        spike_fail <- !is.na(spike_val) & (spike_val > 50)
        qc_flags[spike_fail] <- 4
      }
      
      valid_v <- v[!is.na(v) & !is_argo_missing_val(v)]
      if (length(valid_v) > 1 && length(unique(valid_v)) == 1) {
        qc_flags[!is.na(v)] <- 4
      }
      
      if ("HOOK_FLAG_4" %in% names(df_cyc)) {
        qc_flags[df_cyc$HOOK_FLAG_4] <- 4
      }
      
      qc_flags[is_argo_missing_val(v)] <- 9
      pres_qc_9_mask_doxy <- df_cyc$PRES_QC %in% c("9", 9)
      qc_flags[pres_qc_9_mask_doxy] <- 9
      
      v[qc_flags == 9] <- NA_real_
      df_cyc$DOXY_FINAL <- v
      df_cyc$DOXY_FINAL_QC <- qc_flags
      return(df_cyc)
    }) %>%
    ungroup() %>%
    select(-any_of(c("HOOK_FLAG_4", "is_valid_pres", "max_pres_val", "sec_pres_val", "doxy_deepest", "doxy_second", "deepest_diff")))
  
  # ==============================================================================
  # STEP 6: EXPORT COMBINED CSV RETAINING ALL R, A & D MODES
  # ==============================================================================
  combo_output_path <- file.path(combo_output_dir, sprintf("CHLA_BBP_DOXY_%d.csv", float_id))
  
  df_export_prep <- df_all %>%
    mutate(
      float_num = coalesce(as.character(WMOID), as.character(float_num)),
      across(where(is.logical), as.numeric),
      DOXY_ADJUSTED = coalesce(DOXY_ADJUSTED, DOXY_CALCULATED_ADJ),
      DATE_TIME     = format(DATE_TIME, format = "%Y-%m-%d %H:%M:%S")
    )
  
  df_export_prep <- df_export_prep %>%
    mutate(
      override_condition = (
        (float_num == "4903624") | 
          (float_num == "2904010" & CYCLE_NUMBER >= 49) | 
          (float_num == "2904011" & CYCLE_NUMBER >= 24)
      ),
      CHLA                          = ifelse(override_condition, -99999.0, CHLA),
      CHLA_QC                       = ifelse(override_condition, "9", as.character(CHLA_QC)),
      CHLA_ADJUSTED                 = ifelse(override_condition, -99999.0, CHLA_ADJUSTED),
      CHLA_ADJUSTED_QC              = ifelse(override_condition, "9", as.character(CHLA_ADJUSTED_QC)),
      CHLA_FINAL                    = ifelse(override_condition, -99999.0, CHLA_FINAL),
      CHLA_FINAL_QC                 = ifelse(override_condition, "9", as.character(CHLA_FINAL_QC)),
      
      BBP700                        = ifelse(override_condition, -99999.0, BBP700),
      BBP700_QC                     = ifelse(override_condition, "9", as.character(BBP700_QC)),
      BBP700_ADJUSTED               = ifelse(override_condition, -99999.0, BBP700_ADJUSTED),
      BBP700_ADJUSTED_QC            = ifelse(override_condition, "9", as.character(BBP700_ADJUSTED_QC)),
      BBP700_FINAL                  = ifelse(override_condition, -99999.0, BBP700_FINAL),
      BBP700_FINAL_QC               = ifelse(override_condition, "9", as.character(BBP700_FINAL_QC)),
      
      CHLA_FLUORESCENCE             = ifelse(override_condition, -99999.0, CHLA_FLUORESCENCE),
      CHLA_FLUORESCENCE_QC          = ifelse(override_condition, "9", as.character(CHLA_FLUORESCENCE_QC)),
      CHLA_FLUORESCENCE_FINAL       = ifelse(override_condition, -99999.0, CHLA_FLUORESCENCE_FINAL),
      CHLA_FLUORESCENCE_FINAL_QC    = ifelse(override_condition, "9", as.character(CHLA_FLUORESCENCE_FINAL_QC)),
      
      CHLA_NoLUT                    = ifelse(override_condition, -99999.0, CHLA_NoLUT),
      SCALE_CHLA                    = ifelse(override_condition, -99999.0, SCALE_CHLA),
      DARK_CHLA                     = ifelse(override_condition, -99999.0, DARK_CHLA)
    ) %>%
    select(-override_condition)
  
  df_combo_export <- df_export_prep %>%
    select(-any_of(c(
      "WMOID", "ID", "BBP700_QC_MANUAL",
      "DOXY_CALCULATED_ADJ", "float_ids", "TIME",
      "launch_date", "lauch_date", "LAUNCH_DATE", "sun_test", "SOLAR_ANGLE", 
      "FLUO_CHLA_RANGE", "CHLA_QC_AUTO"
    ))) %>%
    rename_with(toupper) %>%
    relocate(
      FLOAT_NUM, CYCLE_NUMBER, DATE_TIME, LATITUDE, LONGITUDE, PRES, PRES_QC,
      DOXY, DOXY_QC, 
      DOXY_ADJUSTED, DOXY_ADJUSTED_QC, 
      DOXY_FINAL, DOXY_FINAL_QC, DOXY_ADJUSTED_ERROR,
      CHLA, CHLA_QC, 
      CHLA_ADJUSTED, CHLA_ADJUSTED_QC, 
      CHLA_FINAL, CHLA_FINAL_QC,
      BBP700, BBP700_QC, 
      BBP700_ADJUSTED, BBP700_ADJUSTED_QC, 
      BBP700_FINAL, BBP700_FINAL_QC,
      CHLA_FLUORESCENCE, CHLA_FLUORESCENCE_QC,
      CHLA_FLUORESCENCE_FINAL, CHLA_FLUORESCENCE_FINAL_QC
    ) %>%
    filter(!is.na(PRES) & !is_argo_missing_val(PRES))
  
  qc_cols <- grep("_QC$", names(df_combo_export), value = TRUE)
  meta_cols <- c(qc_cols, "FLOAT_NUM", "CYCLE_NUMBER", "DATE_TIME", "LATITUDE", "LONGITUDE", "PRES", "PRES_ADJUSTED")
  val_cols <- setdiff(names(df_combo_export), meta_cols)
  
  df_combo_export <- df_combo_export %>%
    mutate(
      across(all_of(qc_cols), ~ ifelse(is.na(.) | as.character(.) %in% c("NA", "", "99999", "-99999", "0"), "9", as.character(.)))
    )
  
  for (qc_col in qc_cols) {
    param_stem <- sub("_QC$", "", qc_col)
    if (param_stem %in% names(df_combo_export) && !param_stem %in% c("PRES", "PRES_ADJUSTED")) {
      missing_idx <- is_argo_missing_val(df_combo_export[[param_stem]])
      df_combo_export[[param_stem]][missing_idx] <- -99999.0
      df_combo_export[[qc_col]][missing_idx]     <- "9"
    }
  }
  
  df_combo_export <- df_combo_export %>%
    mutate(
      across(all_of(val_cols), ~ {
        num_val <- suppressWarnings(as.numeric(as.character(.)))
        ifelse(is_argo_missing_val(num_val), -99999.0, num_val)
      })
    )
  
  # --------------------------------------------------------------------------
  # FILTER OUT LEVELS CONSISTENTLY MISSING ACROSS ALL BGC PARAMETERS & MODES
  # --------------------------------------------------------------------------
  bgc_pairs <- list(
    list(val = "DOXY",                    qc = "DOXY_QC"),
    list(val = "DOXY_ADJUSTED",           qc = "DOXY_ADJUSTED_QC"),
    list(val = "DOXY_FINAL",              qc = "DOXY_FINAL_QC"),
    list(val = "CHLA",                    qc = "CHLA_QC"),
    list(val = "CHLA_ADJUSTED",           qc = "CHLA_ADJUSTED_QC"),
    list(val = "CHLA_FINAL",              qc = "CHLA_FINAL_QC"),
    list(val = "BBP700",                  qc = "BBP700_QC"),
    list(val = "BBP700_ADJUSTED",         qc = "BBP700_ADJUSTED_QC"),
    list(val = "BBP700_FINAL",            qc = "BBP700_FINAL_QC"),
    list(val = "CHLA_FLUORESCENCE",       qc = "CHLA_FLUORESCENCE_QC"),
    list(val = "CHLA_FLUORESCENCE_FINAL", qc = "CHLA_FLUORESCENCE_FINAL_QC")
  )
  
  existing_pairs <- Filter(function(p) p$val %in% names(df_combo_export), bgc_pairs)
  
  is_missing_pair <- function(p) {
    v_miss <- is_argo_missing_val(df_combo_export[[p$val]])
    if (p$qc %in% names(df_combo_export)) {
      q_miss <- is.na(df_combo_export[[p$qc]]) | as.character(df_combo_export[[p$qc]]) %in% c("9", "-99999", "99999", "", "NA")
      return(v_miss | q_miss)
    }
    return(v_miss)
  }
  
  all_missing_matrix <- sapply(existing_pairs, is_missing_pair)
  row_all_bgc_missing <- apply(all_missing_matrix, 1, all)
  
  n_rows_before <- nrow(df_combo_export)
  df_combo_export <- df_combo_export %>% filter(!row_all_bgc_missing)
  n_rows_removed <- n_rows_before - nrow(df_combo_export)
  
  message(sprintf("[%s] Removed %d depth levels consistently missing across all BGC parameters and modes (%d levels retained).", 
                  active_wmo, n_rows_removed, nrow(df_combo_export)))
  
  write.csv(df_combo_export, combo_output_path, row.names = FALSE)
  message(sprintf("\n>>> SUCCESSFULLY EXPORTED ALL MODES COMBINED CSV -> %s", combo_output_path))
  
  # ==============================================================================
  # STEP 6.5: AUDIT AND REPORT INCONSISTENCIES WITHIN AND ACROSS VARIABLES
  # ==============================================================================
  audit_missing_mismatches(df_combo_export, active_wmo, dir_qc_summary)
  
  # ==============================================================================
  # STEP 7: INDIVIDUAL PER-FLOAT PER-CYCLE PROFILE PLOTS (R, A, D MODES)
  # ==============================================================================
  message(sprintf("Generating per-cycle profile plots for Float %s...", active_wmo))
  
  cycles_to_plot <- sort(unique(df_combo_export$CYCLE_NUMBER[!is.na(df_combo_export$CYCLE_NUMBER)]))
  
  for (cyc in cycles_to_plot) {
    df_cyc <- df_combo_export %>% filter(CYCLE_NUMBER == cyc)
    if (nrow(df_cyc) == 0) next
    
    bottle_cyc_df <- NULL
    if (exists("Bottle_Vali") && nrow(Bottle_Vali) > 0) {
      b_float_col <- intersect(c("Float_num", "FLOAT_NUM", "float_ids", "WMO"), names(Bottle_Vali))
      if (length(b_float_col) > 0) {
        b_match <- Bottle_Vali[as.character(Bottle_Vali[[b_float_col[1]]]) == as.character(active_wmo), ]
        b_cyc_col <- intersect(c("CYCLE_NUMBER", "Cycle", "CYCLE", "cycle"), names(b_match))
        if (length(b_cyc_col) > 0) {
          b_match <- b_match[b_match[[b_cyc_col[1]]] == cyc, ]
        }
        if (nrow(b_match) > 0 && "DOXY_bottle" %in% names(b_match) && "PRES" %in% names(b_match)) {
          bottle_cyc_df <- b_match %>%
            filter(!is_argo_missing_val(DOXY_bottle) & !is_argo_missing_val(PRES)) %>%
            transmute(PRES = as.numeric(PRES), Value = as.numeric(DOXY_bottle))
        }
      }
    }
    
    doxy_plot_file <- file.path(dir_doxy_plots, sprintf("Float_%s_Cycle_%03d_DOXY_Profile.png", active_wmo, cyc))
    generate_cycle_profile_plot(
      df_cyc = df_cyc,
      var_cols = c("DOXY", "DOXY_ADJUSTED", "DOXY_FINAL"),
      var_modes = c("R-mode", "A-mode", "D-mode"),
      var_title = "DOXY",
      unit_expr = expression(Oxygen ~ (mu * "mol/kg")),
      out_file = doxy_plot_file,
      float_id = active_wmo,
      cyc_num = cyc,
      bottle_df = bottle_cyc_df
    )
    
    chla_plot_file <- file.path(dir_chla_plots, sprintf("Float_%s_Cycle_%03d_CHLA_Profile.png", active_wmo, cyc))
    generate_cycle_profile_plot(
      df_cyc = df_cyc,
      var_cols = c("CHLA", "CHLA_ADJUSTED", "CHLA_FINAL"),
      var_modes = c("R-mode", "A-mode", "D-mode"),
      var_title = "CHLA",
      unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)),
      out_file = chla_plot_file,
      float_id = active_wmo,
      cyc_num = cyc
    )
    
    bbp_plot_file <- file.path(dir_bbp_plots, sprintf("Float_%s_Cycle_%03d_BBP700_Profile.png", active_wmo, cyc))
    generate_cycle_profile_plot(
      df_cyc = df_cyc,
      var_cols = c("BBP700", "BBP700_ADJUSTED", "BBP700_FINAL"),
      var_modes = c("R-mode", "A-mode", "D-mode"),
      var_title = "BBP700",
      unit_expr = expression(b[bp](700) ~ (m^{-1})),
      out_file = bbp_plot_file,
      float_id = active_wmo,
      cyc_num = cyc
    )
    
    fluo_plot_file <- file.path(dir_fluo_plots, sprintf("Float_%s_Cycle_%03d_CHLA_FLUO_Profile.png", active_wmo, cyc))
    generate_cycle_profile_plot(
      df_cyc = df_cyc,
      var_cols = c("CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_FINAL"),
      var_modes = c("R-mode", "D-mode"),
      var_title = "CHLA FLUORESCENCE",
      unit_expr = "Fluorescence (RU)",
      out_file = fluo_plot_file,
      float_id = active_wmo,
      cyc_num = cyc
    )
  }
  
  # ==============================================================================
  # STEP 8: SECTION CONTOUR & QC SUMMARY PLOTS (MATCHING PYTHON SCIPY/MATPLOTLIB)
  # ==============================================================================
  df_summary_base <- df_combo_export %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES         = as.numeric(PRES)
    )
  
  create_interpolated_grid_r <- function(clean_df, val_col, grid_nx = 100, grid_ny = 100) {
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
    
    if (nrow(grouped) < 5 ||
        n_distinct(grouped$CYCLE_NUMBER) < 2 ||
        n_distinct(grouped$PRES) < 2) {
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
        nearest_idx <- which.min(dist_sq)
        z_mat[r, c] <- grouped$z_val[nearest_idx]
      }
    }
    
    grid_df <- expand.grid(x = xo_seq, y = yo_seq)
    grid_df$z <- as.vector(z_mat)
    
    if (all(is.na(grid_df$z)) || diff(range(grid_df$z, na.rm = TRUE)) == 0) {
      return(NULL)
    }
    return(grid_df)
  }
  
  generate_section_qc_plot <- function(df_base, var_name, qc_var_name, title_prefix, file_suffix, 
                                       col_palette = "viridis", unit_expr = "", is_chla = FALSE) {
    if (!var_name %in% names(df_base)) return(NULL)
    
    df_plot <- df_base
    
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
    max_pres <- safe_max(pres_all)
    if (is.na(max_pres) || max_pres <= 0) max_pres <- 2000.0
    
    grid_df <- if (nrow(df_valid) > 0) create_interpolated_grid_r(df_valid, var_name) else NULL
    
    # --------------------------------------------------------------------------
    # Upper Plot: Section Contour + Python-Identical Missing Overlay ('x' markers)
    # --------------------------------------------------------------------------
    p_data <- ggplot()
    
    if (!is.null(grid_df)) {
      if (is_chla) {
        z_min <- min(grid_df$z, na.rm = TRUE)
        z_max <- max(grid_df$z, na.rm = TRUE)
        base_breaks <- c(-0.1, 0.01, 0.03, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.2, 2.0, 5.0)
        chla_breaks <- unique(sort(c(min(z_min, -0.1), base_breaks[base_breaks > z_min & base_breaks < z_max], max(z_max, 5.0))))
        
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
                   size = 1.2, alpha = 0.8, na.rm = TRUE) +
        scale_color_viridis_c(option = col_palette, name = unit_expr)
    }
    
    if (nrow(df_missing) > 0) {
      p_data <- p_data + 
        geom_point(data = df_missing, aes(x = CYCLE_NUMBER, y = PRES, shape = "Missing Data (-99999)"), 
                   color = "#7f7f7f", size = 1.2, stroke = 0.7, alpha = 0.35, na.rm = TRUE) +
        scale_shape_manual(name = "", values = c("Missing Data (-99999)" = 4))
    }
    
    p_data <- p_data + 
      scale_y_reverse(limits = c(max_pres, 0), expand = c(0, 0)) +
      labs(
        title = paste0("Float ", float_id, " - ", title_prefix, " Section Profile"),
        x = NULL,
        y = "Pressure (dbar)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.major = element_line(linetype = "dashed", color = "grey80", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA)
      )
    
    # --------------------------------------------------------------------------
    # Lower Plot: QC Flags Scatter Map (Python QC_COLOR_MAP)
    # --------------------------------------------------------------------------
    df_qc_plot <- df_plot %>% 
      filter(!is_argo_missing_val(PRES) & !is_argo_missing_val(CYCLE_NUMBER))
    
    u_qcs <- unique(df_qc_plot$QC_VAL)
    unique_qcs <- u_qcs[order(ifelse(grepl("^[0-9]+$", u_qcs), as.numeric(u_qcs), 99))]
    
    p_qc <- ggplot()
    
    for (q_flag in unique_qcs) {
      sub_df <- df_qc_plot %>% filter(QC_VAL == q_flag)
      p_qc <- p_qc + 
        geom_point(data = sub_df, aes(x = CYCLE_NUMBER, y = PRES, color = QC_VAL), 
                   size = 1.2, alpha = 0.7, shape = 16, na.rm = TRUE)
    }
    
    p_qc <- p_qc + 
      scale_y_reverse(limits = c(max_pres, 0), expand = c(0, 0)) +
      scale_color_manual(values = qc_color_scale, limits = unique_qcs, name = "QC Flag") +
      labs(
        title = paste0(title_prefix, " QC Flags Section Map"),
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
    
    p_combined <- (p_data / p_qc) + 
      plot_layout(heights = c(1, 1), guides = "keep") +
      plot_annotation(
        title = paste0("BGC-Argo Float ", float_id, " - ", title_prefix, " Section DMQC Summary (Verified BD Output)"),
        subtitle = paste0("Generated on: ", today_str),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
    
    out_file <- file.path(dir_qc_summary, sprintf("Float_%d_%s_Section_Summary.png", float_id, file_suffix))
    ggsave(out_file, plot = p_combined, width = 12, height = 8, dpi = 150)
    message(sprintf("[%s] Saved Aligned BD-Verified QC Summary Plot -> %s", active_wmo, out_file))
  }
  
  generate_section_qc_plot(df_summary_base, "DOXY", "DOXY_QC", "DOXY R_mode", "DOXY_R_mode", 
                           col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "DOXY_ADJUSTED", "DOXY_ADJUSTED_QC", "DOXY A_mode", "DOXY_A_mode", 
                           col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "DOXY_FINAL", "DOXY_FINAL_QC", "DOXY D_mode", "DOXY_D_mode", 
                           col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  
  generate_section_qc_plot(df_summary_base, "BBP700", "BBP700_QC", "BBP700 R_mode", "BBP700_R_mode", 
                           col_palette = "turbo", unit_expr = expression(b[bp](700) ~ (m^{-1})))
  generate_section_qc_plot(df_summary_base, "BBP700_ADJUSTED", "BBP700_ADJUSTED_QC", "BBP700 A_mode", "BBP700_A_mode", 
                           col_palette = "turbo", unit_expr = expression(b[bp](700) ~ (m^{-1})))
  generate_section_qc_plot(df_summary_base, "BBP700_FINAL", "BBP700_FINAL_QC", "BBP700 D_mode", "BBP700_D_mode", 
                           col_palette = "turbo", unit_expr = expression(b[bp](700) ~ (m^{-1})))
  
  generate_section_qc_plot(df_summary_base, "CHLA", "CHLA_QC", "CHLA R_mode", "CHLA_R_mode", 
                           col_palette = "turbo", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  generate_section_qc_plot(df_summary_base, "CHLA_ADJUSTED", "CHLA_ADJUSTED_QC", "CHLA A_mode", "CHLA_A_mode", 
                           col_palette = "turbo", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  generate_section_qc_plot(df_summary_base, "CHLA_FINAL", "CHLA_FINAL_QC", "CHLA D_mode", "CHLA_D_mode", 
                           col_palette = "turbo", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  
  generate_section_qc_plot(df_summary_base, "CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_QC", "CHLA_FLUORESCENCE R_mode", "CHLA_FLUORESCENCE_R_mode", 
                           col_palette = "mako", unit_expr = "Fluorescence (RU)")
  generate_section_qc_plot(df_summary_base, "CHLA_FLUORESCENCE_FINAL", "CHLA_FLUORESCENCE_FINAL_QC", "CHLA_FLUORESCENCE D_mode", "CHLA_FLUORESCENCE_D_mode", 
                           col_palette = "mako", unit_expr = "Fluorescence (RU)")
  
  message("\n=== COMPLETE DMQC PROCESSING FOR FLOAT: ", float_id, " ===")
}