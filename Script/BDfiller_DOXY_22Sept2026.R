#!/usr/bin/env Rscript
# ==============================================================================
# SECTION 1: Library Imports
# ==============================================================================
cat("\014")
rm(list = ls())

library(ncdf4)
library(gsw)
library(stringr)

# Set print options equivalent to np.set_printoptions(threshold=np.inf)
options(max.print = .Machine$integer.max)

# Load profile-level CSV file containing oxygen data
csv_dir  <- "C:/Users/Jennifer.McWhorter/Documents/GitHub/DMQC_CHLA_DOXY/Output/COMBO_CHLA_DOXY_DMODE"
csv_file <- "CHLA_DOXY_6999992.csv"  # Appended .csv extension if needed
csv_path <- file.path(csv_dir, csv_file)

if (!file.exists(csv_path)) {
  # Fallback check without extension if file lacks one
  csv_path_alt <- file.path(csv_dir, "CHLA_DOXY_6999992")
  if (file.exists(csv_path_alt)) {
    csv_path <- csv_path_alt
  } else {
    stop(sprintf("CSV file not found at: %s", csv_path))
  }
}

csv_data <- read.csv(csv_path, stringsAsFactors = FALSE)

# Check required columns exist (added DOXY_SLOPE and DOXY_DRIFT)
req_cols <- c("FLOAT_NUM", "CYCLE_NUMBER", "PRES", "DOXY", "DOXY_FINAL", "DOXY_FINAL_QC", "DOXY_SLOPE", "DOXY_DRIFT")
if (!all(req_cols %in% names(csv_data))) {
  stop(sprintf("CSV file must contain the following columns: %s", paste(req_cols, collapse = ", ")))
}

# Force CYCLE_NUMBER to integer to avoid type mismatching during matching
csv_data$CYCLE_NUMBER <- as.integer(csv_data$CYCLE_NUMBER)

# Get distinct float list from CSV
unique_floats <- unique(csv_data$FLOAT_NUM)

# ==============================================================================
# SECTION 2: Helper Functions
# ==============================================================================

get_juld <- function(filename, profile_idx) {
  nc <- nc_open(filename)
  on.exit(nc_close(nc))
  juld_var <- ncvar_get(nc, "JULD")
  return(juld_var[profile_idx])
}

get_iprof_phys <- function(pres_raw, pres_bgc, target_iprof) {
  iprof_phys <- -1
  num_cols <- if (is.matrix(pres_raw)) ncol(pres_raw) else 1
  target_pres <- if (is.matrix(pres_bgc)) pres_bgc[, target_iprof] else pres_bgc
  
  for (col in seq_len(num_cols)) {
    current_pres <- if (is.matrix(pres_raw)) pres_raw[, col] else pres_raw
    diff <- max(abs(current_pres - target_pres), na.rm = TRUE)
    if (diff < 0.1) {
      iprof_phys <- col
      break
    }
  }
  
  if (iprof_phys < 0) {
    stop("Matching PRES values not found")
  } else {
    cat(sprintf("Using profile %d of physical file to determine density.\n", iprof_phys))
    return(iprof_phys)
  }
}

create_working_bd_file <- function(filename, dest_dir) {
  base_name <- basename(filename)
  bd_name <- gsub("^BR", "BD", base_name)
  w_filename <- file.path(dest_dir, paste0("w_", bd_name))
  file.copy(filename, w_filename, overwrite = TRUE)
  return(w_filename)
}

get_phys_filename <- function(bgc_filename, base_dir) {
  base_name <- basename(bgc_filename)
  core_phys_name <- sub("^B", "", base_name)
  
  d_name <- sub("^R", "D", core_phys_name)
  phys_filename <- file.path(base_dir, d_name)
  
  if (!file.exists(phys_filename)) {
    phys_filename <- file.path(base_dir, "D", d_name)
  }
  
  if (!file.exists(phys_filename)) {
    r_name <- sub("^D", "R", core_phys_name)
    phys_filename <- file.path(base_dir, r_name)
    
    if (!file.exists(phys_filename)) {
      phys_filename <- file.path(base_dir, "R", r_name)
    }
    
    if (file.exists(phys_filename)) {
      cat(sprintf("Using phys R file: %s\n", phys_filename))
    }
  }
  
  if (!file.exists(phys_filename)) {
    stop(sprintf("No corresponding phys file found for %s (looked for %s and %s in %s)", 
                 bgc_filename, d_name, sub("^D", "R", core_phys_name), base_dir))
  }
  
  return(phys_filename)
}

get_phys_raw_pres <- function(phys_filename) {
  nc <- nc_open(phys_filename)
  on.exit(nc_close(nc))
  return(ncvar_get(nc, "PRES"))
}

get_dens <- function(phys_filename, verbose = FALSE) {
  nc <- nc_open(phys_filename)
  on.exit(nc_close(nc))
  
  var_names <- names(nc$var)
  num_adj <- 0
  
  temp <- if ("TEMP_ADJUSTED" %in% var_names) { num_adj <- num_adj + 1; ncvar_get(nc, "TEMP_ADJUSTED") } 
  else { cat("TEMP_ADJUSTED not found, using TEMP\n"); ncvar_get(nc, "TEMP") }
  
  psal <- if ("PSAL_ADJUSTED" %in% var_names) { num_adj <- num_adj + 1; ncvar_get(nc, "PSAL_ADJUSTED") } 
  else { cat("PSAL_ADJUSTED not found, using PSAL\n"); ncvar_get(nc, "PSAL") }
  
  pres <- if ("PRES" %in% var_names) { num_adj <- num_adj + 1; ncvar_get(nc, "PRES") } 
  else { cat("PRES not found, using PRES\n"); ncvar_get(nc, "PRES") }
  
  if (num_adj == 3 && verbose) {
    cat("Using ADJUSTED values of p,T,S to determine density\n")
  } else if (num_adj == 0) {
    cat("Using RAW values of p,T,S to determine density\n")
  }
  
  dens <- gsw_rho_t_exact(SA = psal, t = temp, p = pres)
  return(list(dens = dens, psal = psal, temp = temp))
}

update_history <- function(nc, history_list, profile_idx) {
  hix <- nc$dim$N_HISTORY$len
  
  for (var_name in names(history_list)) {
    val <- history_list[[var_name]]
    if (var_name %in% names(nc$var)) {
      ncvar_put(nc, var_name, val, start = c(1, profile_idx, hix), count = c(nchar(val), 1, 1))
    }
  }
}

write_history <- function(nc, profile_idx, inst, ref, comment_op) {
  utc_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  ncatt_put(nc, 0, "history", paste(utc_time, "creation"), prec = "text")
  ncatt_put(nc, 0, "comment_dmqc_operator", comment_op, prec = "text")
  
  hist_data <- list(
    HISTORY_INSTITUTION = inst,
    HISTORY_STEP = "ARSQ",
    HISTORY_SOFTWARE = "SAGE",
    HISTORY_SOFTWARE_RELEASE = "2024",
    HISTORY_REFERENCE = ref,
    HISTORY_DATE = format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC"),
    HISTORY_ACTION = "IP",
    HISTORY_PARAMETER = "DOXY"
  )
  
  update_history(nc, hist_data, profile_idx)
  
  utc_compact <- format(Sys.time(), "%Y%m%d%H%M%S", tz = "UTC")
  if ("DATE_CREATION" %in% names(nc$var)) ncvar_put(nc, "DATE_CREATION", utc_compact)
  if ("DATE_UPDATE" %in% names(nc$var)) ncvar_put(nc, "DATE_UPDATE", utc_compact)
}

write_parameter_data_mode <- function(nc, profile_idx, mode_char = "D") {
  params <- ncvar_get(nc, "STATION_PARAMETERS", collapse_degen = FALSE)
  pdm_array <- ncvar_get(nc, "PARAMETER_DATA_MODE", collapse_degen = FALSE)
  
  if (is.array(params) && length(dim(params)) == 3) {
    prof_params <- params[, , profile_idx]
    
    if (!is.matrix(prof_params)) {
      prof_params <- matrix(prof_params, nrow = dim(params)[1], ncol = dim(params)[2])
    }
    
    param_strings <- apply(prof_params, 2, function(x) paste(x, collapse = ""))
    DOXY_idx <- which(startsWith(trimws(param_strings), "DOXY"))
    
    if (length(DOXY_idx) > 0) {
      if (length(dim(pdm_array)) == 2) {
        pdm_array[DOXY_idx, profile_idx] <- mode_char
      } else {
        pdm_array[DOXY_idx] <- mode_char
      }
    }
  }
  
  ncvar_put(nc, "PARAMETER_DATA_MODE", pdm_array)
  
  data_mode <- ncvar_get(nc, "DATA_MODE", collapse_degen = FALSE)
  data_mode[profile_idx] <- "D"
  ncvar_put(nc, "DATA_MODE", data_mode)
}

# Update Slope and Drift calibration coefficients from CSV
write_DOXY_slope_drift <- function(nc, profile_idx, float_df, target_cycle) {
  cycle_df <- float_df[as.integer(float_df$CYCLE_NUMBER) == as.integer(target_cycle), ]
  
  if (nrow(cycle_df) == 0) return()
  
  # Extract slope and drift for this cycle (taking first row entry)
  slope_val <- cycle_df$DOXY_SLOPE[1]
  drift_val <- cycle_df$DOXY_DRIFT[1]
  
  if (is.na(slope_val) && is.na(drift_val)) return()
  
  # Format calibration string for Argo NetCDF standards
  calib_str <- sprintf("m=%s, d=%s", 
                       if (is.na(slope_val)) "1" else as.character(slope_val), 
                       if (is.na(drift_val)) "0" else as.character(drift_val))
  
  var_names <- names(nc$var)
  
  # 1. Update STANDARD Argo SCIENTIFIC_CALIB_COEFFICIENT if present
  if ("SCIENTIFIC_CALIB_COEFFICIENT" %in% var_names) {
    params <- ncvar_get(nc, "STATION_PARAMETERS", collapse_degen = FALSE)
    if (is.array(params) && length(dim(params)) == 3) {
      prof_params <- params[, , profile_idx]
      if (!is.matrix(prof_params)) {
        prof_params <- matrix(prof_params, nrow = dim(params)[1], ncol = dim(params)[2])
      }
      param_strings <- apply(prof_params, 2, function(x) paste(x, collapse = ""))
      DOXY_idx <- which(startsWith(trimws(param_strings), "DOXY"))
      
      if (length(DOXY_idx) > 0) {
        calib_arr <- ncvar_get(nc, "SCIENTIFIC_CALIB_COEFFICIENT", collapse_degen = FALSE)
        # Pad string to required NetCDF character dimension length
        char_len <- dim(calib_arr)[1]
        padded_str <- str_pad(calib_str, char_len, side = "right")
        
        if (length(dim(calib_arr)) == 3) {
          ncvar_put(nc, "SCIENTIFIC_CALIB_COEFFICIENT", padded_str, 
                    start = c(1, DOXY_idx, profile_idx), count = c(char_len, 1, 1))
        }
      }
    }
  }
  
  # 2. Direct write if custom numeric/character DOXY_SLOPE / DOXY_DRIFT variables exist
  if ("DOXY_SLOPE" %in% var_names && !is.na(slope_val)) {
    ncvar_put(nc, "DOXY_SLOPE", slope_val)
  }
  if ("DOXY_DRIFT" %in% var_names && !is.na(drift_val)) {
    ncvar_put(nc, "DOXY_DRIFT", drift_val)
  }
  
  cat(sprintf("Updated DOXY Slope/Drift for Cycle %d: %s\n", target_cycle, calib_str))
}

# Direct write of DOXY, DOXY_ADJUSTED, and DOXY_ADJUSTED_QC using matched CSV pressure levels
write_DOXY_from_csv <- function(nc, profile_idx, float_df, target_cycle) {
  var_names <- names(nc$var)
  
  pres_nc_full <- if ("PRES" %in% var_names) {
    ncvar_get(nc, "PRES")
  } else {
    ncvar_get(nc, "PRES")
  }
  
  pres_nc <- if (is.matrix(pres_nc_full)) pres_nc_full[, profile_idx] else pres_nc_full
  
  # Ensure target_cycle matching is type safe
  cycle_df <- float_df[as.integer(float_df$CYCLE_NUMBER) == as.integer(target_cycle), ]
  
  if (nrow(cycle_df) == 0) {
    warning(sprintf("No matching CSV rows found for cycle %d", target_cycle))
    DOXY_adj <- ncvar_get(nc, "DOXY_ADJUSTED")
    return(if (is.matrix(DOXY_adj)) DOXY_adj[, profile_idx] else DOXY_adj)
  }
  
  DOXY_full <- ncvar_get(nc, "DOXY")
  DOXY_vector <- if (is.matrix(DOXY_full)) DOXY_full[, profile_idx] else DOXY_full
  
  DOXY_adj_full <- if ("DOXY_ADJUSTED" %in% var_names) ncvar_get(nc, "DOXY_ADJUSTED") else DOXY_full
  DOXY_adj_vector <- if (is.matrix(DOXY_adj_full)) DOXY_adj_full[, profile_idx] else DOXY_adj_full
  
  # Prepare QC vector buffer initialized to spaces
  DOXY_adj_qc <- rep(" ", length(DOXY_adj_vector))
  
  # Match levels by PRES (0.01 dbar threshold tolerance)
  for (i in seq_along(pres_nc)) {
    if (is.na(pres_nc[i])) next
    
    match_idx <- which(abs(cycle_df$PRES - pres_nc[i]) < 0.01)
    if (length(match_idx) > 0) {
      idx <- match_idx[1]
      DOXY_vector[i]     <- cycle_df$DOXY[idx]
      DOXY_adj_vector[i] <- cycle_df$DOXY_FINAL[idx]
      
      # Populate DOXY_ADJUSTED_QC using DOXY_FINAL_QC from CSV (cast to character)
      DOXY_adj_qc[i]     <- as.character(cycle_df$DOXY_FINAL_QC[idx])
    }
  }
  
  # Generate raw DOXY QC flags
  DOXY_qc <- rep(" ", length(DOXY_vector))
  DOXY_qc[!is.na(DOXY_vector)] <- "1"
  
  # Write updated values back into NetCDF
  if (is.matrix(DOXY_full)) {
    ncvar_put(nc, "DOXY", DOXY_vector, start = c(1, profile_idx), count = c(length(DOXY_vector), 1))
    ncvar_put(nc, "DOXY_ADJUSTED", DOXY_adj_vector, start = c(1, profile_idx), count = c(length(DOXY_adj_vector), 1))
    ncvar_put(nc, "DOXY_QC", DOXY_qc, start = c(1, profile_idx), count = c(length(DOXY_qc), 1))
    ncvar_put(nc, "DOXY_ADJUSTED_QC", DOXY_adj_qc, start = c(1, profile_idx), count = c(length(DOXY_adj_qc), 1))
  } else {
    ncvar_put(nc, "DOXY", DOXY_vector)
    ncvar_put(nc, "DOXY_ADJUSTED", DOXY_adj_vector)
    ncvar_put(nc, "DOXY_QC", DOXY_qc)
    ncvar_put(nc, "DOXY_ADJUSTED_QC", DOXY_adj_qc)
  }
  
  return(DOXY_adj_vector)
}

write_DOXY_adjusted_error <- function(nc, profile_idx, err_mbar, psal, temp, pres, dens, DOXY_adj) {
  valid_idx <- !is.na(psal) & !is.na(DOXY_adj)
  DOXY_adj_error <- rep(NA, length(psal))
  
  err_umol_L <- err_mbar * 1.00
  DOXY_adj_error[valid_idx] <- err_umol_L
  
  DOXY_adj_error_umol_kg <- (DOXY_adj_error * 1000) / dens
  
  DOXY_adj_full <- ncvar_get(nc, "DOXY_ADJUSTED")
  if (is.matrix(DOXY_adj_full)) {
    ncvar_put(nc, "DOXY_ADJUSTED_ERROR", DOXY_adj_error_umol_kg, 
              start = c(1, profile_idx), count = c(length(DOXY_adj_error_umol_kg), 1))
  } else {
    ncvar_put(nc, "DOXY_ADJUSTED_ERROR", DOXY_adj_error_umol_kg)
  }
}

get_profile <- function(filename) {
  prof_str <- sub(".*_([0-9]+)[D]?\\.nc$", "\\1", basename(filename))
  prof_num <- as.integer(prof_str)
  
  if (is.na(prof_num)) {
    warning(sprintf("Could not parse profile number from: %s", filename))
    return(NA)
  }
  return(prof_num)
}

# Safely rename file with garbage collection and copy fallback for Windows file locks
safe_rename <- function(from_file, to_file) {
  # Force R to purge lingering C-level handles from ncdf4 library
  gc(verbose = FALSE)
  
  renamed <- suppressWarnings(file.rename(from_file, to_file))
  
  if (!renamed) {
    # Fallback if Windows file lock persists
    copied <- file.copy(from_file, to_file, overwrite = TRUE)
    if (copied) {
      unlink(from_file)
      cat(sprintf("Copied and replaced (lock fallback): %s -> %s\n", basename(from_file), basename(to_file)))
    } else {
      warning(sprintf("Failed to rename or copy %s to %s", from_file, to_file))
    }
  } else {
    cat(sprintf("Renamed %s to %s\n", basename(from_file), basename(to_file)))
  }
}

# ==============================================================================
# SECTION 3: Global Processing Loop
# ==============================================================================

for (i in seq_along(unique_floats)) {
  
  floatid <- as.integer(unique_floats[i])
  float_df <- csv_data[csv_data$FLOAT_NUM == floatid, ]
  
  cat(sprintf("\n==================================================\n"))
  cat(sprintf(" Processing Float ID: %d (%d of %d)\n", floatid, i, length(unique_floats)))
  cat(sprintf("==================================================\n"))
  
  tryCatch({
    inst_float <- "aoml_apex"
    profile_DOXY_qc <- "A"
    
    if (inst_float == "aoml_apex") {
      float_dir <- sprintf("C:/Users/Jennifer.McWhorter/Documents/Data/2026/For_BD_DOXY/%d", floatid)
      comment_dmqc_operator <- "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer McWhorter, NOAA/AOML"
      history_institution <- "AO"
      history_reference <- "WOA2023"
      DOXY_adj_err <- 2
      iprof <- 2 # Python index 1 -> R index 2
    } else if (inst_float == "aoml_navis") {
      float_dir <- sprintf("C:/Users/Jennifer.McWhorter/Documents/Data/2026/For_BD_DOXY/%d", floatid)
      comment_dmqc_operator <- "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer McWhorter, NOAA/AOML"
      history_institution <- "AO"
      history_reference <- "WOA2023"
      DOXY_adj_err <- 5
      iprof <- 1 # Python index 0 -> R index 1
    }
    
    float_dir <- sub("[/\\\\]+$", "", float_dir)
    today_str <- format(Sys.Date(), "%Y-%m-%d")
    output_dir <- file.path(float_dir, today_str)
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    parameter_data_mode <- "D"
    float_id_str <- as.character(floatid)
    
    # Extract unique cycles ONLY from the CSV dataframe
    csv_cycles <- sort(unique(float_df$CYCLE_NUMBER))
    
    if (length(csv_cycles) == 0) {
      warning(sprintf("No valid cycle numbers found in CSV for Float ID '%d'. Skipping...", floatid))
      next
    }
    
    # Process each cycle and immediately finalize file renaming to avoid accumulation
    for (target_cycle in csv_cycles) {
      
      # Search specifically for BD or BR files that match this cycle
      pattern_bd <- sprintf("BD%s_%03d.nc", float_id_str, target_cycle)
      pattern_br <- sprintf("BR%s_%03d.nc", float_id_str, target_cycle)
      
      matched_files <- Sys.glob(file.path(float_dir, c(pattern_bd, pattern_br)))
      
      if (length(matched_files) == 0) {
        # Fallback for filenames without 3-digit zero padding
        pattern_bd_raw <- sprintf("BD%s_%d.nc", float_id_str, target_cycle)
        pattern_br_raw <- sprintf("BR%s_%d.nc", float_id_str, target_cycle)
        matched_files <- Sys.glob(file.path(float_dir, c(pattern_bd_raw, pattern_br_raw)))
      }
      
      if (length(matched_files) == 0) {
        warning(sprintf("Cycle %d present in CSV, but no matching BD/BR file found in %s", target_cycle, float_dir))
        next
      }
      
      # Prefer BD file if present, otherwise use BR file
      bd_match <- matched_files[grepl("^BD", basename(matched_files))]
      bgc_filename <- if (length(bd_match) > 0) bd_match[1] else matched_files[1]
      
      cat(sprintf("Processing Cycle %d: %s\n", target_cycle, bgc_filename))
      
      w_bgc_filename <- create_working_bd_file(bgc_filename, output_dir)
      
      # Execute file write inside tryCatch to ensure connection closes reliably
      tryCatch({
        nc <- nc_open(w_bgc_filename, write = TRUE)
        
        # 1. Update history and global metadata
        write_history(nc, iprof, history_institution, history_reference, comment_dmqc_operator)
        
        # 2. Update parameter data modes
        write_parameter_data_mode(nc, iprof, parameter_data_mode)
        
        # 3. Write DOXY slope and drift values
        write_DOXY_slope_drift(nc, iprof, float_df, target_cycle)
        
        # 4. Populate DOXY, DOXY_ADJUSTED, and DOXY_ADJUSTED_QC directly from CSV profile matches
        DOXY_adjusted <- write_DOXY_from_csv(nc, iprof, float_df, target_cycle)
        
        # 5. Extract density profile matching physical parameters
        phys_filename <- get_phys_filename(bgc_filename, float_dir)
        pres_phys_raw <- get_phys_raw_pres(phys_filename)
        phys_data <- get_dens(phys_filename)
        
        pres_bgc <- ncvar_get(nc, "PRES")
        iprof_phys <- get_iprof_phys(pres_phys_raw, pres_bgc, iprof)
        
        # 6. Calculate & write DOXY adjusted error
        psal_col <- if (is.matrix(phys_data$psal)) phys_data$psal[, iprof_phys] else phys_data$psal
        temp_col <- if (is.matrix(phys_data$temp)) phys_data$temp[, iprof_phys] else phys_data$temp
        pres_col <- if (is.matrix(pres_bgc)) pres_bgc[, iprof] else pres_bgc
        dens_col <- if (is.matrix(phys_data$dens)) phys_data$dens[, iprof_phys] else phys_data$dens
        
        write_DOXY_adjusted_error(nc, iprof, DOXY_adj_err, 
                                  psal_col, 
                                  temp_col, 
                                  pres_col, 
                                  dens_col, 
                                  DOXY_adjusted)
        
        # 7. Update overall profile quality flag
        ncvar_put(nc, "PROFILE_DOXY_QC", profile_DOXY_qc, start = c(iprof), count = c(1))
        
      }, finally = {
        if (exists("nc")) {
          nc_close(nc)
          rm(nc)
        }
      })
      
      # Rename immediately within the cycle step while the context is clear
      base_name <- basename(w_bgc_filename)
      new_name  <- gsub("^w_BD", "BD", base_name)
      new_path  <- file.path(output_dir, new_name)
      
      safe_rename(w_bgc_filename, new_path)
    }
    
  }, error = function(e) {
    message(sprintf("Error processing Float ID %d: %s", floatid, e$message))
  })
}