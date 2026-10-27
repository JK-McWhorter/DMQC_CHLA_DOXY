### BGC-Argo CHLA, BBP700 & DOXY Corrections DMQC Pipeline
### Integrated DMQC processing for Float 4903622 (or batch execution)
### Created by Jen McWhorter and Marin Cornec
### Last updated Oct-2026: Sequential evaluation of cycles until 5 valid 900-1000m minimums are found,
###                       Dark offset derived as median of 5 cycle minimums,
###                       DOXY_QC ingestion & robust plot column protection,
###                       Pre-NPQ QC despiking, 0.9*MLD upward NPQ, FLUO_CHLA retention bug fix,
###                       and standardized per-cycle multi-mode plots (CHLA, FLUO, BBP, DOXY).

# DOXY drift was not applied to floats that were deployed at the end of May 2026, limited profiles for analysis 
## flags 
# 0 - NO QC performed
# 1 - good data 
# 2 - probably good data 
# 3 - probably bad data 
# 4 - bad data 
# 5 - value changed 
# 8 - estimated value 
# 9 - missing value 

# =========================================================================
# 1. ENVIRONMENT CLEANUP & PACKAGE INITIALIZATION
# =========================================================================

cat("\014")
rm(list = ls()) # Clean global environment at start of execution

required_pkgs <- c("dplyr", "ggplot2", "lubridate", "gridExtra", "tidyverse", 
                   "ggeffects", "mapdata", "suncalc", "ncdf4", "oce", "terra", 
                   "sf", "stars", "gstat", "rnaturalearth", "rnaturalearthdata", 
                   "viridis", "knitr", "patchwork", "stringr", "interp")

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# =========================================================================
# 2. CONFIGURATION & DIRECTORY MANAGEMENT
# =========================================================================

# --- Main Base Directory (Single Root for Inputs and Outputs) ---
dir_main_root <- "C:/Users/Jennifer.McWhorter/Documents/GitHub/DMQC_CHLA_DOXY"

# --- Structured Output & Intermediate Directories based on Main Path ---
dir_functions    <- file.path(dir_main_root, "Functions")
dir_output       <- file.path(dir_main_root, "Output")
dir_plot_base    <- file.path(dir_main_root, "Plots")
dir_script       <- file.path(dir_main_root, "Script")
dir_tables       <- file.path(dir_main_root, "Tables")
dir_profs        <- file.path(dir_main_root, "Tables/profs")

# Auxiliary & External Input Directories rerouted to Main Root
dir_data_root    <- file.path(dir_main_root, "Data")
combo_output_dir <- file.path(dir_output, "COMBO_CHLA_DOXY_DMODE")

# Create local storage directories if they do not exist
for (d in c(dir_functions, dir_output, dir_plot_base, dir_script, dir_tables, dir_profs, dir_data_root, combo_output_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Source local function scripts
func.sources <- list.files(dir_functions, pattern = "\\.R$", full.names = TRUE)
for (f in func.sources) {
  source(f)
}

# --- IN-MEMORY DICTIONARY MAPPING FOR LAUNCH DATES ---
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

# --- IN-MEMORY DICTIONARY MAPPING FOR FLOAT TYPES ---
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

# --- DOXY ADJUSTED ERROR WRITER FUNCTION ---
write_DOXY_adjusted_error <- function(nc, profile_idx, err_mbar, psal, temp, pres, dens, DOXY_adj) {
  valid_idx <- !is.na(psal) & !is.na(DOXY_adj) & !is.na(dens) & (dens > 0)
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

# --- Load Calibration Data ---
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


float_ids     <- calib_df$float_ids #C(4903622) to target one float 
target_floats <- c(6999992, 4903904, 1902800, 5907206, 4903624)
today_str     <- format(Sys.Date(), "%Y-%m-%d")

# Region Map Boundaries
reg_countries <- c('Mexico', 'USA', 'Cuba', 'Belize', "Guatemala", "Cayman Islands", "Jamaica")
lons          <- c(-98, -59)
lats          <- c(8, 29)

path_lut_nc     <- file.path(dir_data_root, "LUT/127630.nc")
path_bottle_val <- file.path(dir_data_root, "Bottle_samples/Master_Vali_24_July_2025.csv")
url_argo_index  <- "ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt"

qc_color_scale <- c(
  "0" = "gray70", "1" = "#33A02C", "2" = "#B2DF8A", "3" = "#FDBF6F",
  "4" = "#E31A1C", "5" = "#1F78B4", "8" = "#CAB2D6", "9" = "black"
)

# Standardized Palette: R-mode (Red), A-mode (Brown), D-mode (Green), Bottle (Black)
mode_levels      <- c("R-mode", "A-mode", "D-mode", "Bottle")
color_mapping    <- c("R-mode" = "#E31A1C", "A-mode" = "#8C510A", "D-mode" = "#33A02C", "Bottle" = "#000000")
linetype_mapping <- c("R-mode" = "dotted",  "A-mode" = "dashed",  "D-mode" = "solid",   "Bottle" = "blank")
linewidth_map    <- c("R-mode" = 0.7,       "A-mode" = 0.9,       "D-mode" = 1.1,       "Bottle" = 0)
shape_mapping    <- c("R-mode" = NA,        "A-mode" = NA,        "D-mode" = NA,        "Bottle" = 16)

# =========================================================================
# 3. GLOBAL DATA INGESTION & PREPARATION
# =========================================================================

options(timeout = 600)
options(HTTPUserAgent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)")

index_file_path <- file.path(dir_tables, 'index.txt')
if (!file.exists(index_file_path)) {
  download.file(
    url_argo_index, 
    index_file_path, 
    quiet = FALSE, 
    mode = "wb", 
    method = "libcurl", 
    cacheOK = TRUE
  )
}
warning(paste("last download of index file was on :", file.info(index_file_path)$ctime))

index <- read.table(index_file_path, skip = 9, sep = ",")
files <- as.character(index[, 1])
ident <- strsplit(files, "/")
ident <- matrix(unlist(ident), ncol = 4, byrow = TRUE)
dac <- ident[, 1]
wod <- as.character(ident[, 2])
prof_id <- ident[, 4]
lat <- index[, 3]
lon <- index[, 4]
time_index <- index[, 2]

LUT <- rast(path_lut_nc)
df_LUT <- as.data.frame(LUT, xy = TRUE, na.rm = TRUE)
lut_sf <- st_as_sf(df_LUT, coords = c("x", "y"), crs = 4326)

reg <- map_data("world2Hires")
reg <- subset(reg, region %in% reg_countries)
reg$long <- (360 - reg$long) * -1

if (file.exists(path_bottle_val)) {
  Bottle_Vali <- read.csv(path_bottle_val)
} else {
  warning("Bottle validation file not found at path: ", path_bottle_val)
  Bottle_Vali <- data.frame(Float_num = numeric(0), DOXY_bottle = numeric(0), PRES = numeric(0))
}

# =========================================================================
# 4. BATCH EXECUTION LOOP (INDIVIDUAL FLOAT PROCESSING)
# =========================================================================

for (woko in float_ids) {
  
  float_id   <- as.numeric(woko)
  active_wmo <- as.character(woko)
  
  # --- Structured Output Folders: Date -> Float_ID -> Sensor ---
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
  crit <- which(wod == active_wmo)
  
  if (length(crit) == 0) {
    warning(paste("No matching index records found on GDAC for WMO:", active_wmo))
    next
  }
  
  # --- STEP 1: DOWNLOAD & EXTRACT PROFILES (WITH EXPLICIT NA BYPASS) ---
  for (kiwi in crit) {
    print(c(woko, prof_id[kiwi]))
    
    if (substr(prof_id[kiwi], 14, 14) == "D") next 
    if (is.na(time_index[kiwi])) next 
    
    www <- substr(prof_id[kiwi], 3, 13)
    profile_bio <- NULL
    profile_physic <- NULL
    
    meta_path <- file.path(dir_profs, paste0(woko, "_meta.nc"))
    if (!file.exists(meta_path)) {
      try(download.file(paste0('ftp://ftp.ifremer.fr/ifremer/argo/dac/', woko, '/', woko, '_meta.nc'),
                        meta_path, quiet = FALSE, mode = "wb", method = "libcurl", cacheOK = TRUE), silent = TRUE)
    }
    
    bd_file <- file.path(dir_profs, paste0("BD", www, ".nc"))
    br_file <- file.path(dir_profs, paste0("BR", www, ".nc"))
    d_file  <- file.path(dir_profs, paste0("D", www, ".nc"))
    r_file  <- file.path(dir_profs, paste0("R", www, ".nc"))
    
    if (!file.exists(bd_file) && !file.exists(br_file)) {
      try(download.file(paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", files[kiwi]),
                        file.path(dir_profs, prof_id[kiwi]), mode = "wb", method = "libcurl", quiet = TRUE), silent = TRUE)
    }
    if (!file.exists(bd_file) && !file.exists(br_file)) next
    
    if (!file.exists(d_file) && !file.exists(r_file)) {
      try(download.file(paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/D", www, ".nc"),
                        d_file, mode = "wb", method = "libcurl", quiet = TRUE), silent = TRUE)
    }
    if (!file.exists(d_file) && !file.exists(r_file)) {
      try(download.file(paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/R", www, ".nc"),
                        r_file, mode = "wb", method = "libcurl", quiet = TRUE), silent = TRUE)
    }
    if (!file.exists(d_file) && !file.exists(r_file)) next
    
    profile_bio    <- nc_open(if (file.exists(bd_file)) bd_file else br_file, readunlim = FALSE, write = FALSE)
    profile_physic <- nc_open(if (file.exists(d_file)) d_file else r_file, readunlim = FALSE, write = FALSE)
    
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
    
    # Ingest PRES and bypass sentinels immediately to NA_real_
    pres <- as.vector(ncvar_get(profile_physic, "PRES"))
    pres[pres %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
    
    pres_qc_raw <- ncvar_get(profile_physic, "PRES_QC")
    pres_qc_str <- if (is.matrix(pres_qc_raw)) apply(pres_qc_raw, 2, function(x) paste(x, collapse = "")) else paste(pres_qc_raw, collapse = "")
    pres_qc <- as.numeric(unlist(strsplit(pres_qc_str[1], "")))[1:length(pres)]
    
    temp <- rep(NA_real_, length(pres)); temp_qc <- rep(NA_real_, length(pres))
    sal  <- rep(NA_real_, length(pres)); sal_qc  <- rep(NA_real_, length(pres))
    
    if ("TEMP" %in% names(profile_physic$var)) {
      temp <- as.vector(ncvar_get(profile_physic, "TEMP"))
      temp[temp %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      temp_qc_raw <- ncvar_get(profile_physic, "TEMP_QC")
      temp_qc_str <- if (is.matrix(temp_qc_raw)) apply(temp_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(temp_qc_raw, collapse = "")
      temp_qc <- as.numeric(unlist(strsplit(temp_qc_str, "")))[1:length(pres)]
    }
    
    if ("PSAL" %in% names(profile_physic$var)) {
      sal <- as.vector(ncvar_get(profile_physic, "PSAL"))
      sal[sal %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      sal_qc_raw <- ncvar_get(profile_physic, "PSAL_QC")
      sal_qc_str <- if (is.matrix(sal_qc_raw)) apply(sal_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(sal_qc_raw, collapse = "")
      sal_qc <- as.numeric(unlist(strsplit(sal_qc_str, "")))[1:length(pres)]
    }
    
    MLD <- NA_real_
    if (length(which(!is.na(pres))) > 15 && max(pres, na.rm = TRUE) > 250) {
      sigma_all <- swSigmaTheta(sal, temp, pres)
      valid_mld <- !is.na(pres) & !is.na(sigma_all)
      if (sum(valid_mld) > 5) {
        dep_sigma <- swDepth(pres[valid_mld], lat_val)
        MLD <- MLD_calc(sigma_all[valid_mld], dep_sigma)
      }
    }
    
    chl <- rep(NA_real_, length(pres)); chl_qc <- rep(NA_real_, length(pres))
    chl_adj <- rep(NA_real_, length(pres)) 
    bbp <- rep(NA_real_, length(pres)); bbp_qc <- rep(NA_real_, length(pres))
    fluo <- rep(NA_real_, length(pres)); fluo_qc <- rep(NA_real_, length(pres))
    doxy <- rep(NA_real_, length(pres)); doxy_qc <- rep(NA_real_, length(pres))
    doxy_adj <- rep(NA_real_, length(pres)); doxy_adj_qc <- rep(NA_real_, length(pres))
    scale_chla <- NA_real_; dark_chla <- NA_real_
    
    if ("CHLA" %in% names(profile_bio$var)) {
      chl <- as.vector(ncvar_get(profile_bio, "CHLA"))
      chl[chl %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      chl_qc_v <- paste(ncvar_get(profile_bio, "CHLA_QC"), collapse = "")
      chl_qc <- unlist(strsplit(chl_qc_v, ""))[1:length(pres)]
    }
    
    if ("CHLA_ADJUSTED" %in% names(profile_bio$var)) {
      chl_adj <- as.vector(ncvar_get(profile_bio, "CHLA_ADJUSTED"))
      chl_adj[chl_adj %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
    }
    
    if ("FLUORESCENCE_CHLA" %in% names(profile_bio$var)) {
      fluo <- as.vector(ncvar_get(profile_bio, "FLUORESCENCE_CHLA"))
      fluo[fluo %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      fluo_qc_v <- paste(ncvar_get(profile_bio, "FLUORESCENCE_CHLA_QC"), collapse = "")
      fluo_qc <- unlist(strsplit(fluo_qc_v, ""))[1:length(pres)]
      
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
    
    bbp_var_name <- intersect(c("BBP700", "BBP_700"), names(profile_bio$var))
    if (length(bbp_var_name) > 0) {
      bbp_raw <- as.vector(ncvar_get(profile_bio, bbp_var_name[1]))
      bbp_raw[bbp_raw %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      if (length(bbp_raw) == length(pres)) bbp <- bbp_raw
      qc_name <- paste0(bbp_var_name[1], "_QC")
      if (qc_name %in% names(profile_bio$var)) {
        bbp_qc_v <- paste(ncvar_get(profile_bio, qc_name), collapse = "")
        parsed_qc <- unlist(strsplit(bbp_qc_v, ""))[1:length(pres)]
        bbp_qc <- ifelse(is.na(parsed_qc), 1, as.numeric(parsed_qc))
      }
    }
    
    if ("DOXY" %in% names(profile_bio$var)) {
      doxy <- as.vector(ncvar_get(profile_bio, "DOXY"))
      doxy[doxy %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      if ("DOXY_QC" %in% names(profile_bio$var)) {
        doxy_qc_raw <- ncvar_get(profile_bio, "DOXY_QC")
        doxy_qc_str <- if (is.matrix(doxy_qc_raw)) apply(doxy_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(doxy_qc_raw, collapse = "")
        doxy_qc <- as.numeric(unlist(strsplit(doxy_qc_str, "")))[1:length(pres)]
      }
    }
    
    if ("DOXY_ADJUSTED" %in% names(profile_bio$var)) {
      doxy_adj <- as.vector(ncvar_get(profile_bio, "DOXY_ADJUSTED"))
      doxy_adj[doxy_adj %in% c(-99999, 99999, -99999.0, 99999.0)] <- NA_real_
      if ("DOXY_ADJUSTED_QC" %in% names(profile_bio$var)) {
        doxy_adj_qc_raw <- ncvar_get(profile_bio, "DOXY_ADJUSTED_QC")
        doxy_adj_qc_str <- if (is.matrix(doxy_adj_qc_raw)) apply(doxy_adj_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(doxy_adj_qc_raw, collapse = "")
        doxy_adj_qc <- as.numeric(unlist(strsplit(doxy_adj_qc_str, "")))[1:length(pres)]
      }
    }
    
    data_sub <- data.frame(
      PRES = as.numeric(pres), PRES_QC = pres_qc,
      TEMP = as.numeric(temp), PSAL = as.numeric(sal),
      CHLA = chl, CHLA_QC = chl_qc, CHLA_ADJUSTED = chl_adj,
      FLUO_CHLA = fluo, FLUO_CHLA_QC = fluo_qc,
      BBP700 = bbp, BBP700_QC = bbp_qc,
      DOXY = doxy, DOXY_QC = doxy_qc,
      DOXY_ADJUSTED = doxy_adj, DOXY_ADJUSTED_QC = doxy_adj_qc
    ) %>%
      filter(!is.na(CHLA) | !is.na(DOXY) | !is.na(FLUO_CHLA) | !is.na(BBP700)) %>%
      arrange(PRES)
    
    if (nrow(data_sub) > 0) {
      data_sub$TIME <- time_val
      data_sub$LONGITUDE <- lon_val
      data_sub$LATITUDE  <- lat_val
      data_sub$CYCLE_NUMBER <- as.numeric(str_remove(substr(www, 9, 11), "^0+"))
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
  # STEP 1.5: BBP700 MANUAL QUALITY CONTROL TESTS (PRESERVING NA)
  # =========================================================================
  message("\nProcessing BBP700 Manual Quality Control Tests...")
  float_data$BBP700_QC_MANUAL <- 1
  path_profs <- dir_profs
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    data_sub <- float_data[which(float_data$CYCLE_NUMBER == i), ]
    bbp <- data_sub$BBP700
    dep_bbp <- data_sub$PRES
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
      BBP700_QC_NUM   = suppressWarnings(as.numeric(as.character(BBP700_QC))),
      BBP700_FINAL_QC = pmax(coalesce(BBP700_QC_NUM, 1), coalesce(BBP700_QC_MANUAL, 1), na.rm = TRUE),
      BBP700_FINAL_QC = ifelse(is.na(BBP700) | PRES_QC %in% c("9", 9), 9, BBP700_FINAL_QC),
      BBP700_ADJUSTED = BBP700, # A-mode representation
      BBP700_FINAL    = ifelse(is.na(BBP700) | BBP700_FINAL_QC %in% c(4, 9), NA_real_, BBP700)
    ) %>%
    select(-BBP700_QC_NUM)
  
  # --- STEP 2: GEOSPATIAL MAPS & INITIAL QC ---
  float_data_qc_na <- float_data
  
  p_map <- ggplot() +
    geom_path(data = float_data_qc_na, aes(x = LONGITUDE, y = LATITUDE, color = float_num, group = float_num)) +
    geom_point(data = float_data_qc_na, aes(x = LONGITUDE, y = LATITUDE, color = float_num), size = 3) +
    geom_polygon(data = reg, aes(x = long, y = lat, group = group), fill = "darkgrey", color = NA) + 
    coord_map(xlim = lons, ylim = lats) +
    ylab("Lat (deg. N)") + xlab("Lon (deg. E)") +
    ggtitle(paste('WMO', active_wmo)) + theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(size = 15), 
          axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15), 
          axis.title.y = element_text(size = 15), aspect.ratio = 1)
  
  float_data_qc_na$PRES <- as.numeric(float_data_qc_na$PRES)
  float_data_qc_na$CHLA_QC <- as.numeric(float_data_qc_na$CHLA_QC)
  
  p_DMQC <- ggplot(float_data_qc_na) + facet_wrap(~float_num) +
    geom_point(aes(y = PRES, x = CHLA_QC, color = float_num)) +
    theme_bw() + ggtitle('CHLA_QC_flags by depth') + 
    scale_y_reverse(limits = c(2000, 5)) + labs(colour = "WMO", x = "QC flags", y = "Pres [mbar]") +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 8)) +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.text.y = element_text(size = 15), axis.title.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), 
          aspect.ratio = 1)
  
  p1 <- ggplot(data = float_data_qc_na) +
    geom_path(aes(x = CHLA, y = PRES, color = float_num, group = CYCLE_NUMBER), linewidth = .1) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    labs(colour = "CYCLE") + scale_y_reverse(limits = c(2000, 5)) + scale_x_continuous(position = "top") +
    labs(x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
  # =========================================================================
  # STEP 3: ARGO CHLA COOKBOOK IN SITU DARK CORRECTION
  # Sequential extraction of deep minimums across 5 valid cycles (900-1000 dbar)
  # =========================================================================
  float_data$FLUO_CHLA_RANGE <- ifelse(
    !is.na(float_data$FLUO_CHLA) & 
      !(float_data$FLUO_CHLA %in% c(-99999, 99999, -99999.0, 99999.0)) & 
      float_data$FLUO_CHLA >= -0.2 & 
      float_data$FLUO_CHLA <= 100,
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
    
    max_p <- max(data_sub$PRES, na.rm = TRUE)
    if (is.infinite(max_p) || max_p < 900) {
      cat(sprintf("--> Cycle %d: cast max depth < 900 dbar (%.1f dbar) - Skipping profile\n", cyc_i, max_p))
      next 
    }
    
    # Strictly extract valid observations in the 900-1000 dbar window
    deep_layer <- data_sub %>% 
      filter(!is.na(PRES) & 
               PRES >= 900 & 
               PRES <= 1000 & 
               !is.na(FLUO_CHLA_RANGE) & 
               !(FLUO_CHLA_RANGE %in% c(-99999, 99999, -99999.0, 99999.0)))
    
    if (nrow(deep_layer) < 3) {
      cat(sprintf("--> Cycle %d: no valid data in 900-1000 dbar window (n = %d) - Skipping profile\n", cyc_i, nrow(deep_layer)))
      next
    }
    
    min_val <- min(deep_layer$FLUO_CHLA_RANGE, na.rm = TRUE)
    if (!is.na(min_val) && is.finite(min_val)) {
      cycle_deep_mins   <- c(cycle_deep_mins, min_val)
      qualifying_cycles <- c(qualifying_cycles, cyc_i)
      cat(sprintf("--> Cycle %d: valid 900-1000 dbar minimum = %.4f RU (n = %d) [Total valid: %d/5]\n", 
                  cyc_i, min_val, nrow(deep_layer), length(cycle_deep_mins)))
    }
  }
  
  MED <- if (length(cycle_deep_mins) > 0) median(cycle_deep_mins, na.rm = TRUE) else NA_real_
  
  if (is.na(MED) || is.infinite(MED) || MED < 0) {
    MED <- coalesce(float_data$DARK_FLUO[1], 0)
    cat(sprintf("\n[DARK CORRECTION] Float %s - No qualifying deep data found. Falling back to Factory Dark: %.4f\n", active_wmo, MED))
  } else {
    cat(sprintf("\n[DARK CORRECTION] Float %s - Derived MED from 5 profile minimums (900-1000 dbar): %.4f (from cycles: %s)\n", 
                active_wmo, MED, paste(qualifying_cycles, collapse = ", ")))
  }
  
  scale_factor <- coalesce(float_data$SCALE_FLUO[1], 1)
  
  # Standard dark offset subtraction across water column (natural sensor noise floor preserved)
  float_data$CHLA_DARK                  <- (float_data$FLUO_CHLA_RANGE - MED) * scale_factor
  float_data$CHLA_NoLUT                 <- float_data$CHLA_DARK
  float_data$CHLA_ADJUSTED              <- coalesce(float_data$CHLA_ADJUSTED, float_data$CHLA_DARK)
  float_data$CHLA_FLUORESCENCE          <- (float_data$FLUO_CHLA - coalesce(float_data$DARK_FLUO, 0)) * scale_factor
  float_data$CHLA_FLUORESCENCE_ADJUSTED <- (float_data$FLUO_CHLA - MED) * scale_factor
  
  first_cycle_num <- ifelse(1 %in% float_data$CYCLE_NUMBER, 1, min(as.numeric(float_data$CYCLE_NUMBER), na.rm = TRUE))
  float_data_firstPROF <- float_data %>% filter(CYCLE_NUMBER == first_cycle_num)
  
  p2 <- ggplot(data = float_data_firstPROF) +
    geom_path(aes(x = CHLA, y = PRES, color = "Raw"), linewidth = .3) +
    geom_path(aes(x = CHLA_NoLUT, y = PRES, color = "Adjusted"), linewidth = .3) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    scale_color_manual(name = "Profile Type", breaks = c("Raw", "Adjusted"), 
                       values = c("Raw" = "darkgreen", "Adjusted" = "blue")) +
    scale_y_reverse(limits = c(2000, 5)) + scale_x_continuous(position = "top") +
    labs(title = paste0("WMO ", active_wmo, " - First Profile"), x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
  # =========================================================================
  # STEP 4: PRE-NPQ AUTOMATED QC & DESPIKING (COOKBOOK ORDER OF OPERATIONS)
  # =========================================================================
  float_data$CHLA_QC_AUTO <- 1
  for (i in unique(float_data$CYCLE_NUMBER)) {
    idx <- which(float_data$CYCLE_NUMBER == i)
    chla_v <- float_data$CHLA_NoLUT[idx]
    n_obs  <- length(chla_v)
    qc_vec <- rep(1, n_obs)
    
    # 1. Gross Range Check: [-0.1, 50.0] mg/m3 per Argo Cookbook
    range_fail <- !is.na(chla_v) & (chla_v < -0.1 | chla_v > 50.0)
    qc_vec[range_fail] <- 4
    
    # 2. BGC-Argo Official 5-point Running Median Residual Spike Test
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
    
    # 3. Stuck Value Test
    valid_v <- chla_v[!is.na(chla_v)]
    if (length(valid_v) > 5 && length(unique(valid_v)) == 1) {
      qc_vec[!is.na(chla_v)] <- 4
    }
    
    qc_vec[is.na(chla_v)] <- 9
    qc_vec[float_data$PRES_QC[idx] %in% c("9", 9)] <- 9
    float_data$CHLA_QC_AUTO[idx] <- qc_vec
  }
  
  # =========================================================================
  # STEP 4.5: SOLAR ANGLE & COOKBOOK NPQ (0.9*MLD, SURFACE-PEAK PROTECTED)
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
                             float_data$CHLA_QC_AUTO[idx] < 4)]
      
      if (length(mld_sub) > 1) {
        max_chla_idx <- mld_sub[which.max(float_data$CHLA_NoLUT[mld_sub])]
        max_chla_val <- float_data$CHLA_NoLUT[max_chla_idx]
        pres_at_max  <- float_data$PRES[max_chla_idx]
        shallowest_p <- min(float_data$PRES[mld_sub], na.rm = TRUE)
        
        # Surface-Peak Bypass: Do not extrapolate if maximum occurs at shallowest bin
        if (pres_at_max > shallowest_p && is.finite(max_chla_val)) {
          upward_idx <- idx[which(float_data$PRES[idx] < pres_at_max & float_data$CHLA_QC_AUTO[idx] < 4)]
          if (length(upward_idx) > 0) {
            float_data$CHLA_NoLUT[upward_idx] <- max_chla_val
            float_data$NPQ_APPLIED[upward_idx] <- TRUE
          }
        }
      }
    }
  }
  
  p3a <- ggplot(data = filter(float_data, sun_test == "day")) +
    geom_point(aes(x = CHLA_NoLUT, y = PRES, color = CYCLE_NUMBER), size = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    labs(colour = "CYCLE") + scale_y_reverse(limits = c(100, 0)) + scale_x_continuous(position = "top") +
    labs(x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
  p3 <- ggplot(data = float_data) +
    geom_path(aes(x = CHLA_NoLUT, y = PRES, color = sun_test), linewidth = .1) +
    geom_point(aes(x = CHLA_NoLUT, y = PRES, color = sun_test), size = .5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    labs(colour = "Sun Test") + scale_y_reverse(limits = c(300, 5)) + scale_x_continuous(position = "top") +
    labs(x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
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
  
  # ==============================================================================
  # BGC-ARGO FINAL QC ASSIGNMENT (UPDATING FLAG 5 ON MODIFIED POINTS)
  # ==============================================================================
  matched_result_v2 <- matched_result %>%
    mutate(
      SCALE_CHLA   = scale_chla,
      DARK_CHLA    = dark_chla,
      physio_ratio = ifelse(float_num %in% target_floats, 1, physio_ratio),
      CHLA_FINAL   = ifelse(float_num %in% target_floats, CHLA_NoLUT, CHLA_LUT),
      base_qc      = case_when(
        float_num %in% target_floats               ~ 2,
        !is.na(physio_ratio) & physio_ratio != 2 ~ 1,
        TRUE                                      ~ 2
      )
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
      qc_flags_chla[is.na(chla_v)]                  <- 9
      pres_qc_9_mask <- df_cyc$PRES_QC %in% c("9", 9)
      qc_flags_chla[pres_qc_9_mask] <- 9
      
      qc_flags_fluo <- rep(1, n_obs)
      range_fail_fluo <- !is.na(fluo_v) & (fluo_v < -0.2 | fluo_v > 100)
      qc_flags_fluo[range_fail_fluo] <- 4
      qc_flags_fluo[is.na(fluo_v)]   <- 9
      qc_flags_fluo[pres_qc_9_mask]  <- 9
      
      qc_flags_fluo_adj <- ifelse(length(cycle_deep_mins) >= 5, 1, 2)
      qc_flags_fluo_adj <- rep(qc_flags_fluo_adj, n_obs)
      qc_flags_fluo_adj[qc_flags_fluo >= 4] <- qc_flags_fluo[qc_flags_fluo >= 4]
      qc_flags_fluo_adj[is.na(df_cyc$CHLA_FLUORESCENCE_ADJUSTED)] <- 9
      qc_flags_fluo_adj[pres_qc_9_mask] <- 9
      
      df_cyc$CHLA_FINAL                    <- chla_v
      df_cyc$CHLA_FINAL_QC                 <- qc_flags_chla
      df_cyc$CHLA_FLUORESCENCE_QC          <- qc_flags_fluo
      df_cyc$CHLA_FLUORESCENCE_ADJUSTED_QC <- qc_flags_fluo_adj
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
  iprof        <- if (current_float_type == "aoml_navis") 1 else 2
  
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
      DOXY_QC       = if ("DOXY_QC" %in% names(.)) DOXY_QC else 1,
      DOXY_ADJUSTED = if ("DOXY_ADJUSTED" %in% names(.)) DOXY_ADJUSTED else NA_real_,
      DOXY_ADJUSTED_QC = if ("DOXY_ADJUSTED_QC" %in% names(.)) DOXY_ADJUSTED_QC else 1,
      CHLA_ADJUSTED = coalesce(CHLA_ADJUSTED, CHLA_NoLUT)
    ) %>%
    left_join(calib_df, by = "float_ids")
  
  active_slope       <- df_all$doxy_slope[1]
  active_drift       <- df_all$doxy_drift[1]
  active_launch_date <- df_all$launch_date[1]
  
  delta_days <- as.numeric(difftime(df_all$DATE_TIME, df_all$launch_date, units = "days"))
  gain       <- coalesce(df_all$doxy_slope, 1) * (1 + (coalesce(df_all$doxy_drift, 0) / 100) * (delta_days / 365.25))
  df_all$DOXY_CALCULATED_ADJ <- df_all$DOXY * gain
  
  df_all <- df_all %>%
    mutate(
      DENSITY = oce::swRho(salinity = PSAL, temperature = TEMP, pressure = PRES_ADJUSTED),
      DOXY_ADJUSTED_ERROR = write_DOXY_adjusted_error(
        nc = NULL, profile_idx = iprof, err_mbar = DOXY_adj_err, 
        psal = PSAL, temp = TEMP, pres = PRES_ADJUSTED, 
        dens = DENSITY, DOXY_adj = DOXY_CALCULATED_ADJ
      )
    )
  
  # Hook Flagging & DOXY QC Tests
  df_all <- df_all %>%
    group_by(CYCLE_NUMBER) %>%
    mutate(
      is_valid_pres = !is.na(PRES_ADJUSTED) & !is.na(DOXY_CALCULATED_ADJ),
      max_pres_val  = ifelse(any(is_valid_pres), max(PRES_ADJUSTED[is_valid_pres]), NA_real_),
      sec_pres_val  = ifelse(sum(is_valid_pres) >= 2, max(PRES_ADJUSTED[is_valid_pres & PRES_ADJUSTED < max_pres_val]), NA_real_),
      doxy_deepest  = ifelse(!is.na(max_pres_val) & PRES_ADJUSTED == max_pres_val, DOXY_CALCULATED_ADJ, NA_real_),
      doxy_second   = ifelse(!is.na(sec_pres_val) & PRES_ADJUSTED == sec_pres_val, DOXY_CALCULATED_ADJ, NA_real_),
      deepest_diff  = abs(max(doxy_deepest, na.rm = TRUE) - max(doxy_second, na.rm = TRUE)),
      HOOK_FLAG_4   = !is.na(PRES_ADJUSTED) & !is.na(max_pres_val) & (PRES_ADJUSTED == max_pres_val) & (!is.na(deepest_diff) & deepest_diff > 2)
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
      
      valid_v <- v[!is.na(v)]
      if (length(valid_v) > 1 && length(unique(valid_v)) == 1) {
        qc_flags[!is.na(v)] <- 4
      }
      
      if ("HOOK_FLAG_4" %in% names(df_cyc)) {
        qc_flags[df_cyc$HOOK_FLAG_4] <- 4
      }
      
      qc_flags[is.na(v)] <- 9
      pres_qc_9_mask_doxy <- df_cyc$PRES_QC %in% c("9", 9)
      qc_flags[pres_qc_9_mask_doxy] <- 9
      
      df_cyc$DOXY_FINAL_QC <- qc_flags
      return(df_cyc)
    }) %>%
    ungroup() %>%
    select(-any_of(c("HOOK_FLAG_4", "is_valid_pres", "max_pres_val", "sec_pres_val", "doxy_deepest", "doxy_second", "deepest_diff")))
  
  # Bottle samples alignment
  df_bottle_clean <- Bottle_Vali %>%
    filter(Float_num == float_id & !is.na(DOXY_bottle) & !is.na(PRES)) %>%
    select(PRES_ADJUSTED = PRES, Oxygen = DOXY_bottle) %>%
    mutate(Mode = factor("Bottle", levels = mode_levels))
  
  # ==============================================================================
  # STEP 5: PER-CYCLE PROFILES (CHLA, FLUORESCENCE, BBP700, DOXY)
  # Standard Legend: R-mode (red), A-mode (brown), D-mode (green), Bottle (black)
  # ==============================================================================
  all_cycles <- sort(unique(df_all$CYCLE_NUMBER[!is.na(df_all$CYCLE_NUMBER)]))
  message("\nProcessing standardized per-cycle mode plots for Float ", float_id, " across ", length(all_cycles), " cycles...")
  
  for (cyc in all_cycles) {
    df_sub <- df_all %>% filter(CYCLE_NUMBER == cyc)
    if (nrow(df_sub) == 0) next
    
    cyc_lat      <- df_sub$LATITUDE[1]
    cyc_lon      <- df_sub$LONGITUDE[1]
    cyc_date_str <- format(df_sub$DATE_TIME[1], "%Y-%m-%d")
    
    # -------------------------------------------------------------------------
    # 1. PER-CYCLE CHLA PROFILE PLOT
    # -------------------------------------------------------------------------
    df_chla_cycle_clean <- df_sub %>%
      select(PRES, CHLA, CHLA_ADJUSTED, CHLA_FINAL) %>%
      pivot_longer(
        cols      = c("CHLA", "CHLA_ADJUSTED", "CHLA_FINAL"),
        names_to  = "Mode_Raw",
        values_to = "Chlorophyll"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "CHLA"          ~ "R-mode",
          Mode_Raw == "CHLA_ADJUSTED" ~ "A-mode",
          Mode_Raw == "CHLA_FINAL"    ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(Chlorophyll) & !is.na(PRES)) %>%
      arrange(Mode, PRES) # D-mode plotted last on top
    
    if (nrow(df_chla_cycle_clean) > 0) {
      active_chla_modes <- levels(droplevels(df_chla_cycle_clean$Mode))
      
      p_chla_cycle <- ggplot(df_chla_cycle_clean, aes(x = Chlorophyll, y = PRES, color = Mode, linetype = Mode, linewidth = Mode)) +
        geom_path(alpha = 0.9) +
        scale_color_manual(values = color_mapping, limits = active_chla_modes, drop = TRUE) +
        scale_linetype_manual(values = linetype_mapping, limits = active_chla_modes, drop = TRUE) +
        scale_linewidth_manual(values = linewidth_map, limits = active_chla_modes, drop = TRUE) +
        scale_y_reverse(limits = c(max(df_chla_cycle_clean$PRES, na.rm = TRUE), 0)) +
        scale_x_continuous(position = "top") +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " CHLA Profiles (R-mode, A-mode, D-mode)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = expression(Chlorophyll ~ a ~ (mg/m^3)),
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              linetype  = linetype_mapping[active_chla_modes],
              linewidth = linewidth_map[active_chla_modes]
            )
          ),
          linetype  = "none",
          linewidth = "none"
        )
      
      ggsave(file.path(dir_chla_plots, sprintf("Float_%d_Cycle_%03d_CHLA_Profile.png", float_id, cyc)), 
             plot = p_chla_cycle, width = 9.5, height = 7, dpi = 300)
    }
    
    # -------------------------------------------------------------------------
    # 2. PER-CYCLE CHLA_FLUORESCENCE PROFILE PLOT
    # -------------------------------------------------------------------------
    df_fluo_cycle_clean <- df_sub %>%
      select(PRES, CHLA_FLUORESCENCE, CHLA_FLUORESCENCE_ADJUSTED) %>%
      pivot_longer(
        cols      = c("CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_ADJUSTED"),
        names_to  = "Mode_Raw",
        values_to = "Fluorescence"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "CHLA_FLUORESCENCE"          ~ "R-mode",
          Mode_Raw == "CHLA_FLUORESCENCE_ADJUSTED" ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(Fluorescence) & !is.na(PRES)) %>%
      arrange(Mode, PRES)
    
    if (nrow(df_fluo_cycle_clean) > 0) {
      active_fluo_modes <- levels(droplevels(df_fluo_cycle_clean$Mode))
      
      p_fluo_cycle <- ggplot(df_fluo_cycle_clean, aes(x = Fluorescence, y = PRES, color = Mode, linetype = Mode, linewidth = Mode)) +
        geom_path(alpha = 0.9) +
        scale_color_manual(values = color_mapping, limits = active_fluo_modes, drop = TRUE) +
        scale_linetype_manual(values = linetype_mapping, limits = active_fluo_modes, drop = TRUE) +
        scale_linewidth_manual(values = linewidth_map, limits = active_fluo_modes, drop = TRUE) +
        scale_y_reverse(limits = c(max(df_fluo_cycle_clean$PRES, na.rm = TRUE), 0)) +
        scale_x_continuous(position = "top") +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " CHLA Fluorescence Profiles (R-mode, D-mode)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = "Chlorophyll Fluorescence (RU)",
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              linetype  = linetype_mapping[active_fluo_modes],
              linewidth = linewidth_map[active_fluo_modes]
            )
          ),
          linetype  = "none",
          linewidth = "none"
        )
      
      ggsave(file.path(dir_fluo_plots, sprintf("Float_%d_Cycle_%03d_CHLA_FLUORESCENCE_Profile.png", float_id, cyc)), 
             plot = p_fluo_cycle, width = 9.5, height = 7, dpi = 300)
    }
    
    # -------------------------------------------------------------------------
    # 3. PER-CYCLE BBP700 PROFILE PLOT
    # -------------------------------------------------------------------------
    df_bbp_cycle_clean <- df_sub %>%
      select(PRES, BBP700, BBP700_ADJUSTED, BBP700_FINAL) %>%
      pivot_longer(
        cols      = c("BBP700", "BBP700_ADJUSTED", "BBP700_FINAL"),
        names_to  = "Mode_Raw",
        values_to = "Backscattering"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "BBP700"          ~ "R-mode",
          Mode_Raw == "BBP700_ADJUSTED" ~ "A-mode",
          Mode_Raw == "BBP700_FINAL"    ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(Backscattering) & !is.na(PRES)) %>%
      arrange(Mode, PRES)
    
    if (nrow(df_bbp_cycle_clean) > 0) {
      active_bbp_modes <- levels(droplevels(df_bbp_cycle_clean$Mode))
      
      p_bbp_cycle <- ggplot(df_bbp_cycle_clean, aes(x = Backscattering, y = PRES, color = Mode, linetype = Mode, linewidth = Mode)) +
        geom_path(alpha = 0.9) +
        scale_color_manual(values = color_mapping, limits = active_bbp_modes, drop = TRUE) +
        scale_linetype_manual(values = linetype_mapping, limits = active_bbp_modes, drop = TRUE) +
        scale_linewidth_manual(values = linewidth_map, limits = active_bbp_modes, drop = TRUE) +
        scale_y_reverse(limits = c(max(df_bbp_cycle_clean$PRES, na.rm = TRUE), 0)) +
        scale_x_continuous(position = "top") +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " BBP700 Profiles (R-mode, A-mode, D-mode)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = expression(b[bp](700) ~ (m^{-1})),
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              linetype  = linetype_mapping[active_bbp_modes],
              linewidth = linewidth_map[active_bbp_modes]
            )
          ),
          linetype  = "none",
          linewidth = "none"
        )
      
      ggsave(file.path(dir_bbp_plots, sprintf("Float_%d_Cycle_%03d_BBP700_Profile.png", float_id, cyc)), 
             plot = p_bbp_cycle, width = 9.5, height = 7, dpi = 300)
    }
    
    # -------------------------------------------------------------------------
    # 4. PER-CYCLE DOXY PROFILE PLOT
    # -------------------------------------------------------------------------
    df_doxy_cycle_clean <- df_sub %>%
      select(PRES_ADJUSTED, DOXY, DOXY_ADJUSTED, DOXY_CALCULATED_ADJ) %>%
      pivot_longer(
        cols      = c("DOXY", "DOXY_ADJUSTED", "DOXY_CALCULATED_ADJ"),
        names_to  = "Mode_Raw",
        values_to = "Oxygen"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "DOXY"                ~ "R-mode",
          Mode_Raw == "DOXY_ADJUSTED"       ~ "A-mode",
          Mode_Raw == "DOXY_CALCULATED_ADJ" ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = mode_levels)
      ) %>%
      filter(!is.na(Oxygen) & !is.na(PRES_ADJUSTED))
    
    has_bottle <- (cyc == 1) && (nrow(df_bottle_clean) > 0)
    df_plot_combined <- if (has_bottle) {
      bind_rows(df_doxy_cycle_clean, df_bottle_clean)
    } else {
      df_doxy_cycle_clean
    }
    
    df_plot_combined <- df_plot_combined %>% arrange(Mode, PRES_ADJUSTED)
    
    if (nrow(df_plot_combined) > 0) {
      active_doxy_modes <- levels(droplevels(df_plot_combined$Mode))
      
      p_doxy_cycle <- ggplot() +
        geom_path(
          data = filter(df_plot_combined, Mode != "Bottle"),
          aes(x = Oxygen, y = PRES_ADJUSTED, color = Mode, linetype = Mode, linewidth = Mode),
          alpha = 0.9
        ) +
        geom_point(
          data = filter(df_plot_combined, Mode == "Bottle"),
          aes(x = Oxygen, y = PRES_ADJUSTED, color = Mode, shape = Mode),
          size = 2.5, alpha = 0.95
        ) +
        scale_color_manual(values = color_mapping, limits = active_doxy_modes, drop = TRUE) +
        scale_linetype_manual(values = linetype_mapping, limits = active_doxy_modes, drop = TRUE) +
        scale_linewidth_manual(values = linewidth_map, limits = active_doxy_modes, drop = TRUE) +
        scale_shape_manual(values = shape_mapping, limits = active_doxy_modes, drop = TRUE) +
        scale_y_reverse(limits = c(2000, 0)) +
        scale_x_continuous(position = "top", limits = c(90, 250)) +
        labs(
          title    = paste0("Float ", float_id, " Cycle ", cyc, " Oxygen Profiles (Slope: ", active_slope, ", Drift: ", active_drift, "%/yr)"),
          subtitle = paste0("Lat: ", round(cyc_lat, 2), "°, Lon: ", round(cyc_lon, 2), "° | Cycle Date: ", cyc_date_str),
          caption  = paste0("Plot generated on: ", today_str),
          x        = expression(Oxygen ~ (mu * "mol/kg")),
          y        = "Pressure (dbar)",
          color    = "Data Source"
        ) +
        theme_bw(base_size = 13) +
        theme(
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.title      = element_text(face = "bold", size = 11),
          legend.text       = element_text(size = 10),
          legend.key.width  = unit(1.4, "cm"),
          legend.key.height = unit(0.5, "cm")
        ) +
        guides(
          color = guide_legend(
            override.aes = list(
              shape     = shape_mapping[active_doxy_modes],
              linetype  = linetype_mapping[active_doxy_modes],
              linewidth = linewidth_map[active_doxy_modes],
              size      = ifelse(active_doxy_modes == "Bottle", 2.5, NA)
            )
          ),
          linetype  = "none",
          shape     = "none",
          linewidth = "none"
        )
      
      ggsave(file.path(dir_doxy_plots, sprintf("Float_%d_Cycle_%03d_Profile.png", float_id, cyc)), 
             plot = p_doxy_cycle, width = 9.5, height = 7, dpi = 300)
    }
  }
  
  # ==============================================================================
  # STEP 6: SECTION CONTOUR & QC SUMMARY PLOTS (SAFE QC_VAL PROTECTION)
  # ==============================================================================
  df_summary_base <- df_all %>%
    filter(!is.na(PRES)) %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES         = as.numeric(PRES)
    )
  
  create_interpolated_grid <- function(df, val_col, grid_nx = 100, grid_ny = 100) {
    clean_df <- df %>% 
      filter(!is.na(!!sym(val_col))) %>% 
      select(x = CYCLE_NUMBER, y = PRES, z = !!sym(val_col)) %>% 
      group_by(x, y) %>% 
      summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
    
    if (nrow(clean_df) < 5 || length(unique(clean_df$x)) < 2 || length(unique(clean_df$y)) < 2) return(NULL)
    if (diff(range(clean_df$z, na.rm = TRUE)) == 0) return(NULL)
    
    interp_res <- suppressWarnings(tryCatch({
      interp::interp(
        x = clean_df$x, y = clean_df$y, z = clean_df$z,
        xo = seq(min(clean_df$x), max(clean_df$x), length.out = grid_nx),
        yo = seq(min(clean_df$y), max(clean_df$y), length.out = grid_ny),
        linear = TRUE, extrap = FALSE, duplicate = "mean"
      )
    }, error = function(e) return(NULL)))
    
    if (is.null(interp_res) || is.null(interp_res$z) || !is.matrix(interp_res$z)) return(NULL)
    if (nrow(interp_res$z) != grid_nx || ncol(interp_res$z) != grid_ny) return(NULL)
    
    z_vector <- as.vector(interp_res$z)
    grid_df  <- expand.grid(x = interp_res$x, y = interp_res$y)
    if (length(z_vector) != nrow(grid_df) || all(is.na(z_vector))) return(NULL)
    grid_df$z <- z_vector
    return(grid_df %>% filter(!is.na(z)))
  }
  
  generate_section_qc_plot <- function(df_base, var_name, qc_var_name, title_prefix, file_suffix, 
                                       col_palette = "viridis", unit_expr = "", is_chla = FALSE) {
    if (!var_name %in% names(df_base) || !any(!is.na(df_base[[var_name]]))) return(NULL)
    grid_df <- create_interpolated_grid(df_base, var_name)
    
    df_qc <- df_base %>% filter(!is.na(!!sym(var_name)))
    # Defensive lookup: if qc_var_name is not present, default gracefully to "1"
    if (qc_var_name %in% names(df_qc)) {
      df_qc$QC_VAL <- as.character(coalesce(as.character(df_qc[[qc_var_name]]), "1"))
    } else {
      df_qc$QC_VAL <- "1"
    }
    
    if (!is.null(grid_df)) {
      if (is_chla) {
        chla_breaks <- c(-Inf, 0.01, 0.03, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.2, 2.0, Inf)
        p_data <- ggplot() +
          geom_contour_filled(data = grid_df, aes(x = x, y = y, z = z), breaks = chla_breaks) +
          scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
          scale_fill_viridis_d(option = "turbo", name = unit_expr, drop = FALSE) +
          labs(title = paste0("Float ", float_id, " - Interpolated ", title_prefix, " Section Profile"),
               x = "Cycle Number", y = "Pressure (dbar)") +
          theme_bw(base_size = 12) + theme(legend.position = "right")
      } else {
        p_data <- ggplot() +
          geom_contour_filled(data = grid_df, aes(x = x, y = y, z = z), bins = 15) +
          scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
          scale_fill_viridis_d(option = col_palette, name = unit_expr) +
          labs(title = paste0("Float ", float_id, " - Interpolated ", title_prefix, " Section Profile"),
               x = "Cycle Number", y = "Pressure (dbar)") +
          theme_bw(base_size = 12) + theme(legend.position = "right")
      }
    } else {
      p_data <- ggplot(df_qc, aes(x = CYCLE_NUMBER, y = PRES)) +
        geom_tile(aes(fill = !!sym(var_name)), width = 1.2, height = 10) +
        scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
        scale_fill_viridis_c(option = col_palette, name = unit_expr) +
        labs(title = paste0("Float ", float_id, " - ", title_prefix, " Section Profile"),
             x = "Cycle Number", y = "Pressure (dbar)") +
        theme_bw(base_size = 12) + theme(legend.position = "right")
    }
    
    p_qc <- ggplot(df_qc, aes(x = CYCLE_NUMBER, y = PRES, color = QC_VAL)) +
      geom_point(alpha = 0.7, size = 1.2) +
      scale_y_reverse(limits = c(max(df_base$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(title = paste0(title_prefix, " QC Flags Section Map"), x = "Cycle Number", y = "Pressure (dbar)") +
      theme_bw(base_size = 12) + theme(legend.position = "right")
    
    p_combined <- (p_data / p_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - ", title_prefix, " Section DMQC Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    ggsave(file.path(dir_qc_summary, sprintf("Float_%d_%s_Section_Summary.png", float_id, file_suffix)), plot = p_combined, width = 12, height = 8, dpi = 150)
  }
  
  generate_section_qc_plot(df_summary_base, "DOXY", "DOXY_QC", "DOXY RAW (R-mode)", "DOXY_RAW", col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "DOXY_ADJUSTED", "DOXY_ADJUSTED_QC", "DOXY_ADJUSTED (A-mode)", "DOXY_ADJUSTED", col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "DOXY_FINAL", "DOXY_FINAL_QC", "DOXY_FINAL (D-mode)", "DOXY_FINAL", col_palette = "plasma", unit_expr = expression(Oxygen ~ (mu * "mol/kg")))
  generate_section_qc_plot(df_summary_base, "CHLA", "CHLA_QC", "CHLA RAW (R-mode)", "CHLA_RAW", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  generate_section_qc_plot(df_summary_base, "CHLA_ADJUSTED", "CHLA_FINAL_QC", "CHLA_ADJUSTED (A-mode)", "CHLA_ADJUSTED", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  generate_section_qc_plot(df_summary_base, "CHLA_FINAL", "CHLA_FINAL_QC", "CHLA_FINAL (D-mode)", "CHLA_FINAL", unit_expr = expression(Chlorophyll ~ a ~ (mg/m^3)), is_chla = TRUE)
  generate_section_qc_plot(df_summary_base, "CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_QC", "CHLA_FLUORESCENCE RAW (R-mode)", "CHLA_FLUORESCENCE_RAW", col_palette = "mako", unit_expr = "Fluorescence (RU)")
  generate_section_qc_plot(df_summary_base, "CHLA_FLUORESCENCE_ADJUSTED", "CHLA_FLUORESCENCE_ADJUSTED_QC", "CHLA_FLUORESCENCE_ADJUSTED (A-mode)", "CHLA_FLUORESCENCE_ADJUSTED", col_palette = "mako", unit_expr = "Fluorescence (RU)")
  generate_section_qc_plot(df_summary_base, "BBP700", "BBP700_QC", "BBP700 RAW (R-mode)", "BBP700_RAW", col_palette = "turbo", unit_expr = expression(b[bp](700) ~ (m^{-1})))
  generate_section_qc_plot(df_summary_base, "BBP700_FINAL", "BBP700_FINAL_QC", "BBP700_FINAL (D-mode)", "BBP700_FINAL", col_palette = "turbo", unit_expr = expression(b[bp](700) ~ (m^{-1})))
  
  # ==============================================================================
  # STEP 7: EXPORT COMBINED CSV (REPLACING NA WITH -99999.0 AT THE VERY END)
  # ==============================================================================
  combo_output_path <- file.path(combo_output_dir, sprintf("CHLA_BBP_DOXY_%d.csv", float_id))
  
  df_export_prep <- df_all %>%
    mutate(
      float_num = coalesce(as.character(WMOID), as.character(float_num)),
      across(where(is.logical), as.numeric)
    )
  
  df_export_prep <- df_export_prep %>%
    mutate(
      override_condition = (
        (float_num == "4903624") |
          (float_num == "2904010" & CYCLE_NUMBER >= 49) |
          (float_num == "2904011" & CYCLE_NUMBER >= 24)
      ),
      CHLA_FINAL                    = ifelse(override_condition, NA_real_, CHLA_FINAL),
      BBP700_FINAL                  = ifelse(override_condition, NA_real_, BBP700_FINAL),
      CHLA_FINAL_QC                 = ifelse(override_condition, 9, CHLA_FINAL_QC),
      BBP700_FINAL_QC               = ifelse(override_condition, 9, BBP700_FINAL_QC),
      CHLA                          = ifelse(override_condition, NA_real_, CHLA),
      CHLA_ADJUSTED                 = ifelse(override_condition, NA_real_, CHLA_ADJUSTED),
      CHLA_FLUORESCENCE             = ifelse(override_condition, NA_real_, CHLA_FLUORESCENCE),
      CHLA_FLUORESCENCE_QC          = ifelse(override_condition, 9, CHLA_FLUORESCENCE_QC),
      CHLA_FLUORESCENCE_ADJUSTED    = ifelse(override_condition, NA_real_, CHLA_FLUORESCENCE_ADJUSTED),
      CHLA_FLUORESCENCE_ADJUSTED_QC = ifelse(override_condition, 9, CHLA_FLUORESCENCE_ADJUSTED_QC),
      CHLA_NoLUT                    = ifelse(override_condition, NA_real_, CHLA_NoLUT),
      BBP700                        = ifelse(override_condition, NA_real_, BBP700),
      BBP700_QC                     = ifelse(override_condition, 9, BBP700_QC),
      SCALE_CHLA                    = ifelse(override_condition, NA_real_, SCALE_CHLA),
      DARK_CHLA                     = ifelse(override_condition, NA_real_, DARK_CHLA)
    ) %>%
    select(-override_condition)
  
  df_export_prep <- df_export_prep %>%
    mutate(
      CHLA_FINAL_QC   = ifelse(is.na(CHLA_FINAL), 9, CHLA_FINAL_QC),
      BBP700_FINAL_QC = ifelse(is.na(BBP700_FINAL), 9, BBP700_FINAL_QC),
      DOXY_FINAL_QC   = ifelse(is.na(DOXY_FINAL), 9, DOXY_FINAL_QC)
    )
  
  # Select final columns, clean names, and swap NAs to -99999.0 ONLY HERE
  df_combo_export <- df_export_prep %>%
    select(-any_of(c(
      "WMOID", "ID", "DOXY_QC", "DOXY_ADJUSTED_QC", "BBP700_QC_MANUAL",
      "DOXY_ADJUSTED", "DOXY_CALCULATED_ADJ", "float_ids", "TIME",
      "launch_date", "lauch_date", "LAUNCH_DATE", "sun_test", "SOLAR_ANGLE", "FLUO_CHLA_RANGE", "CHLA_QC_AUTO"
    ))) %>%
    rename_with(toupper) %>%
    relocate(
      FLOAT_NUM, CYCLE_NUMBER,
      CHLA_FINAL, CHLA_FINAL_QC,
      CHLA_FLUORESCENCE, CHLA_FLUORESCENCE_QC,
      CHLA_FLUORESCENCE_ADJUSTED, CHLA_FLUORESCENCE_ADJUSTED_QC,
      BBP700_FINAL, BBP700_FINAL_QC,
      DOXY_FINAL, DOXY_FINAL_QC, DOXY_ADJUSTED_ERROR
    ) %>%
    filter(!is.na(PRES)) %>%
    mutate(across(where(is.logical), ~ as.numeric(.))) %>%
    mutate(across(where(is.numeric), ~ tidyr::replace_na(., -99999.0))) %>%
    mutate(across(where(is.character), ~ tidyr::replace_na(., "-99999")))
  
  write.csv(df_combo_export, combo_output_path, row.names = FALSE)
  message(sprintf("\n>>> SUCCESSFULLY EXPORTED SINGLE COMBINED CSV -> %s", combo_output_path))
  message("\n=== COMPLETE DMQC PROCESSING FOR FLOAT: ", float_id, " ===")
}