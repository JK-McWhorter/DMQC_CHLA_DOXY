### BGC-Argo CHLA, BBP700 & DOXY Corrections DMQC Pipeline
### Integrated DMQC processing for all listed float IDs
### Created by Jen McWhorter and Marin Cornec
### Last updated 01-Oct-2026

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
                   "viridis", "knitr", "patchwork", "stringr")

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
  DOXY_adj_error <- rep(NA, length(psal))
  
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

# --- DIAGNOSTIC PRINT: VERIFY PARSED LAUNCH DATES PER FLOAT ID ---
cat("\n==============================================================================\n")
cat("   CALIBRATION LAUNCH DATE INGESTION CHECK\n")
cat("==============================================================================\n")
print(knitr::kable(calib_df %>% select(float_ids, launch_date, doxy_slope, doxy_drift, float_type), format = "simple"))
cat("==============================================================================\n\n")

# Extract processing targets and dates
float_ids     <- calib_df$float_ids
target_floats <- c(6999992, 4903904, 1902800, 5907206, 4903624) # Restricted strictly for CHLA No_LUT
today_str     <- format(Sys.Date(), "%Y-%m-%d")

# Region Map Boundaries
reg_countries <- c('Mexico', 'USA', 'Cuba', 'Belize', "Guatemala", "Cayman Islands", "Jamaica")
lons          <- c(-98, -59)
lats          <- c(8, 29)

# --- Specific File Inputs ---
path_lut_nc     <- file.path(dir_data_root, "LUT/127630.nc")
path_bottle_val <- file.path(dir_data_root, "Bottle_samples/Master_Vali_24_July_2025.csv")

# --- Web URLs for Data Ingestion ---
url_argo_index <- "ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt"

# Standard Argo QC Flag Color Mapping
qc_color_scale <- c(
  "0" = "gray70",  # No QC
  "1" = "#33A02C", # Good
  "2" = "#B2DF8A", # Probably Good
  "3" = "#FDBF6F", # Probably Bad
  "4" = "#E31A1C", # Bad
  "5" = "#1F78B4", # Value Changed
  "8" = "#CAB2D6", # Estimated
  "9" = "black"    # Missing
)

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
variables <- as.character(index[, 8])
var_mode  <- as.character(index[, 9])
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
  dir_bbp_plots   <- file.path(dir_float_plots, "BBP")
  dir_lut_plots   <- file.path(dir_float_plots, "CHLA", "LUT")
  dir_doxy_plots  <- file.path(dir_float_plots, "DOXY")
  dir_qc_summary  <- file.path(dir_float_plots, "QC_SUMMARIES")
  
  for (p_dir in c(dir_chla_plots, dir_bbp_plots, dir_lut_plots, dir_doxy_plots, dir_qc_summary)) {
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
  
  # --- STEP 1: DOWNLOAD & EXTRACT PROFILES ---
  for (kiwi in crit) {
    print(c(woko, prof_id[kiwi]))
    
    if (substr(prof_id[kiwi], 14, 14) == "D") { print("Descent Profile"); next } 
    if (is.na(time_index[kiwi])) { print("Time=NA"); next } 
    
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
    
    if (!file.exists(bd_file) && !file.exists(br_file)) {
      print("no bio file"); next
    }
    
    if (!file.exists(d_file) && !file.exists(r_file)) {
      try(download.file(paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/D", www, ".nc"),
                        d_file, mode = "wb", method = "libcurl", quiet = TRUE), silent = TRUE)
    }
    
    if (!file.exists(d_file) && !file.exists(r_file)) {
      try(download.file(paste0("ftp://ftp.ifremer.fr/ifremer/argo/dac/", dac[kiwi], "/", wod[kiwi], "/profiles/R", www, ".nc"),
                        r_file, mode = "wb", method = "libcurl", quiet = TRUE), silent = TRUE)
    }
    
    if (!file.exists(d_file) && !file.exists(r_file)) {
      print("no phy file"); next
    }
    
    if (www == "4903624_094") { print("weird profile"); next }
    
    profile_bio    <- nc_open(if (file.exists(bd_file)) bd_file else br_file, readunlim = FALSE, write = FALSE)
    profile_physic <- nc_open(if (file.exists(d_file)) d_file else r_file, readunlim = FALSE, write = FALSE)
    
    position_qc <- substr(ncvar_get(profile_physic, "POSITION_QC"), 1, 1)
    if (position_qc == 3 || position_qc == 4) print("bad position")
    
    lat_val <- ncvar_get(profile_physic, "LATITUDE")[1]
    lon_val <- ncvar_get(profile_physic, "LONGITUDE")[1]
    
    if (is.na(lat_val) || is.na(lon_val)) {
      print("no geoloc")
      nc_close(profile_bio); nc_close(profile_physic)
      next
    }
    
    jd <- ncvar_get(profile_physic, "JULD")[1]
    origin <- as.POSIXct("1950-01-01 00:00:00", tz = "UTC") 
    time_val <- origin + jd * 3600 * 24
    
    if (is.na(time_val)) {
      print("bad date")
      nc_close(profile_bio); nc_close(profile_physic)
      next
    }
    
    pres <- as.vector(ncvar_get(profile_physic, "PRES"))
    pres_qc_raw <- ncvar_get(profile_physic, "PRES_QC")
    if (is.matrix(pres_qc_raw)) {
      pres_qc_str <- apply(pres_qc_raw, 2, function(x) paste(x, collapse = ""))
    } else {
      pres_qc_str <- paste(pres_qc_raw, collapse = "")
    }
    pres_qc <- as.numeric(unlist(strsplit(pres_qc_str[1], "")))[1:length(pres)]
    
    temp <- rep(NA, length(pres)); temp_qc <- rep(NA, length(pres))
    sal  <- rep(NA, length(pres)); sal_qc  <- rep(NA, length(pres))
    
    if ("TEMP" %in% names(profile_physic$var)) {
      temp <- as.vector(ncvar_get(profile_physic, "TEMP"))
      temp_qc_raw <- ncvar_get(profile_physic, "TEMP_QC")
      temp_qc_str <- if (is.matrix(temp_qc_raw)) apply(temp_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(temp_qc_raw, collapse = "")
      temp_qc <- as.numeric(unlist(strsplit(temp_qc_str, "")))[1:length(pres)]
    }
    
    if ("PSAL" %in% names(profile_physic$var)) {
      sal <- as.vector(ncvar_get(profile_physic, "PSAL"))
      sal_qc_raw <- ncvar_get(profile_physic, "PSAL_QC")
      sal_qc_str <- if (is.matrix(sal_qc_raw)) apply(sal_qc_raw, 2, function(x) paste(x, collapse = ""))[1] else paste(sal_qc_raw, collapse = "")
      sal_qc <- as.numeric(unlist(strsplit(sal_qc_str, "")))[1:length(pres)]
    }
    
    MLD <- NA
    if (length(which(!is.na(pres))) > 15 && max(pres, na.rm = TRUE) > 250) {
      sigma_all <- swSigmaTheta(sal, temp, pres)
      valid_mld <- !is.na(pres) & !is.na(sigma_all)
      if (sum(valid_mld) > 5) {
        dep_sigma <- swDepth(pres[valid_mld], lat_val)
        MLD <- MLD_calc(sigma_all[valid_mld], dep_sigma)
      }
    }
    
    chl <- rep(NA, length(pres)); chl_qc <- rep(NA, length(pres))
    chl_adj <- rep(NA, length(pres)) 
    bbp <- rep(NA, length(pres)); bbp_qc <- rep(NA, length(pres))
    fluo <- rep(NA, length(pres)); fluo_qc <- rep(NA, length(pres))
    doxy <- rep(NA, length(pres)); doxy_qc <- rep(NA, length(pres))
    doxy_adj <- rep(NA, length(pres)); doxy_adj_qc <- rep(NA, length(pres))
    scale_chla <- NA; dark_chla <- NA
    
    if ("CHLA" %in% names(profile_bio$var)) {
      chl <- as.vector(ncvar_get(profile_bio, "CHLA"))
      chl_qc_v <- paste(ncvar_get(profile_bio, "CHLA_QC"), collapse = "")
      chl_qc <- unlist(strsplit(chl_qc_v, ""))[1:length(pres)]
    }
    
    if ("CHLA_ADJUSTED" %in% names(profile_bio$var)) {
      chl_adj <- as.vector(ncvar_get(profile_bio, "CHLA_ADJUSTED"))
    }
    
    if ("FLUORESCENCE_CHLA" %in% names(profile_bio$var)) {
      fluo <- as.vector(ncvar_get(profile_bio, "FLUORESCENCE_CHLA"))
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
    
    if ("BBP700" %in% names(profile_bio$var)) {
      bbp <- as.vector(ncvar_get(profile_bio, "BBP700"))
      bbp_qc_v <- paste(ncvar_get(profile_bio, "BBP700_QC"), collapse = "")
      bbp_qc <- unlist(strsplit(bbp_qc_v, ""))[1:length(pres)]
    }
    
    if ("DOXY" %in% names(profile_bio$var)) {
      doxy <- as.vector(ncvar_get(profile_bio, "DOXY"))
      doxy_qc_v <- paste(ncvar_get(profile_bio, "DOXY_QC"), collapse = "")
      doxy_qc <- unlist(strsplit(doxy_qc_v, ""))[1:length(pres)]
    }
    
    if ("DOXY_ADJUSTED" %in% names(profile_bio$var)) {
      doxy_adj <- as.vector(ncvar_get(profile_bio, "DOXY_ADJUSTED"))
      doxy_adj_qc_v <- paste(ncvar_get(profile_bio, "DOXY_ADJUSTED_QC"), collapse = "")
      doxy_adj_qc <- unlist(strsplit(doxy_adj_qc_v, ""))[1:length(pres)]
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
  # STEP 1.5: BBP700 MANUAL QUALITY CONTROL TESTS
  # =========================================================================
  message("\nProcessing BBP700 Manual Quality Control Tests...")
  float_data$BBP700_QC_MANUAL <- 1 # Initiate new QC column (default = good data)
  
  path_profs <- dir_profs # Alias for local nc lookup
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    data_sub <- float_data[which(float_data$CYCLE_NUMBER == i), ]
    bbp <- data_sub$BBP700
    dep_bbp <- data_sub$PRES
    
    if (all(is.na(bbp))) next
    
    ### Offset test
    try({
      bbp_qc_test <- bbp_offset(bbp)
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    ### Hook test
    www <- paste0(data_sub$float_num[1], "_", str_pad(data_sub$CYCLE_NUMBER[1], 3, pad = "0"))
    bd_path <- file.path(path_profs, paste0("BD", www, ".nc"))
    br_path <- file.path(path_profs, paste0("BR", www, ".nc"))
    meta_path <- file.path(path_profs, paste0(data_sub$float_num[1], "_meta.nc"))
    
    if (file.exists(bd_path) || file.exists(br_path)) {
      try({
        profile_bio <- nc_open(if (file.exists(bd_path)) bd_path else br_path, readunlim = FALSE, write = FALSE)
        meta_file <- nc_open(meta_path, readunlim = FALSE, write = FALSE)
        
        config_mission <- ncvar_get(profile_bio, "CONFIG_MISSION_NUMBER")[1]
        conf_par_name  <- ncvar_get(meta_file, "CONFIG_PARAMETER_NAME")
        conf_value     <- which(trimws(conf_par_name) == "CONFIG_ParkPressure_dbar")
        
        if (length(conf_value) > 0) {
          park_depth <- ncvar_get(meta_file, "CONFIG_PARAMETER_VALUE")[conf_value[1], config_mission]
          
          if (!is.na(park_depth)) {
            bbp_qc_test <- bbp_hook(bbp, dep_bbp, park_depth)
            for (j in 1:length(bbp_qc_test)) {
              if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
                float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
              }
            }
          }
        }
        nc_close(profile_bio)
        nc_close(meta_file)
      }, silent = TRUE)
    }
    
    ### Range test
    try({
      bbp_qc_test <- Range_test_bbp(bbp)$bbp_qc
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    ### High deep values
    try({
      bbp_qc_test <- bbp_high_deep_val(bbp, dep_bbp)
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    ### Negative values
    try({
      bbp_qc_test <- bbp_neg_val(bbp, dep_bbp)
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    ### Noisy values
    try({
      bbp_qc_test <- bbp_noisy(bbp, dep_bbp)
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
    
    ### Missing values
    try({
      bbp_qc_test <- bbp_missing(bbp, dep_bbp)
      for (j in 1:length(bbp_qc_test)) {
        if (!is.na(bbp_qc_test[j]) && bbp_qc_test[j] != 1) {
          float_data$BBP700_QC_MANUAL[which(float_data$CYCLE_NUMBER == i)][j] <- bbp_qc_test[j]
        }
      }
    }, silent = TRUE)
  }
  
  # Derive BBP700_FINAL and BBP700_FINAL_QC
  float_data <- float_data %>%
    mutate(
      BBP700_QC_NUM = suppressWarnings(as.numeric(as.character(BBP700_QC))),
      BBP700_FINAL_QC = pmax(coalesce(BBP700_QC_NUM, 1), coalesce(BBP700_QC_MANUAL, 1), na.rm = TRUE),
      BBP700_FINAL_QC = ifelse(is.na(BBP700) | PRES_QC %in% c("9", 9), 9, BBP700_FINAL_QC),
      BBP700_FINAL    = ifelse(is.na(BBP700) | BBP700_FINAL_QC %in% c(4, 9), 99999.0, BBP700)
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
  
  # --- STEP 3: DARK OFFSET DERIVATION & DEEP CHLA ZEROING ---
  float_data_qc_na_chla <- float_data_qc_na %>% filter(!is.na(CHLA) & CHLA > -0.2 & CHLA < 100)
  
  min_first_five <- NULL
  if (nrow(float_data_qc_na_chla) > 0) {
    first_five_cycles <- float_data_qc_na_chla %>%
      group_by(CYCLE_NUMBER) %>%
      summarise(min_time = min(TIME, na.rm = TRUE), max_pres = max(PRES, na.rm = TRUE)) %>%
      arrange(min_time) %>%
      pull(CYCLE_NUMBER) %>%
      head(5)
    
    for (i in first_five_cycles) {
      data_sub <- float_data_qc_na_chla %>% filter(CYCLE_NUMBER == i)
      deep_data <- data_sub %>% filter(PRES >= 600)
      if (nrow(deep_data) == 0) deep_data <- data_sub
      
      if (nrow(deep_data) > 0 && any(!is.na(deep_data$FLUO_CHLA))) {
        min_prof <- min(deep_data$FLUO_CHLA, na.rm = TRUE)
        min_first_five <- c(min_first_five, min_prof)
      }
    }
  }
  
  MED <- median(min_first_five, na.rm = TRUE)
  if (is.na(MED) || is.infinite(MED)) MED <- 0 
  
  float_data$min_FLUOCHLA_cycle1 <- ifelse(length(min_first_five) >= 1, min_first_five[1], NA)
  float_data$min_FLUOCHLA_cycle2 <- ifelse(length(min_first_five) >= 2, min_first_five[2], NA)
  float_data$min_FLUOCHLA_cycle3 <- ifelse(length(min_first_five) >= 3, min_first_five[3], NA)
  float_data$min_FLUOCHLA_cycle4 <- ifelse(length(min_first_five) >= 4, min_first_five[4], NA)
  float_data$min_FLUOCHLA_cycle5 <- ifelse(length(min_first_five) >= 5, min_first_five[5], NA)
  float_data$Median_first5 <- MED
  
  float_data$CHLA_NoLUT <- (float_data$FLUO_CHLA - MED) * coalesce(float_data$SCALE_FLUO, 1)
  float_data$CHLA_FLUORESCENCE <- (float_data$FLUO_CHLA - coalesce(float_data$DARK_FLUO, 0)) * coalesce(float_data$SCALE_FLUO, 1)
  float_data$CHLA_FLUORESCENCE_ADJUSTED <- (float_data$FLUO_CHLA - MED) * coalesce(float_data$SCALE_FLUO, 1)
  
  if (!"CHLA_ADJUSTED" %in% names(float_data) || all(is.na(float_data$CHLA_ADJUSTED))) {
    float_data$CHLA_ADJUSTED <- float_data$CHLA_NoLUT
  }
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    data_sub <- float_data[which(float_data$CYCLE_NUMBER == i),]
    
    if (nrow(data_sub) >= 3 && sum(!is.na(data_sub$FLUO_CHLA)) >= 3) {
      filt_val <- RunningFilter(2, data_sub$FLUO_CHLA, na.fill = TRUE, ends.fill = TRUE, Method = "Median")
      dep_min_prof <- data_sub$PRES[which.min(filt_val)]
    } else if (nrow(data_sub) > 0 && sum(!is.na(data_sub$FLUO_CHLA)) > 0) {
      dep_min_prof <- data_sub$PRES[which.min(data_sub$FLUO_CHLA)]
    } else {
      dep_min_prof <- NA
    }
    
    if (!is.na(dep_min_prof)) {
      float_data$CHLA_NoLUT[which(float_data$CYCLE_NUMBER == i & float_data$PRES >= dep_min_prof)] <- 0
    }
  }
  
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
  
  # --- STEP 4: SOLAR ANGLE & NPQ CORRECTION ---
  float_data <- tibble::rowid_to_column(float_data, "ID")
  float_data$SOLAR_ANGLE <- NA
  
  for (i in unique(float_data$CYCLE_NUMBER)) {
    idx <- which(float_data$CYCLE_NUMBER == i)
    float_data$SOLAR_ANGLE[idx] <- getSunlightPosition(
      unique(float_data$TIME[idx])[1],
      unique(float_data$LATITUDE[idx])[1],
      unique(float_data$LONGITUDE[idx])[1])$altitude * (180 / pi)
  }
  
  df_day   <- float_data %>% filter(SOLAR_ANGLE >= 0)
  df_night <- float_data %>% filter(SOLAR_ANGLE < 0)
  
  df_day$sun_test <- 'day'; df_night$sun_test <- 'night'
  df_suntest <- bind_rows(df_day, df_night)
  
  float_data_MLD <- select(df_suntest, float_num, CYCLE_NUMBER, LONGITUDE, LATITUDE, TIME, PRES, CHLA_NoLUT, CHLA, MLD, sun_test)
  float_data_ML_up <- float_data_MLD %>% group_by(CYCLE_NUMBER) %>% filter(sun_test == 'day' & PRES <= MLD)
  float_data_ML_down <- float_data_MLD %>% group_by(CYCLE_NUMBER) %>% filter(sun_test == 'day' & PRES > MLD)
  
  if (nrow(float_data_ML_up) > 0) {
    float_data_ML_up_CHLA <- float_data_ML_up %>% group_by(CYCLE_NUMBER) %>% filter(PRES == max(PRES))
    float_data_v2 <- select(float_data_ML_up_CHLA, CYCLE_NUMBER, CHLA_NoLUT) %>% rename(CHLA_new = CHLA_NoLUT)
    
    test <- merge(float_data_v2, float_data_ML_up)
    test$CHLA_NoLUT <- test$CHLA_new
    test2 <- full_join(test, float_data_ML_down, by = c("CYCLE_NUMBER", "float_num", "LONGITUDE", "LATITUDE", "TIME", "PRES", "CHLA_NoLUT", "CHLA", "MLD", "sun_test"))
    data_new <- test2
    data_new$CHLA_new[is.na(data_new$CHLA_new)] <- data_new$CHLA_NoLUT[is.na(data_new$CHLA_new)]
    day_adjusted <- select(data_new, float_num, CYCLE_NUMBER, TIME, PRES, CHLA, CHLA_NoLUT, MLD, LATITUDE, LONGITUDE)
  } else {
    day_adjusted <- select(df_day, float_num, CYCLE_NUMBER, TIME, PRES, CHLA, CHLA_NoLUT, MLD, LATITUDE, LONGITUDE)
  }
  
  night_selected <- select(df_night, float_num, CYCLE_NUMBER, TIME, PRES, CHLA, CHLA_NoLUT, MLD, LATITUDE, LONGITUDE)
  day_adjusted$sun_test <- 'day'; night_selected$sun_test <- 'night'
  df_all_22 <- bind_rows(day_adjusted, night_selected)
  
  p3a <- ggplot(data = day_adjusted) +
    geom_point(aes(x = CHLA_NoLUT, y = PRES, color = CYCLE_NUMBER), linewidth = .1) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    labs(colour = "CYCLE") + scale_y_reverse(limits = c(100, 0)) + scale_x_continuous(position = "top") +
    labs(x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
  p3 <- ggplot(data = df_all_22) +
    geom_path(aes(x = CHLA_NoLUT, y = PRES, color = sun_test), linewidth = .1) +
    geom_point(aes(x = CHLA_NoLUT, y = PRES, color = sun_test), size = .5) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15), 
          axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 15), 
          legend.text = element_text(size = 15), legend.title = element_text(size = 15), aspect.ratio = 1) + 
    labs(colour = "Sun Test") + scale_y_reverse(limits = c(300, 5)) + scale_x_continuous(position = "top") +
    labs(x = expression(paste("Chlorophyll a [", mg, "/m"^"3", "]")), y = "Pressure [dbar]")
  
  df_all_22$CHLA_NPQ_V1 <- df_all_22$CHLA_NoLUT
  for (i in unique(df_all_22$CYCLE_NUMBER)) {
    idx <- which(df_all_22$CYCLE_NUMBER == i)
    if (length(idx) > 0 && unique(df_all_22$sun_test[idx])[1] != "night") {
      data_sub <- df_all_22[idx, ]
      df_all_22$CHLA_NPQ_V1[idx] <- NPQ_MLD(data_sub$CHLA_NoLUT, data_sub$PRES, unique(data_sub$MLD)[1])
    }
  }
  
  npq_lookup <- df_all_22 %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES = as.numeric(PRES)
    ) %>%
    distinct(CYCLE_NUMBER, PRES, .keep_all = TRUE) %>%
    select(CYCLE_NUMBER, PRES, CHLA_NPQ_V1)
  
  float_data <- float_data %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES = as.numeric(PRES)
    ) %>%
    left_join(npq_lookup, by = c("CYCLE_NUMBER", "PRES")) %>%
    mutate(
      CHLA_NoLUT = ifelse(!is.na(CHLA_NPQ_V1), CHLA_NPQ_V1, CHLA_NoLUT)
    ) %>%
    select(-CHLA_NPQ_V1)
  
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
      mutate(
        CHLA_LUT     = CHLA_NoLUT, 
        physio_ratio = 1
      )
  }
  
  matched_result$match_dist_m <- as.numeric(st_distance(float_sf, lut_sf[nearest_indices, ], by_element = TRUE))
  DMQC_first_prof <- subset(matched_result, CYCLE_NUMBER == first_cycle_num)
  
  p4_LUT <- ggplot() + 
    scale_y_reverse(limits = c(500, 0)) + scale_x_continuous(position = "top") + theme_grey(base_size = 15) + 
    geom_path(data = DMQC_first_prof, aes(x = CHLA_LUT, y = PRES, color = 'CHLA LUT')) + 
    geom_point(data = DMQC_first_prof, aes(x = CHLA_LUT, y = PRES, color = 'CHLA LUT')) + 
    geom_path(data = DMQC_first_prof, aes(x = CHLA_NoLUT, y = PRES, color = 'CHLA No LUT')) + 
    geom_point(data = DMQC_first_prof, aes(x = CHLA_NoLUT, y = PRES, color = 'CHLA No LUT')) + 
    geom_path(data = DMQC_first_prof, aes(x = CHLA, y = PRES, color = 'CHLA Raw')) + 
    geom_point(data = DMQC_first_prof, aes(x = CHLA, y = PRES, color = 'CHLA Raw')) + 
    scale_color_manual(name = "First Profile", breaks = c('CHLA LUT', 'CHLA No LUT', 'CHLA Raw'), 
                       values = c('CHLA LUT' = 'darkgreen', 'CHLA No LUT' = 'green', 'CHLA Raw' = 'blue')) +
    ylab('Pressure [dbar]') + xlab('Chlorophyll [mg/m3]') + labs(title = paste0('WMO# ', active_wmo, ', LUT')) + 
    theme(aspect.ratio = 1, plot.title = element_text(size = 15, face = "bold", vjust = 2), legend.position = "right")
  
  plot_data <- rbind(
    data.frame(CHLA = matched_result$CHLA_LUT, PRES = matched_result$PRES, CYCLE_NUMBER = as.numeric(matched_result$CYCLE_NUMBER), DataType = "CHLA LUT"),
    data.frame(CHLA = matched_result$CHLA_NoLUT, PRES = matched_result$PRES, CYCLE_NUMBER = as.numeric(matched_result$CYCLE_NUMBER), DataType = "CHLA No LUT"),
    data.frame(CHLA = matched_result$CHLA, PRES = matched_result$PRES, CYCLE_NUMBER = as.numeric(matched_result$CYCLE_NUMBER), DataType = "CHLA Raw")
  )
  
  p5_side_by_side <- ggplot(plot_data, aes(x = CHLA, y = PRES, color = CYCLE_NUMBER, group = CYCLE_NUMBER)) + 
    scale_y_reverse(limits = c(500, 0)) + scale_x_continuous(position = "top") + theme_grey(base_size = 15) + 
    geom_path(linewidth = 0.2) + geom_point(size = 0.2) + 
    facet_wrap(~ DataType, nrow = 1) + scale_color_gradient(name = "Cycle Number", low = "blue", high = "red") +
    ylab('Pressure [dbar]') + xlab('Chlorophyll [mg/m3]') + labs(title = paste0('WMO# ', active_wmo, ', LUT Comparison')) + 
    theme(aspect.ratio = 1, plot.title = element_text(size = 15, face = "bold", vjust = 2),
          legend.position = "right", strip.text = element_text(face = "bold", size = 12))
  
  # --- STEP 6: SAVE CHLA & LUT PLOTS TO STRUCTURED DIRECTORIES ---
  ggsave(file.path(dir_chla_plots, paste0("CHLA_trajectory_map_", active_wmo, ".png")), plot = p_map, width = 8, height = 8, dpi = 300)
  ggsave(file.path(dir_chla_plots, paste0("CHLA_QC_flags_", active_wmo, ".png")), plot = p_DMQC, width = 8, height = 8, dpi = 300)
  ggsave(file.path(dir_chla_plots, paste0("CHLA_raw_profiles_", active_wmo, ".png")), plot = p1, width = 8, height = 8, dpi = 300)
  ggsave(file.path(dir_chla_plots, paste0("CHLA_first_profile_comparison_", active_wmo, ".png")), plot = p2, width = 8, height = 8, dpi = 300)
  ggsave(file.path(dir_chla_plots, paste0("CHLA_day_adjustments_mld_", active_wmo, ".png")), plot = p3a, width = 8, height = 8, dpi = 300)
  ggsave(file.path(dir_chla_plots, paste0("CHLA_diurnal_adjustments_", active_wmo, ".png")), plot = p3, width = 8, height = 8, dpi = 300)
  
  g4 <- arrangeGrob(p4_LUT, ncol = 1, nrow = 1)
  ggsave(file.path(dir_lut_plots, paste0("CHLA_LUT_", active_wmo, ".png")), g4, width = 10, height = 8, dpi = 300)
  
  g5 <- arrangeGrob(p5_side_by_side, ncol = 1, nrow = 1)
  ggsave(file.path(dir_lut_plots, paste0("CHLA_LUT_", active_wmo, "_allcycles.png")), g5, width = 12, height = 4, dpi = 300)
  
  # ==============================================================================
  # BGC-ARGO COOKBOOK CHLA & CHLA_FLUORESCENCE QUALITY CONTROL TESTS
  # ==============================================================================
  
  matched_result_v2 <- matched_result %>%
    mutate(
      SCALE_CHLA = scale_chla,
      DARK_CHLA  = dark_chla,
      physio_ratio = ifelse(float_num %in% target_floats, 1, physio_ratio),
      CHLA_FINAL = ifelse(float_num %in% target_floats, CHLA_NoLUT, CHLA_LUT),
      
      base_qc = case_when(
        float_num %in% target_floats               ~ 2,
        !is.na(physio_ratio) & physio_ratio != 2 ~ 1,
        TRUE                                      ~ 2
      ),
      
      is_daytime = SOLAR_ANGLE >= 0
    ) %>%
    group_by(CYCLE_NUMBER) %>%
    group_modify(~ {
      df_cyc <- .x
      chla_v  <- df_cyc$CHLA_FINAL
      fluo_v  <- df_cyc$CHLA_FLUORESCENCE
      n_obs   <- length(chla_v)
      
      qc_flags_chla <- df_cyc$base_qc
      npq_mask <- df_cyc$is_daytime & !is.na(df_cyc$MLD) & (df_cyc$PRES <= df_cyc$MLD)
      qc_flags_chla[npq_mask] <- 5
      
      range_fail_chla <- !is.na(chla_v) & (chla_v < -0.2 | chla_v > 100)
      qc_flags_chla[range_fail_chla] <- 4
      
      spike_fail_chla <- rep(FALSE, n_obs)
      if (n_obs >= 3) {
        v_prev <- c(NA, chla_v[1:(n_obs - 1)])
        v_next <- c(chla_v[2:n_obs], NA)
        
        spike_val <- abs(chla_v - 0.5 * (v_prev + v_next)) - 0.5 * abs(v_next - v_prev)
        spike_fail_chla <- !is.na(spike_val) & (spike_val > 5.0)
        qc_flags_chla[spike_fail_chla] <- 4
      }
      
      valid_v_chla <- chla_v[!is.na(chla_v)]
      if (length(valid_v_chla) > 1 && length(unique(valid_v_chla)) == 1) {
        qc_flags_chla[!is.na(chla_v)] <- 4
      }
      
      qc_flags_chla[df_cyc$CHLA_QC %in% c("4", 4)] <- 4
      qc_flags_chla[is.na(chla_v)]                  <- 9
      
      pres_qc_9_mask <- df_cyc$PRES_QC %in% c("9", 9)
      qc_flags_chla[pres_qc_9_mask] <- df_cyc$PRES_QC[pres_qc_9_mask]
      
      qc_flags_fluo <- rep(1, n_obs)
      range_fail_fluo <- !is.na(fluo_v) & (fluo_v < -0.2 | fluo_v > 100)
      qc_flags_fluo[range_fail_fluo] <- 4
      
      valid_v_fluo <- fluo_v[!is.na(fluo_v)]
      if (length(valid_v_fluo) > 1 && length(unique(valid_v_fluo)) == 1) {
        qc_flags_fluo[!is.na(fluo_v)] <- 4
      }
      
      qc_flags_fluo[df_cyc$PRES_QC %in% c("4", 4)] <- 4
      qc_flags_fluo[is.na(fluo_v)]                  <- 9
      qc_flags_fluo[pres_qc_9_mask] <- df_cyc$PRES_QC[pres_qc_9_mask]
      
      qc_flags_fluo_adj <- ifelse(length(min_first_five) >= 5, 1, 2)
      qc_flags_fluo_adj <- rep(qc_flags_fluo_adj, n_obs)
      
      qc_flags_fluo_adj[qc_flags_fluo >= 4] <- qc_flags_fluo[qc_flags_fluo >= 4]
      qc_flags_fluo_adj[is.na(df_cyc$CHLA_FLUORESCENCE_ADJUSTED)] <- 9
      qc_flags_fluo_adj[pres_qc_9_mask] <- df_cyc$PRES_QC[pres_qc_9_mask]
      
      df_cyc$CHLA_FINAL_QC                  <- qc_flags_chla
      df_cyc$CHLA_FLUORESCENCE_QC           <- qc_flags_fluo
      df_cyc$CHLA_FLUORESCENCE_ADJUSTED_QC  <- qc_flags_fluo_adj
      
      return(df_cyc)
    }) %>%
    ungroup() %>%
    select(-any_of(c(
      "SOLAR_ANGLE", "FLUO_CHLA", "FLUO_CHLA_QC", 
      "match_dist_m", "DARK_FLUO", "SCALE_FLUO", 
      "CHLA_LUT", "base_qc", "is_daytime"
    )))
  
  # ==============================================================================
  # MASTER BGC-ARGO DMQC OXYGEN ADJUSTMENT PIPELINE
  # ==============================================================================
  
  # --- Instrument Profile & Error Configuration ---
  current_float_type <- float_type_dict$float_type[float_type_dict$float_ids == float_id]
  if (length(current_float_type) == 0) current_float_type <- "aoml_apex"
  
  if (current_float_type == "aoml_navis") {
    inst_float   <- "aoml_navis"
    DOXY_adj_err <- 5      # 5 umol/kg error configuration for Navis
    iprof        <- 1
  } else {
    inst_float   <- "aoml_apex"
    DOXY_adj_err <- 2      # Default 2 umol/kg error configuration for Apex
    iprof        <- 2
  }
  
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
      DOXY_ADJUSTED = if ("DOXY_ADJUSTED" %in% names(.)) DOXY_ADJUSTED else NA_real_
    ) %>%
    left_join(calib_df, by = "float_ids")
  
  active_slope       <- df_all$doxy_slope[1]
  active_drift       <- df_all$doxy_drift[1]
  active_launch_date <- df_all$launch_date[1]
  
  message("\n==============================================================================")
  message(">>> PROCESSING OXYGEN CORRECTIONS FOR FLOAT: ", float_id)
  message("    Type: ", inst_float, " | Error Config: ", DOXY_adj_err, " umol/kg")
  message("    Using Slope: ", active_slope, " | Drift: ", active_drift, "%/yr | Launch Date: ", active_launch_date)
  message("==============================================================================\n")
  
  if (all(is.na(df_all$DOXY)) && all(is.na(df_all$DOXY_ADJUSTED))) {
    warning("No valid DOXY or DOXY_ADJUSTED data found in matched_result_v2 for WMO: ", active_wmo)
    next
  }
  
  # 1. Compute relative days since launch & gain per row
  delta_days <- as.numeric(difftime(df_all$DATE_TIME, df_all$launch_date, units = "days"))
  gain       <- coalesce(df_all$doxy_slope, 1) * (1 + (coalesce(df_all$doxy_drift, 0) / 100) * (delta_days / 365.25))
  df_all$DOXY_CALCULATED_ADJ <- df_all$DOXY * gain
  
  # Calculate Seawater Density (kg/m3) and DOXY_ADJUSTED_ERROR (umol/kg)
  df_all <- df_all %>%
    mutate(
      DENSITY = oce::swRho(salinity = PSAL, temperature = TEMP, pressure = PRES_ADJUSTED),
      DOXY_ADJUSTED_ERROR = write_DOXY_adjusted_error(
        nc = NULL, 
        profile_idx = iprof, 
        err_mbar = DOXY_adj_err, 
        psal = PSAL, 
        temp = TEMP, 
        pres = PRES_ADJUSTED, 
        dens = DENSITY, 
        DOXY_adj = DOXY_CALCULATED_ADJ
      )
    )
  
  # 2. Hook Flagging for deepest pressure levels (Flag QC = 4 if diff > 2)
  all_cycles <- sort(unique(df_all$CYCLE_NUMBER[!is.na(df_all$CYCLE_NUMBER)]))
  
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
  
  hook_corrections_log <- df_all %>%
    filter(HOOK_FLAG_4) %>%
    select(
      Float_ID = float_ids,
      Cycle = CYCLE_NUMBER,
      Pres_dbar = PRES_ADJUSTED,
      Original_DOXY = DOXY_CALCULATED_ADJ,
      Difference = deepest_diff
    )
  
  if (nrow(hook_corrections_log) > 0) {
    cat("\n------------------------------------------------------------------------------\n")
    cat("   HOOK FLAGGED AS BAD (QC = 4) FOR DEEPEST PRES LEVEL (> 2 THRESHOLD)\n")
    cat("------------------------------------------------------------------------------\n")
    print(hook_corrections_log)
  } else {
    cat("\n--> No deepest level hook flagging (> 2 threshold) required.\n")
  }
  
  # ==============================================================================
  # ARGO COOKBOOK DOXY QUALITY CONTROL TESTS (RANGE, SPIKE, STUCK VALUE, HOOK)
  # ==============================================================================
  
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
      
      spike_fail <- rep(FALSE, n_obs)
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
      qc_flags[pres_qc_9_mask_doxy] <- df_cyc$PRES_QC[pres_qc_9_mask_doxy]
      
      df_cyc$DOXY_FINAL_QC <- qc_flags
      return(df_cyc)
    }) %>%
    ungroup() %>%
    select(-any_of(c("HOOK_FLAG_4", "is_valid_pres", "max_pres_val", "sec_pres_val", "doxy_deepest", "doxy_second", "deepest_diff")))
  
  mode_levels <- c("R-mode", "A-mode", "D-mode", "Bottle")
  color_mapping    <- c("R-mode" = "#5C4033", "A-mode" = "#FF7F00", "D-mode" = "#33A02C", "Bottle" = "#000000")
  linetype_mapping <- c("R-mode" = "solid",   "A-mode" = "solid",   "D-mode" = "solid",   "Bottle" = "blank")
  shape_mapping    <- c("R-mode" = NA,        "A-mode" = NA,        "D-mode" = NA,        "Bottle" = 16)
  
  df_bottle_clean <- Bottle_Vali %>%
    filter(Float_num == float_id & !is.na(DOXY_bottle) & !is.na(PRES)) %>%
    select(PRES_ADJUSTED = PRES, Oxygen = DOXY_bottle) %>%
    mutate(Mode = factor("Bottle", levels = mode_levels))
  
  # ==============================================================================
  # PER-CYCLE PROFILES: CHLA & DOXY MODE COMPARISON PLOTS
  # ==============================================================================
  
  message("\nProcessing CHLA mode plots for Float ", float_id, " across ", length(all_cycles), " cycles...")
  
  for (cyc in all_cycles) {
    df_sub <- df_all %>% filter(CYCLE_NUMBER == cyc)
    if (nrow(df_sub) == 0) next
    
    cyc_lat      <- df_sub$LATITUDE[1]
    cyc_lon      <- df_sub$LONGITUDE[1]
    cyc_date_str <- format(df_sub$DATE_TIME[1], "%Y-%m-%d")
    
    chla_target_cols <- c("CHLA", "CHLA_ADJUSTED", "CHLA_FINAL")
    
    df_chla_cycle_clean <- df_sub %>%
      pivot_longer(
        cols      = any_of(chla_target_cols),
        names_to  = "Mode_Raw",
        values_to = "Chlorophyll"
      ) %>%
      mutate(
        Mode_Name = case_when(
          Mode_Raw == "CHLA"          ~ "R-mode",
          Mode_Raw == "CHLA_ADJUSTED" ~ "A-mode",
          Mode_Raw == "CHLA_FINAL"    ~ "D-mode"
        ),
        Mode = factor(Mode_Name, levels = c("R-mode", "A-mode", "D-mode"))
      ) %>%
      select(PRES, Chlorophyll, Mode) %>%
      filter(!is.na(Chlorophyll) & !is.na(PRES)) %>%
      arrange(Mode, PRES)
    
    if (nrow(df_chla_cycle_clean) == 0) next
    
    p_chla_cycle <- ggplot(df_chla_cycle_clean, aes(x = Chlorophyll, y = PRES, color = Mode, linetype = Mode)) +
      geom_path(linewidth = 0.6, alpha = 0.9) +
      scale_color_manual(values = color_mapping[c("R-mode", "A-mode", "D-mode")], drop = FALSE) +
      scale_linetype_manual(values = linetype_mapping[c("R-mode", "A-mode", "D-mode")], drop = FALSE) +
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
        legend.key.height = unit(0.5, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.8), color = "gray80"),
        plot.caption      = element_text(size = 9, color = "gray40", hjust = 1)
      ) +
      guides(
        color = guide_legend(
          override.aes = list(
            linetype  = c("solid", "solid", "solid"),
            linewidth = c(0.8, 0.8, 0.8)
          )
        ),
        linetype = "none"
      )
    
    chla_file_name <- file.path(dir_chla_plots, sprintf("Float_%d_Cycle_%03d_CHLA_Profile.png", float_id, cyc))
    ggsave(filename = chla_file_name, plot = p_chla_cycle, width = 9.5, height = 7, dpi = 300)
  }
  
  message("\nProcessing DOXY plots for Float ", float_id, " across ", length(all_cycles), " cycles...")
  
  for (cyc in all_cycles) {
    
    df_sub <- df_all %>% filter(CYCLE_NUMBER == cyc)
    if (nrow(df_sub) == 0) next
    
    cyc_lat      <- df_sub$LATITUDE[1]
    cyc_lon      <- df_sub$LONGITUDE[1]
    cyc_date_str <- format(df_sub$DATE_TIME[1], "%Y-%m-%d")
    
    all_target_cols <- c("DOXY", "DOXY_ADJUSTED", "DOXY_CALCULATED_ADJ")
    
    df_cycle_clean <- df_sub %>%
      pivot_longer(
        cols      = all_of(all_target_cols),
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
      select(PRES_ADJUSTED, Oxygen, Mode) %>%
      filter(!is.na(Oxygen) & !is.na(PRES_ADJUSTED))
    
    df_plot_combined <- if (cyc == 1) {
      bind_rows(df_cycle_clean, df_bottle_clean)
    } else {
      df_cycle_clean
    }
    
    df_plot_combined <- df_plot_combined %>% arrange(Mode, PRES_ADJUSTED)
    
    p_cycle <- ggplot() +
      geom_path(
        data = filter(df_plot_combined, Mode != "Bottle"),
        aes(x = Oxygen, y = PRES_ADJUSTED, color = Mode, linetype = Mode),
        linewidth = 0.6, alpha = 0.9
      ) +
      geom_point(
        data = filter(df_plot_combined, Mode == "Bottle"),
        aes(x = Oxygen, y = PRES_ADJUSTED, color = Mode, shape = Mode),
        size = 2.5, alpha = 0.95
      ) +
      scale_color_manual(values = color_mapping, drop = FALSE) +
      scale_linetype_manual(values = linetype_mapping, drop = FALSE) +
      scale_shape_manual(values = shape_mapping, drop = FALSE) +
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
        legend.key.height = unit(0.5, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.8), color = "gray80"),
        plot.caption      = element_text(size = 9, color = "gray40", hjust = 1)
      ) +
      guides(
        color = guide_legend(
          override.aes = list(
            shape     = shape_mapping,
            linetype  = linetype_mapping,
            linewidth = c(0.8, 0.8, 0.8, 0),
            size      = c(NA, NA, NA, 2.5)
          )
        ),
        linetype = "none",
        shape    = "none"
      )
    
    file_name <- file.path(dir_doxy_plots, sprintf("Float_%d_Cycle_%03d_Profile.png", float_id, cyc))
    
    ggsave(
      filename = file_name,
      plot     = p_cycle,
      width    = 9.5,
      height   = 7,
      dpi      = 300
    )
  }
  
  # ==============================================================================
  # DUAL-PANEL DMQC SUMMARY PLOTS (PROFILE ON LEFT vs. QC FLAGS ON RIGHT)
  # ==============================================================================
  
  message(sprintf("Generating side-by-side QC flag summary plots for Float: %d...", float_id))
  
  df_summary_base <- df_all %>%
    filter(!is.na(PRES)) %>%
    mutate(
      CYCLE_NUMBER = as.numeric(CYCLE_NUMBER),
      PRES         = as.numeric(PRES)
    )
  
  # 1. SUMMARY PLOT FOR CHLA_FINAL
  if ("CHLA_FINAL" %in% names(df_summary_base) && any(!is.na(df_summary_base$CHLA_FINAL))) {
    
    df_chla_sum <- df_summary_base %>%
      filter(!is.na(CHLA_FINAL)) %>%
      mutate(CHLA_FINAL_QC = as.character(CHLA_FINAL_QC))
    
    p_chla_prof <- ggplot(df_chla_sum, aes(x = CHLA_FINAL, y = PRES, color = CYCLE_NUMBER, group = CYCLE_NUMBER)) +
      geom_path(linewidth = 0.25, alpha = 0.7) +
      scale_y_reverse(limits = c(max(df_chla_sum$PRES, na.rm = TRUE), 0)) +
      scale_x_continuous(position = "top") +
      scale_color_viridis_c(name = "Cycle") +
      labs(
        title = paste0("Float ", float_id, " - CHLA_FINAL Profiles"),
        x     = expression(Chlorophyll ~ a ~ (mg/m^3)),
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_chla_qc <- ggplot(df_chla_sum, aes(x = CHLA_FINAL_QC, y = PRES, color = CHLA_FINAL_QC)) +
      geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width = 0.15, height = 0)) +
      scale_y_reverse(limits = c(max(df_chla_sum$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(
        title = "CHLA QC Flags vs Depth",
        x     = "Argo QC Flag",
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_chla_combined <- (p_chla_prof + p_chla_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - CHLA DMQC Flag & Profile Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    chla_qc_summary_file <- file.path(dir_qc_summary, sprintf("Float_%d_CHLA_QC_Summary.png", float_id))
    ggsave(chla_qc_summary_file, plot = p_chla_combined, width = 12, height = 7, dpi = 300)
  }
  
  # 2. SUMMARY PLOT FOR CHLA_FLUORESCENCE
  if ("CHLA_FLUORESCENCE" %in% names(df_summary_base) && any(!is.na(df_summary_base$CHLA_FLUORESCENCE))) {
    
    df_fluo_sum <- df_summary_base %>%
      filter(!is.na(CHLA_FLUORESCENCE)) %>%
      mutate(CHLA_FLUORESCENCE_QC = as.character(CHLA_FLUORESCENCE_QC))
    
    p_fluo_prof <- ggplot(df_fluo_sum, aes(x = CHLA_FLUORESCENCE, y = PRES, color = CYCLE_NUMBER, group = CYCLE_NUMBER)) +
      geom_path(linewidth = 0.25, alpha = 0.7) +
      scale_y_reverse(limits = c(max(df_fluo_sum$PRES, na.rm = TRUE), 0)) +
      scale_x_continuous(position = "top") +
      scale_color_viridis_c(option = "mako", name = "Cycle") +
      labs(
        title = paste0("Float ", float_id, " - CHLA_FLUORESCENCE Profiles"),
        x     = "Chlorophyll Fluorescence (RU)",
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_fluo_qc <- ggplot(df_fluo_sum, aes(x = CHLA_FLUORESCENCE_QC, y = PRES, color = CHLA_FLUORESCENCE_QC)) +
      geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width = 0.15, height = 0)) +
      scale_y_reverse(limits = c(max(df_fluo_sum$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(
        title = "CHLA_FLUORESCENCE QC Flags vs Depth",
        x     = "Argo QC Flag",
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_fluo_combined <- (p_fluo_prof + p_fluo_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - CHLA_FLUORESCENCE DMQC Flag & Profile Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    fluo_qc_summary_file <- file.path(dir_qc_summary, sprintf("Float_%d_CHLA_FLUORESCENCE_QC_Summary.png", float_id))
    ggsave(fluo_qc_summary_file, plot = p_fluo_combined, width = 12, height = 7, dpi = 300)
  }
  
  # 3. SUMMARY PLOT FOR BBP700_FINAL
  if ("BBP700_FINAL" %in% names(df_summary_base) && any(!is.na(df_summary_base$BBP700_FINAL) & df_summary_base$BBP700_FINAL != 99999.0)) {
    
    df_bbp_sum <- df_summary_base %>%
      filter(!is.na(BBP700_FINAL) & BBP700_FINAL != 99999.0) %>%
      mutate(BBP700_FINAL_QC = as.character(BBP700_FINAL_QC))
    
    p_bbp_prof <- ggplot(df_bbp_sum, aes(x = BBP700_FINAL, y = PRES, color = CYCLE_NUMBER, group = CYCLE_NUMBER)) +
      geom_path(linewidth = 0.25, alpha = 0.7) +
      scale_y_reverse(limits = c(max(df_bbp_sum$PRES, na.rm = TRUE), 0)) +
      scale_x_continuous(position = "top") +
      scale_color_viridis_c(option = "turbo", name = "Cycle") +
      labs(
        title = paste0("Float ", float_id, " - BBP700_FINAL Profiles"),
        x     = expression(b[bp](700) ~ (m^{-1})),
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_bbp_qc <- ggplot(df_bbp_sum, aes(x = BBP700_FINAL_QC, y = PRES, color = BBP700_FINAL_QC)) +
      geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width = 0.15, height = 0)) +
      scale_y_reverse(limits = c(max(df_bbp_sum$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(
        title = "BBP700 QC Flags vs Depth",
        x     = "Argo QC Flag",
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_bbp_combined <- (p_bbp_prof + p_bbp_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - BBP700 DMQC Flag & Profile Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    bbp_qc_summary_file <- file.path(dir_qc_summary, sprintf("Float_%d_BBP700_QC_Summary.png", float_id))
    ggsave(bbp_qc_summary_file, plot = p_bbp_combined, width = 12, height = 7, dpi = 300)
  }
  
  # 4. SUMMARY PLOT FOR DOXY_FINAL
  if ("DOXY_FINAL" %in% names(df_summary_base) && any(!is.na(df_summary_base$DOXY_FINAL))) {
    
    df_doxy_sum <- df_summary_base %>%
      filter(!is.na(DOXY_FINAL)) %>%
      mutate(DOXY_FINAL_QC = as.character(DOXY_FINAL_QC))
    
    p_doxy_prof <- ggplot(df_doxy_sum, aes(x = DOXY_FINAL, y = PRES, color = CYCLE_NUMBER, group = CYCLE_NUMBER)) +
      geom_path(linewidth = 0.25, alpha = 0.7) +
      scale_y_reverse(limits = c(max(df_doxy_sum$PRES, na.rm = TRUE), 0)) +
      scale_x_continuous(position = "top") +
      scale_color_viridis_c(option = "plasma", name = "Cycle") +
      labs(
        title = paste0("Float ", float_id, " - DOXY_FINAL Profiles"),
        x     = expression(Dissolved ~ Oxygen ~ (mu * "mol/kg")),
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_doxy_qc <- ggplot(df_doxy_sum, aes(x = DOXY_FINAL_QC, y = PRES, color = DOXY_FINAL_QC)) +
      geom_point(alpha = 0.6, size = 1.5, position = position_jitter(width = 0.15, height = 0)) +
      scale_y_reverse(limits = c(max(df_doxy_sum$PRES, na.rm = TRUE), 0)) +
      scale_color_manual(values = qc_color_scale, drop = FALSE, name = "QC Flag") +
      labs(
        title = "DOXY QC Flags vs Depth",
        x     = "Argo QC Flag",
        y     = "Depth / Pressure (dbar)"
      ) +
      theme_bw(base_size = 12) +
      theme(aspect.ratio = 1.2, legend.position = "bottom")
    
    p_doxy_combined <- (p_doxy_prof + p_doxy_qc) +
      plot_annotation(
        title    = paste0("BGC-Argo Float ", float_id, " - DOXY DMQC Flag & Profile Summary"),
        subtitle = paste0("Generated on: ", today_str),
        theme    = theme(plot.title = element_text(face = "bold", size = 14))
      )
    
    doxy_qc_summary_file <- file.path(dir_qc_summary, sprintf("Float_%d_DOXY_QC_Summary.png", float_id))
    ggsave(doxy_qc_summary_file, plot = p_doxy_combined, width = 12, height = 7, dpi = 300)
  }
  
  # ==============================================================================
  # EXPORT SINGLE COMBINED CSV FOR CHLA, BBP700, FLUORESCENCE & DOXY DMODE DATA
  # ==============================================================================
  
  combo_output_path <- file.path(combo_output_dir, sprintf("CHLA_BBP_DOXY_%d.csv", float_id))
  
  # 1. Combine WMOID and float_num into float_num & coerce logical columns
  df_export_prep <- df_all %>%
    mutate(
      float_num = coalesce(as.character(WMOID), as.character(float_num)),
      across(where(is.logical), as.numeric)
    )
  
  # --- Apply Float ID and Cycle Specific Overrides ---
  df_export_prep <- df_export_prep %>%
    mutate(
      override_condition = (
        (float_num == "4903624") |
          (float_num == "2904010" & CYCLE_NUMBER >= 49) |
          (float_num == "2904011" & CYCLE_NUMBER >= 24)
      ),
      # Target CHLA_FINAL, BBP700_FINAL, and associated QC flags
      CHLA_FINAL      = ifelse(override_condition, 99999.0, CHLA_FINAL),
      BBP700_FINAL    = ifelse(override_condition, 99999.0, BBP700_FINAL),
      CHLA_FINAL_QC   = ifelse(override_condition, 9, CHLA_FINAL_QC),
      BBP700_FINAL_QC = ifelse(override_condition, 9, BBP700_FINAL_QC),
      
      # Additional requested parameters set to 99999.0 on override condition
      CHLA                       = ifelse(override_condition, 99999.0, CHLA),
      CHLA_ADJUSTED              = ifelse(override_condition, 99999.0, CHLA_ADJUSTED),
      CHLA_FLUORESCENCE          = ifelse(override_condition, 99999.0, CHLA_FLUORESCENCE),
      CHLA_FLUORESCENCE_QC = ifelse(override_condition, 9, CHLA_FLUORESCENCE_QC),
      CHLA_FLUORESCENCE_ADJUSTED = ifelse(override_condition, 99999.0, CHLA_FLUORESCENCE_ADJUSTED),
      CHLA_FLUORESCENCE_ADJUSTED_QC = ifelse(override_condition, 9, CHLA_FLUORESCENCE_ADJUSTED_QC),
      CHLA_NoLUT                 = ifelse(override_condition, 99999.0, CHLA_NoLUT),
      BBP700                     = ifelse(override_condition, 99999.0, BBP700),
      BBP700_QC = ifelse(override_condition, 9, BBP700_QC),
      SCALE_CHLA                 = ifelse(override_condition, 99999.0, SCALE_CHLA),
      DARK_CHLA                  = ifelse(override_condition, 99999.0, DARK_CHLA)
    ) %>%
    select(-override_condition)
  
  # 2. Apply rules for CHLA_FINAL and CHLA_FINAL_QC based on CHLA_ADJUSTED / CHLA_ADJUSTED_QC
  if ("CHLA_ADJUSTED" %in% names(df_export_prep)) {
    df_export_prep <- df_export_prep %>%
      mutate(
        CHLA_FINAL = ifelse(is.na(CHLA_ADJUSTED) | CHLA_ADJUSTED == 99999.0, 99999.0, CHLA_FINAL),
        CHLA_FINAL_QC = ifelse(is.na(CHLA_ADJUSTED) | CHLA_ADJUSTED == 99999.0, 9, CHLA_FINAL_QC)
      )
  } else {
    df_export_prep <- df_export_prep %>%
      mutate(
        CHLA_FINAL = ifelse(is.na(CHLA_FINAL), 99999.0, CHLA_FINAL),
        CHLA_FINAL_QC = ifelse(is.na(CHLA_FINAL), 9, CHLA_FINAL_QC)
      )
  }
  
  if ("CHLA_ADJUSTED_QC" %in% names(df_export_prep)) {
    df_export_prep <- df_export_prep %>%
      mutate(
        CHLA_FINAL_QC = ifelse(is.na(CHLA_ADJUSTED_QC) | CHLA_ADJUSTED_QC %in% c(9, "9"), 9, CHLA_FINAL_QC)
      )
  }
  
  # 3. Clean up columns, rename to uppercase, reorder, filter out PRES == 99999, and fill missing values with 99999.0
  df_combo_export <- df_export_prep %>%
    select(-any_of(c(
      "WMOID", "ID", "DOXY_QC", "BBP700_QC_MANUAL",
      "DOXY_ADJUSTED", "DOXY_ADJUSTED_QC",
      "DOXY_CALCULATED_ADJ", "float_ids", "TIME",
      "launch_date", "lauch_date", "LAUNCH_DATE"
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
    filter(PRES != 99999) %>%
    mutate(across(everything(), ~ tidyr::replace_na(., 99999.0)))
  
  write.csv(df_combo_export, combo_output_path, row.names = FALSE)
  message(sprintf("\n>>> SUCCESSFULLY EXPORTED SINGLE COMBINED CSV -> %s", combo_output_path))
  message("\n=== COMPLETE DMQC PROCESSING FOR FLOAT: ", float_id, " ===")
}