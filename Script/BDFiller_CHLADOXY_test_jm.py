#!/usr/bin/env python3
# coding: utf-8

# ==============================================================================
# SECTION 1: Library Imports & Helper Functions
# ==============================================================================
import datetime
import gc
import glob
import os
import re
import shutil
import sys
from datetime import datetime as dt, timezone

import gsw
import netCDF4
import netCDF4 as nc
import numpy as np
import pandas as pd


FLOAT_TYPES = {
    4903904: "aoml_navis",
}
# Derived list of WMO Float IDs to process
WMO_FLOAT_IDS = list(FLOAT_TYPES.keys())

# Common Configuration Defaults
comment_dmqc_operator_chla = (
    "PRIMARY | https://orcid.org/0009-0006-7862-6267 | Brandon Navarro,"
    " NOAA/AOML;"
)
history_parameter_chla = "CHLA"
history_institution = "AO"  # AO for AOML
history_reference_chla = "MULT"
history_software_chla = "RBIO"
history_software_release_chla = "2025"

scientific_calibration_comment_CHLA = (
    "CHLA adjustment (Slope: specified in http://dx.doi.org/10.13155/35385 and"
    " computed with MLD_LIMIT = 0.03, and following recommendations of Sauzede"
    " et al., 2025 (https://doi.org/10.17882/105732), Quenching: Xing et al.,"
    " 2018, Terrats et al., 2020)"
)
scientific_calibration_comment_CHLA_FLU = (
    "CHLA_FLUORESCENCE  (specified in http://dx.doi.org/10.13155/35385 and"
    " computed with MLD_LIMIT = 0.03)"
)
scientific_calibration_equation_CHLA = (
    "CHLA_ADJUSTED = CHLA_NPQ for PRES in [0, ZMaxFluo ], CHLA_ADJUSTED ="
    " ((FLUORESCENCE_CHLA-MEDIAN(PRELIM_DARK_CHLA)*SCALE_CHLA)/PHYSIO_RATIO"
)
scientific_calibration_coefficient_CHLA = "PHYSIO_RATIO=1.0"
CHLA_Adjusted_ERROR_est = 0.07

iprof_chla = 0
data_state_indicator = ["2", "C", "", ""]
parameter_data_mode = "D"


# ==============================================================================
# SECTION 2: Helper Functions - CHLA Processing
# ==============================================================================


def create_working_bd_file(filename):
  """Create a working copy of the B file with a leading 'w_' in the name."""
  path, name = os.path.split(filename)
  bd_name = name.replace("BR", "BD")
  w_filename = os.path.join(path, f"w_{bd_name}")
  shutil.copyfile(filename, w_filename)
  return w_filename


def get_profile_chla(filename):
  """Extract the profile index from filename, return as int."""
  profile = filename[-6:-3]
  return int(profile)


def organize_b_files(bd_files, br_files):
  """Sort B*nc files by profile. If both BR and BD exist, prefer BD."""
  ptr_bd = 0
  ptr_br = 0

  max_br = get_profile_chla(br_files[-1]) if br_files else 0
  max_bd = get_profile_chla(bd_files[-1]) if bd_files else 0
  max_prof = max(max_bd, max_br)

  sorted_b_files = []
  for idx in range(1, max_prof + 1):
    file_found = False
    for ptr in range(ptr_bd, len(bd_files)):
      if get_profile_chla(bd_files[ptr]) == idx:
        sorted_b_files.append(bd_files[ptr])
        ptr_bd = ptr + 1
        file_found = True
        break
    if not file_found:
      for ptr in range(ptr_br, len(br_files)):
        if get_profile_chla(br_files[ptr]) == idx:
          sorted_b_files.append(br_files[ptr])
          ptr_br = ptr + 1
          break
  return sorted_b_files


def update_history_chla(nc_ds, dct, iprof_idx):
  """Update HISTORY array entries natively using netCDF4."""
  hix = nc_ds.dimensions["N_HISTORY"].size
  for name, value in dct.items():
    if name in nc_ds.variables:
      char_len = nc_ds.dimensions[nc_ds[name].dimensions[-1]].size
      padded_val = str(value).ljust(char_len)[:char_len]
      nc_ds[name][hix, iprof_idx, :] = nc.stringtochar(
          np.array(padded_val, dtype=f"S{char_len}")
      )


def write_history_chla(bgc_file, iprof_idx):
  """Write global attributes and HISTORY variables for CHLA."""
  bgc_file.history = datetime.datetime.utcnow().strftime(
      "%Y-%m-%dT%H:%M:%SZ creation"
  )
  bgc_file.setncattr("comment_dmqc_operator", comment_dmqc_operator_chla)

  history_step = "ARSQ"
  history_action = "IP"
  UTCcurrent = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")

  update_history_chla(
      bgc_file,
      {
          "HISTORY_INSTITUTION": history_institution,
          "HISTORY_STEP": history_step,
          "HISTORY_SOFTWARE": history_software_chla,
          "HISTORY_SOFTWARE_RELEASE": history_software_release_chla,
          "HISTORY_REFERENCE": history_reference_chla,
          "HISTORY_DATE": UTCcurrent,
          "HISTORY_ACTION": history_action,
          "HISTORY_PARAMETER": history_parameter_chla,
      },
      iprof_idx,
  )

  bgc_file.variables["DATE_UPDATE"][:] = nc.stringtochar(
      np.array(UTCcurrent, dtype="S14")
  )


def write_parameter_data_mode_chla(bgc_file, iprof_idx=0):
  """Set PARAMETER_DATA_MODE and DATA_MODE synchronously across all N_PROF profiles."""
  n_prof = bgc_file.dimensions["N_PROF"].size
  n_param = bgc_file.dimensions["N_PARAM"].size

  pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
  data_mode = bgc_file.variables["DATA_MODE"][:]

  for iprof in range(n_prof):
    param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
    if isinstance(param_mat, np.ma.MaskedArray):
      param_mat = param_mat.filled(b" ")

    if iprof == iprof_idx:
      data_mode[iprof] = "D"
      for j in range(n_param):
        param_str = "".join([
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in param_mat[j]
        ]).strip()
        if param_str.startswith("CHLA"):
          pdm[iprof, j] = "D"
        elif param_str and pdm[iprof, j] == "R":
          pdm[iprof, j] = "A"
    else:
      if data_mode[iprof] == "D":
        for j in range(n_param):
          param_str = "".join([
              c.decode("utf-8", errors="ignore")
              if isinstance(c, bytes)
              else str(c)
              for c in param_mat[j]
          ]).strip()
          if param_str and pdm[iprof, j] == "R":
            pdm[iprof, j] = "A"

  bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
  bgc_file.variables["DATA_MODE"][:] = data_mode


def get_profile_qc_grade(qc_masked_array):
  """Calculate Argo profile QC letter grade (A-F, Z)."""
  if hasattr(qc_masked_array, "compressed"):
    unmasked_vals = qc_masked_array.compressed()
  else:
    unmasked_vals = qc_masked_array

  valid_qcs = []
  for q in unmasked_vals:
    if isinstance(q, np.ma.core.MaskedConstant) or q is np.ma.masked:
      continue
    s = q.decode("utf-8").strip() if isinstance(q, bytes) else str(q).strip()
    if s != "" and s != "9":
      valid_qcs.append(s)

  total_points = len(valid_qcs)
  if total_points == 0:
    return "Z"

  bad_count = sum(1 for q in valid_qcs if q in ["3", "4"])
  good_count = total_points - bad_count
  pct_good = (good_count / total_points) * 100.0

  if pct_good == 100.0:
    return "A"
  elif 75.0 <= pct_good < 100.0:
    return "B"
  elif 50.0 <= pct_good < 75.0:
    return "C"
  elif 25.0 <= pct_good < 50.0:
    return "D"
  elif 0.0 < pct_good < 25.0:
    return "E"
  else:
    return "F"


def write_scientific_calib_chla(bgc_file, idx_profile, bio_dmqc_csv_path):
  """Write SCIENTIFIC_CALIB_* variables for CHLA safely against missing early cycles."""
  df_bio = pd.read_csv(bio_dmqc_csv_path)

  dark_cols = sorted(
      [col for col in df_bio.columns if col.startswith("MIN_FLUOCHLA_CYCLE")]
  )
  if dark_cols:
    dark_vals = df_bio[dark_cols].iloc[0].dropna().astype(int).tolist()
    dark_str = " ".join(map(str, dark_vals)) if dark_vals else "NA"
  else:
    dark_str = "NA"

  scale_val = (
      df_bio["SCALE_CHLA"].iloc[0] if "SCALE_CHLA" in df_bio.columns else "NA"
  )
  cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

  if not cycle_df.empty and "PHYSIO_RATIO" in cycle_df.columns:
    physio_val = cycle_df["PHYSIO_RATIO"].iloc[0]
  else:
    physio_val = "1"

  calib_coefficient = (
      f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val},"
      f" PHYSIO_RATIO = {physio_val}"
  )
  calib_coefficient_flu = (
      f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}"
  )
  UTCcurrent = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")

  str256_len = bgc_file.dimensions["STRING256"].size

  SciCalComArray_CHLA = np.ma.empty(shape=(str256_len), dtype="|S1")
  SciCalComArray_CHLA[:] = ""
  SciCalComArray_CHLA.mask = True
  SciCalComArray_CHLA[: len(scientific_calibration_comment_CHLA)] = list(
      scientific_calibration_comment_CHLA
  )

  SciCalComArray_CHLA_FLU = np.ma.empty(shape=(str256_len), dtype="|S1")
  SciCalComArray_CHLA_FLU[:] = ""
  SciCalComArray_CHLA_FLU.mask = True
  SciCalComArray_CHLA_FLU[: len(scientific_calibration_comment_CHLA_FLU)] = (
      list(scientific_calibration_comment_CHLA_FLU)
  )

  SciCalEquArray_CHLA = np.ma.empty(shape=(str256_len), dtype="|S1")
  SciCalEquArray_CHLA[:] = ""
  SciCalEquArray_CHLA.mask = True
  SciCalEquArray_CHLA[: len(scientific_calibration_equation_CHLA)] = list(
      scientific_calibration_equation_CHLA
  )

  SciCalCoeArray_CHLA = np.ma.empty(shape=(str256_len), dtype="|S1")
  SciCalCoeArray_CHLA[:] = ""
  SciCalCoeArray_CHLA.mask = True
  SciCalCoeArray_CHLA[: len(calib_coefficient)] = list(calib_coefficient)

  SciCalCoeArray_CHLA_FLU = np.ma.empty(shape=(str256_len), dtype="|S1")
  SciCalCoeArray_CHLA_FLU[:] = ""
  SciCalCoeArray_CHLA_FLU.mask = True
  SciCalCoeArray_CHLA_FLU[: len(calib_coefficient_flu)] = list(
      calib_coefficient_flu
  )

  SciCalDateArray = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))

  n_prof_size = bgc_file.dimensions["N_PROF"].size
  n_param_size = bgc_file.dimensions["N_PARAM"].size
  chla_found = False

  for iprof_idx in range(n_prof_size):
    param_list = bgc_file.variables["STATION_PARAMETERS"][iprof_idx].data.astype(
        str
    )
    for j in range(n_param_size):
      param_str = "".join(param_list[j]).strip()
      if param_str.startswith("CHLA_F"):
        bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = (
            SciCalComArray_CHLA_FLU
        )
        bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][
            iprof_idx, 0, j, :
        ] = SciCalCoeArray_CHLA_FLU
        bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = (
            SciCalDateArray
        )
      elif param_str.startswith("CHLA"):
        bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = (
            SciCalComArray_CHLA
        )
        bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = (
            SciCalEquArray_CHLA
        )
        bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][
            iprof_idx, 0, j, :
        ] = SciCalCoeArray_CHLA
        bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = (
            SciCalDateArray
        )
        chla_found = True

  if not chla_found:
    print(
        "Warning: CHLA parameter was not found in file for profile cycle"
        f" {idx_profile}."
    )


def write_chla_BBP_adjusted(
    bgc_file, idx_profile, bio_dmqc_csv_path, iprof_idx=iprof_chla
):
  """Populate CHLA_ADJUSTED and CHLA_FLUORESCENCE_ADJUSTED without value/error mismatches."""
  df_bio = pd.read_csv(bio_dmqc_csv_path)
  df_bio = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

  n_levels = bgc_file.dimensions["N_LEVELS"].size

  CHLA_Adjusted_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  CHLA_Adjusted_Array[:] = 99999.0
  CHLA_Adjusted_Array.mask = True

  CHLA_AdjustedQC_Array = np.ma.empty(shape=(n_levels,), dtype="|S1")
  CHLA_AdjustedQC_Array[:] = b"9"
  CHLA_AdjustedQC_Array.mask = True

  CHLA_Adjusted_ERROR_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  CHLA_Adjusted_ERROR_Array[:] = 99999.0
  CHLA_Adjusted_ERROR_Array.mask = True

  CHLA_FLUORESCENCE_Adjusted_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  CHLA_FLUORESCENCE_Adjusted_Array[:] = 99999.0
  CHLA_FLUORESCENCE_Adjusted_Array.mask = True

  CHLA_FLUORESCENCE_AdjustedQC_Array = np.ma.empty(shape=(n_levels,), dtype="|S1")
  CHLA_FLUORESCENCE_AdjustedQC_Array[:] = b"9"
  CHLA_FLUORESCENCE_AdjustedQC_Array.mask = True

  CHLA_FLUORESCENCE_Adjusted_ERROR_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  CHLA_FLUORESCENCE_Adjusted_ERROR_Array[:] = 99999.0
  CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask = True

  assigned_nc_pres_vals = set()

  for row in range(len(df_bio)):
    row_data = df_bio.iloc[row]
    csv_pres = np.float32(row_data["PRES"])

    for i in range(n_levels):
      nc_pres = np.float32(bgc_file.variables["PRES"][iprof_idx, i])

      if nc_pres in assigned_nc_pres_vals:
        continue

      if np.isclose(csv_pres, nc_pres, atol=0.05):
        raw_qc = row_data["CHLA_FINAL_QC"]
        qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"

        if qc_str in ["4", "9"]:
          CHLA_Adjusted_Array[i] = 99999.0
          CHLA_Adjusted_ERROR_Array[i] = 99999.0
          CHLA_Adjusted_Array.mask[i] = True
          CHLA_Adjusted_ERROR_Array.mask[i] = True
        else:
          CHLA_Adjusted_Array[i] = np.float32(row_data["CHLA_FINAL"])
          CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
          CHLA_Adjusted_Array.mask[i] = False
          CHLA_Adjusted_ERROR_Array.mask[i] = False

        CHLA_AdjustedQC_Array[i] = qc_str.encode("utf-8")
        CHLA_AdjustedQC_Array.mask[i] = False

        fluo_val = (
            row_data["CHLA_FLUORESCENCE"]
            if "CHLA_FLUORESCENCE" in row_data
            else np.nan
        )
        if pd.isna(fluo_val) or fluo_val == 99999.0 or qc_str in ["4", "9"]:
          CHLA_FLUORESCENCE_Adjusted_Array[i] = 99999.0
          CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = 99999.0
          CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = True
          CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = True
        else:
          CHLA_FLUORESCENCE_Adjusted_Array[i] = np.float32(fluo_val)
          CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = np.float32(
              CHLA_Adjusted_ERROR_est
          )
          CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = False
          CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = False

        CHLA_FLUORESCENCE_AdjustedQC_Array[i] = qc_str.encode("utf-8")
        CHLA_FLUORESCENCE_AdjustedQC_Array.mask[i] = False

        assigned_nc_pres_vals.add(nc_pres)
        break

  bgc_file.variables["CHLA_ADJUSTED"][iprof_idx] = CHLA_Adjusted_Array
  bgc_file.variables["CHLA_ADJUSTED_QC"][iprof_idx] = CHLA_AdjustedQC_Array
  bgc_file.variables["CHLA_ADJUSTED_ERROR"][iprof_idx] = (
      CHLA_Adjusted_ERROR_Array
  )

  if "CHLA_FLUORESCENCE_ADJUSTED" in bgc_file.variables:
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED"][iprof_idx] = (
        CHLA_FLUORESCENCE_Adjusted_Array
    )
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_QC"][iprof_idx] = (
        CHLA_FLUORESCENCE_AdjustedQC_Array
    )
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_ERROR"][iprof_idx] = (
        CHLA_FLUORESCENCE_Adjusted_ERROR_Array
    )

  profile_chla_qc = get_profile_qc_grade(CHLA_AdjustedQC_Array)
  if "PROFILE_CHLA_QC" in bgc_file.variables:
    bgc_file.variables["PROFILE_CHLA_QC"][iprof_idx] = profile_chla_qc

  if "PROFILE_CHLA_FLUORESCENCE_QC" in bgc_file.variables:
    profile_fluo_qc = get_profile_qc_grade(CHLA_FLUORESCENCE_AdjustedQC_Array)
    bgc_file.variables["PROFILE_CHLA_FLUORESCENCE_QC"][iprof_idx] = (
        profile_fluo_qc
    )


# ==============================================================================
# SECTION 3: Helper Functions - DOXY Processing
# ==============================================================================


def get_iprof_phys(pres_raw, pres_bgc, target_iprof):
  """Find column index in physical PRES matching target BGC PRES profile."""
  iprof_phys = -1
  target_pres = pres_bgc[:, target_iprof] if pres_bgc.ndim > 1 else pres_bgc
  num_cols = pres_raw.shape[1] if pres_raw.ndim > 1 else 1

  for col in range(num_cols):
    current_pres = pres_raw[:, col] if pres_raw.ndim > 1 else pres_raw
    valid_mask = ~np.isnan(current_pres) & ~np.isnan(target_pres)
    if not np.any(valid_mask):
      continue

    diff = np.max(np.abs(current_pres[valid_mask] - target_pres[valid_mask]))
    if diff < 0.1:
      iprof_phys = col
      break

  if iprof_phys < 0:
    raise ValueError("Matching PRES values not found")

  print(f"Using profile {iprof_phys} of physical file to determine density.")
  return iprof_phys


def create_working_doxy_bd_file(filename, dest_dir):
  """Copy file to destination directory as 'w_BD...' working file."""
  base_name = os.path.basename(filename)
  bd_name = re.sub(r"^BR", "BD", base_name)
  w_filename = os.path.join(dest_dir, f"w_{bd_name}")
  shutil.copyfile(filename, w_filename)
  return w_filename


def get_phys_filename(bgc_filename, base_dir):
  """Locate matching core physical NetCDF file."""
  base_name = os.path.basename(bgc_filename)
  core_phys_name = re.sub(r"^B", "", base_name)

  d_name = re.sub(r"^R", "D", core_phys_name)
  phys_filename = os.path.join(base_dir, d_name)

  if not os.path.exists(phys_filename):
    phys_filename = os.path.join(base_dir, "D", d_name)

  if not os.path.exists(phys_filename):
    r_name = re.sub(r"^D", "R", core_phys_name)
    phys_filename = os.path.join(base_dir, r_name)
    if not os.path.exists(phys_filename):
      phys_filename = os.path.join(base_dir, "R", r_name)

  if not os.path.exists(phys_filename):
    raise FileNotFoundError(
        f"No corresponding phys file found for {bgc_filename} in {base_dir}"
    )

  return phys_filename


def get_phys_raw_pres(phys_filename):
  """Read PRES variable from physical file."""
  with nc.Dataset(phys_filename, "r") as ds:
    return ds.variables["PRES"][:]


def get_dens(phys_filename, verbose=False):
  """Calculate potential density (rho) using TEOS-10 GSW."""
  with nc.Dataset(phys_filename, "r") as ds:
    var_names = ds.variables.keys()

    temp = (
        ds.variables["TEMP_ADJUSTED"][:]
        if "TEMP_ADJUSTED" in var_names
        else ds.variables["TEMP"][:]
    )
    psal = (
        ds.variables["PSAL_ADJUSTED"][:]
        if "PSAL_ADJUSTED" in var_names
        else ds.variables["PSAL"][:]
    )
    pres = ds.variables["PRES"][:]

    dens = gsw.rho_t_exact(SA=psal, t=temp, p=pres)
    return {"dens": dens, "psal": psal, "temp": temp}


def update_history_doxy(nc_ds, dct, iprof_idx):
  """Update HISTORY array entries natively using netCDF4 (Identical to CHLA implementation)."""
  hix = nc_ds.dimensions["N_HISTORY"].size
  for name, value in dct.items():
    if name in nc_ds.variables:
      char_len = nc_ds.dimensions[nc_ds[name].dimensions[-1]].size
      padded_val = str(value).ljust(char_len)[:char_len]
      nc_ds[name][hix, iprof_idx, :] = nc.stringtochar(
          np.array(padded_val, dtype=f"S{char_len}")
      )


def write_history_doxy(ds, profile_idx, inst, ref, comment_op):
  """Update global history attributes for DOXY processing (Identical structure to CHLA)."""
  ds.history = datetime.datetime.utcnow().strftime(
      "%Y-%m-%dT%H:%M:%SZ creation"
  )
  ds.setncattr("comment_dmqc_operator", comment_op)

  history_step = "ARSQ"
  history_action = "IP"
  history_software = "SAGE"
  history_software_release = "2024"
  history_parameter = "DOXY"
  UTCcurrent = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")

  update_history_doxy(
      ds,
      {
          "HISTORY_INSTITUTION": inst,
          "HISTORY_STEP": history_step,
          "HISTORY_SOFTWARE": history_software,
          "HISTORY_SOFTWARE_RELEASE": history_software_release,
          "HISTORY_REFERENCE": ref,
          "HISTORY_DATE": UTCcurrent,
          "HISTORY_ACTION": history_action,
          "HISTORY_PARAMETER": history_parameter,
      },
      profile_idx,
  )

  if "DATE_UPDATE" in ds.variables:
    ds.variables["DATE_UPDATE"][:] = nc.stringtochar(
        np.array(UTCcurrent, dtype="S14")
    )


def write_parameter_data_mode_doxy(bgc_file, iprof_idx=0):
  """Set PARAMETER_DATA_MODE and DATA_MODE synchronously across all N_PROF profiles (using CHLA method)."""
  n_prof = bgc_file.dimensions["N_PROF"].size
  n_param = bgc_file.dimensions["N_PARAM"].size

  pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
  data_mode = bgc_file.variables["DATA_MODE"][:]

  for iprof in range(n_prof):
    param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
    if isinstance(param_mat, np.ma.MaskedArray):
      param_mat = param_mat.filled(b" ")

    if iprof == iprof_idx:
      data_mode[iprof] = "D"
      for j in range(n_param):
        param_str = "".join([
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in param_mat[j]
        ]).strip()
        if param_str.startswith("DOXY"):
          pdm[iprof, j] = "D"
        elif param_str and pdm[iprof, j] == "R":
          pdm[iprof, j] = "A"
    else:
      if data_mode[iprof] == "D":
        for j in range(n_param):
          param_str = "".join([
              c.decode("utf-8", errors="ignore")
              if isinstance(c, bytes)
              else str(c)
              for c in param_mat[j]
          ]).strip()
          if param_str and pdm[iprof, j] == "R":
            pdm[iprof, j] = "A"

  bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
  bgc_file.variables["DATA_MODE"][:] = data_mode

  # DIAGNOSTIC PRINT STATEMENT
  print(
      f"[VERIFICATION] write_parameter_data_mode_doxy matched CHLA method.\n"
      f"  DATA_MODE: {data_mode.tolist()}\n"
      f"  PARAMETER_DATA_MODE: {pdm.tolist()}"
  )


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle):
  """Write DOXY slope and drift calibration coefficients."""
  cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)]
  if cycle_df.empty:
    return

  slope_val = cycle_df["DOXY_SLOPE"].iloc[0]
  drift_val = cycle_df["DOXY_DRIFT"].iloc[0]

  if pd.isna(slope_val) and pd.isna(drift_val):
    return

  s_str = "1" if pd.isna(slope_val) else str(slope_val)
  d_str = "0" if pd.isna(drift_val) else str(drift_val)
  calib_str = f"m={s_str}, d={d_str}"

  var_names = ds.variables.keys()

  if "SCIENTIFIC_CALIB_COEFFICIENT" in var_names:
    params = ds.variables["STATION_PARAMETERS"][:]
    if params.ndim == 3:
      prof_params = params[profile_idx, :, :]
      if isinstance(prof_params, np.ma.MaskedArray):
        prof_params = prof_params.filled(b" ")

      param_strings = []
      for row in prof_params:
        row_chars = [
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in row
        ]
        param_strings.append("".join(row_chars).strip())

      doxy_indices = [
          idx for idx, s in enumerate(param_strings) if s.startswith("DOXY")
      ]

      if doxy_indices:
        calib_var = ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"]
        doxy_idx = doxy_indices[0]
        char_len = calib_var.shape[-1]
        padded_str = calib_str.ljust(char_len)

        if calib_var.ndim == 3:
          calib_var[profile_idx, doxy_idx, :] = nc.stringtochar(
              np.array(padded_str, dtype=f"S{char_len}")
          )

  if "DOXY_SLOPE" in var_names and not pd.isna(slope_val):
    ds.variables["DOXY_SLOPE"][:] = slope_val
  if "DOXY_DRIFT" in var_names and not pd.isna(drift_val):
    ds.variables["DOXY_DRIFT"][:] = drift_val

  print(f"Updated DOXY Slope/Drift for Cycle {target_cycle}: {calib_str}")


def write_DOXY_from_csv(ds, profile_idx, float_df, target_cycle):
  """Populate DOXY, DOXY_ADJUSTED, and QC variables following CHLA profile target methods."""
  var_names = ds.variables.keys()
  n_prof = ds.dimensions["N_PROF"].size
  n_levels = ds.dimensions["N_LEVELS"].size

  pres_nc_full = ds.variables["PRES"][:]
  pres_nc = (
      pres_nc_full[:, profile_idx] if pres_nc_full.ndim > 1 else pres_nc_full
  )

  cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()

  DOXY_Adjusted_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  DOXY_Adjusted_Array[:] = 99999.0
  DOXY_Adjusted_Array.mask = True

  DOXY_AdjustedQC_Array = np.ma.empty(shape=(n_levels,), dtype="|S1")
  DOXY_AdjustedQC_Array[:] = b"9"
  DOXY_AdjustedQC_Array.mask = False

  DOXY_Array = np.ma.empty(
      shape=(n_levels,), fill_value=99999.0, dtype="float32"
  )
  DOXY_Array[:] = 99999.0
  DOXY_Array.mask = True

  DOXY_QC_Array = np.ma.empty(shape=(n_levels,), dtype="|S1")
  DOXY_QC_Array[:] = b"9"
  DOXY_QC_Array.mask = False

  assigned_nc_pres_vals = set()

  if not cycle_df.empty:
    for row in range(len(cycle_df)):
      row_data = cycle_df.iloc[row]
      csv_pres = np.float32(row_data["PRES"])

      for i in range(n_levels):
        nc_pres = np.float32(pres_nc[i])
        if nc_pres in assigned_nc_pres_vals or np.isnan(nc_pres):
          continue

        if np.isclose(csv_pres, nc_pres, atol=0.05):
          raw_qc = row_data.get("DOXY_FINAL_QC")
          qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"

          raw_doxy = row_data.get("DOXY")
          if pd.notna(raw_doxy) and raw_doxy != 99999.0:
            DOXY_Array[i] = np.float32(raw_doxy)
            DOXY_Array.mask[i] = False
            DOXY_QC_Array[i] = b"1"

          raw_doxy_final = row_data.get("DOXY_FINAL")
          if (
              pd.isna(raw_doxy_final)
              or raw_doxy_final == 99999.0
              or qc_str in ["4", "9"]
          ):
            DOXY_Adjusted_Array[i] = 99999.0
            DOXY_Adjusted_Array.mask[i] = True
          else:
            DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
            DOXY_Adjusted_Array.mask[i] = False

          DOXY_AdjustedQC_Array[i] = qc_str.encode("utf-8")
          assigned_nc_pres_vals.add(nc_pres)
          break

  # Write target profile
  if "DOXY" in var_names:
    if ds.variables["DOXY"].ndim > 1:
      ds.variables["DOXY"][:, profile_idx] = DOXY_Array
      ds.variables["DOXY_QC"][:, profile_idx] = DOXY_QC_Array
    else:
      ds.variables["DOXY"][:] = DOXY_Array
      ds.variables["DOXY_QC"][:] = DOXY_QC_Array

  if "DOXY_ADJUSTED" in var_names:
    if ds.variables["DOXY_ADJUSTED"].ndim > 1:
      ds.variables["DOXY_ADJUSTED"][:, profile_idx] = DOXY_Adjusted_Array
      ds.variables["DOXY_ADJUSTED_QC"][:, profile_idx] = DOXY_AdjustedQC_Array
    else:
      ds.variables["DOXY_ADJUSTED"][:] = DOXY_Adjusted_Array
      ds.variables["DOXY_ADJUSTED_QC"][:] = DOXY_AdjustedQC_Array

  profile_doxy_qc = get_profile_qc_grade(DOXY_AdjustedQC_Array)
  if "PROFILE_DOXY_QC" in ds.variables:
    ds.variables["PROFILE_DOXY_QC"][profile_idx] = profile_doxy_qc

  # DIAGNOSTIC PRINT STATEMENT FOR SECONDARY PROFILE HANDLING
  print("\n" + "=" * 60)
  print(
      f"[VERIFICATION] write_DOXY_from_csv matched CHLA secondary profile method for Cycle {target_cycle}:"
  )
  for iprof in range(n_prof):
    adj_qc_raw = (
        ds.variables["DOXY_ADJUSTED_QC"][:, iprof]
        if "DOXY_ADJUSTED_QC" in var_names
        else []
    )
    blank_count = sum(
        1
        for q in adj_qc_raw
        if (q.decode("utf-8") if isinstance(q, bytes) else str(q)) == " "
    )
    nine_count = sum(
        1
        for q in adj_qc_raw
        if (q.decode("utf-8") if isinstance(q, bytes) else str(q)) == "9"
    )
    val_count = sum(
        1
        for q in adj_qc_raw
        if (q.decode("utf-8") if isinstance(q, bytes) else str(q))
        not in [" ", "9"]
    )
    target_flag = " (Target Profile Modified)" if iprof == profile_idx else " (Secondary Profile Untouched)"
    print(
        f"  Profile [{iprof}]{target_flag}: Blank (' ') = {blank_count},"
        f" Fill ('9') = {nine_count}, Valid QC = {val_count}"
    )
  print("=" * 60 + "\n")

  return DOXY_Adjusted_Array.filled(99999.0)


def write_DOXY_adjusted_error(
    ds, profile_idx, err_mbar, psal, temp, pres, dens, doxy_adj
):
  """Calculate and assign DOXY_ADJUSTED_ERROR in µmol/kg for target Delayed-Mode profile."""
  valid_idx = ~np.isnan(psal) & ~np.isnan(doxy_adj) & (doxy_adj != 99999.0)
  doxy_adj_error = np.full(len(psal), 99999.0, dtype="float32")

  err_umol_L = err_mbar * 1.00
  doxy_adj_error[valid_idx] = err_umol_L
  doxy_adj_error_umol_kg = (doxy_adj_error * 1000.0) / dens
  doxy_adj_error_ma = np.ma.masked_values(doxy_adj_error_umol_kg, 99999.0)

  doxy_adj_full = ds.variables["DOXY_ADJUSTED"][:]
  if doxy_adj_full.ndim > 1:
    ds.variables["DOXY_ADJUSTED_ERROR"][:, profile_idx] = doxy_adj_error_ma
  else:
    ds.variables["DOXY_ADJUSTED_ERROR"][:] = doxy_adj_error_ma


def safe_rename(from_file, to_file):
  """Safely rename file handling Windows/Linux file locks."""
  gc.collect()
  try:
    os.replace(from_file, to_file)
    print(
        f"Renamed {os.path.basename(from_file)} to {os.path.basename(to_file)}"
    )
  except OSError:
    try:
      shutil.copyfile(from_file, to_file)
      os.remove(from_file)
      print(
          "Copied and replaced (lock fallback):"
          f" {os.path.basename(from_file)} -> {os.path.basename(to_file)}"
      )
    except Exception as e:
      print(
          f"Warning: Failed to rename or copy {from_file} to {to_file}: {e}"
      )


# ==============================================================================
# SECTION 4: Processing Loop Across All Float IDs
# ==============================================================================

total_floats = len(WMO_FLOAT_IDS)

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
  print("\n" + "=" * 70)
  print(f" PROCESSING WMO FLOAT ID: {WMOfloatid} ({idx_f} of {total_floats})")
  print("=" * 70)

  float_dir = f"/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
  bio_dmqc_csv_path = (
      f"/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_DOXY_{WMOfloatid}.csv"
  )
  output_lut_dir = os.path.join(float_dir, "LUT")

  if not os.path.exists(bio_dmqc_csv_path):
    print(
        f"ERROR: CSV file not found at {bio_dmqc_csv_path}. Skipping Float"
        f" {WMOfloatid}..."
    )
    continue

  if not os.path.exists(float_dir):
    print(
        f"ERROR: Float directory not found at {float_dir}. Skipping Float"
        f" {WMOfloatid}..."
    )
    continue

  os.makedirs(output_lut_dir, exist_ok=True)

  # --------------------------------------------------------------------------
  # STEP 1: CHLA BD Filler
  # --------------------------------------------------------------------------
  print(f"\n--- [Float {WMOfloatid}] Step 1: Running CHLA BD Filler ---")

  all_bd_files = sorted(
      glob.glob(os.path.join(float_dir, f"BD*{WMOfloatid}_*.nc"))
  )
  all_br_files = sorted(
      glob.glob(os.path.join(float_dir, f"AOML_BR*{WMOfloatid}_*.nc"))
  )
  if not all_br_files:
    all_br_files = sorted(
        glob.glob(os.path.join(float_dir, f"BR*{WMOfloatid}_*.nc"))
    )

  sorted_b_files = organize_b_files(all_bd_files, all_br_files)
  print(f"{len(sorted_b_files)} relevant B files found in {float_dir}")

  new_bd_files = []

  for bgc_filename in sorted_b_files:
    print(f"Processing CHLA for {os.path.basename(bgc_filename)}")
    idx_profile = int(bgc_filename[-6:-3])
    w_bgc_filename = create_working_bd_file(bgc_filename)

    try:
      bgc_file = netCDF4.Dataset(w_bgc_filename, "a")

      write_history_chla(bgc_file, iprof_chla)
      write_parameter_data_mode_chla(bgc_file, iprof_chla)
      write_scientific_calib_chla(bgc_file, idx_profile, bio_dmqc_csv_path)

      for i in range(bgc_file.dimensions["N_PROF"].size):
        if i != iprof_chla:
          continue
        bgc_file.variables["DATA_STATE_INDICATOR"][i] = np.ma.array(
            data_state_indicator, mask=[False, False, True, True], dtype="|S1"
        )

      write_chla_BBP_adjusted(
          bgc_file, idx_profile, bio_dmqc_csv_path, iprof_chla
      )

      bgc_file.close()
      new_bd_files.append(w_bgc_filename)
    except Exception as e:
      print(f"Error processing CHLA for {bgc_filename}: {e}")

  # Move output files into the LUT folder for Step 2
  for file in new_bd_files:
    path, name = os.path.split(file)
    new_name = name.replace("w_AOML_BD", "BD").replace("w_BD", "BD")
    new_path = os.path.join(output_lut_dir, new_name)
    safe_rename(file, new_path)

  # --------------------------------------------------------------------------
  # STEP 2: DOXY BD Filler (Using output from Step 1)
  # --------------------------------------------------------------------------
  print(f"\n--- [Float {WMOfloatid}] Step 2: Running DOXY BD Filler ---")

  csv_data = pd.read_csv(bio_dmqc_csv_path)

  req_cols = [
      "FLOAT_NUM",
      "CYCLE_NUMBER",
      "PRES",
      "DOXY",
      "DOXY_FINAL",
      "DOXY_FINAL_QC",
      "DOXY_SLOPE",
      "DOXY_DRIFT",
  ]
  if not all(col in csv_data.columns for col in req_cols):
    print(
        f"Warning: CSV file for float {WMOfloatid} missing required DOXY"
        " columns. Skipping Step 2..."
    )
    continue

  csv_data["CYCLE_NUMBER"] = csv_data["CYCLE_NUMBER"].astype(int)
  unique_floats = csv_data["FLOAT_NUM"].unique()

  for float_id_raw in unique_floats:
    floatid = int(float_id_raw)
    float_df = csv_data[csv_data["FLOAT_NUM"] == floatid]

    try:
      inst_float = FLOAT_TYPES.get(floatid, "aoml_apex")

      if inst_float == "aoml_apex":
        comment_dmqc_operator = (
            "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer"
            " McWhorter, NOAA/AOML"
        )
        history_institution = "AO"
        history_reference = "WOA2023"
        DOXY_adj_err = 2
        iprof_doxy = 0
      elif inst_float == "aoml_navis":
        comment_dmqc_operator = (
            "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer"
            " McWhorter, NOAA/AOML"
        )
        history_institution = "AO"
        history_reference = "WOA2023"
        DOXY_adj_err = 5
        iprof_doxy = 0

      today_str = dt.now().strftime("%Y-%m-%d")
      final_doxy_out_dir = os.path.join(output_lut_dir, today_str)
      os.makedirs(final_doxy_out_dir, exist_ok=True)

      csv_cycles = sorted(float_df["CYCLE_NUMBER"].unique())

      for target_cycle in csv_cycles:
        pattern_bd = f"BD*{floatid}_{target_cycle:03d}.nc"
        pattern_bd_raw = f"BD*{floatid}_{target_cycle}.nc"

        matched_files = glob.glob(
            os.path.join(output_lut_dir, pattern_bd)
        ) or glob.glob(os.path.join(output_lut_dir, pattern_bd_raw))

        if not matched_files:
          print(
              f"Warning: Cycle {target_cycle} present in CSV, but no matching"
              f" BD file found in {output_lut_dir}"
          )
          continue

        bgc_filename = matched_files[0]
        print(
            f"Processing DOXY for Cycle {target_cycle}"
            f" ({inst_float.upper()}): {os.path.basename(bgc_filename)}"
        )

        w_bgc_filename = create_working_doxy_bd_file(
            bgc_filename, final_doxy_out_dir
        )

        try:
          ds = nc.Dataset(w_bgc_filename, "r+")

          # 1. History Metadata
          write_history_doxy(
              ds,
              iprof_doxy,
              history_institution,
              history_reference,
              comment_dmqc_operator,
          )

          # 2. Synchronize Parameter Data Mode and Station Parameters
          write_parameter_data_mode_doxy(ds, iprof_idx=iprof_doxy)

          # 3. Write Slope & Drift
          write_DOXY_slope_drift(ds, iprof_doxy, float_df, target_cycle)

          # 4. Write DOXY values from CSV & Update PROFILE_DOXY_QC
          doxy_adjusted = write_DOXY_from_csv(
              ds, iprof_doxy, float_df, target_cycle
          )

          # 5. Get physical profile density
          phys_filename = get_phys_filename(bgc_filename, float_dir)
          pres_phys_raw = get_phys_raw_pres(phys_filename)
          phys_data = get_dens(phys_filename)

          pres_bgc = ds.variables["PRES"][:]
          iprof_phys = get_iprof_phys(pres_phys_raw, pres_bgc, iprof_doxy)

          # 6. DOXY Error Calculation
          psal_col = (
              phys_data["psal"][:, iprof_phys]
              if phys_data["psal"].ndim > 1
              else phys_data["psal"]
          )
          temp_col = (
              phys_data["temp"][:, iprof_phys]
              if phys_data["temp"].ndim > 1
              else phys_data["temp"]
          )
          pres_col = (
              pres_bgc[:, iprof_doxy] if pres_bgc.ndim > 1 else pres_bgc
          )
          dens_col = (
              phys_data["dens"][:, iprof_phys]
              if phys_data["dens"].ndim > 1
              else phys_data["dens"]
          )

          write_DOXY_adjusted_error(
              ds,
              iprof_doxy,
              DOXY_adj_err,
              psal_col,
              temp_col,
              pres_col,
              dens_col,
              doxy_adjusted,
          )

          ds.close()

        except Exception as e:
          print(f"Error updating DOXY in file {w_bgc_filename}: {e}")
          if "ds" in locals() and ds.isopen():
            ds.close()

        # Final rename in output directory
        base_name = os.path.basename(w_bgc_filename)
        new_name = re.sub(r"^w_BD", "BD", base_name)
        new_path = os.path.join(final_doxy_out_dir, new_name)

        safe_rename(w_bgc_filename, new_path)

    except Exception as e:
      print(f"Error processing DOXY for Float ID {floatid}: {e}")

print("\nProcessing complete for all WMO Float IDs.")