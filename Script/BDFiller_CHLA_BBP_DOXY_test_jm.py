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
import traceback
from datetime import datetime as dt, timezone

import netCDF4
import netCDF4 as nc
import numpy as np
import pandas as pd


# Dynamic mapping for Float IDs and Float Types
FLOAT_TYPES = {
    4903622: "aoml_apex",
    # 2904010: "aoml_apex",
    # 2904011: "aoml_apex",
    # 4903624: "aoml_apex",
    # 4903625: "aoml_apex",
    # 4903904: "aoml_navis",
    # 6999992: "aoml_navis",
    # 7902327: "aoml_navis",
    # 3902693: "aoml_apex",
    # 1902800: "aoml_apex",
    # 7901009: "aoml_navis",
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
BBP700_Adjusted_ERROR_est = 0.0005  # Default estimated error for BBP700

data_state_indicator = ["2", "C", "", ""]
parameter_data_mode = "D"

# Standard Argo valid_min and valid_max ranges for physical & coordinate parameters
ARGO_VALID_RANGES = {
    "LATITUDE": (-90.0, 90.0),
    "LONGITUDE": (-180.0, 180.0),
    "PRES": (0.0, 12000.0),
    "TEMP_DOXY": (-2.5, 40.0),
    "CHLA_FLUORESCENCE": (0.0, 50.0),
    "CHLA_FLUORESCENCE_ADJUSTED": (0.0, 50.0),
    "RPHASE_DOXY": (0.0, 500.0),
    "HUMIDITY_NITRATE": (0.0, 100.0),
    "TPHASE_DOXY": (0.0, 500.0),
    "PHASE_DELAY_DOXY": (0.0, 500.0),
    "FREQUENCY_DOXY": (0.0, 100000.0),
    "DOXY": (0.0, 600.0),
    "DOXY_ADJUSTED": (0.0, 600.0),
}


# ==============================================================================
# SECTION 2: Sanitation & Correction Helper Functions
# ==============================================================================


def clean_and_fill_qc_variables(bgc_file):
    """Safely inspect and clean ALL QC variables in the NetCDF file.
    Replaces blank or null QC flags with '9' (missing data). Strictly avoids modifying PH variables.
    """
    n_prof = bgc_file.dimensions["N_PROF"].size

    for var_name, var in bgc_file.variables.items():
        if not var_name.endswith("_QC"):
            continue

        # Preserve PH portion completely
        if "PH" in var_name:
            continue

        for iprof in range(n_prof):
            try:
                if var.dtype.kind in ["S", "U", "O"]:
                    raw_var_slice = var[iprof] if var.ndim == 1 else var[iprof, :]
                    slice_data = np.array(raw_var_slice, copy=True)

                    if hasattr(slice_data, "filled"):
                        slice_data = slice_data.filled(b" ")

                    char_arr = slice_data.astype("|S1")
                    blank_mask = (
                        (char_arr == b" ")
                        | (char_arr == b"")
                        | (char_arr == b"\x00")
                    )

                    blank_count = int(np.sum(blank_mask))

                    if blank_count > 0:
                        fill_code = b"F" if var_name.startswith("PROFILE_") else b"9"
                        char_arr[blank_mask] = fill_code

                        if var.ndim == 1:
                            bgc_file.variables[var_name][iprof] = char_arr[0]
                        else:
                            bgc_file.variables[var_name][iprof, :] = char_arr

            except Exception as e:
                print(f"Error cleaning QC var {var_name}: {e}")


def add_missing_valid_range_attributes(bgc_file):
    """Ensure required Argo valid_min and valid_max attributes exist for essential variables."""
    for var_name, (vmin, vmax) in ARGO_VALID_RANGES.items():
        if var_name in bgc_file.variables:
            var = bgc_file.variables[var_name]

            if "valid_min" not in var.ncattrs():
                var.setncattr("valid_min", np.float32(vmin))

            if "valid_max" not in var.ncattrs():
                var.setncattr("valid_max", np.float32(vmax))


def remove_forbidden_attributes(bgc_file):
    """Remove valid_min and valid_max attributes from variables where forbidden by Argo NetCDF standards.
    Retains attributes for coordinate and primary measurement parameters listed in ARGO_VALID_RANGES.
    """
    forbidden_attrs = ["valid_min", "valid_max"]
    for var_name, var in bgc_file.variables.items():
        if var_name not in ARGO_VALID_RANGES:
            for attr in forbidden_attrs:
                if attr in var.ncattrs():
                    var.delncattr(attr)


# ==============================================================================
# SECTION 3: Helper Functions - CHLA & BBP700 Processing
# ==============================================================================


def detect_parameter_profile(ds, param_prefix):
    """Dynamically locate the profile index (N_PROF) containing the target parameter prefix."""
    if "STATION_PARAMETERS" not in ds.variables:
        return 0

    n_prof = ds.dimensions["N_PROF"].size
    n_param = ds.dimensions["N_PARAM"].size

    for iprof in range(n_prof):
        param_mat = ds.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        for j in range(n_param):
            param_str = "".join([
                c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
                for c in param_mat[j]
            ]).strip()
            if param_str.startswith(param_prefix):
                return iprof
    return 0


def create_working_bd_file(filename, dest_dir=None):
    """Create a working copy of the B file with a leading 'w_' in the name."""
    path, name = os.path.split(filename)
    bd_name = name.replace("BR", "BD")
    out_dir = dest_dir if dest_dir else path
    w_filename = os.path.join(out_dir, f"w_{bd_name}")
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
    bgc_file.history = datetime.datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ creation"
    )
    bgc_file.setncattr("comment_dmqc_operator", comment_dmqc_operator_chla)

    history_step = "ARSQ"
    history_action = "IP"
    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

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
    """Set PARAMETER_DATA_MODE and DATA_MODE strictly for CHLA and BBP700 on target profile iprof_idx."""
    n_param = bgc_file.dimensions["N_PARAM"].size

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof_idx]
    if isinstance(param_mat, np.ma.MaskedArray):
        param_mat = param_mat.filled(b" ")

    # Only set target profile data_mode to 'D' if it contains CHLA/BBP700
    for j in range(n_param):
        param_str = "".join([
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in param_mat[j]
        ]).strip()
        if param_str == "PRES":
            pdm[iprof_idx, j] = "R"  # Force PRES to 'R'
        elif param_str.startswith("CHLA") or param_str.startswith("BBP"):
            pdm[iprof_idx, j] = "D"
            data_mode[iprof_idx] = "D"

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
    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

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
    bgc_file, idx_profile, bio_dmqc_csv_path, iprof_idx=0
):
    """Populate CHLA_ADJUSTED and CHLA_FLUORESCENCE_ADJUSTED on target profile iprof_idx."""
    df_bio = pd.read_csv(bio_dmqc_csv_path)
    df_bio = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    n_levels = bgc_file.dimensions["N_LEVELS"].size

    CHLA_Adjusted_Array = np.ma.empty(
        shape=(n_levels,), fill_value=99999.0, dtype="float32"
    )
    CHLA_Adjusted_Array[:] = 99999.0
    CHLA_Adjusted_Array.mask = True

    CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

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

    CHLA_FLUORESCENCE_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

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
                raw_qc = row_data["CHLA_FINAL_QC"] if "CHLA_FINAL_QC" in row_data else "9"
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

                assigned_nc_pres_vals.add(nc_pres)
                break

    bgc_file.variables["CHLA_ADJUSTED"][iprof_idx] = CHLA_Adjusted_Array
    bgc_file.variables["CHLA_ADJUSTED_QC"][iprof_idx, :] = CHLA_AdjustedQC_Array
    bgc_file.variables["CHLA_ADJUSTED_ERROR"][iprof_idx] = (
        CHLA_Adjusted_ERROR_Array
    )

    if "CHLA_FLUORESCENCE_ADJUSTED" in bgc_file.variables:
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED"][iprof_idx] = (
            CHLA_FLUORESCENCE_Adjusted_Array
        )
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_QC"][iprof_idx, :] = (
            CHLA_FLUORESCENCE_AdjustedQC_Array
        )
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_ERROR"][iprof_idx] = (
            CHLA_FLUORESCENCE_Adjusted_ERROR_Array
        )

    # Compute profile QC grade specifically for CHLA on target profile iprof_idx
    prof_qc = get_profile_qc_grade(CHLA_AdjustedQC_Array)
    prof_fluo_qc = get_profile_qc_grade(CHLA_FLUORESCENCE_AdjustedQC_Array)

    if "PROFILE_CHLA_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_CHLA_QC"][iprof_idx] = np.array([prof_qc], dtype="|S1")
    if "PROFILE_CHLA_FLUORESCENCE_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_CHLA_FLUORESCENCE_QC"][iprof_idx] = np.array([prof_fluo_qc], dtype="|S1")


def write_BBP700_adjusted(
    bgc_file, idx_profile, bio_dmqc_csv_path, iprof_idx=0
):
    """Populate BBP700_ADJUSTED, BBP700_ADJUSTED_QC, BBP700_ADJUSTED_ERROR, and PROFILE_BBP700_QC."""
    if "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    df_bio = pd.read_csv(bio_dmqc_csv_path)
    df_bio = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    if "BBP700_FINAL" not in df_bio.columns or "BBP700_FINAL_QC" not in df_bio.columns:
        print(f"Warning: BBP700_FINAL/QC columns missing in CSV for cycle {idx_profile}.")
        return

    n_levels = bgc_file.dimensions["N_LEVELS"].size

    BBP700_Adjusted_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    BBP700_Adjusted_Array[:] = 99999.0
    BBP700_Adjusted_Array.mask = True

    BBP700_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

    BBP700_Adjusted_ERROR_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    BBP700_Adjusted_ERROR_Array[:] = 99999.0
    BBP700_Adjusted_ERROR_Array.mask = True

    assigned_nc_pres_vals = set()

    for row in range(len(df_bio)):
        row_data = df_bio.iloc[row]
        csv_pres = np.float32(row_data["PRES"])

        for i in range(n_levels):
            nc_pres = np.float32(bgc_file.variables["PRES"][iprof_idx, i])

            if nc_pres in assigned_nc_pres_vals:
                continue

            if np.isclose(csv_pres, nc_pres, atol=0.05):
                raw_qc = row_data["BBP700_FINAL_QC"]
                qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"

                raw_bbp = row_data["BBP700_FINAL"]

                if pd.isna(raw_bbp) or raw_bbp == 99999.0 or qc_str in ["4", "9"]:
                    BBP700_Adjusted_Array[i] = 99999.0
                    BBP700_Adjusted_ERROR_Array[i] = 99999.0
                    BBP700_Adjusted_Array.mask[i] = True
                    BBP700_Adjusted_ERROR_Array.mask[i] = True
                else:
                    BBP700_Adjusted_Array[i] = np.float32(raw_bbp)
                    BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                    BBP700_Adjusted_Array.mask[i] = False
                    BBP700_Adjusted_ERROR_Array.mask[i] = False

                BBP700_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                assigned_nc_pres_vals.add(nc_pres)
                break

    bgc_file.variables["BBP700_ADJUSTED"][iprof_idx] = BBP700_Adjusted_Array
    if "BBP700_ADJUSTED_QC" in bgc_file.variables:
        bgc_file.variables["BBP700_ADJUSTED_QC"][iprof_idx, :] = BBP700_AdjustedQC_Array
    if "BBP700_ADJUSTED_ERROR" in bgc_file.variables:
        bgc_file.variables["BBP700_ADJUSTED_ERROR"][iprof_idx] = BBP700_Adjusted_ERROR_Array

    # Compute and set profile QC grade specifically for BBP700
    prof_bbp_qc = get_profile_qc_grade(BBP700_AdjustedQC_Array)
    if "PROFILE_BBP700_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_BBP700_QC"][iprof_idx] = np.array([prof_bbp_qc], dtype="|S1")


# ==============================================================================
# SECTION 4: Helper Functions - DOXY Processing
# ==============================================================================


def create_working_doxy_bd_file(filename, dest_dir):
    """Copy file to destination directory as 'w_BD...' working file."""
    base_name = os.path.basename(filename)
    bd_name = re.sub(r"^BR", "BD", base_name)
    w_filename = os.path.join(dest_dir, f"w_{bd_name}")
    shutil.copyfile(filename, w_filename)
    return w_filename


def update_history_doxy(nc_ds, dct, iprof_idx):
    """Update HISTORY array entries natively using netCDF4."""
    hix = nc_ds.dimensions["N_HISTORY"].size
    for name, value in dct.items():
        if name in nc_ds.variables:
            char_len = nc_ds.dimensions[nc_ds[name].dimensions[-1]].size
            padded_val = str(value).ljust(char_len)[:char_len]
            nc_ds[name][hix, iprof_idx, :] = nc.stringtochar(
                np.array(padded_val, dtype=f"S{char_len}")
            )


def write_history_doxy(ds, profile_idx, inst, ref, comment_op):
    """Update global history attributes for DOXY processing."""
    print(f"[DEBUG] Substep 1: Writing HISTORY metadata for DOXY...")
    ds.history = datetime.datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ creation"
    )
    ds.setncattr("comment_dmqc_operator", comment_op)

    history_step = "ARSQ"
    history_action = "IP"
    history_software = "BITTIG"
    history_software_release = "2024"
    history_parameter = "DOXY"
    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

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
    """Set PARAMETER_DATA_MODE and DATA_MODE strictly for DOXY on target profile iprof_idx.
    Leaves non-target parameters (PH_IN_SITU_TOTAL, NITRATE, etc.) untouched.
    """
    n_param = bgc_file.dimensions["N_PARAM"].size

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof_idx]
    if isinstance(param_mat, np.ma.MaskedArray):
        param_mat = param_mat.filled(b" ")

    for j in range(n_param):
        param_str = "".join([
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in param_mat[j]
        ]).strip()
        if param_str == "PRES":
            pdm[iprof_idx, j] = "R"  # Force PRES to 'R'
        elif param_str.startswith("DOXY"):
            pdm[iprof_idx, j] = "D"
            data_mode[iprof_idx] = "D"

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle):
    """Write DOXY slope and drift calibration coefficients safely handling missing data and array dimensions."""
    print(f"[DEBUG] Substep 3: Writing DOXY slope & drift coefficients...")
    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)]
    if cycle_df.empty:
        print(f"  [Warning] No cycle data found in CSV for cycle {target_cycle}")
        return

    slope_val = cycle_df["DOXY_SLOPE"].iloc[0] if "DOXY_SLOPE" in cycle_df.columns else np.nan
    drift_val = cycle_df["DOXY_DRIFT"].iloc[0] if "DOXY_DRIFT" in cycle_df.columns else np.nan

    s_str = "1" if pd.isna(slope_val) else str(slope_val)
    d_str = "0" if pd.isna(drift_val) else str(drift_val)
    calib_str = f"m={s_str}, d={d_str}"

    if "SCIENTIFIC_CALIB_COEFFICIENT" in ds.variables:
        calib_var = ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"]
        n_prof = ds.dimensions["N_PROF"].size
        n_param = ds.dimensions["N_PARAM"].size
        char_len = calib_var.shape[-1]

        padded_str = calib_str.ljust(char_len)[:char_len]
        char_arr = np.array(list(padded_str), dtype="|S1")

        for iprof in range(n_prof):
            station_params = ds.variables["STATION_PARAMETERS"][iprof]
            if isinstance(station_params, np.ma.MaskedArray):
                station_params = station_params.filled(b" ")

            for j in range(n_param):
                param_str = "".join([
                    c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
                    for c in station_params[j]
                ]).strip()

                if param_str.startswith("DOXY"):
                    if calib_var.ndim == 4:
                        calib_var[iprof, 0, j, :] = char_arr
                    elif calib_var.ndim == 3:
                        calib_var[iprof, j, :] = char_arr

    if "DOXY_SLOPE" in ds.variables and not pd.isna(slope_val):
        ds.variables["DOXY_SLOPE"][:] = slope_val
    if "DOXY_DRIFT" in ds.variables and not pd.isna(drift_val):
        ds.variables["DOXY_DRIFT"][:] = drift_val

    print(f"[DEBUG] Updated DOXY Slope/Drift for Cycle {target_cycle}: {calib_str}")


def write_DOXY_from_csv(ds, profile_idx, float_df, target_cycle):
    """Populate DOXY, DOXY_ADJUSTED, DOXY_ADJUSTED_ERROR, and QC variables from CSV."""
    print(f"[DEBUG] Substep 4: Running write_DOXY_from_csv...")
    var_names = ds.variables.keys()
    n_prof = ds.dimensions["N_PROF"].size
    n_levels = ds.dimensions["N_LEVELS"].size

    pres_nc_full = ds.variables["PRES"][:]

    if pres_nc_full.ndim > 1:
        if pres_nc_full.shape[0] == n_prof:
            pres_nc = pres_nc_full[profile_idx, :]
        else:
            pres_nc = pres_nc_full[:, profile_idx]
    else:
        pres_nc = pres_nc_full

    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()

    # Numeric Arrays (Masked)
    DOXY_Adjusted_Array = np.ma.empty(
        shape=(n_levels,), fill_value=99999.0, dtype="float32"
    )
    DOXY_Adjusted_Array[:] = 99999.0
    DOXY_Adjusted_Array.mask = True

    DOXY_Adjusted_Error_Array = np.ma.empty(
        shape=(n_levels,), fill_value=99999.0, dtype="float32"
    )
    DOXY_Adjusted_Error_Array[:] = 99999.0
    DOXY_Adjusted_Error_Array.mask = True

    DOXY_Array = np.ma.empty(
        shape=(n_levels,), fill_value=99999.0, dtype="float32"
    )
    DOXY_Array[:] = 99999.0
    DOXY_Array.mask = True

    DOXY_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")
    DOXY_QC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

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

                    raw_doxy_error = row_data.get("DOXY_ADJUSTED_ERROR")
                    if (
                        pd.notna(raw_doxy_error)
                        and raw_doxy_error != 99999.0
                        and qc_str not in ["4", "9"]
                    ):
                        DOXY_Adjusted_Error_Array[i] = np.float32(raw_doxy_error)
                        DOXY_Adjusted_Error_Array.mask[i] = False
                    else:
                        DOXY_Adjusted_Error_Array[i] = 99999.0
                        DOXY_Adjusted_Error_Array.mask[i] = True

                    DOXY_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                    assigned_nc_pres_vals.add(nc_pres)
                    break

    if "DOXY" in var_names:
        doxy_var = ds.variables["DOXY"]
        if doxy_var.ndim > 1:
            if doxy_var.shape[0] == n_prof:
                doxy_var[profile_idx, :] = DOXY_Array
                ds.variables["DOXY_QC"][profile_idx, :] = DOXY_QC_Array
            else:
                doxy_var[:, profile_idx] = DOXY_Array
                ds.variables["DOXY_QC"][:, profile_idx] = DOXY_QC_Array
        else:
            doxy_var[:] = DOXY_Array
            ds.variables["DOXY_QC"][:] = DOXY_QC_Array

    if "DOXY_ADJUSTED" in var_names:
        doxy_adj_var = ds.variables["DOXY_ADJUSTED"]
        if doxy_adj_var.ndim > 1:
            if doxy_adj_var.shape[0] == n_prof:
                doxy_adj_var[profile_idx, :] = DOXY_Adjusted_Array
                ds.variables["DOXY_ADJUSTED_QC"][profile_idx, :] = DOXY_AdjustedQC_Array
            else:
                doxy_adj_var[:, profile_idx] = DOXY_Adjusted_Array
                ds.variables["DOXY_ADJUSTED_QC"][:, profile_idx] = DOXY_AdjustedQC_Array
        else:
            doxy_adj_var[:] = DOXY_Adjusted_Array
            ds.variables["DOXY_ADJUSTED_QC"][:] = DOXY_AdjustedQC_Array

    if "DOXY_ADJUSTED_ERROR" in var_names:
        doxy_adj_err_var = ds.variables["DOXY_ADJUSTED_ERROR"]
        if doxy_adj_err_var.ndim > 1:
            if doxy_adj_err_var.shape[0] == n_prof:
                doxy_adj_err_var[profile_idx, :] = DOXY_Adjusted_Error_Array
            else:
                doxy_adj_err_var[:, profile_idx] = DOXY_Adjusted_Error_Array
        else:
            doxy_adj_err_var[:] = DOXY_Adjusted_Error_Array

    # Compute profile QC grade specifically for DOXY on target profile profile_idx
    prof_qc = get_profile_qc_grade(DOXY_AdjustedQC_Array)
    if "PROFILE_DOXY_QC" in ds.variables:
        prof_qc_var = ds.variables["PROFILE_DOXY_QC"]
        qc_char = np.array([prof_qc], dtype="|S1")
        if prof_qc_var.ndim == 1:
            prof_qc_var[profile_idx] = qc_char
        elif prof_qc_var.ndim == 2:
            prof_qc_var[profile_idx, 0] = qc_char

    print(f"[DEBUG] Substep 4 complete for write_DOXY_from_csv.")
    return (
        DOXY_Adjusted_Array.filled(99999.0)
        if hasattr(DOXY_Adjusted_Array, "filled")
        else DOXY_Adjusted_Array
    )


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
# SECTION 5: Processing Loop Across All Float IDs
# ==============================================================================

total_floats = len(WMO_FLOAT_IDS)

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
    print("\n" + "=" * 70)
    print(f" PROCESSING WMO FLOAT ID: {WMOfloatid} ({idx_f} of {total_floats})")
    print("=" * 70)

    main_float_dir = f"/data/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
    bio_dmqc_csv_path = (
        f"/data/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_BBP_DOXY_{WMOfloatid}.csv"
    )

    output_lut_dir = os.path.join(main_float_dir, "LUT")
    today_str = dt.now().strftime("%Y-%m-%d")
    final_doxy_out_dir = os.path.join(output_lut_dir, today_str)

    if not os.path.exists(bio_dmqc_csv_path):
        print(
            f"ERROR: CSV file not found at {bio_dmqc_csv_path}. Skipping Float"
            f" {WMOfloatid}..."
        )
        continue

    if not os.path.exists(main_float_dir):
        print(
            f"ERROR: Main float directory not found at {main_float_dir}. Skipping Float"
            f" {WMOfloatid}..."
        )
        continue

    os.makedirs(output_lut_dir, exist_ok=True)
    os.makedirs(final_doxy_out_dir, exist_ok=True)

    # --------------------------------------------------------------------------
    # STEP 1: CHLA & BBP700 BD Filler (Reads from main directory -> Outputs into LUT)
    # --------------------------------------------------------------------------
    print(f"\n--- [Float {WMOfloatid}] Step 1: Running CHLA & BBP700 BD Filler ---")

    all_bd_files = sorted(
        glob.glob(os.path.join(main_float_dir, f"BD*{WMOfloatid}_*.nc"))
    )
    all_br_files = sorted(
        glob.glob(os.path.join(main_float_dir, f"AOML_BR*{WMOfloatid}_*.nc"))
    )
    if not all_br_files:
        all_br_files = sorted(
            glob.glob(os.path.join(main_float_dir, f"BR*{WMOfloatid}_*.nc"))
        )

    sorted_b_files = organize_b_files(all_bd_files, all_br_files)
    print(f"{len(sorted_b_files)} relevant B files found in main directory: {main_float_dir}")

    new_bd_files = []

    for bgc_filename in sorted_b_files:
        print(f"Processing CHLA and BBP700 for {os.path.basename(bgc_filename)}")
        idx_profile = int(bgc_filename[-6:-3])

        w_bgc_filename = create_working_bd_file(bgc_filename, dest_dir=output_lut_dir)

        try:
            bgc_file = netCDF4.Dataset(w_bgc_filename, "a")

            iprof_chla = detect_parameter_profile(bgc_file, "CHLA")
            iprof_bbp = detect_parameter_profile(bgc_file, "BBP700")

            write_history_chla(bgc_file, iprof_chla)
            write_parameter_data_mode_chla(bgc_file, iprof_chla)
            write_scientific_calib_chla(bgc_file, idx_profile, bio_dmqc_csv_path)

            for i in range(bgc_file.dimensions["N_PROF"].size):
                if i not in [iprof_chla, iprof_bbp]:
                    continue
                bgc_file.variables["DATA_STATE_INDICATOR"][i] = np.ma.array(
                    data_state_indicator, mask=[False, False, True, True], dtype="|S1"
                )

            # Write adjusted CHLA and BBP700 data
            write_chla_BBP_adjusted(
                bgc_file, idx_profile, bio_dmqc_csv_path, iprof_chla
            )
            write_BBP700_adjusted(
                bgc_file, idx_profile, bio_dmqc_csv_path, iprof_bbp
            )

            # Clean ALL QC variables (skipping PH) and fix missing/forbidden attributes
            clean_and_fill_qc_variables(bgc_file)
            add_missing_valid_range_attributes(bgc_file)
            remove_forbidden_attributes(bgc_file)

            bgc_file.close()
            new_bd_files.append(w_bgc_filename)
        except Exception as e:
            print(f"Error processing CHLA/BBP700 for {bgc_filename}: {e}")

    for file in new_bd_files:
        path, name = os.path.split(file)
        new_name = re.sub(r"^w_AOML_BD", "BD", name)
        new_name = re.sub(r"^w_BD", "BD", new_name)
        new_path = os.path.join(output_lut_dir, new_name)
        safe_rename(file, new_path)

    # --------------------------------------------------------------------------
    # STEP 2: DOXY BD Filler (Reads from LUT directory -> Outputs into LUT/YYYY-MM-DD)
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
        "DOXY_ADJUSTED_ERROR",
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

            comment_dmqc_operator = (
                "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer"
                " McWhorter, NOAA/AOML"
            )
            history_institution = "AO"
            history_reference = "WOA2023"

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
                        f" BD file found in LUT directory: {output_lut_dir}"
                    )
                    continue

                bgc_filename = matched_files[0]
                print(
                    f"\nProcessing DOXY for Cycle {target_cycle}"
                    f" ({inst_float.upper()}): {os.path.basename(bgc_filename)}"
                )

                w_bgc_filename = create_working_doxy_bd_file(
                    bgc_filename, final_doxy_out_dir
                )

                try:
                    ds = nc.Dataset(w_bgc_filename, "r+")

                    iprof_doxy = detect_parameter_profile(ds, "DOXY")

                    write_history_doxy(
                        ds,
                        iprof_doxy,
                        history_institution,
                        history_reference,
                        comment_dmqc_operator,
                    )

                    write_parameter_data_mode_doxy(ds, iprof_idx=iprof_doxy)
                    write_DOXY_slope_drift(ds, iprof_doxy, float_df, target_cycle)

                    doxy_adjusted = write_DOXY_from_csv(
                        ds, iprof_doxy, float_df, target_cycle
                    )

                    # Clean ALL QC variables (skipping PH) and fix missing/forbidden attributes
                    clean_and_fill_qc_variables(ds)
                    add_missing_valid_range_attributes(ds)
                    remove_forbidden_attributes(ds)

                    ds.close()

                except Exception as e:
                    print(f"Error updating DOXY in file {w_bgc_filename}: {e}")
                    traceback.print_exc()
                    if "ds" in locals() and ds.isopen():
                        ds.close()

                base_name = os.path.basename(w_bgc_filename)
                new_name = re.sub(r"^w_BD", "BD", base_name)
                new_path = os.path.join(final_doxy_out_dir, new_name)

                safe_rename(w_bgc_filename, new_path)

        except Exception as e:
            print(f"Error processing DOXY for Float ID {floatid}: {e}")

print("\nProcessing complete for all WMO Float IDs.")