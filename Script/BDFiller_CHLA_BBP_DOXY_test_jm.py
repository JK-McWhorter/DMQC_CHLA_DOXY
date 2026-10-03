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
    2904010: "aoml_apex",
    2904011: "aoml_apex",
    4903624: "aoml_apex",
    4903625: "aoml_apex",
    4903904: "aoml_navis",
    6999992: "aoml_navis",
    7902327: "aoml_navis",
    3902693: "aoml_apex",
    1902800: "aoml_apex",
    7901009: "aoml_navis",
}

# Derived list of WMO Float IDs to process
WMO_FLOAT_IDS = list(FLOAT_TYPES.keys())

# Common Configuration Defaults
comment_dmqc_operator_chla = (
    "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer"
    " McWhorter, NOAA/AOML"
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

# Argo Cookbook DOXY standard strings
scientific_calibration_equation_DOXY = "DOXY_ADJUSTED = DOXY * slope + drift"
scientific_calibration_comment_DOXY = (
    "DOXY calibration computed against WOA climatology following Bittig et al. (2018) methodology"
)

# Argo Cookbook Broken Sensor strings
scientific_calibration_comment_BROKEN = "Sensor failure / broken sensor - data flagged bad or missing (Argo BGC QC Manual)"
scientific_calibration_equation_BROKEN = "Not applicable"
scientific_calibration_coefficient_BROKEN = "Not applicable"

data_state_indicator = ["2", "C", "", ""]

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


def apply_float_override_conditions(df, WMOfloatid):
    """Applies Float ID and Cycle-Specific overrides mirroring R dplyr mutate block for broken sensors."""
    if "CYCLE_NUMBER" not in df.columns:
        return df

    float_str = str(WMOfloatid)
    cycle_num = df["CYCLE_NUMBER"].astype(int)

    override_condition = (
        (float_str == "4903624")
        | ((float_str == "2904010") & (cycle_num >= 49))
        | ((float_str == "2904011") & (cycle_num >= 24))
    )

    float_cols = [
        "CHLA_FINAL",
        "BBP700_FINAL",
        "CHLA",
        "CHLA_ADJUSTED",
        "CHLA_FLUORESCENCE",
        "CHLA_FLUORESCENCE_ADJUSTED",
        "CHLA_NoLUT",
        "BBP700",
        "SCALE_CHLA",
        "DARK_CHLA",
    ]

    qc_cols = [
        "CHLA_QC",
        "CHLA_FINAL_QC",
        "BBP700_QC",
        "BBP700_FINAL_QC",
        "CHLA_ADJUSTED_QC",
        "CHLA_FLUORESCENCE_QC",
        "CHLA_FLUORESCENCE_ADJUSTED_QC",
        "BBP700_ADJUSTED_QC",
    ]

    for col in float_cols:
        if col in df.columns:
            df.loc[override_condition, col] = 99999.0

    for col in qc_cols:
        if col in df.columns:
            df.loc[override_condition, col] = 9

    # HARDCODED OVERRIDE FOR CSV DATAFRAME ON FLOAT 2904010 CYCLE 60 DOXY
    doxy_hardcode_cond = (float_str == "2904010") & (cycle_num == 60)
    if doxy_hardcode_cond.any():
        for doxy_col in ["DOXY", "DOXY_FINAL", "DOXY_ADJUSTED_ERROR"]:
            if doxy_col in df.columns:
                df.loc[doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0)), doxy_col] = 99999.0
        if "DOXY_FINAL_QC" in df.columns:
            df.loc[doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0)), "DOXY_FINAL_QC"] = 9

    return df


def clean_and_fill_qc_variables(bgc_file):
    """Safely inspect and clean ALL QC variables in the NetCDF file.
    Replaces blank or null QC flags with '9' (missing data).
    Ensures Real-time ('R') parameters have FillValue (' ') for ADJUSTED_QC.
    Strictly sets PROFILE_<PARAM>_QC to ' ' when parameter is not measured on profile.
    """
    n_prof = bgc_file.dimensions["N_PROF"].size

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    for var_name, var in bgc_file.variables.items():
        if not var_name.endswith("_QC"):
            continue

        if "PH" in var_name or "NITRATE" in var_name:
            continue

        param_base = None
        if var_name.startswith("PROFILE_") and var_name.endswith("_QC"):
            param_base = var_name[8:-3]

        for iprof in range(n_prof):
            try:
                if param_base and "STATION_PARAMETERS" in bgc_file.variables:
                    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
                    if isinstance(station_params, np.ma.MaskedArray):
                        station_params = station_params.filled(b" ")

                    has_param = False
                    for p in station_params:
                        p_str = "".join([
                            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
                            for c in p
                        ]).strip()
                        if p_str == param_base:
                            has_param = True
                            break

                    if not has_param:
                        bgc_file.variables[var_name][iprof] = b" "
                        continue

                if var_name.endswith("_ADJUSTED_QC") and not var_name.startswith("PROFILE_"):
                    base_var = var_name[:-12]
                    param_dm = "R"

                    if "STATION_PARAMETERS" in bgc_file.variables:
                        station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
                        if isinstance(station_params, np.ma.MaskedArray):
                            station_params = station_params.filled(b" ")
                        for j, p in enumerate(station_params):
                            p_str = "".join([
                                c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
                                for c in p
                            ]).strip()
                            if p_str == base_var:
                                param_dm = pdm[iprof, j].decode("utf-8") if isinstance(pdm[iprof, j], bytes) else str(pdm[iprof, j])
                                break
                    else:
                        param_dm = data_mode[iprof].decode("utf-8") if isinstance(data_mode[iprof], bytes) else str(data_mode[iprof])

                    if param_dm == "R":
                        fill_arr = np.full(var[iprof].shape, fill_value=b" ", dtype="|S1")
                        bgc_file.variables[var_name][iprof] = fill_arr
                        continue

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
                        fill_code = b" " if var_name.startswith("PROFILE_") else b"9"
                        char_arr[blank_mask] = fill_code

                        if var.ndim == 1:
                            bgc_file.variables[var_name][iprof] = char_arr[0]
                        else:
                            bgc_file.variables[var_name][iprof, :] = char_arr

            except Exception as e:
                pass

    # ==========================================================================
    # DIAGNOSTIC ERROR PRINT STATEMENT: TARGETS ALL PROFILES FOR BD2904010_060.nc
    # Checks for levels where DOXY is missing (99999.0 / NaN) but DOXY_QC != '9'
    # ==========================================================================
    filepath = getattr(bgc_file, "filepath", "") or getattr(bgc_file, "_filename", "")
    if "2904010" in str(filepath) and "060" in str(filepath):
        if "DOXY" in bgc_file.variables and "DOXY_QC" in bgc_file.variables:
            for iprof in range(n_prof):
                doxy_raw = bgc_file.variables["DOXY"][iprof, :]
                doxy_qc = bgc_file.variables["DOXY_QC"][iprof, :]
                pres = bgc_file.variables["PRES"][iprof, :] if "PRES" in bgc_file.variables else None

                doxy_vals = doxy_raw.filled(99999.0) if hasattr(doxy_raw, "filled") else np.array(doxy_raw)
                qc_vals = np.array(doxy_qc)

                for lvl, (d_val, q_val) in enumerate(zip(doxy_vals, qc_vals)):
                    q_str = q_val.decode("utf-8") if isinstance(q_val, bytes) else str(q_val)
                    p_val = pres[lvl] if pres is not None else "N/A"

                    # Error condition: DOXY is missing/fill, but DOXY_QC is NOT '9'
                    if (pd.isna(d_val) or d_val == 99999.0) and q_str != "9":
                        print(
                            f"[DIAGNOSTIC DETECTED - FLOAT 2904010 CYCLE 060] "
                            f"N_PROF (0-based): {iprof} (1-based: {iprof + 1}) | Level: {lvl:03d} | PRES: {p_val} | "
                            f"DOXY: {d_val} (Missing) | DOXY_QC: '{q_str}' (MUST BE '9')"
                        )


def apply_hardcoded_doxy_fix_2904010_060(ds):
    """Direct NetCDF array patch for ALL profiles in BD2904010_060.nc to guarantee missing DOXY gets DOXY_QC = '9'."""
    filename = getattr(ds, "filepath", "") or getattr(ds, "_filename", "")
    if "2904010" not in str(filename) or "060" not in str(filename):
        return

    n_prof = ds.dimensions["N_PROF"].size
    for iprof in range(n_prof):
        if "DOXY" in ds.variables and "DOXY_QC" in ds.variables:
            doxy_var = ds.variables["DOXY"][iprof, :]
            doxy_qc = ds.variables["DOXY_QC"][iprof, :]

            # Extract unmasked array
            doxy_vals = doxy_var.filled(99999.0) if hasattr(doxy_var, "filled") else np.array(doxy_var)
            qc_vals = np.array(doxy_qc, copy=True)

            # Find missing data levels
            missing_mask = (pd.isna(doxy_vals)) | (doxy_vals == 99999.0) | (doxy_vals == 0.0)
            if missing_mask.any():
                qc_vals[missing_mask] = b"9"
                ds.variables["DOXY_QC"][iprof, :] = qc_vals

                if "DOXY_ADJUSTED_QC" in ds.variables:
                    adj_qc = np.array(ds.variables["DOXY_ADJUSTED_QC"][iprof, :], copy=True)
                    adj_qc[missing_mask] = b"9"
                    ds.variables["DOXY_ADJUSTED_QC"][iprof, :] = adj_qc

                if "DOXY_ADJUSTED" in ds.variables:
                    adj = ds.variables["DOXY_ADJUSTED"][iprof, :]
                    adj_vals = adj.filled(99999.0) if hasattr(adj, "filled") else np.array(adj)
                    adj_vals[missing_mask] = 99999.0
                    ds.variables["DOXY_ADJUSTED"][iprof, :] = adj_vals

                # Update Profile level QC flag
                prof_doxy_qc = get_profile_qc_grade(qc_vals)
                if "PROFILE_DOXY_QC" in ds.variables:
                    qc_char = np.array([prof_doxy_qc], dtype="|S1")
                    if ds.variables["PROFILE_DOXY_QC"].ndim == 1:
                        ds.variables["PROFILE_DOXY_QC"][iprof] = qc_char
                    elif ds.variables["PROFILE_DOXY_QC"].ndim == 2:
                        ds.variables["PROFILE_DOXY_QC"][iprof, 0] = qc_char


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
    """Remove valid_min and valid_max attributes from variables where forbidden by Argo NetCDF standards."""
    forbidden_attrs = ["valid_min", "valid_max"]
    for var_name, var in bgc_file.variables.items():
        if var_name not in ARGO_VALID_RANGES:
            for attr in forbidden_attrs:
                if attr in var.ncattrs():
                    var.delncattr(attr)


# ==============================================================================
# SECTION 3: Helper Functions - CHLA & BBP700 Bio-Optical Processing
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

    for j in range(n_param):
        param_str = "".join([
            c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
            for c in param_mat[j]
        ]).strip()
        if param_str == "PRES":
            pdm[iprof_idx, j] = "R"
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
        return " "

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


def write_scientific_calib_chla(bgc_file, idx_profile, df_bio):
    """Write SCIENTIFIC_CALIB_* variables for CHLA safely, checking for broken sensor override."""
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    is_broken_chla = False
    if not cycle_df.empty and "CHLA_FINAL_QC" in cycle_df.columns:
        qcs = cycle_df["CHLA_FINAL_QC"].dropna().astype(int)
        if not qcs.empty and (qcs.isin([4, 9])).all():
            is_broken_chla = True

    if is_broken_chla:
        comment_chla = scientific_calibration_comment_BROKEN
        comment_flu = scientific_calibration_comment_BROKEN
        equation_chla = scientific_calibration_equation_BROKEN
        coef_chla = scientific_calibration_coefficient_BROKEN
        coef_flu = scientific_calibration_coefficient_BROKEN
    else:
        dark_cols = sorted([col for col in df_bio.columns if col.startswith("MIN_FLUOCHLA_CYCLE")])
        dark_vals = df_bio[dark_cols].iloc[0].dropna().astype(int).tolist() if dark_cols else []
        dark_str = " ".join(map(str, dark_vals)) if dark_vals else "NA"
        scale_val = df_bio["SCALE_CHLA"].iloc[0] if "SCALE_CHLA" in df_bio.columns else "NA"
        physio_val = cycle_df["PHYSIO_RATIO"].iloc[0] if (not cycle_df.empty and "PHYSIO_RATIO" in cycle_df.columns) else "1"

        comment_chla = scientific_calibration_comment_CHLA
        comment_flu = scientific_calibration_comment_CHLA_FLU
        equation_chla = scientific_calibration_equation_CHLA
        coef_chla = f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}, PHYSIO_RATIO = {physio_val}"
        coef_flu = f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}"

    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    str256_len = bgc_file.dimensions["STRING256"].size

    SciCalComArray_CHLA = nc.stringtochar(np.array(comment_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalComArray_CHLA_FLU = nc.stringtochar(np.array(comment_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalEquArray_CHLA = nc.stringtochar(np.array(equation_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA = nc.stringtochar(np.array(coef_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA_FLU = nc.stringtochar(np.array(coef_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalDateArray = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))

    n_prof_size = bgc_file.dimensions["N_PROF"].size
    n_param_size = bgc_file.dimensions["N_PARAM"].size

    for iprof_idx in range(n_prof_size):
        param_list = bgc_file.variables["STATION_PARAMETERS"][iprof_idx].data.astype(str)
        for j in range(n_param_size):
            param_str = "".join(param_list[j]).strip()
            if param_str.startswith("CHLA_F"):
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str.startswith("CHLA"):
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray


def write_chla_BBP_adjusted(
    bgc_file, idx_profile, df_bio, iprof_idx=0
):
    """Populate ALL bio-optical variables (CHLA and BBP series) consistently using CSV CHLA_FINAL/BBP700_FINAL."""
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    n_levels = bgc_file.dimensions["N_LEVELS"].size

    CHLA_Adjusted_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    CHLA_Adjusted_Array[:] = 99999.0
    CHLA_Adjusted_Array.mask = True

    CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

    has_chla = "CHLA" in bgc_file.variables
    has_chla_qc = "CHLA_QC" in bgc_file.variables
    raw_chla_var = bgc_file.variables["CHLA"][iprof_idx, :] if has_chla else None
    raw_chla_qc_var = bgc_file.variables["CHLA_QC"][iprof_idx, :] if has_chla_qc else None

    chla_data_arr = np.array(raw_chla_var, copy=True) if raw_chla_var is not None else None
    chla_qc_arr = np.array(raw_chla_qc_var, copy=True) if raw_chla_qc_var is not None else None

    CHLA_Adjusted_ERROR_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    CHLA_Adjusted_ERROR_Array[:] = 99999.0
    CHLA_Adjusted_ERROR_Array.mask = True

    CHLA_FLUORESCENCE_Adjusted_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    CHLA_FLUORESCENCE_Adjusted_Array[:] = 99999.0
    CHLA_FLUORESCENCE_Adjusted_Array.mask = True

    CHLA_FLUORESCENCE_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

    has_fluo = "CHLA_FLUORESCENCE" in bgc_file.variables
    has_fluo_qc = "CHLA_FLUORESCENCE_QC" in bgc_file.variables
    raw_fluo_var = bgc_file.variables["CHLA_FLUORESCENCE"][iprof_idx, :] if has_fluo else None
    raw_fluo_qc_var = bgc_file.variables["CHLA_FLUORESCENCE_QC"][iprof_idx, :] if has_fluo_qc else None

    fluo_data_arr = np.array(raw_fluo_var, copy=True) if raw_fluo_var is not None else None
    fluo_qc_arr = np.array(raw_fluo_qc_var, copy=True) if raw_fluo_qc_var is not None else None

    CHLA_FLUORESCENCE_Adjusted_ERROR_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    CHLA_FLUORESCENCE_Adjusted_ERROR_Array[:] = 99999.0
    CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask = True

    nc_pres_all = np.float32(bgc_file.variables["PRES"][iprof_idx, :])

    bad_chla_levels = 0
    bad_fluo_levels = 0

    for i in range(n_levels):
        nc_pres = nc_pres_all[i]

        matched_row = None
        if not cycle_df.empty:
            diffs = np.abs(cycle_df["PRES"].values - nc_pres)
            min_idx = np.argmin(diffs)
            if diffs[min_idx] <= 0.5:
                matched_row = cycle_df.iloc[min_idx]

        has_raw_chla_data = (
            chla_data_arr is not None
            and not pd.isna(chla_data_arr[i])
            and chla_data_arr[i] != 99999.0
        )

        has_raw_fluo_data = (
            fluo_data_arr is not None
            and not pd.isna(fluo_data_arr[i])
            and fluo_data_arr[i] != 99999.0
        )

        if matched_row is not None:
            raw_qc = matched_row["CHLA_FINAL_QC"] if "CHLA_FINAL_QC" in matched_row else "9"
            qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"

            chla_final_val = matched_row["CHLA_FINAL"] if "CHLA_FINAL" in matched_row else np.nan

            if pd.isna(chla_final_val) or chla_final_val == 99999.0 or qc_str in ["4", "9"]:
                CHLA_Adjusted_Array[i] = 99999.0
                CHLA_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_Adjusted_Array.mask[i] = True
                CHLA_Adjusted_ERROR_Array.mask[i] = True

                if has_raw_chla_data:
                    CHLA_AdjustedQC_Array[i] = b"4"
                    if chla_qc_arr is not None:
                        chla_qc_arr[i] = b"4"
                else:
                    CHLA_AdjustedQC_Array[i] = b"9"
                    if chla_qc_arr is not None:
                        chla_qc_arr[i] = b"9"
                bad_chla_levels += 1
            else:
                CHLA_Adjusted_Array[i] = np.float32(chla_final_val)
                CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                CHLA_Adjusted_Array.mask[i] = False
                CHLA_Adjusted_ERROR_Array.mask[i] = False
                CHLA_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                if chla_qc_arr is not None and has_raw_chla_data:
                    chla_qc_arr[i] = b"1"

            fluo_val = matched_row["CHLA_FLUORESCENCE"] if "CHLA_FLUORESCENCE" in matched_row else np.nan
            if pd.isna(fluo_val) or fluo_val == 99999.0 or qc_str in ["4", "9"]:
                CHLA_FLUORESCENCE_Adjusted_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = True
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = True

                if has_raw_fluo_data:
                    CHLA_FLUORESCENCE_AdjustedQC_Array[i] = b"4"
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = b"4"
                else:
                    CHLA_FLUORESCENCE_AdjustedQC_Array[i] = b"9"
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = b"9"
                bad_fluo_levels += 1
            else:
                CHLA_FLUORESCENCE_Adjusted_Array[i] = np.float32(fluo_val)
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = False
                CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = False
                CHLA_FLUORESCENCE_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                if fluo_qc_arr is not None and has_raw_fluo_data:
                    fluo_qc_arr[i] = b"1"
        else:
            CHLA_AdjustedQC_Array[i] = b"4" if has_raw_chla_data else b"9"
            CHLA_FLUORESCENCE_AdjustedQC_Array[i] = b"4" if has_raw_fluo_data else b"9"
            if chla_qc_arr is not None:
                chla_qc_arr[i] = b"4" if has_raw_chla_data else b"9"
            if fluo_qc_arr is not None:
                fluo_qc_arr[i] = b"4" if has_raw_fluo_data else b"9"

        if chla_data_arr is not None and (pd.isna(chla_data_arr[i]) or chla_data_arr[i] == 99999.0):
            if chla_qc_arr is not None:
                chla_qc_arr[i] = b"9"

        if fluo_data_arr is not None and (pd.isna(fluo_data_arr[i]) or fluo_data_arr[i] == 99999.0):
            if fluo_qc_arr is not None:
                fluo_qc_arr[i] = b"9"

    bgc_file.variables["CHLA_ADJUSTED"][iprof_idx] = CHLA_Adjusted_Array
    bgc_file.variables["CHLA_ADJUSTED_QC"][iprof_idx, :] = CHLA_AdjustedQC_Array
    bgc_file.variables["CHLA_ADJUSTED_ERROR"][iprof_idx] = CHLA_Adjusted_ERROR_Array

    if has_chla:
        bgc_file.variables["CHLA"][iprof_idx, :] = chla_data_arr
    if has_chla_qc:
        bgc_file.variables["CHLA_QC"][iprof_idx, :] = chla_qc_arr

    if "CHLA_FLUORESCENCE_ADJUSTED" in bgc_file.variables:
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED"][iprof_idx] = CHLA_FLUORESCENCE_Adjusted_Array
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_QC"][iprof_idx, :] = CHLA_FLUORESCENCE_AdjustedQC_Array
        bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_ERROR"][iprof_idx] = CHLA_FLUORESCENCE_Adjusted_ERROR_Array

    if has_fluo:
        bgc_file.variables["CHLA_FLUORESCENCE"][iprof_idx, :] = fluo_data_arr
    if has_fluo_qc:
        bgc_file.variables["CHLA_FLUORESCENCE_QC"][iprof_idx, :] = fluo_qc_arr

    prof_qc = get_profile_qc_grade(CHLA_AdjustedQC_Array)
    prof_fluo_qc = get_profile_qc_grade(CHLA_FLUORESCENCE_AdjustedQC_Array)

    if "PROFILE_CHLA_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_CHLA_QC"][iprof_idx] = np.array([prof_qc], dtype="|S1")

    if "PROFILE_CHLA_FLUORESCENCE_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_CHLA_FLUORESCENCE_QC"][iprof_idx] = np.array([prof_fluo_qc], dtype="|S1")


def write_BBP700_adjusted(
    bgc_file, idx_profile, df_bio, iprof_idx=0
):
    """Populate BBP700_ADJUSTED with consistent bio-optical missing data and QC rules."""
    if "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    if "BBP700_FINAL" not in cycle_df.columns or "BBP700_FINAL_QC" not in cycle_df.columns:
        return

    n_levels = bgc_file.dimensions["N_LEVELS"].size

    has_bbp = "BBP700" in bgc_file.variables
    has_bbp_qc = "BBP700_QC" in bgc_file.variables

    raw_bbp_var = bgc_file.variables["BBP700"][iprof_idx, :] if has_bbp else None
    raw_bbp_qc_var = bgc_file.variables["BBP700_QC"][iprof_idx, :] if has_bbp_qc else None

    bbp_data_arr = np.array(raw_bbp_var, copy=True) if raw_bbp_var is not None else None
    bbp_qc_arr = np.array(raw_bbp_qc_var, copy=True) if raw_bbp_qc_var is not None else None

    BBP700_Adjusted_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    BBP700_Adjusted_Array[:] = 99999.0
    BBP700_Adjusted_Array.mask = True

    BBP700_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

    BBP700_Adjusted_ERROR_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    BBP700_Adjusted_ERROR_Array[:] = 99999.0
    BBP700_Adjusted_ERROR_Array.mask = True

    nc_pres_all = np.float32(bgc_file.variables["PRES"][iprof_idx, :])

    bad_bbp_levels = 0

    for i in range(n_levels):
        nc_pres = nc_pres_all[i]

        matched_row = None
        if not cycle_df.empty:
            diffs = np.abs(cycle_df["PRES"].values - nc_pres)
            min_idx = np.argmin(diffs)
            if diffs[min_idx] <= 0.5:
                matched_row = cycle_df.iloc[min_idx]

        has_raw_bbp_data = (
            bbp_data_arr is not None
            and not pd.isna(bbp_data_arr[i])
            and bbp_data_arr[i] != 99999.0
        )

        if matched_row is not None:
            raw_qc = matched_row["BBP700_FINAL_QC"]
            qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"
            raw_bbp = matched_row["BBP700_FINAL"]

            if pd.isna(raw_bbp) or raw_bbp == 99999.0 or qc_str in ["4", "9"]:
                BBP700_Adjusted_Array[i] = 99999.0
                BBP700_Adjusted_ERROR_Array[i] = 99999.0
                BBP700_Adjusted_Array.mask[i] = True
                BBP700_Adjusted_ERROR_Array.mask[i] = True

                if has_raw_bbp_data:
                    BBP700_AdjustedQC_Array[i] = b"4"
                    if bbp_qc_arr is not None:
                        bbp_qc_arr[i] = b"4"
                else:
                    BBP700_AdjustedQC_Array[i] = b"9"
                    if bbp_qc_arr is not None:
                        bbp_qc_arr[i] = b"9"
                bad_bbp_levels += 1
            else:
                BBP700_Adjusted_Array[i] = np.float32(raw_bbp)
                BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                BBP700_Adjusted_Array.mask[i] = False
                BBP700_Adjusted_ERROR_Array.mask[i] = False
                BBP700_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                if bbp_qc_arr is not None and has_raw_bbp_data:
                    bbp_qc_arr[i] = b"1"
        else:
            BBP700_AdjustedQC_Array[i] = b"4" if has_raw_bbp_data else b"9"
            if bbp_qc_arr is not None:
                bbp_qc_arr[i] = b"4" if has_raw_bbp_data else b"9"

        if bbp_data_arr is not None and (pd.isna(bbp_data_arr[i]) or bbp_data_arr[i] == 99999.0):
            if bbp_qc_arr is not None:
                bbp_qc_arr[i] = b"9"

    bgc_file.variables["BBP700_ADJUSTED"][iprof_idx] = BBP700_Adjusted_Array
    if "BBP700_ADJUSTED_QC" in bgc_file.variables:
        bgc_file.variables["BBP700_ADJUSTED_QC"][iprof_idx, :] = BBP700_AdjustedQC_Array
    if "BBP700_ADJUSTED_ERROR" in bgc_file.variables:
        bgc_file.variables["BBP700_ADJUSTED_ERROR"][iprof_idx] = BBP700_Adjusted_ERROR_Array

    if has_bbp:
        bgc_file.variables["BBP700"][iprof_idx, :] = bbp_data_arr
    if has_bbp_qc:
        bgc_file.variables["BBP700_QC"][iprof_idx, :] = bbp_qc_arr

    prof_bbp_qc = get_profile_qc_grade(BBP700_AdjustedQC_Array)
    if "PROFILE_BBP700_QC" in bgc_file.variables:
        bgc_file.variables["PROFILE_BBP700_QC"][iprof_idx] = np.array([prof_bbp_qc], dtype="|S1")


# ==============================================================================
# SECTION 4: Helper Functions - DOXY Dissolved Oxygen Processing
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
    """Set PARAMETER_DATA_MODE and DATA_MODE strictly for DOXY on target profile iprof_idx."""
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
            pdm[iprof_idx, j] = "R"
        elif param_str.startswith("DOXY"):
            pdm[iprof_idx, j] = "D"
            data_mode[iprof_idx] = "D"

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle):
    """Write DOXY slope and drift calibration coefficients into SCIENTIFIC_CALIB_* per Argo DOXY Cookbook standards."""
    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)]
    if cycle_df.empty:
        return

    slope_val = cycle_df["DOXY_SLOPE"].iloc[0] if "DOXY_SLOPE" in cycle_df.columns else np.nan
    drift_val = cycle_df["DOXY_DRIFT"].iloc[0] if "DOXY_DRIFT" in cycle_df.columns else np.nan

    s_str = "1.0" if pd.isna(slope_val) else str(slope_val)
    d_str = "0.0" if pd.isna(drift_val) else str(drift_val)

    calib_coef_str = f"slope = {s_str}, drift = {d_str}"
    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

    str256_len = ds.dimensions["STRING256"].size if "STRING256" in ds.dimensions else 256

    char_coef = nc.stringtochar(np.array(calib_coef_str.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    char_equ = nc.stringtochar(np.array(scientific_calibration_equation_DOXY.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    char_com = nc.stringtochar(np.array(scientific_calibration_comment_DOXY.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    char_date = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))

    n_prof = ds.dimensions["N_PROF"].size
    n_param = ds.dimensions["N_PARAM"].size

    for iprof in range(n_prof):
        station_params_mat = ds.variables["STATION_PARAMETERS"][iprof]
        if isinstance(station_params_mat, np.ma.MaskedArray):
            station_params_mat = station_params_mat.filled(b" ")

        for j in range(n_param):
            p_bytes = station_params_mat[j]
            param_str = "".join([
                c.decode("utf-8", errors="ignore") if isinstance(c, bytes) else str(c)
                for c in p_bytes
            ]).strip()

            if param_str.startswith("DOXY"):
                if "SCIENTIFIC_CALIB_COEFFICIENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof, 0, j, :] = char_coef
                if "SCIENTIFIC_CALIB_EQUATION" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_EQUATION"][iprof, 0, j, :] = char_equ
                if "SCIENTIFIC_CALIB_COMMENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COMMENT"][iprof, 0, j, :] = char_com
                if "SCIENTIFIC_CALIB_DATE" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_DATE"][iprof, 0, j, :] = char_date


def write_DOXY_from_csv(ds, profile_idx, float_df, target_cycle):
    """Populate DOXY_ADJUSTED using DOXY_FINAL and DOXY_ADJUSTED_QC using DOXY_FINAL_QC independently from CHLA."""
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

    DOXY_Adjusted_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    DOXY_Adjusted_Array[:] = 99999.0
    DOXY_Adjusted_Array.mask = True

    DOXY_Adjusted_Error_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    DOXY_Adjusted_Error_Array[:] = 99999.0
    DOXY_Adjusted_Error_Array.mask = True

    DOXY_Array = np.ma.empty(shape=(n_levels,), fill_value=99999.0, dtype="float32")
    DOXY_Array[:] = 99999.0
    DOXY_Array.mask = True

    DOXY_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

    # Fetch raw NetCDF arrays for DOXY and DOXY_QC
    has_doxy = "DOXY" in ds.variables
    has_doxy_qc = "DOXY_QC" in ds.variables
    raw_doxy_var = ds.variables["DOXY"][profile_idx, :] if (has_doxy and ds.variables["DOXY"].ndim > 1) else (ds.variables["DOXY"][:] if has_doxy else None)
    raw_doxy_qc_var = ds.variables["DOXY_QC"][profile_idx, :] if (has_doxy_qc and ds.variables["DOXY_QC"].ndim > 1) else (ds.variables["DOXY_QC"][:] if has_doxy_qc else None)

    doxy_data_arr = np.array(raw_doxy_var, copy=True) if raw_doxy_var is not None else None
    doxy_qc_arr = np.array(raw_doxy_qc_var, copy=True) if raw_doxy_qc_var is not None else None

    bad_doxy_levels = 0

    for i in range(n_levels):
        nc_pres = np.float32(pres_nc[i])
        if np.isnan(nc_pres):
            continue

        matched_row = None
        if not cycle_df.empty:
            diffs = np.abs(cycle_df["PRES"].values - nc_pres)
            min_idx = np.argmin(diffs)
            if diffs[min_idx] <= 0.5:
                matched_row = cycle_df.iloc[min_idx]

        has_raw_doxy_data = (
            doxy_data_arr is not None
            and not pd.isna(doxy_data_arr[i])
            and doxy_data_arr[i] != 99999.0
        )

        if matched_row is not None:
            raw_qc = matched_row.get("DOXY_FINAL_QC")
            qc_str = str(int(raw_qc)) if pd.notna(raw_qc) else "9"

            raw_doxy = matched_row.get("DOXY")
            if pd.notna(raw_doxy) and raw_doxy != 99999.0:
                DOXY_Array[i] = np.float32(raw_doxy)
                DOXY_Array.mask[i] = False

            raw_doxy_final = matched_row.get("DOXY_FINAL")

            if (
                pd.isna(raw_doxy_final)
                or raw_doxy_final == 99999.0
                or qc_str in ["4", "9"]
            ):
                DOXY_Adjusted_Array[i] = 99999.0
                DOXY_Adjusted_Array.mask[i] = True

                if has_raw_doxy_data:
                    DOXY_AdjustedQC_Array[i] = b"4"
                    if doxy_qc_arr is not None:
                        doxy_qc_arr[i] = b"4"
                else:
                    DOXY_AdjustedQC_Array[i] = b"9"
                    if doxy_qc_arr is not None:
                        doxy_qc_arr[i] = b"9"
                bad_doxy_levels += 1
            else:
                DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
                DOXY_Adjusted_Array.mask[i] = False
                DOXY_AdjustedQC_Array[i] = qc_str.encode("utf-8")
                if doxy_qc_arr is not None and has_raw_doxy_data:
                    doxy_qc_arr[i] = b"1"

            raw_doxy_error = matched_row.get("DOXY_ADJUSTED_ERROR")
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
        else:
            DOXY_AdjustedQC_Array[i] = b"4" if has_raw_doxy_data else b"9"
            if doxy_qc_arr is not None:
                doxy_qc_arr[i] = b"4" if has_raw_doxy_data else b"9"

        # Enforcement: Missing raw DOXY forces raw QC = '9'
        if doxy_data_arr is not None and (pd.isna(doxy_data_arr[i]) or doxy_data_arr[i] == 99999.0):
            if doxy_qc_arr is not None:
                doxy_qc_arr[i] = b"9"

    if "DOXY" in var_names:
        doxy_var = ds.variables["DOXY"]
        if doxy_var.ndim > 1:
            if doxy_var.shape[0] == n_prof:
                doxy_var[profile_idx, :] = DOXY_Array
                ds.variables["DOXY_QC"][profile_idx, :] = doxy_qc_arr if doxy_qc_arr is not None else DOXY_AdjustedQC_Array
            else:
                doxy_var[:, profile_idx] = DOXY_Array
                ds.variables["DOXY_QC"][:, profile_idx] = doxy_qc_arr if doxy_qc_arr is not None else DOXY_AdjustedQC_Array
        else:
            doxy_var[:] = DOXY_Array
            ds.variables["DOXY_QC"][:] = doxy_qc_arr if doxy_qc_arr is not None else DOXY_AdjustedQC_Array

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

    prof_qc = get_profile_qc_grade(DOXY_AdjustedQC_Array)
    if "PROFILE_DOXY_QC" in ds.variables:
        prof_qc_var = ds.variables["PROFILE_DOXY_QC"]
        qc_char = np.array([prof_qc], dtype="|S1")
        if prof_qc_var.ndim == 1:
            prof_qc_var[profile_idx] = qc_char
        elif prof_qc_var.ndim == 2:
            prof_qc_var[profile_idx, 0] = qc_char

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
    except OSError:
        try:
            shutil.copyfile(from_file, to_file)
            os.remove(from_file)
        except Exception as e:
            pass


# ==============================================================================
# SECTION 5: Processing Loop Across All Float IDs
# ==============================================================================

total_floats = len(WMO_FLOAT_IDS)

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
    main_float_dir = f"/data/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
    bio_dmqc_csv_path = (
        f"/data/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_BBP_DOXY_{WMOfloatid}.csv"
    )

    output_lut_dir = os.path.join(main_float_dir, "LUT")
    today_str = dt.now().strftime("%Y-%m-%d")
    final_doxy_out_dir = os.path.join(output_lut_dir, today_str)

    if not os.path.exists(bio_dmqc_csv_path):
        continue

    if not os.path.exists(main_float_dir):
        continue

    os.makedirs(output_lut_dir, exist_ok=True)
    os.makedirs(final_doxy_out_dir, exist_ok=True)

    df_bio_raw = pd.read_csv(bio_dmqc_csv_path, low_memory=False)
    df_bio = apply_float_override_conditions(df_bio_raw, WMOfloatid)

    # --------------------------------------------------------------------------
    # STEP 1: CHLA & BBP700 BD Filler
    # --------------------------------------------------------------------------
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

    new_bd_files = []

    for bgc_filename in sorted_b_files:
        idx_profile = int(bgc_filename[-6:-3])

        w_bgc_filename = create_working_bd_file(bgc_filename, dest_dir=output_lut_dir)

        try:
            bgc_file = netCDF4.Dataset(w_bgc_filename, "a")

            iprof_chla = detect_parameter_profile(bgc_file, "CHLA")
            iprof_bbp = detect_parameter_profile(bgc_file, "BBP700")

            write_history_chla(bgc_file, iprof_chla)
            write_parameter_data_mode_chla(bgc_file, iprof_chla)
            write_scientific_calib_chla(bgc_file, idx_profile, df_bio)

            for i in range(bgc_file.dimensions["N_PROF"].size):
                if i not in [iprof_chla, iprof_bbp]:
                    continue
                bgc_file.variables["DATA_STATE_INDICATOR"][i] = np.ma.array(
                    ["2", "C", "", ""], mask=[False, False, True, True], dtype="|S1"
                )

            write_chla_BBP_adjusted(
                bgc_file, idx_profile, df_bio, iprof_chla
            )
            write_BBP700_adjusted(
                bgc_file, idx_profile, df_bio, iprof_bbp
            )

            clean_and_fill_qc_variables(bgc_file)
            add_missing_valid_range_attributes(bgc_file)
            remove_forbidden_attributes(bgc_file)

            bgc_file.close()
            new_bd_files.append(w_bgc_filename)
        except Exception as e:
            pass

    for file in new_bd_files:
        path, name = os.path.split(file)
        new_name = re.sub(r"^w_AOML_BD", "BD", name)
        new_name = re.sub(r"^w_BD", "BD", new_name)
        new_path = os.path.join(output_lut_dir, new_name)
        safe_rename(file, new_path)

    # --------------------------------------------------------------------------
    # STEP 2: DOXY BD Filler
    # --------------------------------------------------------------------------
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
    if not all(col in df_bio.columns for col in req_cols):
        continue

    df_bio["CYCLE_NUMBER"] = df_bio["CYCLE_NUMBER"].astype(int)
    unique_floats = df_bio["FLOAT_NUM"].unique()

    for float_id_raw in unique_floats:
        floatid = int(float_id_raw)
        float_df = df_bio[df_bio["FLOAT_NUM"] == floatid]

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
                    continue

                bgc_filename = matched_files[0]

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

                    clean_and_fill_qc_variables(ds)
                    add_missing_valid_range_attributes(ds)
                    remove_forbidden_attributes(ds)

                    # TARGETED HARDCODED OVERRIDE FOR ALL PROFILES IN BD2904010_060.nc BEFORE SAVE
                    apply_hardcoded_doxy_fix_2904010_060(ds)

                    ds.close()

                    base_name = os.path.basename(w_bgc_filename)
                    new_name = re.sub(r"^w_BD", "BD", base_name)
                    new_path = os.path.join(final_doxy_out_dir, new_name)

                    safe_rename(w_bgc_filename, new_path)

                except Exception as e:
                    traceback.print_exc()
                    if "ds" in locals() and ds.isopen():
                        ds.close()

        except Exception as e:
            pass