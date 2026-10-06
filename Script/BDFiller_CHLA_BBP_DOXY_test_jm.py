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
from collections import Counter
from datetime import datetime as dt, timezone

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
    "CHLA_FLUORESCENCE (specified in http://dx.doi.org/10.13155/35385 and"
    " computed with MLD_LIMIT = 0.03)"
)
scientific_calibration_equation_CHLA = (
    "CHLA_ADJUSTED = CHLA_NPQ for PRES in [0, ZMaxFluo], where "
    "CHLA = ((FLUORESCENCE_CHLA - MEDIAN(PRELIM_DARK_CHLA)) * SCALE_CHLA) /"
    " PHYSIO_RATIO"
)
scientific_calibration_coefficient_CHLA = "PHYSIO_RATIO=1.0"
CHLA_Adjusted_ERROR_est = 0.07
BBP700_Adjusted_ERROR_est = 0.0005  # Default estimated error for BBP700

# Argo Cookbook DOXY standard strings
scientific_calibration_equation_DOXY = "DOXY_ADJUSTED = DOXY * slope + drift"
scientific_calibration_comment_DOXY = (
    "DOXY calibration computed against WOA climatology following Bittig et"
    " al. (2018) methodology"
)

# Argo Cookbook Broken Sensor strings
scientific_calibration_comment_BROKEN = (
    "Sensor failure / broken sensor - data flagged bad or missing (Argo BGC QC"
    " Manual)"
)
scientific_calibration_equation_BROKEN = "Not applicable"
scientific_calibration_coefficient_BROKEN = "Not applicable"

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


def is_broken_sensor_float(WMOfloatid, cycle_num):
    """Evaluates if the float and cycle represent a known broken bio-optical sensor."""
    float_str = str(WMOfloatid)
    cycle_num = int(cycle_num)
    return (
        (float_str == "4903624")
        | ((float_str == "2904010") & (cycle_num >= 49))
        | ((float_str == "2904011") & (cycle_num >= 24))
    )


def apply_float_override_conditions(df, WMOfloatid):
    """Applies Float ID and Cycle-Specific overrides for broken bio-optical sensors in CSV data."""
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
            df.loc[override_condition, col] = 4  # Flag bad data per Argo manual

    # HARDCODED OVERRIDE FOR CSV DATAFRAME ON FLOAT 2904010 CYCLE 60 DOXY
    doxy_hardcode_cond = (float_str == "2904010") & (cycle_num == 60)
    if doxy_hardcode_cond.any():
        for doxy_col in ["DOXY", "DOXY_FINAL", "DOXY_ADJUSTED_ERROR"]:
            if doxy_col in df.columns:
                df.loc[
                    doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0)),
                    doxy_col,
                ] = 99999.0

        doxy_qc_target_cols = [
            col for col in ["DOXY_FINAL_QC", "DOXY_QC"] if col in df.columns
        ]
        if doxy_qc_target_cols:
            df.loc[
                doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0)),
                doxy_qc_target_cols,
            ] = 9

    return df


def clean_and_fill_qc_variables(bgc_file):
    """Safely inspect and clean ALL QC variables across ALL profiles in the NetCDF file."""
    n_prof = bgc_file.dimensions["N_PROF"].size
    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:] if "PARAMETER_DATA_MODE" in bgc_file.variables else None

    for var_name, var in list(bgc_file.variables.items()):
        if not var_name.endswith("_QC"):
            continue

        if "PH" in var_name or "NITRATE" in var_name or "TEMP_CPU_CHLA" in var_name:
            continue

        # 1. PROFILE_<PARAM>_QC Sanitation
        if var_name.startswith("PROFILE_"):
            param_base = var_name[8:-3]
            for iprof in range(n_prof):
                try:
                    if "STATION_PARAMETERS" in bgc_file.variables:
                        station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
                        if isinstance(station_params, np.ma.MaskedArray):
                            station_params = station_params.filled(b" ")

                        has_param = False
                        for p in station_params:
                            p_str = "".join([
                                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                                for c in p
                            ]).strip()
                            if p_str == param_base:
                                has_param = True
                                break

                        if not has_param:
                            bgc_file.variables[var_name][iprof] = b" "
                except Exception:
                    traceback.print_exc()
            continue

        # 2. _ADJUSTED_QC and Raw _QC Synchronization
        is_adjusted_qc = var_name.endswith("_ADJUSTED_QC")
        base_param = var_name[:-12] if is_adjusted_qc else var_name[:-3]

        for iprof in range(n_prof):
            try:
                param_dm = "R"
                if "STATION_PARAMETERS" in bgc_file.variables and pdm is not None:
                    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
                    if isinstance(station_params, np.ma.MaskedArray):
                        station_params = station_params.filled(b" ")
                    for j, p in enumerate(station_params):
                        p_str = "".join([
                            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                            for c in p
                        ]).strip()
                        if p_str == base_param:
                            raw_pdm = pdm[iprof, j]
                            if isinstance(raw_pdm, (bytes, np.bytes_)):
                                param_dm = raw_pdm.tobytes().decode("utf-8").strip() if hasattr(raw_pdm, "tobytes") else raw_pdm.decode("utf-8").strip()
                            else:
                                param_dm = str(raw_pdm).strip()
                            break

                if is_adjusted_qc and param_dm == "R":
                    target_shape = var[iprof].shape if var.ndim > 1 else var.shape
                    fill_arr = np.full(target_shape, fill_value=b" ", dtype="|S1")
                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = fill_arr
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = fill_arr if fill_arr.ndim == 0 else fill_arr[0]
                    else:
                        bgc_file.variables[var_name][:] = fill_arr
                    continue

                p_var_name = base_param if not is_adjusted_qc else f"{base_param}_ADJUSTED"
                if p_var_name in bgc_file.variables:
                    p_var = bgc_file.variables[p_var_name]
                    p_raw = p_var[iprof, :] if p_var.ndim > 1 else p_var[:]
                    p_vals = p_raw.filled(99999.0) if hasattr(p_raw, "filled") else np.array(p_raw)

                    q_raw = var[iprof, :] if var.ndim > 1 else var[iprof]
                    q_vals = np.array(q_raw, copy=True)
                    if hasattr(q_vals, "filled"):
                        q_vals = q_vals.filled(b" ")

                    q_char = q_vals.astype("|S1")

                    missing_mask = (pd.isna(p_vals)) | (p_vals == 99999.0)
                    valid_data_mask = ~missing_mask

                    if missing_mask.any():
                        q_char[missing_mask] = b"9"

                    invalid_qc_mask = valid_data_mask & (
                        (q_char == b" ")
                        | (q_char == b"")
                        | (q_char == b"0")
                        | (q_char == b"\x00")
                    )

                    if invalid_qc_mask.any():
                        q_char[invalid_qc_mask] = b"1"

                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = q_char
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = q_char if q_char.ndim == 0 else q_char[0]
                    else:
                        bgc_file.variables[var_name][:] = q_char

            except Exception:
                traceback.print_exc()


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
                (
                    c.decode("utf-8", errors="ignore")
                    if isinstance(c, (bytes, np.bytes_))
                    else str(c)
                )
                for c in param_mat[j]
            ]).strip()
            if param_str.startswith(param_prefix):
                return iprof
    return 0


def create_working_bd_file(filename, dest_dir=None):
    """Creates a working copy of the input B file with a leading 'w_' in dest_dir or same directory."""
    path, name = os.path.split(filename)
    bd_name = re.sub(r"^(AOML_)?BR", "BD", name)
    out_dir = dest_dir if dest_dir else path
    w_filename = os.path.join(out_dir, f"w_{bd_name}")
    shutil.copyfile(filename, w_filename)
    return w_filename


def get_profile_chla(filename):
    """Extract profile index from filename, return as int."""
    profile = filename[-6:-3]
    return int(profile)


def organize_b_files(bd_files, br_files):
    """Sort B*nc files by profile index."""
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


def update_history_entry(nc_ds, dct, iprof_idx):
    """Update HISTORY array entries natively using netCDF4."""
    hix = nc_ds.dimensions["N_HISTORY"].size
    target_hix = max(0, hix - 1)
    for name, value in dct.items():
        if name in nc_ds.variables:
            char_len = nc_ds.dimensions[nc_ds[name].dimensions[-1]].size
            padded_val = str(value).ljust(char_len)[:char_len]
            char_arr = np.array(padded_val, dtype=f"S{char_len}")
            nc_ds[name][target_hix, iprof_idx, :] = nc.stringtochar(char_arr)


def write_history_metadata(bgc_file, iprof_chla, iprof_doxy):
    """Write global attributes and HISTORY variables for CHLA and DOXY."""
    bgc_file.history = datetime.datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ creation"
    )
    bgc_file.setncattr("comment_dmqc_operator", comment_dmqc_operator_chla)

    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

    update_history_entry(
        bgc_file,
        {
            "HISTORY_INSTITUTION": history_institution,
            "HISTORY_STEP": "ARSQ",
            "HISTORY_SOFTWARE": history_software_chla,
            "HISTORY_SOFTWARE_RELEASE": history_software_release_chla,
            "HISTORY_REFERENCE": history_reference_chla,
            "HISTORY_DATE": UTCcurrent,
            "HISTORY_ACTION": "IP  ",
            "HISTORY_PARAMETER": history_parameter_chla,
        },
        iprof_chla,
    )

    update_history_entry(
        bgc_file,
        {
            "HISTORY_INSTITUTION": "AO",
            "HISTORY_STEP": "ARSQ",
            "HISTORY_SOFTWARE": "BITTIG",
            "HISTORY_SOFTWARE_RELEASE": "2024",
            "HISTORY_REFERENCE": "WOA2023",
            "HISTORY_DATE": UTCcurrent,
            "HISTORY_ACTION": "IP  ",
            "HISTORY_PARAMETER": "DOXY",
        },
        iprof_doxy,
    )

    if "DATE_UPDATE" in bgc_file.variables:
        bgc_file.variables["DATE_UPDATE"][:] = nc.stringtochar(
            np.array(UTCcurrent, dtype="S14")
        )


def write_parameter_data_modes(bgc_file):
    """Set PARAMETER_DATA_MODE and DATA_MODE strictly for PRES, CHLA, CHLA_FLUORESCENCE, BBP700, and DOXY."""
    n_param = bgc_file.dimensions["N_PARAM"].size
    n_prof = bgc_file.dimensions["N_PROF"].size
    str_param_len = bgc_file.variables["STATION_PARAMETERS"].shape[2]

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    for iprof in range(n_prof):
        param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        existing_params = []
        for j in range(n_param):
            p_bytes = param_mat[j]
            p_str = "".join([
                (
                    c.decode("utf-8", errors="ignore")
                    if isinstance(c, (bytes, np.bytes_))
                    else str(c)
                )
                for c in p_bytes
            ]).strip()

            existing_params.append(p_str)

            if p_str == "PRES":
                pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"
            elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "BBP700", "DOXY"]:
                pdm[iprof, j] = b"D" if pdm.dtype.kind in ["S", "U", "O"] else "D"
                data_mode[iprof] = (
                    b"D" if data_mode.dtype.kind in ["S", "U", "O"] else "D"
                )

        for required_param in ["CHLA", "CHLA_FLUORESCENCE", "BBP700", "DOXY"]:
            if (
                required_param in bgc_file.variables
                and required_param not in existing_params
            ):
                try:
                    empty_slot_idx = existing_params.index("")
                    padded_param = required_param.ljust(str_param_len)[:str_param_len]
                    char_param = nc.stringtochar(
                        np.array(padded_param, dtype=f"S{str_param_len}")
                    )
                    bgc_file.variables["STATION_PARAMETERS"][
                        iprof, empty_slot_idx, :
                    ] = char_param
                    pdm[iprof, empty_slot_idx] = (
                        b"D" if pdm.dtype.kind in ["S", "U", "O"] else "D"
                    )
                    data_mode[iprof] = (
                        b"D" if data_mode.dtype.kind in ["S", "U", "O"] else "D"
                    )
                    existing_params[empty_slot_idx] = required_param
                except ValueError:
                    pass

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def get_profile_qc_grade(qc_masked_array):
    """Calculate Argo profile QC letter grade (A-F, ' ')."""
    if hasattr(qc_masked_array, "compressed"):
        unmasked_vals = qc_masked_array.compressed()
    else:
        unmasked_vals = qc_masked_array

    raw_qcs = []
    for q in unmasked_vals:
        if isinstance(q, np.ma.core.MaskedConstant) or q is np.ma.masked:
            continue
        s = q.decode("utf-8").strip() if isinstance(q, bytes) else str(q).strip()
        if s != "":
            raw_qcs.append(s)

    if len(raw_qcs) == 0 or all(q in ["4", "9"] for q in raw_qcs):
        return " "

    valid_qcs = [q for q in raw_qcs if q not in ["4", "9"]]
    total_points = len(valid_qcs)

    if total_points == 0:
        return " "

    bad_count = sum(1 for q in valid_qcs if q in ["3"])
    good_count = total_points - bad_count
    pct_good = (good_count / total_points) * 100.0

    if bad_count == 0:
        return "A"
    elif 75.0 <= pct_good < 100.0:
        return "B"
    elif 50.0 <= pct_good < 75.0:
        return "C"
    elif 25.0 <= pct_good < 75.0:
        return "D"
    elif 0.0 < pct_good < 25.0:
        return "E"
    else:
        return " "


def write_scientific_calib_chla(bgc_file, WMOfloatid, idx_profile, df_bio):
    """Write SCIENTIFIC_CALIB_* variables for CHLA & BBP700 per Argo Cookbook standard for broken sensors."""
    is_broken_chla = is_broken_sensor_float(WMOfloatid, idx_profile)

    if not is_broken_chla:
        cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
        if not cycle_df.empty and "CHLA_FINAL_QC" in cycle_df.columns:
            qcs = cycle_df["CHLA_FINAL_QC"].dropna().astype(int)
            if not qcs.empty and (qcs.isin([4, 9])).all():
                is_broken_chla = True

    if is_broken_chla:
        comment_chla = scientific_calibration_comment_BROKEN
        comment_flu = scientific_calibration_comment_BROKEN
        comment_bbp = scientific_calibration_comment_BROKEN
        equation_chla = scientific_calibration_equation_BROKEN
        equation_bbp = scientific_calibration_equation_BROKEN
        coef_chla = scientific_calibration_coefficient_BROKEN
        coef_flu = scientific_calibration_coefficient_BROKEN
        coef_bbp = scientific_calibration_coefficient_BROKEN
    else:
        dark_cols = sorted(
            [col for col in df_bio.columns if col.startswith("MIN_FLUOCHLA_CYCLE")]
        )
        dark_vals = (
            df_bio[dark_cols].iloc[0].dropna().astype(int).tolist()
            if dark_cols
            else []
        )
        dark_str = " ".join(map(str, dark_vals)) if dark_vals else "NA"
        scale_val = (
            df_bio["SCALE_CHLA"].iloc[0] if "SCALE_CHLA" in df_bio.columns else "NA"
        )
        physio_val = (
            cycle_df["PHYSIO_RATIO"].iloc[0]
            if (not cycle_df.empty and "PHYSIO_RATIO" in cycle_df.columns)
            else "1"
        )

        comment_chla = scientific_calibration_comment_CHLA
        comment_flu = scientific_calibration_comment_CHLA_FLU
        comment_bbp = "BBP700 spike test and regional transformation"
        equation_chla = scientific_calibration_equation_CHLA
        equation_bbp = "BBP700_ADJUSTED = (BBP700 - DARK_BBP700) * SCALE_BBP700"
        coef_chla = (
            f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val},"
            f" PHYSIO_RATIO = {physio_val}"
        )
        coef_flu = f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}"
        coef_bbp = "Not applicable"

    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    str256_len = bgc_file.dimensions["STRING256"].size

    SciCalComArray_CHLA = nc.stringtochar(
        np.array(comment_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalComArray_CHLA_FLU = nc.stringtochar(
        np.array(comment_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalComArray_BBP = nc.stringtochar(
        np.array(comment_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalEquArray_CHLA = nc.stringtochar(
        np.array(
            equation_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"
        )
    )
    SciCalEquArray_BBP = nc.stringtochar(
        np.array(
            equation_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"
        )
    )
    SciCalCoeArray_CHLA = nc.stringtochar(
        np.array(coef_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalCoeArray_CHLA_FLU = nc.stringtochar(
        np.array(coef_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalCoeArray_BBP = nc.stringtochar(
        np.array(coef_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}")
    )
    SciCalDateArray = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))

    n_prof_size = bgc_file.dimensions["N_PROF"].size
    n_param_size = bgc_file.dimensions["N_PARAM"].size

    for iprof_idx in range(n_prof_size):
        station_params = bgc_file.variables["STATION_PARAMETERS"][iprof_idx]
        if isinstance(station_params, np.ma.MaskedArray):
            station_params = station_params.filled(b" ")

        for j in range(n_param_size):
            p_bytes = station_params[j]
            param_str = "".join([
                (
                    c.decode("utf-8", errors="ignore")
                    if isinstance(c, (bytes, np.bytes_))
                    else str(c)
                )
                for c in p_bytes
            ]).strip()

            if param_str == "CHLA_FLUORESCENCE":
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = (
                    SciCalComArray_CHLA_FLU
                )
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][
                    iprof_idx, 0, j, :
                ] = SciCalCoeArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = (
                    SciCalDateArray
                )
            elif param_str == "CHLA":
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
            elif param_str == "BBP700":
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = (
                    SciCalComArray_BBP
                )
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = (
                    SciCalEquArray_BBP
                )
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][
                    iprof_idx, 0, j, :
                ] = SciCalCoeArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = (
                    SciCalDateArray
                )


def write_chla_BBP_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio, iprof_idx=0):
    """Populate CHLA_ADJUSTED and CHLA_FLUORESCENCE_ADJUSTED across ALL N_PROF."""
    n_prof = bgc_file.dimensions["N_PROF"].size
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    for prof in range(n_prof):
        CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_Array[:] = 99999.0

        CHLA_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9", dtype="|S1"
        )

        has_chla = "CHLA" in bgc_file.variables
        has_chla_qc = "CHLA_QC" in bgc_file.variables
        raw_chla_var = (
            bgc_file.variables["CHLA"][prof, :]
            if (has_chla and bgc_file.variables["CHLA"].ndim > 1)
            else (bgc_file.variables["CHLA"][:] if has_chla else None)
        )
        raw_chla_qc_var = (
            bgc_file.variables["CHLA_QC"][prof, :]
            if (has_chla_qc and bgc_file.variables["CHLA_QC"].ndim > 1)
            else (bgc_file.variables["CHLA_QC"][:] if has_chla_qc else None)
        )

        chla_data_arr = (
            np.array(raw_chla_var, copy=True) if raw_chla_var is not None else None
        )
        chla_qc_arr = (
            np.array(raw_chla_qc_var, copy=True)
            if raw_chla_qc_var is not None
            else None
        )

        CHLA_Adjusted_ERROR_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        CHLA_Adjusted_ERROR_Array[:] = 99999.0

        CHLA_FLUORESCENCE_Adjusted_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        CHLA_FLUORESCENCE_Adjusted_Array[:] = 99999.0

        CHLA_FLUORESCENCE_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9", dtype="|S1"
        )

        has_fluo = "CHLA_FLUORESCENCE" in bgc_file.variables
        has_fluo_qc = "CHLA_FLUORESCENCE_QC" in bgc_file.variables
        raw_fluo_var = (
            bgc_file.variables["CHLA_FLUORESCENCE"][prof, :]
            if (has_fluo and bgc_file.variables["CHLA_FLUORESCENCE"].ndim > 1)
            else (bgc_file.variables["CHLA_FLUORESCENCE"][:] if has_fluo else None)
        )
        raw_fluo_qc_var = (
            bgc_file.variables["CHLA_FLUORESCENCE_QC"][prof, :]
            if (has_fluo_qc and bgc_file.variables["CHLA_FLUORESCENCE_QC"].ndim > 1)
            else (
                bgc_file.variables["CHLA_FLUORESCENCE_QC"][:]
                if has_fluo_qc
                else None
            )
        )

        fluo_data_arr = (
            np.array(raw_fluo_var, copy=True) if raw_fluo_var is not None else None
        )
        fluo_qc_arr = (
            np.array(raw_fluo_qc_var, copy=True)
            if raw_fluo_qc_var is not None
            else None
        )

        CHLA_FLUORESCENCE_Adjusted_ERROR_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        CHLA_FLUORESCENCE_Adjusted_ERROR_Array[:] = 99999.0

        nc_pres_all = np.float32(
            bgc_file.variables["PRES"][prof, :]
            if bgc_file.variables["PRES"].ndim > 1
            else bgc_file.variables["PRES"][:]
        )
        has_pres_qc = "PRES_QC" in bgc_file.variables
        raw_pres_qc_var = (
            bgc_file.variables["PRES_QC"][prof, :]
            if (has_pres_qc and bgc_file.variables["PRES_QC"].ndim > 1)
            else (bgc_file.variables["PRES_QC"][:] if has_pres_qc else None)
        )

        for i in range(n_levels):
            nc_pres = nc_pres_all[i]
            pres_qc_char = (
                raw_pres_qc_var[i].decode("utf-8")
                if (raw_pres_qc_var is not None and isinstance(raw_pres_qc_var[i], bytes))
                else (
                    str(raw_pres_qc_var[i])
                    if raw_pres_qc_var is not None
                    else "9"
                )
            )

            is_pres_missing_or_bad = (
                np.isnan(nc_pres) or nc_pres == 99999.0 or pres_qc_char in ["4", "9"]
            )

            is_raw_chla_missing = (
                chla_data_arr is None
                or pd.isna(chla_data_arr[i])
                or chla_data_arr[i] == 99999.0
            )

            is_raw_fluo_missing = (
                fluo_data_arr is None
                or pd.isna(fluo_data_arr[i])
                or fluo_data_arr[i] == 99999.0
            )

            matched_row = None
            if not cycle_df.empty and not is_broken and not is_pres_missing_or_bad:
                diffs = np.abs(cycle_df["PRES"].values - nc_pres)
                min_idx = np.argmin(diffs)
                if diffs[min_idx] <= 0.5:
                    matched_row = cycle_df.iloc[min_idx]

            raw_chla_qc_char = (
                chla_qc_arr[i].decode("utf-8")
                if (chla_qc_arr is not None and isinstance(chla_qc_arr[i], bytes))
                else str(chla_qc_arr[i]) if chla_qc_arr is not None else "9"
            )
            raw_fluo_qc_char = (
                fluo_qc_arr[i].decode("utf-8")
                if (fluo_qc_arr is not None and isinstance(fluo_qc_arr[i], bytes))
                else str(fluo_qc_arr[i]) if fluo_qc_arr is not None else "9"
            )

            # Strict Missingness Alignment for CHLA_ADJUSTED
            if is_pres_missing_or_bad or is_raw_chla_missing or is_broken:
                CHLA_Adjusted_Array[i] = 99999.0
                CHLA_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_Adjusted_Array.mask[i] = True
                CHLA_Adjusted_ERROR_Array.mask[i] = True
                target_pres_qc = (
                    b"4" if (is_broken or pres_qc_char == "4") else b"9"
                )
                CHLA_AdjustedQC_Array[i] = target_pres_qc

                # TWO-WAY SYNCHRONIZATION BACK TO RAW CHLA AND CHLA_QC
                if chla_data_arr is not None:
                    chla_data_arr[i] = 99999.0
                if chla_qc_arr is not None:
                    chla_qc_arr[i] = target_pres_qc
            elif matched_row is not None:
                csv_chla_qc = (
                    matched_row["CHLA_FINAL_QC"]
                    if "CHLA_FINAL_QC" in matched_row
                    else "9"
                )
                chla_qc_str = str(int(csv_chla_qc)) if pd.notna(csv_chla_qc) else "9"
                chla_final_val = (
                    matched_row["CHLA_FINAL"] if "CHLA_FINAL" in matched_row else np.nan
                )

                if (
                    raw_chla_qc_char in ["4", "9"]
                    or chla_qc_str in ["4", "9"]
                    or pd.isna(chla_final_val)
                    or chla_final_val == 99999.0
                ):
                    CHLA_Adjusted_Array[i] = 99999.0
                    CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    CHLA_Adjusted_Array.mask[i] = True
                    CHLA_Adjusted_ERROR_Array.mask[i] = True
                    target_flag = (
                        raw_chla_qc_char
                        if raw_chla_qc_char in ["4", "9"]
                        else (chla_qc_str if chla_qc_str in ["4", "9"] else "9")
                    )
                    CHLA_AdjustedQC_Array[i] = target_flag.encode("utf-8")

                    # TWO-WAY SYNCHRONIZATION BACK TO RAW CHLA AND CHLA_QC
                    if chla_data_arr is not None:
                        chla_data_arr[i] = 99999.0
                    if chla_qc_arr is not None:
                        chla_qc_arr[i] = target_flag.encode("utf-8")
                else:
                    CHLA_Adjusted_Array[i] = np.float32(chla_final_val)
                    CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                    CHLA_Adjusted_Array.mask[i] = False
                    CHLA_Adjusted_ERROR_Array.mask[i] = False
                    CHLA_AdjustedQC_Array[i] = (
                        chla_qc_str.encode("utf-8")
                        if chla_qc_str in ["1", "2", "3"]
                        else b"1"
                    )
            else:
                CHLA_Adjusted_Array[i] = 99999.0
                CHLA_Adjusted_Array.mask[i] = True
                CHLA_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_Adjusted_ERROR_Array.mask[i] = True
                CHLA_AdjustedQC_Array[i] = b"9"
                if chla_data_arr is not None:
                    chla_data_arr[i] = 99999.0
                if chla_qc_arr is not None:
                    chla_qc_arr[i] = b"9"

            # Strict Missingness Alignment for FLUORESCENCE_ADJUSTED
            if is_pres_missing_or_bad or is_raw_fluo_missing or is_broken:
                CHLA_FLUORESCENCE_Adjusted_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = True
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = True
                target_pres_qc = (
                    b"4" if (is_broken or pres_qc_char == "4") else b"9"
                )
                CHLA_FLUORESCENCE_AdjustedQC_Array[i] = target_pres_qc

                # TWO-WAY SYNCHRONIZATION BACK TO RAW CHLA_FLUORESCENCE AND CHLA_FLUORESCENCE_QC
                if fluo_data_arr is not None:
                    fluo_data_arr[i] = 99999.0
                if fluo_qc_arr is not None:
                    fluo_qc_arr[i] = target_pres_qc
            elif matched_row is not None:
                fluo_adj_qc_str = (
                    str(int(matched_row["CHLA_FLUORESCENCE_ADJUSTED_QC"]))
                    if (
                        "CHLA_FLUORESCENCE_ADJUSTED_QC" in matched_row
                        and pd.notna(matched_row["CHLA_FLUORESCENCE_ADJUSTED_QC"])
                    )
                    else raw_chla_qc_char
                )
                fluo_val = (
                    matched_row["CHLA_FLUORESCENCE"]
                    if "CHLA_FLUORESCENCE" in matched_row
                    else np.nan
                )

                if (
                    raw_fluo_qc_char in ["4", "9"]
                    or fluo_adj_qc_str in ["4", "9"]
                    or pd.isna(fluo_val)
                    or fluo_val == 99999.0
                ):
                    CHLA_FLUORESCENCE_Adjusted_Array[i] = 99999.0
                    CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = 99999.0
                    CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = True
                    CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = True
                    target_fluo_flag = (
                        raw_fluo_qc_char
                        if raw_fluo_qc_char in ["4", "9"]
                        else (
                            fluo_adj_qc_str if fluo_adj_qc_str in ["4", "9"] else "9"
                        )
                    )
                    CHLA_FLUORESCENCE_AdjustedQC_Array[i] = target_fluo_flag.encode("utf-8")

                    # TWO-WAY SYNCHRONIZATION BACK TO RAW CHLA_FLUORESCENCE AND CHLA_FLUORESCENCE_QC
                    if fluo_data_arr is not None:
                        fluo_data_arr[i] = 99999.0
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = target_fluo_flag.encode("utf-8")
                else:
                    CHLA_FLUORESCENCE_Adjusted_Array[i] = np.float32(fluo_val)
                    CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = np.float32(
                        CHLA_Adjusted_ERROR_est
                    )
                    CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = False
                    CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = False
                    CHLA_FLUORESCENCE_AdjustedQC_Array[i] = (
                        fluo_adj_qc_str.encode("utf-8")
                        if fluo_adj_qc_str in ["1", "2", "3"]
                        else b"1"
                    )
            else:
                CHLA_FLUORESCENCE_Adjusted_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_Array.mask[i] = True
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i] = 99999.0
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask[i] = True
                CHLA_FLUORESCENCE_AdjustedQC_Array[i] = b"9"
                if fluo_data_arr is not None:
                    fluo_data_arr[i] = 99999.0
                if fluo_qc_arr is not None:
                    fluo_qc_arr[i] = b"9"

        # Explicit 2D Array Slicing Assignment
        bgc_file.variables["CHLA_ADJUSTED"][prof, :] = CHLA_Adjusted_Array
        bgc_file.variables["CHLA_ADJUSTED_QC"][prof, :] = CHLA_AdjustedQC_Array
        bgc_file.variables["CHLA_ADJUSTED_ERROR"][prof, :] = CHLA_Adjusted_ERROR_Array

        if has_chla:
            bgc_file.variables["CHLA"][prof, :] = chla_data_arr
        if has_chla_qc:
            bgc_file.variables["CHLA_QC"][prof, :] = chla_qc_arr

        if "CHLA_FLUORESCENCE_ADJUSTED" in bgc_file.variables:
            bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED"][prof, :] = (
                CHLA_FLUORESCENCE_Adjusted_Array
            )
            bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_QC"][prof, :] = (
                CHLA_FLUORESCENCE_AdjustedQC_Array
            )
            bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_ERROR"][prof, :] = (
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array
            )

        if has_fluo:
            bgc_file.variables["CHLA_FLUORESCENCE"][prof, :] = fluo_data_arr
        if has_fluo_qc:
            bgc_file.variables["CHLA_FLUORESCENCE_QC"][prof, :] = fluo_qc_arr

        # Profile QC Grade Assignment
        prof_qc = get_profile_qc_grade(CHLA_AdjustedQC_Array)
        prof_fluo_qc = get_profile_qc_grade(CHLA_FLUORESCENCE_AdjustedQC_Array)

        if "PROFILE_CHLA_QC" in bgc_file.variables:
            bgc_file.variables["PROFILE_CHLA_QC"][prof] = np.array(
                [prof_qc], dtype="|S1"
            )

        if "PROFILE_CHLA_FLUORESCENCE_QC" in bgc_file.variables:
            bgc_file.variables["PROFILE_CHLA_FLUORESCENCE_QC"][prof] = np.array(
                [prof_fluo_qc], dtype="|S1"
            )


def write_BBP700_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio, iprof_idx=0):
    """Populate BBP700_ADJUSTED enforcing strict raw missingness alignment across ALL N_PROF."""
    if "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    n_prof = bgc_file.dimensions["N_PROF"].size
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    for prof in range(n_prof):
        has_bbp = "BBP700" in bgc_file.variables
        has_bbp_qc = "BBP700_QC" in bgc_file.variables

        raw_bbp_var = (
            bgc_file.variables["BBP700"][prof, :]
            if (has_bbp and bgc_file.variables["BBP700"].ndim > 1)
            else (bgc_file.variables["BBP700"][:] if has_bbp else None)
        )
        raw_bbp_qc_var = (
            bgc_file.variables["BBP700_QC"][prof, :]
            if (has_bbp_qc and bgc_file.variables["BBP700_QC"].ndim > 1)
            else (bgc_file.variables["BBP700_QC"][:] if has_bbp_qc else None)
        )

        bbp_data_arr = (
            np.array(raw_bbp_var, copy=True) if raw_bbp_var is not None else None
        )
        bbp_qc_arr = (
            np.array(raw_bbp_qc_var, copy=True) if raw_bbp_qc_var is not None else None
        )

        BBP700_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        BBP700_Adjusted_Array[:] = 99999.0

        BBP700_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9", dtype="|S1"
        )

        BBP700_Adjusted_ERROR_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        BBP700_Adjusted_ERROR_Array[:] = 99999.0

        nc_pres_all = np.float32(
            bgc_file.variables["PRES"][prof, :]
            if bgc_file.variables["PRES"].ndim > 1
            else bgc_file.variables["PRES"][:]
        )
        has_pres_qc = "PRES_QC" in bgc_file.variables
        raw_pres_qc_var = (
            bgc_file.variables["PRES_QC"][prof, :]
            if (has_pres_qc and bgc_file.variables["PRES_QC"].ndim > 1)
            else (bgc_file.variables["PRES_QC"][:] if has_pres_qc else None)
        )

        for i in range(n_levels):
            nc_pres = nc_pres_all[i]
            pres_qc_char = (
                raw_pres_qc_var[i].decode("utf-8")
                if (raw_pres_qc_var is not None and isinstance(raw_pres_qc_var[i], bytes))
                else (
                    str(raw_pres_qc_var[i])
                    if raw_pres_qc_var is not None
                    else "9"
                )
            )

            is_pres_missing_or_bad = (
                np.isnan(nc_pres) or nc_pres == 99999.0 or pres_qc_char in ["4", "9"]
            )

            is_raw_bbp_missing = (
                bbp_data_arr is None
                or pd.isna(bbp_data_arr[i])
                or bbp_data_arr[i] == 99999.0
            )

            matched_row = None
            if not cycle_df.empty and not is_broken and not is_pres_missing_or_bad:
                diffs = np.abs(cycle_df["PRES"].values - nc_pres)
                min_idx = np.argmin(diffs)
                if diffs[min_idx] <= 0.5:
                    matched_row = cycle_df.iloc[min_idx]

            raw_qc_char = (
                bbp_qc_arr[i].decode("utf-8")
                if (bbp_qc_arr is not None and isinstance(bbp_qc_arr[i], bytes))
                else str(bbp_qc_arr[i]) if bbp_qc_arr is not None else "9"
            )

            if is_pres_missing_or_bad or is_raw_bbp_missing or is_broken:
                BBP700_Adjusted_Array[i] = 99999.0
                BBP700_Adjusted_ERROR_Array[i] = 99999.0
                BBP700_Adjusted_Array.mask[i] = True
                BBP700_Adjusted_ERROR_Array.mask[i] = True
                target_pres_qc = (
                    b"4" if (is_broken or pres_qc_char == "4") else b"9"
                )
                BBP700_AdjustedQC_Array[i] = target_pres_qc

                # TWO-WAY SYNCHRONIZATION BACK TO RAW BBP700 AND BBP700_QC
                if bbp_data_arr is not None:
                    bbp_data_arr[i] = 99999.0
                if bbp_qc_arr is not None:
                    bbp_qc_arr[i] = target_pres_qc
            elif matched_row is not None:
                csv_qc = (
                    matched_row["BBP700_FINAL_QC"]
                    if "BBP700_FINAL_QC" in matched_row
                    else "9"
                )
                qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "9"
                raw_bbp = (
                    matched_row["BBP700_FINAL"]
                    if "BBP700_FINAL" in matched_row
                    else np.nan
                )

                if (
                    raw_qc_char in ["4", "9"]
                    or qc_str in ["4", "9"]
                    or pd.isna(raw_bbp)
                    or raw_bbp == 99999.0
                ):
                    BBP700_Adjusted_Array[i] = 99999.0
                    BBP700_Adjusted_ERROR_Array[i] = 99999.0
                    BBP700_Adjusted_Array.mask[i] = True
                    BBP700_Adjusted_ERROR_Array.mask[i] = True
                    target_flag = (
                        raw_qc_char
                        if raw_qc_char in ["4", "9"]
                        else (qc_str if qc_str in ["4", "9"] else "9")
                    )
                    BBP700_AdjustedQC_Array[i] = target_flag.encode("utf-8")

                    # TWO-WAY SYNCHRONIZATION BACK TO RAW BBP700 AND BBP700_QC
                    if bbp_data_arr is not None:
                        bbp_data_arr[i] = 99999.0
                    if bbp_qc_arr is not None:
                        bbp_qc_arr[i] = target_flag.encode("utf-8")
                else:
                    BBP700_Adjusted_Array[i] = np.float32(raw_bbp)
                    BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                    BBP700_Adjusted_Array.mask[i] = False
                    BBP700_Adjusted_ERROR_Array.mask[i] = False
                    BBP700_AdjustedQC_Array[i] = (
                        qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1"
                    )

        # Explicit 2D Array Slicing Assignment
        bgc_file.variables["BBP700_ADJUSTED"][prof, :] = BBP700_Adjusted_Array
        if "BBP700_ADJUSTED_QC" in bgc_file.variables:
            bgc_file.variables["BBP700_ADJUSTED_QC"][prof, :] = (
                BBP700_AdjustedQC_Array
            )
        if "BBP700_ADJUSTED_ERROR" in bgc_file.variables:
            bgc_file.variables["BBP700_ADJUSTED_ERROR"][prof, :] = (
                BBP700_Adjusted_ERROR_Array
            )

        if has_bbp:
            bgc_file.variables["BBP700"][prof, :] = bbp_data_arr
        if has_bbp_qc:
            bgc_file.variables["BBP700_QC"][prof, :] = bbp_qc_arr

        prof_bbp_qc = get_profile_qc_grade(BBP700_AdjustedQC_Array)
        if "PROFILE_BBP700_QC" in bgc_file.variables:
            bgc_file.variables["PROFILE_BBP700_QC"][prof] = np.array(
                [prof_bbp_qc], dtype="|S1"
            )


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle):
    """Write DOXY slope and drift calibration coefficients into SCIENTIFIC_CALIB_* per Argo DOXY Cookbook standards."""
    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)]
    if cycle_df.empty:
        return

    slope_val = (
        cycle_df["DOXY_SLOPE"].iloc[0]
        if "DOXY_SLOPE" in cycle_df.columns
        else np.nan
    )
    drift_val = (
        cycle_df["DOXY_DRIFT"].iloc[0]
        if "DOXY_DRIFT" in cycle_df.columns
        else np.nan
    )

    s_str = "1.0" if pd.isna(slope_val) else str(slope_val)
    d_str = "0.0" if pd.isna(drift_val) else str(drift_val)

    calib_coef_str = f"slope = {s_str}, drift = {d_str}"
    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

    str256_len = (
        ds.dimensions["STRING256"].size if "STRING256" in ds.dimensions else 256
    )

    char_coef = nc.stringtochar(
        np.array(
            calib_coef_str.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"
        )
    )
    char_equ = nc.stringtochar(
        np.array(
            scientific_calibration_equation_DOXY.ljust(str256_len)[:str256_len],
            dtype=f"S{str256_len}",
        )
    )
    char_com = nc.stringtochar(
        np.array(
            scientific_calibration_comment_DOXY.ljust(str256_len)[:str256_len],
            dtype=f"S{str256_len}",
        )
    )
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
                (
                    c.decode("utf-8", errors="ignore")
                    if isinstance(c, (bytes, np.bytes_))
                    else str(c)
                )
                for c in p_bytes
            ]).strip()

            if param_str == "DOXY":
                if "SCIENTIFIC_CALIB_COEFFICIENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof, 0, j, :] = (
                        char_coef
                    )
                if "SCIENTIFIC_CALIB_EQUATION" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_EQUATION"][iprof, 0, j, :] = char_equ
                if "SCIENTIFIC_CALIB_COMMENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COMMENT"][iprof, 0, j, :] = char_com
                if "SCIENTIFIC_CALIB_DATE" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_DATE"][iprof, 0, j, :] = char_date


def write_DOXY_from_csv(ds, float_df, target_cycle):
    """Populate DOXY_ADJUSTED ensuring DOXY_ADJUSTED_QC strictly matches raw DOXY_QC across ALL profiles."""
    var_names = ds.variables.keys()
    n_prof = ds.dimensions["N_PROF"].size
    n_levels = ds.dimensions["N_LEVELS"].size

    pres_nc_full = ds.variables["PRES"][:]
    has_pres_qc = "PRES_QC" in ds.variables
    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()

    for iprof in range(n_prof):
        if pres_nc_full.ndim > 1:
            pres_nc = (
                pres_nc_full[iprof, :]
                if pres_nc_full.shape[0] == n_prof
                else pres_nc_full[:, iprof]
            )
        else:
            pres_nc = pres_nc_full

        raw_pres_qc_var = (
            ds.variables["PRES_QC"][iprof, :]
            if (has_pres_qc and ds.variables["PRES_QC"].ndim > 1)
            else (ds.variables["PRES_QC"][:] if has_pres_qc else None)
        )

        DOXY_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        DOXY_Adjusted_Array[:] = 99999.0

        DOXY_Adjusted_Error_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        DOXY_Adjusted_Error_Array[:] = 99999.0

        DOXY_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9", dtype="|S1"
        )

        has_doxy = "DOXY" in ds.variables
        has_doxy_qc = "DOXY_QC" in ds.variables

        raw_doxy_var = (
            ds.variables["DOXY"][iprof, :]
            if (has_doxy and ds.variables["DOXY"].ndim > 1)
            else (ds.variables["DOXY"][:] if has_doxy else None)
        )
        raw_doxy_qc_var = (
            ds.variables["DOXY_QC"][iprof, :]
            if (has_doxy_qc and ds.variables["DOXY_QC"].ndim > 1)
            else (ds.variables["DOXY_QC"][:] if has_doxy_qc else None)
        )

        doxy_data_arr = (
            np.array(raw_doxy_var, copy=True) if raw_doxy_var is not None else None
        )
        doxy_qc_arr = (
            np.array(raw_doxy_qc_var, copy=True)
            if raw_doxy_qc_var is not None
            else None
        )

        for i in range(n_levels):
            nc_pres = np.float32(pres_nc[i])
            pres_qc_char = (
                raw_pres_qc_var[i].decode("utf-8")
                if (raw_pres_qc_var is not None and isinstance(raw_pres_qc_var[i], bytes))
                else (
                    str(raw_pres_qc_var[i])
                    if raw_pres_qc_var is not None
                    else "9"
                )
            )

            # PRES DOMINANCE RULE
            is_pres_missing_or_bad = (
                np.isnan(nc_pres) or nc_pres == 99999.0 or pres_qc_char in ["4", "9"]
            )

            has_valid_raw_data = (
                doxy_data_arr is not None
                and not pd.isna(doxy_data_arr[i])
                and doxy_data_arr[i] != 99999.0
            )

            raw_qc_char = (
                doxy_qc_arr[i].decode("utf-8")
                if (doxy_qc_arr is not None and isinstance(doxy_qc_arr[i], bytes))
                else str(doxy_qc_arr[i]) if doxy_qc_arr is not None else "9"
            )

            if is_pres_missing_or_bad or not has_valid_raw_data:
                DOXY_Adjusted_Array[i] = 99999.0
                DOXY_Adjusted_Array.mask[i] = True
                DOXY_Adjusted_Error_Array[i] = 99999.0
                DOXY_Adjusted_Error_Array.mask[i] = True
                target_pres_qc = (
                    b"4" if pres_qc_char == "4" else b"9"
                )
                DOXY_AdjustedQC_Array[i] = target_pres_qc
                # TWO-WAY SYNCHRONIZATION BACK TO RAW DOXY_QC
                if doxy_qc_arr is not None:
                    doxy_qc_arr[i] = target_pres_qc
                continue

            matched_row = None
            if not cycle_df.empty:
                diffs = np.abs(cycle_df["PRES"].values - nc_pres)
                min_idx = np.argmin(diffs)
                if diffs[min_idx] <= 0.5:
                    matched_row = cycle_df.iloc[min_idx]

            if matched_row is not None:
                csv_qc = matched_row.get("DOXY_FINAL_QC")
                qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "9"

                raw_doxy_final = matched_row.get("DOXY_FINAL")

                if (
                    raw_qc_char in ["4", "9"]
                    or qc_str in ["4", "9"]
                    or pd.isna(raw_doxy_final)
                    or raw_doxy_final == 99999.0
                ):
                    DOXY_Adjusted_Array[i] = 99999.0
                    DOXY_Adjusted_Array.mask[i] = True
                    DOXY_Adjusted_Error_Array[i] = 99999.0
                    DOXY_Adjusted_Error_Array.mask[i] = True
                    target_flag = (
                        raw_qc_char
                        if raw_qc_char in ["4", "9"]
                        else (qc_str if qc_str in ["4", "9"] else "9")
                    )
                    DOXY_AdjustedQC_Array[i] = target_flag.encode("utf-8")
                    # TWO-WAY SYNCHRONIZATION BACK TO RAW DOXY_QC
                    if doxy_qc_arr is not None:
                        doxy_qc_arr[i] = target_flag.encode("utf-8")
                else:
                    DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
                    DOXY_Adjusted_Array.mask[i] = False
                    # Carry DOXY_QC flag directly over to DOXY_ADJUSTED_QC
                    DOXY_AdjustedQC_Array[i] = (
                        raw_qc_char.encode("utf-8")
                        if raw_qc_char in ["1", "2", "3"]
                        else (
                            qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1"
                        )
                    )

                    raw_doxy_error = matched_row.get("DOXY_ADJUSTED_ERROR")
                    if pd.notna(raw_doxy_error) and raw_doxy_error != 99999.0:
                        DOXY_Adjusted_Error_Array[i] = np.float32(raw_doxy_error)
                        DOXY_Adjusted_Error_Array.mask[i] = False
                    else:
                        DOXY_Adjusted_Error_Array[i] = 99999.0
                        DOXY_Adjusted_Error_Array.mask[i] = True
            else:
                DOXY_Adjusted_Array[i] = 99999.0
                DOXY_Adjusted_Array.mask[i] = True
                DOXY_Adjusted_Error_Array[i] = 99999.0
                DOXY_Adjusted_Error_Array.mask[i] = True
                DOXY_AdjustedQC_Array[i] = b"9"
                if doxy_qc_arr is not None:
                    doxy_qc_arr[i] = b"9"

        if "DOXY_QC" in var_names and doxy_qc_arr is not None:
            if ds.variables["DOXY_QC"].ndim > 1:
                ds.variables["DOXY_QC"][iprof, :] = doxy_qc_arr
            else:
                ds.variables["DOXY_QC"][:] = doxy_qc_arr

        if "DOXY_ADJUSTED" in var_names:
            doxy_adj_var = ds.variables["DOXY_ADJUSTED"]
            if doxy_adj_var.ndim > 1:
                if doxy_adj_var.shape[0] == n_prof:
                    doxy_adj_var[iprof, :] = DOXY_Adjusted_Array
                    ds.variables["DOXY_ADJUSTED_QC"][iprof, :] = DOXY_AdjustedQC_Array
                else:
                    doxy_adj_var[:, iprof] = DOXY_Adjusted_Array
                    ds.variables["DOXY_ADJUSTED_QC"][:, iprof] = DOXY_AdjustedQC_Array
            else:
                doxy_adj_var[:] = DOXY_Adjusted_Array
                ds.variables["DOXY_ADJUSTED_QC"][:] = DOXY_AdjustedQC_Array

        if "DOXY_ADJUSTED_ERROR" in var_names:
            doxy_adj_err_var = ds.variables["DOXY_ADJUSTED_ERROR"]
            if doxy_adj_err_var.ndim > 1:
                if doxy_adj_err_var.shape[0] == n_prof:
                    doxy_adj_err_var[iprof, :] = DOXY_Adjusted_Error_Array
                else:
                    doxy_adj_err_var[:, iprof] = DOXY_Adjusted_Error_Array
            else:
                doxy_adj_err_var[:] = DOXY_Adjusted_Error_Array

        prof_qc = get_profile_qc_grade(DOXY_AdjustedQC_Array)
        if "PROFILE_DOXY_QC" in ds.variables:
            prof_qc_var = ds.variables["PROFILE_DOXY_QC"]
            qc_char = np.array([prof_qc], dtype="|S1")
            if prof_qc_var.ndim == 1:
                prof_qc_var[iprof] = qc_char
            elif prof_qc_var.ndim == 2:
                prof_qc_var[iprof, 0] = qc_char


def safe_rename(from_file, to_file):
    """Safely rename file handling Windows/Linux file locks."""
    gc.collect()
    try:
        os.replace(from_file, to_file)
    except OSError:
        try:
            shutil.copyfile(from_file, to_file)
            os.remove(from_file)
        except Exception:
            pass


def extract_parameter_data_mode(ds, target_param, iprof):
    """Helper to cleanly resolve parameter data mode string for a given parameter name and profile index."""
    if (
        "STATION_PARAMETERS" not in ds.variables
        or "PARAMETER_DATA_MODE" not in ds.variables
    ):
        return "UNKNOWN"

    station_params_mat = ds.variables["STATION_PARAMETERS"][iprof]
    if isinstance(station_params_mat, np.ma.MaskedArray):
        station_params_mat = station_params_mat.filled(b" ")

    pdm_raw = ds.variables["PARAMETER_DATA_MODE"][iprof]

    for j in range(ds.dimensions["N_PARAM"].size):
        param_chars = station_params_mat[j]
        p_str = "".join([
            (
                c.decode("utf-8", errors="ignore")
                if isinstance(c, (bytes, np.bytes_))
                else str(c)
            )
            for c in param_chars
        ]).strip()

        if p_str == target_param:
            val = pdm_raw[j]
            if isinstance(val, (bytes, np.bytes_)):
                return (
                    val.tobytes().decode("utf-8").strip()
                    if hasattr(val, "tobytes")
                    else val.decode("utf-8").strip()
                )
            return str(val).strip()

    return "NOT_FOUND"


def inspect_output_files(created_files):
    """Automated inspection function to verify GDAC compliance of output NetCDF files."""
    print("\n" + "=" * 80)
    print("      AUTOMATED OUTPUT INSPECTION & GDAC COMPLIANCE VERIFICATION")
    print("=" * 80)

    total_files_checked = 0
    total_violations_found = 0

    for nc_file_path in created_files:
        if not os.path.exists(nc_file_path):
            continue

        total_files_checked += 1
        file_name = os.path.basename(nc_file_path)
        file_violations = 0

        try:
            ds = nc.Dataset(nc_file_path, "r")
            n_prof = ds.dimensions["N_PROF"].size if "N_PROF" in ds.dimensions else 1

            for prof in range(n_prof):
                for param in ["CHLA", "CHLA_FLUORESCENCE", "BBP700", "DOXY"]:
                    if (
                        param not in ds.variables
                        or f"{param}_ADJUSTED" not in ds.variables
                    ):
                        continue

                    raw_var = ds.variables[param]
                    adj_var = ds.variables[f"{param}_ADJUSTED"]
                    qc_raw_var = ds.variables[f"{param}_QC"]
                    qc_adj_var = ds.variables[f"{param}_ADJUSTED_QC"]
                    err_adj_var = (
                        ds.variables[f"{param}_ADJUSTED_ERROR"]
                        if f"{param}_ADJUSTED_ERROR" in ds.variables
                        else None
                    )

                    raw_data = (
                        raw_var[prof, :].filled(99999.0)
                        if raw_var.ndim > 1
                        else raw_var[:].filled(99999.0)
                    )
                    adj_data = (
                        adj_var[prof, :].filled(99999.0)
                        if adj_var.ndim > 1
                        else adj_var[:].filled(99999.0)
                    )
                    qc_raw = (
                        qc_raw_var[prof, :].astype(str)
                        if qc_raw_var.ndim > 1
                        else qc_raw_var[:].astype(str)
                    )
                    qc_adj = (
                        qc_adj_var[prof, :].astype(str)
                        if qc_adj_var.ndim > 1
                        else qc_adj_var[:].astype(str)
                    )
                    err_data = (
                        err_adj_var[prof, :].filled(99999.0)
                        if (err_adj_var is not None and err_adj_var.ndim > 1)
                        else (
                            err_adj_var[:].filled(99999.0)
                            if err_adj_var is not None
                            else None
                        )
                    )

                    raw_missing_mask = (raw_data == 99999.0) | pd.isna(raw_data)

                    # Check 1: Raw missing but Adjusted not missing
                    viol_data = int(np.sum(raw_missing_mask & (adj_data != 99999.0)))
                    # Check 2: Raw missing but Adjusted QC not '9' or '4'
                    viol_qc = int(
                        np.sum(raw_missing_mask & (~np.isin(qc_adj, ["9", "4"])))
                    )
                    # Check 3: Raw missing but Adjusted Error not missing
                    viol_err = (
                        int(np.sum(raw_missing_mask & (err_data != 99999.0)))
                        if err_data is not None
                        else 0
                    )
                    # Check 4: Adjusted QC is missing ('9') but Raw QC is valid ('1'-'4')
                    viol_qc_sync = int(
                        np.sum((qc_adj == "9") & np.isin(qc_raw, ["1", "2", "3", "4"]))
                    )

                    param_viols = viol_data + viol_qc + viol_err + viol_qc_sync
                    if param_viols > 0:
                        file_violations += param_viols
                        print(
                            f"  [FAIL] File: {file_name} | Profile {prof+1} | Parameter:"
                            f" {param} (mismatches: {param_viols})"
                        )

            ds.close()

            if file_violations == 0:
                print(
                    f"  [PASS] File: {file_name} -> 100% GDAC Compliant across all"
                    " profiles."
                )
            else:
                total_violations_found += file_violations

        except Exception as e:
            print(f"  [ERROR] Inspection failed for file {file_name}: {e}")

    print("-" * 80)
    print(
        f"INSPECTION SUMMARY: Checked {total_files_checked} files | Total GDAC"
        f" Violations: {total_violations_found}"
    )
    if total_violations_found == 0 and total_files_checked > 0:
        print(
            "RESULT: SUCCESS - All generated BD files pass FillValue & QC"
            " consistency tests."
        )
    print("=" * 80 + "\n")


# ==============================================================================
# SECTION 5: Efficient Single-Pass Processing Loop
# ==============================================================================

processed_output_paths = []

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
    main_float_dir = f"/data/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
    bio_dmqc_csv_path = (
        f"/data/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_BBP_DOXY_{WMOfloatid}.csv"
    )

    output_lut_dir = os.path.join(main_float_dir, "LUT")
    today_str = dt.now().strftime("%Y-%m-%d")
    final_out_dir = os.path.join(output_lut_dir, today_str)

    if not os.path.exists(bio_dmqc_csv_path) or not os.path.exists(main_float_dir):
        continue

    os.makedirs(output_lut_dir, exist_ok=True)
    os.makedirs(final_out_dir, exist_ok=True)

    df_bio_raw = pd.read_csv(bio_dmqc_csv_path, low_memory=False)
    df_bio = apply_float_override_conditions(df_bio_raw, WMOfloatid)

    # Locate input B files directly in main_float_dir
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

    # Process each B file by creating a working 'w_' file
    for bgc_filename in sorted_b_files:
        idx_profile = int(bgc_filename[-6:-3])

        # 1. READ IN: Create working 'w_' copy in main_float_dir
        w_bgc_filename = create_working_bd_file(
            bgc_filename, dest_dir=main_float_dir
        )

        try:
            ds = nc.Dataset(w_bgc_filename, "r+")

            iprof_chla = detect_parameter_profile(ds, "CHLA")
            iprof_bbp = detect_parameter_profile(ds, "BBP700")
            iprof_doxy = detect_parameter_profile(ds, "DOXY")

            # 1. Update PARAMETER_DATA_MODE and DATA_MODE strictly per Argo standards FIRST
            write_parameter_data_modes(ds)

            # 2. Update CHLA & BBP700 Delayed-Mode Data
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_chla)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_bbp)

            # 3. Update DOXY Delayed-Mode Data
            write_DOXY_slope_drift(ds, 0, df_bio, idx_profile)
            write_DOXY_from_csv(ds, df_bio, idx_profile)

            # 4. Update History Metadata for both CHLA and DOXY
            write_history_metadata(ds, iprof_chla, iprof_doxy)

            # 5. Clean QC flags and valid range attributes across ALL variables (EXCLUDING TEMP_CPU_CHLA)
            clean_and_fill_qc_variables(ds)
            add_missing_valid_range_attributes(ds)
            remove_forbidden_attributes(ds)

            for i in range(ds.dimensions["N_PROF"].size):
                if i in [iprof_chla, iprof_bbp, iprof_doxy]:
                    ds.variables["DATA_STATE_INDICATOR"][i, :] = nc.stringtochar(
                        np.array("2C  ", dtype="S4")
                    )

            ds.close()
            new_bd_files.append(w_bgc_filename)

        except Exception:
            traceback.print_exc()
            if "ds" in locals() and ds.isopen():
                ds.close()

    # 2. WRITE OUT: Strip w_ / AOML_ prefixes and move filled BD files to LUT output dir
    for w_file in new_bd_files:
        path, name = os.path.split(w_file)

        new_name = name.replace("w_AOML_BD", "BD")
        new_name = new_name.replace("w_AOML_BR", "BD")
        new_name = new_name.replace("w_BD", "BD")
        new_name = new_name.replace("w_BR", "BD")

        new_path = os.path.join(final_out_dir, new_name)
        safe_rename(w_file, new_path)
        processed_output_paths.append(new_path)

# ==============================================================================
# SECTION 6: Output File Inspection & Compliance Verification
# ==============================================================================
inspect_output_files(processed_output_paths)