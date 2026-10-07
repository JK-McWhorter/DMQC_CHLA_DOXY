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

WMO_FLOAT_IDS = list(FLOAT_TYPES.keys())

# Common Configuration Defaults
comment_dmqc_operator_chla = (
    "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer"
    " McWhorter, NOAA/AOML"
)
history_parameter_chla = "CHLA"
history_institution = "AO"
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

# Standard error estimates for fallback values
CHLA_Adjusted_ERROR_est = 0.07
BBP700_Adjusted_ERROR_est = 0.0005
DOXY_Adjusted_ERROR_est_default = 10.0

scientific_calibration_equation_DOXY = "DOXY_ADJUSTED = DOXY * slope + drift"
scientific_calibration_comment_DOXY = (
    "DOXY calibration computed against WOA climatology following Bittig et"
    " al. (2018) methodology"
)

scientific_calibration_comment_BROKEN = (
    "Sensor failure / broken sensor - data flagged bad or missing (Argo BGC QC"
    " Manual)"
)
scientific_calibration_equation_BROKEN = "Not applicable"
scientific_calibration_coefficient_BROKEN = "Not applicable"

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
    float_str = str(WMOfloatid)
    cycle_num = int(cycle_num)
    return (
        (float_str == "4903624")
        | ((float_str == "2904010") & (cycle_num >= 49))
        | ((float_str == "2904011") & (cycle_num >= 24))
    )


def apply_float_override_conditions(df, WMOfloatid):
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
        "FLUORESCENCE_CHLA",
        "FLUORESCENCE_CHLA_ADJUSTED",
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
        "FLUORESCENCE_CHLA_QC",
        "FLUORESCENCE_CHLA_ADJUSTED_QC",
        "BBP700_ADJUSTED_QC",
    ]

    for col in float_cols:
        if col in df.columns:
            df.loc[override_condition, col] = 99999.0

    for col in qc_cols:
        if col in df.columns:
            df.loc[override_condition, col] = 4

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


def has_valid_measurements(bgc_file, param, iprof):
    """Check if a parameter has actual numerical measurements at profile iprof."""
    if param not in bgc_file.variables:
        return False
    var = bgc_file.variables[param]
    p_slice = var[iprof, :] if var.ndim > 1 else var[:]
    p_vals = p_slice.filled(99999.0) if hasattr(p_slice, "filled") else np.array(p_slice)
    valid_mask = (~pd.isna(p_vals)) & (p_vals != 99999.0)
    return valid_mask.any()


def write_parameter_data_modes(bgc_file):
    """Set PARAMETER_DATA_MODE and DATA_MODE based on actual parameter presence per profile."""
    file_name = os.path.basename(bgc_file.filepath()) if hasattr(bgc_file, "filepath") else "File"
    n_param = bgc_file.dimensions["N_PARAM"].size
    n_prof = bgc_file.dimensions["N_PROF"].size

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    for iprof in range(n_prof):
        param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        pdm_mapping_log = []

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

            if p_str:
                if p_str == "PRES":
                    pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"
                elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY"]:
                    if has_valid_measurements(bgc_file, p_str, iprof):
                        pdm[iprof, j] = b"D" if pdm.dtype.kind in ["S", "U", "O"] else "D"
                        data_mode[iprof] = b"D" if data_mode.dtype.kind in ["S", "U", "O"] else "D"
                    else:
                        pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"

                m_val = pdm[iprof, j].decode("utf-8") if isinstance(pdm[iprof, j], bytes) else str(pdm[iprof, j])
                pdm_mapping_log.append(f"{p_str}: '{m_val}'")

        dm_str = data_mode[iprof].decode("utf-8") if isinstance(data_mode[iprof], bytes) else str(data_mode[iprof])
        print(f"[{file_name}] Profile [{iprof}] DATA_MODE: '{dm_str}' | PARAMETER_DATA_MODE => {', '.join(pdm_mapping_log)}")

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def get_profile_qc_grade(qc_masked_array, parameter_present=True):
    """
    Compute overall profile QC grade letter ('A', 'B', 'C', 'D', 'E', 'F', or ' ').
    """
    if not parameter_present:
        return " "

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
        return "F"

    valid_qcs = [q for q in raw_qcs if q not in ["4", "9"]]
    total_points = len(valid_qcs)

    if total_points == 0:
        return "F"

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
        return "F"


def get_param_data_mode(bgc_file, iprof, param_base):
    """Get the specific parameter data mode ('R', 'A', 'D') for profile iprof and parameter param_base."""
    if "STATION_PARAMETERS" not in bgc_file.variables or "PARAMETER_DATA_MODE" not in bgc_file.variables:
        dm_var = bgc_file.variables.get("DATA_MODE")
        if dm_var is not None:
            m = dm_var[iprof]
            return m.decode("utf-8").strip() if isinstance(m, bytes) else str(m).strip()
        return "R"

    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
    pdm_array = bgc_file.variables["PARAMETER_DATA_MODE"][iprof]

    if isinstance(station_params, np.ma.MaskedArray):
        station_params = station_params.filled(b" ")

    for j, p_bytes in enumerate(station_params):
        p_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
            for c in p_bytes
        ]).strip()

        if p_str == param_base:
            m_val = pdm_array[j]
            return m_val.decode("utf-8").strip() if isinstance(m_val, bytes) else str(m_val).strip()

    # Default to main DATA_MODE if parameter not found explicitly in STATION_PARAMETERS
    dm_var = bgc_file.variables.get("DATA_MODE")
    if dm_var is not None:
        m = dm_var[iprof]
        return m.decode("utf-8").strip() if isinstance(m, bytes) else str(m).strip()
    return "R"


def clean_and_fill_qc_variables(bgc_file):
    """
    Safely inspect, harmonize, carry over raw data/QC, and fix level-QC 
    and PROFILE_<PARAM>_QC across ALL profiles (iprof = 0, 1, 2, ...),
    strictly enforcing FillValue on _ADJUSTED_QC when in 'R' mode.
    """
    n_prof = bgc_file.dimensions["N_PROF"].size
    file_name = os.path.basename(bgc_file.filepath()) if hasattr(bgc_file, "filepath") else "File"

    # =========================================================================
    # STEP A: Carry over raw data and raw QC to ADJUSTED variables if missing
    # =========================================================================
    base_params = [
        "DOXY", "CHLA", "CHLA_FLUORESCENCE", "BBP700", "FLUORESCENCE_CHLA",
        "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE443", "DOWN_IRRADIANCE490", "DOWN_IRRADIANCE555"
    ]
    for param in base_params:
        adj_param = f"{param}_ADJUSTED"
        adj_qc_param = f"{param}_ADJUSTED_QC"
        raw_qc_param = f"{param}_QC"
        adj_err_param = f"{param}_ADJUSTED_ERROR"

        if param in bgc_file.variables and adj_param in bgc_file.variables:
            raw_var = bgc_file.variables[param]
            adj_var = bgc_file.variables[adj_param]
            raw_qc_var = bgc_file.variables[raw_qc_param] if raw_qc_param in bgc_file.variables else None
            adj_qc_var = bgc_file.variables[adj_qc_param] if adj_qc_param in bgc_file.variables else None
            adj_err_var = bgc_file.variables[adj_err_param] if adj_err_param in bgc_file.variables else None

            err_default = 10.0 if "DOXY" in param else (0.0005 if "BBP" in param else 0.07)

            for iprof in range(n_prof):
                mode = get_param_data_mode(bgc_file, iprof, param)
                if mode == "R":
                    continue  # Do not populate adjusted variables in Real-time mode

                r_slice = raw_var[iprof, :] if raw_var.ndim > 1 else raw_var[:]
                a_slice = adj_var[iprof, :] if adj_var.ndim > 1 else adj_var[:]

                r_vals = r_slice.filled(99999.0) if hasattr(r_slice, "filled") else np.array(r_slice)
                a_vals = a_slice.filled(99999.0) if hasattr(a_slice, "filled") else np.array(a_slice)

                raw_has_data = (~pd.isna(r_vals)) & (r_vals != 99999.0)
                adj_is_missing = pd.isna(a_vals) | (a_vals == 99999.0)
                carry_over_mask = raw_has_data & adj_is_missing

                if carry_over_mask.any():
                    a_slice_updated = np.array(a_slice, copy=True)
                    a_slice_updated[carry_over_mask] = r_vals[carry_over_mask]

                    if adj_var.ndim > 1:
                        bgc_file.variables[adj_param][iprof, :] = a_slice_updated
                    else:
                        bgc_file.variables[adj_param][:] = a_slice_updated

                    if raw_qc_var is not None and adj_qc_var is not None:
                        q_slice = raw_qc_var[iprof, :] if raw_qc_var.ndim > 1 else raw_qc_var[:]
                        q_vals = np.array(q_slice, copy=True)
                        if hasattr(q_vals, "filled"):
                            q_vals = q_vals.filled(b" ")
                        q_char = q_vals.astype("|S1")

                        adj_q_slice = adj_qc_var[iprof, :] if adj_qc_var.ndim > 1 else adj_qc_var[:]
                        adj_q_char = np.array(adj_q_slice, copy=True)
                        if hasattr(adj_q_char, "filled"):
                            adj_q_char = adj_q_char.filled(b" ")
                        adj_q_char = adj_q_char.astype("|S1")

                        adj_q_char[carry_over_mask] = q_char[carry_over_mask]

                        if adj_qc_var.ndim > 1:
                            bgc_file.variables[adj_qc_param][iprof, :] = adj_q_char
                        else:
                            bgc_file.variables[adj_qc_param][:] = adj_q_char

                    if adj_err_var is not None:
                        err_slice = adj_err_var[iprof, :] if adj_err_var.ndim > 1 else adj_err_var[:]
                        err_vals = np.array(err_slice, copy=True)
                        err_vals[carry_over_mask] = np.float32(err_default)

                        if adj_err_var.ndim > 1:
                            bgc_file.variables[adj_err_param][iprof, :] = err_vals
                        else:
                            bgc_file.variables[adj_err_param][:] = err_vals

                    print(
                        f"[{file_name}] Profile [{iprof}] {param}: Carried over {np.sum(carry_over_mask)} levels"
                        f" from {param} & {raw_qc_param} to {adj_param} & {adj_qc_param}"
                    )

    # =========================================================================
    # STEP B: Clean Level-by-Level QC Arrays Across ALL Profile Slots
    # =========================================================================
    for var_name, var in list(bgc_file.variables.items()):
        if not var_name.endswith("_QC") or var_name.startswith("PROFILE_"):
            continue

        if "PH" in var_name or "NITRATE" in var_name or "TEMP_CPU_CHLA" in var_name:
            continue

        is_adjusted_qc = var_name.endswith("_ADJUSTED_QC")
        base_param = var_name[:-12] if is_adjusted_qc else var_name[:-3]

        for iprof in range(n_prof):
            try:
                mode = get_param_data_mode(bgc_file, iprof, base_param)

                # STRICT RULE: In Real-Time ('R') Mode, _ADJUSTED_QC MUST be FillValue (' ')
                if is_adjusted_qc and mode == "R":
                    n_lev = bgc_file.dimensions["N_LEVELS"].size
                    fill_qc = np.full(shape=(n_lev,), fill_value=b" ", dtype="|S1")
                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = fill_qc
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = b" "
                    else:
                        bgc_file.variables[var_name][:] = fill_qc

                    adj_var_name = f"{base_param}_ADJUSTED"
                    adj_err_name = f"{base_param}_ADJUSTED_ERROR"

                    if adj_var_name in bgc_file.variables:
                        adj_fill = np.full(shape=(n_lev,), fill_value=99999.0, dtype="float32")
                        bgc_file.variables[adj_var_name][iprof, :] = np.ma.masked_equal(adj_fill, 99999.0)

                    if adj_err_name in bgc_file.variables:
                        err_fill = np.full(shape=(n_lev,), fill_value=99999.0, dtype="float32")
                        bgc_file.variables[adj_err_name][iprof, :] = np.ma.masked_equal(err_fill, 99999.0)

                    continue

                p_var_name = (
                    f"{base_param}_ADJUSTED"
                    if is_adjusted_qc and f"{base_param}_ADJUSTED" in bgc_file.variables
                    else base_param
                )

                if p_var_name in bgc_file.variables:
                    p_var = bgc_file.variables[p_var_name]
                    p_raw = p_var[iprof, :] if p_var.ndim > 1 else p_var[:]
                    p_vals = p_raw.filled(99999.0) if hasattr(p_raw, "filled") else np.array(p_raw)

                    q_raw = var[iprof, :] if var.ndim > 1 else var[iprof]
                    q_vals = np.array(q_raw, copy=True)
                    if hasattr(q_vals, "filled"):
                        q_vals = q_vals.filled(b" ")

                    q_char = q_vals.astype("|S1")

                    valid_data_mask = (~pd.isna(p_vals)) & (p_vals != 99999.0)
                    missing_mask = ~valid_data_mask

                    # Assign '9' to missing levels
                    if missing_mask.any():
                        q_char[missing_mask] = b"9"

                    # If ADJUSTED_QC is missing but raw _QC exists, copy raw QC flag
                    if is_adjusted_qc and f"{base_param}_QC" in bgc_file.variables:
                        raw_qc_var = bgc_file.variables[f"{base_param}_QC"]
                        raw_q_slice = raw_qc_var[iprof, :] if raw_qc_var.ndim > 1 else raw_qc_var[:]
                        raw_q_char = np.array(raw_q_slice, copy=True)
                        if hasattr(raw_q_char, "filled"):
                            raw_q_char = raw_q_char.filled(b" ")
                        raw_q_char = raw_q_char.astype("|S1")

                        can_copy_mask = valid_data_mask & np.isin(raw_q_char, [b"1", b"2", b"3", b"4"])
                        if can_copy_mask.any():
                            q_char[can_copy_mask] = raw_q_char[can_copy_mask]

                    # Force '1' to ANY non-'1'..'4' flags at levels with valid numeric data
                    invalid_qc_mask = valid_data_mask & (~np.isin(q_char, [b"1", b"2", b"3", b"4"]))

                    if invalid_qc_mask.any():
                        q_char[invalid_qc_mask] = b"1"

                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = q_char
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = q_char if q_char.ndim == 0 else q_char[0]
                    else:
                        bgc_file.variables[var_name][:] = q_char

                    if var_name in [
                        "DOXY_QC",
                        "DOXY_ADJUSTED_QC",
                        "FLUORESCENCE_CHLA_QC",
                        "FLUORESCENCE_CHLA_ADJUSTED_QC",
                        "CHLA_QC",
                        "CHLA_ADJUSTED_QC",
                        "CHLA_FLUORESCENCE_QC",
                        "CHLA_FLUORESCENCE_ADJUSTED_QC",
                        "BBP700_QC",
                        "BBP700_ADJUSTED_QC",
                    ]:
                        q_str_vals = [c.decode("utf-8", errors="ignore") for c in q_char[valid_data_mask]]
                        unique_qcs, counts = np.unique(q_str_vals, return_counts=True)
                        qc_summary = (
                            ", ".join([f"'{k}': {v}" for k, v in zip(unique_qcs, counts)])
                            if len(unique_qcs) > 0
                            else "None"
                        )
                        print(
                            f"[{file_name}] Profile [{iprof}] {var_name} values at valid data levels "
                            f"({np.sum(valid_data_mask)} total levels): {qc_summary}"
                        )

            except Exception:
                traceback.print_exc()

    # =========================================================================
    # STEP C: Compute PROFILE_<PARAM>_QC Grade Letters ('A'-'F' or ' ') for ALL Profiles
    # =========================================================================
    for var_name in list(bgc_file.variables.keys()):
        if not var_name.startswith("PROFILE_") or not var_name.endswith("_QC"):
            continue

        param_base = var_name[8:-3]
        adj_qc_var_name = (
            f"{param_base}_ADJUSTED_QC"
            if f"{param_base}_ADJUSTED_QC" in bgc_file.variables
            else f"{param_base}_QC"
        )

        for iprof in range(n_prof):
            param_is_present = False

            if "STATION_PARAMETERS" in bgc_file.variables:
                station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
                if isinstance(station_params, np.ma.MaskedArray):
                    station_params = station_params.filled(b" ")

                for p in station_params:
                    p_str = "".join([
                        (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                        for c in p
                    ]).strip()
                    if p_str == param_base and has_valid_measurements(bgc_file, param_base, iprof):
                        param_is_present = True
                        break

            if adj_qc_var_name in bgc_file.variables:
                target_qc_var = bgc_file.variables[adj_qc_var_name]
                target_qc_array = (
                    target_qc_var[iprof, :] if target_qc_var.ndim > 1 else target_qc_var[:]
                )
                prof_grade = get_profile_qc_grade(target_qc_array, parameter_present=param_is_present)
                grade_bytes = prof_grade.encode("utf-8")

                if bgc_file.variables[var_name].ndim == 1:
                    bgc_file.variables[var_name][iprof] = grade_bytes
                elif bgc_file.variables[var_name].ndim == 2:
                    bgc_file.variables[var_name][iprof, 0] = grade_bytes

                print(
                    f"[{file_name}] Profile [{iprof}] {var_name} grade calculated: '{prof_grade}'"
                )


def add_missing_valid_range_attributes(bgc_file):
    """Assign valid_min and valid_max only to allowed variables in ARGO_VALID_RANGES."""
    for var_name, (vmin, vmax) in ARGO_VALID_RANGES.items():
        if var_name in bgc_file.variables:
            var = bgc_file.variables[var_name]

            if "valid_min" not in var.ncattrs():
                var.setncattr("valid_min", np.float32(vmin))

            if "valid_max" not in var.ncattrs():
                var.setncattr("valid_max", np.float32(vmax))


def remove_forbidden_attributes(bgc_file):
    """Strictly remove valid_min and valid_max from all variables NOT in ARGO_VALID_RANGES."""
    forbidden_attrs = ["valid_min", "valid_max"]
    for var_name, var in bgc_file.variables.items():
        if var_name not in ARGO_VALID_RANGES:
            for attr in forbidden_attrs:
                if attr in var.ncattrs():
                    var.delncattr(attr)


def detect_parameter_profile(ds, param_prefix):
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
    path, name = os.path.split(filename)
    bd_name = re.sub(r"^(AOML_)?BR", "BD", name)
    out_dir = dest_dir if dest_dir else path
    w_filename = os.path.join(out_dir, f"w_{bd_name}")
    shutil.copyfile(filename, w_filename)
    return w_filename


def get_profile_chla(filename):
    profile = filename[-6:-3]
    return int(profile)


def organize_b_files(bd_files, br_files):
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
    hix = nc_ds.dimensions["N_HISTORY"].size
    target_hix = max(0, hix - 1)
    for name, value in dct.items():
        if name in nc_ds.variables:
            char_len = nc_ds.dimensions[nc_ds[name].dimensions[-1]].size
            padded_val = str(value).ljust(char_len)[:char_len]
            char_arr = np.array(padded_val, dtype=f"S{char_len}")
            nc_ds[name][target_hix, iprof_idx, :] = nc.stringtochar(char_arr)


def write_history_metadata(bgc_file, iprof_chla, iprof_doxy):
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


def write_scientific_calib_chla(bgc_file, WMOfloatid, idx_profile, df_bio):
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

            if param_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
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
    n_prof = bgc_file.dimensions["N_PROF"].size
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    # Detect CSV column names for Chlorophyll Fluorescence
    fluo_csv_col = None
    for candidate in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE_ADJUSTED"]:
        if candidate in df_bio.columns:
            fluo_csv_col = candidate
            break

    fluo_qc_csv_col = None
    for candidate in ["CHLA_FLUORESCENCE_QC", "FLUORESCENCE_CHLA_QC", "CHLA_FLUORESCENCE_ADJUSTED_QC"]:
        if candidate in df_bio.columns:
            fluo_qc_csv_col = candidate
            break

    for prof in range(n_prof):
        mode = get_param_data_mode(bgc_file, prof, "CHLA")

        CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_Array[:] = 99999.0

        CHLA_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1"
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

        FLUORESCENCE_CHLA_Adjusted_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        FLUORESCENCE_CHLA_Adjusted_Array[:] = 99999.0

        FLUORESCENCE_CHLA_AdjustedQC_Array = np.full(
            shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1"
        )

        target_fluo_name = "FLUORESCENCE_CHLA" if "FLUORESCENCE_CHLA" in bgc_file.variables else "CHLA_FLUORESCENCE"
        target_fluo_qc_name = f"{target_fluo_name}_QC"

        has_fluo = target_fluo_name in bgc_file.variables
        has_fluo_qc = target_fluo_qc_name in bgc_file.variables

        raw_fluo_var = (
            bgc_file.variables[target_fluo_name][prof, :]
            if (has_fluo and bgc_file.variables[target_fluo_name].ndim > 1)
            else (bgc_file.variables[target_fluo_name][:] if has_fluo else None)
        )
        raw_fluo_qc_var = (
            bgc_file.variables[target_fluo_qc_name][prof, :]
            if (has_fluo_qc and bgc_file.variables[target_fluo_qc_name].ndim > 1)
            else (
                bgc_file.variables[target_fluo_qc_name][:]
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

        FLUORESCENCE_CHLA_Adjusted_ERROR_Array = np.ma.masked_all(
            shape=(n_levels,), dtype="float32"
        )
        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[:] = 99999.0

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

        if mode in ["D", "A"]:
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
                    else str(chla_qc_arr[i]) if chla_qc_arr is not None else "1"
                )
                raw_fluo_qc_char = (
                    fluo_qc_arr[i].decode("utf-8")
                    if (fluo_qc_arr is not None and isinstance(fluo_qc_arr[i], bytes))
                    else str(fluo_qc_arr[i]) if fluo_qc_arr is not None else "1"
                )

                # CHLA Alignment
                if is_pres_missing_or_bad or is_raw_chla_missing or is_broken:
                    CHLA_Adjusted_Array[i] = 99999.0
                    CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    CHLA_Adjusted_Array.mask[i] = True
                    CHLA_Adjusted_ERROR_Array.mask[i] = True
                    target_pres_qc = (
                        b"4" if (is_broken or pres_qc_char == "4") else b"9"
                    )
                    CHLA_AdjustedQC_Array[i] = target_pres_qc

                    if chla_data_arr is not None:
                        chla_data_arr[i] = 99999.0
                    if chla_qc_arr is not None:
                        chla_qc_arr[i] = target_pres_qc
                elif matched_row is not None:
                    csv_chla_qc = (
                        matched_row["CHLA_FINAL_QC"]
                        if "CHLA_FINAL_QC" in matched_row
                        else "1"
                    )
                    chla_qc_str = str(int(csv_chla_qc)) if pd.notna(csv_chla_qc) else "1"
                    chla_final_val = (
                        matched_row["CHLA_FINAL"] if "CHLA_FINAL" in matched_row else np.nan
                    )

                    if (
                        raw_chla_qc_char in ["4"]
                        or chla_qc_str in ["4"]
                        or pd.isna(chla_final_val)
                        or chla_final_val == 99999.0
                    ):
                        CHLA_Adjusted_Array[i] = 99999.0
                        CHLA_Adjusted_ERROR_Array[i] = 99999.0
                        CHLA_Adjusted_Array.mask[i] = True
                        CHLA_Adjusted_ERROR_Array.mask[i] = True
                        target_flag = b"4" if (raw_chla_qc_char == "4" or chla_qc_str == "4") else b"1"
                        CHLA_AdjustedQC_Array[i] = target_flag

                        if chla_data_arr is not None:
                            chla_data_arr[i] = 99999.0
                        if chla_qc_arr is not None:
                            chla_qc_arr[i] = target_flag
                    else:
                        CHLA_Adjusted_Array[i] = np.float32(chla_final_val)
                        CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                        CHLA_Adjusted_Array.mask[i] = False
                        CHLA_Adjusted_ERROR_Array.mask[i] = False
                        valid_flag = (
                            chla_qc_str.encode("utf-8")
                            if chla_qc_str in ["1", "2", "3"]
                            else b"1"
                        )
                        CHLA_AdjustedQC_Array[i] = valid_flag
                        if chla_qc_arr is not None:
                            chla_qc_arr[i] = valid_flag
                else:
                    CHLA_Adjusted_Array[i] = np.float32(chla_data_arr[i])
                    CHLA_Adjusted_Array.mask[i] = False
                    CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                    CHLA_Adjusted_ERROR_Array.mask[i] = False

                    fallback_flag = (
                        raw_chla_qc_char.encode("utf-8")
                        if raw_chla_qc_char in ["1", "2", "3"]
                        else b"1"
                    )
                    CHLA_AdjustedQC_Array[i] = fallback_flag
                    if chla_qc_arr is not None:
                        chla_qc_arr[i] = fallback_flag

                # FLUORESCENCE Alignment
                if is_pres_missing_or_bad or is_broken:
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                    target_pres_qc = (
                        b"4" if (is_broken or pres_qc_char == "4") else b"9"
                    )
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = target_pres_qc

                    if fluo_data_arr is not None:
                        fluo_data_arr[i] = 99999.0
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = target_pres_qc
                elif matched_row is not None and fluo_csv_col is not None and pd.notna(matched_row.get(fluo_csv_col)):
                    fluo_csv_val = matched_row.get(fluo_csv_col)
                    fluo_adj_qc_val = matched_row.get(fluo_qc_csv_col) if fluo_qc_csv_col else "1"
                    fluo_adj_qc_str = str(int(fluo_adj_qc_val)) if pd.notna(fluo_adj_qc_val) else "1"

                    if (
                        raw_fluo_qc_char in ["4"]
                        or fluo_adj_qc_str in ["4"]
                        or pd.isna(fluo_csv_val)
                        or fluo_csv_val == 99999.0
                    ):
                        FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                        FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                        target_fluo_flag = b"4" if (raw_fluo_qc_char == "4" or fluo_adj_qc_str == "4") else b"1"
                        FLUORESCENCE_CHLA_AdjustedQC_Array[i] = target_fluo_flag

                        if fluo_data_arr is not None:
                            fluo_data_arr[i] = 99999.0
                        if fluo_qc_arr is not None:
                            fluo_qc_arr[i] = target_fluo_flag
                    else:
                        FLUORESCENCE_CHLA_Adjusted_Array[i] = np.float32(fluo_csv_val)
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                        FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = False
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = False
                        valid_fluo_flag = (
                            fluo_adj_qc_str.encode("utf-8")
                            if fluo_adj_qc_str in ["1", "2", "3"]
                            else b"1"
                        )
                        FLUORESCENCE_CHLA_AdjustedQC_Array[i] = valid_fluo_flag

                        if fluo_data_arr is not None:
                            fluo_data_arr[i] = np.float32(fluo_csv_val)
                        if fluo_qc_arr is not None:
                            fluo_qc_arr[i] = valid_fluo_flag
                elif fluo_data_arr is not None and not is_raw_fluo_missing:
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = np.float32(fluo_data_arr[i])
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = False
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = False

                    fallback_flag = (
                        raw_fluo_qc_char.encode("utf-8")
                        if raw_fluo_qc_char in ["1", "2", "3"]
                        else b"1"
                    )
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = fallback_flag
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = fallback_flag
                else:
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"9"
                    if fluo_qc_arr is not None:
                        fluo_qc_arr[i] = b"9"

        bgc_file.variables["CHLA_ADJUSTED"][prof, :] = CHLA_Adjusted_Array
        bgc_file.variables["CHLA_ADJUSTED_QC"][prof, :] = CHLA_AdjustedQC_Array
        bgc_file.variables["CHLA_ADJUSTED_ERROR"][prof, :] = CHLA_Adjusted_ERROR_Array

        if has_chla:
            bgc_file.variables["CHLA"][prof, :] = chla_data_arr
        if has_chla_qc:
            bgc_file.variables["CHLA_QC"][prof, :] = chla_qc_arr

        # Write to both potential variable name variants in NetCDF
        for var_prefix in ["FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]:
            adj_var = f"{var_prefix}_ADJUSTED"
            adj_qc_var = f"{var_prefix}_ADJUSTED_QC"
            adj_err_var = f"{var_prefix}_ADJUSTED_ERROR"
            raw_v = var_prefix
            raw_q = f"{var_prefix}_QC"

            if adj_var in bgc_file.variables:
                bgc_file.variables[adj_var][prof, :] = FLUORESCENCE_CHLA_Adjusted_Array
            if adj_qc_var in bgc_file.variables:
                bgc_file.variables[adj_qc_var][prof, :] = FLUORESCENCE_CHLA_AdjustedQC_Array
            if adj_err_var in bgc_file.variables:
                bgc_file.variables[adj_err_var][prof, :] = FLUORESCENCE_CHLA_Adjusted_ERROR_Array

            if raw_v in bgc_file.variables and fluo_data_arr is not None:
                bgc_file.variables[raw_v][prof, :] = fluo_data_arr
            if raw_q in bgc_file.variables and fluo_qc_arr is not None:
                bgc_file.variables[raw_q][prof, :] = fluo_qc_arr


def write_BBP700_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio, iprof_idx=0):
    if "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    n_prof = bgc_file.dimensions["N_PROF"].size
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    for prof in range(n_prof):
        mode = get_param_data_mode(bgc_file, prof, "BBP700")

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
            shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1"
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

        if mode in ["D", "A"]:
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
                    else str(bbp_qc_arr[i]) if bbp_qc_arr is not None else "1"
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

                    if bbp_data_arr is not None:
                        bbp_data_arr[i] = 99999.0
                    if bbp_qc_arr is not None:
                        bbp_qc_arr[i] = target_pres_qc
                elif matched_row is not None:
                    csv_qc = (
                        matched_row["BBP700_FINAL_QC"]
                        if "BBP700_FINAL_QC" in matched_row
                        else "1"
                    )
                    qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "1"
                    raw_bbp = (
                        matched_row["BBP700_FINAL"]
                        if "BBP700_FINAL" in matched_row
                        else np.nan
                    )

                    if (
                        raw_qc_char in ["4"]
                        or qc_str in ["4"]
                        or pd.isna(raw_bbp)
                        or raw_bbp == 99999.0
                    ):
                        BBP700_Adjusted_Array[i] = 99999.0
                        BBP700_Adjusted_ERROR_Array[i] = 99999.0
                        BBP700_Adjusted_Array.mask[i] = True
                        BBP700_Adjusted_ERROR_Array.mask[i] = True
                        target_flag = b"4" if (raw_qc_char == "4" or qc_str == "4") else b"1"
                        BBP700_AdjustedQC_Array[i] = target_flag

                        if bbp_data_arr is not None:
                            bbp_data_arr[i] = 99999.0
                        if bbp_qc_arr is not None:
                            bbp_qc_arr[i] = target_flag
                    else:
                        BBP700_Adjusted_Array[i] = np.float32(raw_bbp)
                        BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                        BBP700_Adjusted_Array.mask[i] = False
                        BBP700_Adjusted_ERROR_Array.mask[i] = False
                        valid_bbp_flag = (
                            qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1"
                        )
                        BBP700_AdjustedQC_Array[i] = valid_bbp_flag
                        if bbp_qc_arr is not None:
                            bbp_qc_arr[i] = valid_bbp_flag
                else:
                    BBP700_Adjusted_Array[i] = np.float32(bbp_data_arr[i])
                    BBP700_Adjusted_Array.mask[i] = False
                    BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                    BBP700_Adjusted_ERROR_Array.mask[i] = False

                    fallback_flag = (
                        raw_qc_char.encode("utf-8")
                        if raw_qc_char in ["1", "2", "3"]
                        else b"1"
                    )
                    BBP700_AdjustedQC_Array[i] = fallback_flag
                    if bbp_qc_arr is not None:
                        bbp_qc_arr[i] = fallback_flag

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


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle):
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
    """Populate DOXY_ADJUSTED and synchronize DOXY_ADJUSTED_QC across ALL profiles."""
    var_names = ds.variables.keys()
    n_prof = ds.dimensions["N_PROF"].size
    n_levels = ds.dimensions["N_LEVELS"].size

    pres_nc_full = ds.variables["PRES"][:]
    has_pres_qc = "PRES_QC" in ds.variables
    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()

    for iprof in range(n_prof):
        mode = get_param_data_mode(ds, iprof, "DOXY")
        if mode == "R":
            continue

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
                else str(doxy_qc_arr[i]) if doxy_qc_arr is not None else "1"
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
                qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "1"

                raw_doxy_final = matched_row.get("DOXY_FINAL")

                if (
                    raw_qc_char in ["4"]
                    or qc_str in ["4"]
                    or pd.isna(raw_doxy_final)
                    or raw_doxy_final == 99999.0
                ):
                    DOXY_Adjusted_Array[i] = 99999.0
                    DOXY_Adjusted_ERROR_Array[i] = 99999.0
                    DOXY_Adjusted_Array.mask[i] = True
                    DOXY_Adjusted_ERROR_Array.mask[i] = True
                    target_flag = b"4" if (raw_qc_char == "4" or qc_str == "4") else b"1"
                    DOXY_AdjustedQC_Array[i] = target_flag
                    if doxy_qc_arr is not None:
                        doxy_qc_arr[i] = target_flag
                else:
                    DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
                    DOXY_Adjusted_Array.mask[i] = False
                    valid_flag = (
                        raw_qc_char.encode("utf-8")
                        if raw_qc_char in ["1", "2", "3"]
                        else (
                            qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1"
                        )
                    )
                    DOXY_AdjustedQC_Array[i] = valid_flag
                    if doxy_qc_arr is not None:
                        doxy_qc_arr[i] = valid_flag

                    raw_doxy_error = matched_row.get("DOXY_ADJUSTED_ERROR")
                    if pd.notna(raw_doxy_error) and raw_doxy_error != 99999.0:
                        DOXY_Adjusted_Error_Array[i] = np.float32(raw_doxy_error)
                        DOXY_Adjusted_Error_Array.mask[i] = False
                    else:
                        DOXY_Adjusted_Error_Array[i] = np.float32(DOXY_Adjusted_ERROR_est_default)
                        DOXY_Adjusted_Error_Array.mask[i] = False
            else:
                DOXY_Adjusted_Array[i] = np.float32(doxy_data_arr[i])
                DOXY_Adjusted_Array.mask[i] = False
                DOXY_Adjusted_Error_Array[i] = np.float32(DOXY_Adjusted_ERROR_est_default)
                DOXY_Adjusted_Error_Array.mask[i] = False

                fallback_qc = (
                    raw_qc_char.encode("utf-8")
                    if raw_qc_char in ["1", "2", "3"]
                    else b"1"
                )
                DOXY_AdjustedQC_Array[i] = fallback_qc
                if doxy_qc_arr is not None:
                    doxy_qc_arr[i] = fallback_qc

        if "DOXY_QC" in var_names and doxy_qc_arr is not None:
            if ds.variables["DOXY_QC"].ndim > 1:
                ds.variables["DOXY_QC"][iprof, :] = doxy_qc_arr
            else:
                ds.variables["DOXY_QC"][:] = doxy_qc_arr

        if "DOXY_ADJUSTED" in var_names:
            doxy_adj_var = ds.variables["DOXY_ADJUSTED"]
            doxy_adj_qc_var = ds.variables["DOXY_ADJUSTED_QC"] if "DOXY_ADJUSTED_QC" in var_names else None

            if doxy_adj_var.ndim > 1:
                if doxy_adj_var.shape[0] == n_prof:
                    doxy_adj_var[iprof, :] = DOXY_Adjusted_Array
                    if doxy_adj_qc_var is not None:
                        doxy_adj_qc_var[iprof, :] = DOXY_AdjustedQC_Array
                else:
                    doxy_adj_var[:, iprof] = DOXY_Adjusted_Array
                    if doxy_adj_qc_var is not None:
                        doxy_adj_qc_var[:, iprof] = DOXY_AdjustedQC_Array
            else:
                doxy_adj_var[:] = DOXY_Adjusted_Array
                if doxy_adj_qc_var is not None:
                    doxy_adj_qc_var[:] = DOXY_AdjustedQC_Array

        if "DOXY_ADJUSTED_ERROR" in var_names:
            doxy_adj_err_var = ds.variables["DOXY_ADJUSTED_ERROR"]
            if doxy_adj_err_var.ndim > 1:
                if doxy_adj_err_var.shape[0] == n_prof:
                    doxy_adj_err_var[iprof, :] = DOXY_Adjusted_Error_Array
                else:
                    doxy_adj_err_var[:, iprof] = DOXY_Adjusted_Error_Array
            else:
                doxy_adj_err_var[:] = DOXY_Adjusted_Error_Array


def safe_rename(from_file, to_file):
    gc.collect()
    try:
        os.replace(from_file, to_file)
    except OSError:
        try:
            shutil.copyfile(from_file, to_file)
            os.remove(from_file)
        except Exception:
            pass


# ==============================================================================
# SECTION 3: Processing Loop
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

        w_bgc_filename = create_working_bd_file(
            bgc_filename, dest_dir=main_float_dir
        )

        try:
            ds = nc.Dataset(w_bgc_filename, "r+")

            iprof_chla = detect_parameter_profile(ds, "CHLA")
            iprof_bbp = detect_parameter_profile(ds, "BBP700")
            iprof_doxy = detect_parameter_profile(ds, "DOXY")

            # 1. Update PARAMETER_DATA_MODE and DATA_MODE based on actual measurements
            write_parameter_data_modes(ds)

            # 2. Process CHLA, BBP700, DOXY, and FLUORESCENCE_CHLA delayed-mode parameters
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_chla)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_bbp)
            write_DOXY_slope_drift(ds, 0, df_bio, idx_profile)
            write_DOXY_from_csv(ds, df_bio, idx_profile)

            # 3. Carry over missing adjusted levels and clean QC & PROFILE_<PARAM>_QC variables across ALL profiles
            clean_and_fill_qc_variables(ds)

            # 4. Metadata and attributes
            write_history_metadata(ds, iprof_chla, iprof_doxy)
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

    for w_file in new_bd_files:
        path, name = os.path.split(w_file)

        new_name = name.replace("w_AOML_BD", "BD")
        new_name = new_name.replace("w_AOML_BR", "BD")
        new_name = new_name.replace("w_BD", "BD")
        new_name = new_name.replace("w_BR", "BD")

        new_path = os.path.join(final_out_dir, new_name)
        safe_rename(w_file, new_path)
        processed_output_paths.append(new_path)