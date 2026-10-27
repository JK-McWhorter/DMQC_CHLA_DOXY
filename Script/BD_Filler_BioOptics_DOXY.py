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

import matplotlib
matplotlib.use("Agg")  # Headless cluster/server backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4 as nc
import numpy as np
import pandas as pd
from PIL import Image
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import griddata

try:
    import config
    TEXT_COLOR = getattr(config, "text_color", "#2c3e50")
    BG_COLOR = getattr(config, "background_color", "#ffffff")
except ImportError:
    TEXT_COLOR = "#2c3e50"
    BG_COLOR = "#ffffff"

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
DOXY_EXCLUDE_FLOAT_IDS = {4903910, 4903911, 7902327, 3902693, 1902800}

comment_dmqc_operator_chla = (
    "PRIMARY | https://orcid.org/0000-0003-1297-6599 | Jennifer McWhorter, NOAA/AOML"
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
    "CHLA = ((FLUORESCENCE_CHLA - MEDIAN(PRELIM_DARK_CHLA)) * SCALE_CHLA) / PHYSIO_RATIO"
)
scientific_calibration_coefficient_CHLA = "PHYSIO_RATIO=1.0"

CHLA_Adjusted_ERROR_est = 0.07
BBP700_Adjusted_ERROR_est = 0.0005
DOXY_Adjusted_ERROR_est_default = 10.0

scientific_calibration_equation_DOXY = "DOXY_ADJUSTED = DOXY * slope + drift"
scientific_calibration_comment_DOXY = (
    "DOXY calibration computed against WOA climatology following Bittig et al. (2018) methodology"
)

scientific_calibration_comment_BROKEN = (
    "Sensor failure / broken sensor - data flagged bad or missing (Argo BGC QC Manual)"
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
# SECTION 2: Architecture & Sanitation Helper Functions
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
        "CHLA_FINAL", "BBP700_FINAL", "CHLA", "CHLA_ADJUSTED",
        "CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_ADJUSTED",
        "FLUORESCENCE_CHLA", "FLUORESCENCE_CHLA_ADJUSTED",
        "CHLA_NoLUT", "BBP700", "SCALE_CHLA", "DARK_CHLA",
    ]

    qc_cols = [
        "CHLA_QC", "CHLA_FINAL_QC", "BBP700_QC", "BBP700_FINAL_QC",
        "CHLA_ADJUSTED_QC", "CHLA_FLUORESCENCE_QC",
        "CHLA_FLUORESCENCE_ADJUSTED_QC", "FLUORESCENCE_CHLA_QC",
        "FLUORESCENCE_CHLA_ADJUSTED_QC", "BBP700_ADJUSTED_QC",
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


def detect_parameter_profile(ds, param_prefix):
    """
    Identifies the profile index (N_PROF) where a given parameter resides.
    Accounts for decoupled APEX sub-profiles (CTD vs. Optics vs. Oxygen)
    and resolves legacy naming aliases (FLUORESCENCE_CHLA <-> CHLA_FLUORESCENCE).
    """
    if "STATION_PARAMETERS" not in ds.variables:
        return 0

    n_prof = ds.dimensions["N_PROF"].size
    n_param = ds.dimensions["N_PARAM"].size

    target_prefixes = [param_prefix]
    if param_prefix in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_prefixes = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
    elif param_prefix in ["PH", "PH_IN_SITU_TOTAL"]:
        target_prefixes = ["PH", "PH_IN_SITU_TOTAL"]

    for iprof in range(n_prof):
        param_mat = ds.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        for j in range(n_param):
            p_bytes = param_mat[j]
            param_str = "".join([
                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                for c in p_bytes
            ]).strip()

            for tp in target_prefixes:
                if param_str == tp or param_str.startswith(tp):
                    return iprof
    return 0


def has_valid_measurements(bgc_file, param, iprof=None):
    if iprof is None:
        iprof = detect_parameter_profile(bgc_file, param)

    param_candidates = [param]
    if param == "PRES":
        param_candidates.extend(["PRES_ADJUSTED"])
    elif param == "CHLA":
        param_candidates.extend(["CHLA_ADJUSTED", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"])
    elif param in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        param_candidates.extend(["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED", "FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"])
    elif param == "BBP700":
        param_candidates.extend(["BBP700_ADJUSTED"])
    elif param == "DOXY":
        param_candidates.extend(["DOXY_ADJUSTED"])
    elif param in ["PH", "PH_IN_SITU_TOTAL"]:
        param_candidates.extend(["PH_IN_SITU_TOTAL", "PH_IN_SITU_TOTAL_ADJUSTED"])

    for p in param_candidates:
        if p in bgc_file.variables:
            var = bgc_file.variables[p]
            p_slice = var[iprof, :] if var.ndim > 1 else var[:]
            p_vals = p_slice.filled(99999.0) if hasattr(p_slice, "filled") else np.array(p_slice)
            valid_mask = (~pd.isna(p_vals)) & (p_vals != 99999.0)
            if valid_mask.any():
                return True
    return False


def get_param_data_mode(bgc_file, iprof, param_base):
    if param_base == "PRES":
        return "R"

    if iprof is None:
        iprof = detect_parameter_profile(bgc_file, param_base)

    if "STATION_PARAMETERS" not in bgc_file.variables or "PARAMETER_DATA_MODE" not in bgc_file.variables:
        dm_var = bgc_file.variables.get("DATA_MODE")
        if dm_var is not None:
            m = dm_var[iprof] if dm_var.ndim > 0 else dm_var[:]
            return m.decode("utf-8").strip() if isinstance(m, bytes) else str(m).strip()
        return "R"

    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
    pdm_array = bgc_file.variables["PARAMETER_DATA_MODE"][iprof]

    if isinstance(station_params, np.ma.MaskedArray):
        station_params = station_params.filled(b" ")

    target_names = [param_base]
    if param_base in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_names = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
    elif param_base in ["PH", "PH_IN_SITU_TOTAL"]:
        target_names = ["PH", "PH_IN_SITU_TOTAL"]

    for j, p_bytes in enumerate(station_params):
        p_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
            for c in p_bytes
        ]).strip()

        if p_str in target_names or any(p_str.startswith(tn) for tn in target_names):
            m_val = pdm_array[j]
            return m_val.decode("utf-8").strip() if isinstance(m_val, bytes) else str(m_val).strip()

    dm_var = bgc_file.variables.get("DATA_MODE")
    if dm_var is not None:
        m = dm_var[iprof] if dm_var.ndim > 0 else dm_var[:]
        return m.decode("utf-8").strip() if isinstance(m, bytes) else str(m).strip()
    return "R"


def write_parameter_data_modes(bgc_file):
    n_param = bgc_file.dimensions["N_PARAM"].size
    n_prof = bgc_file.dimensions["N_PROF"].size

    pdm = bgc_file.variables["PARAMETER_DATA_MODE"][:]
    data_mode = bgc_file.variables["DATA_MODE"][:]

    for iprof in range(n_prof):
        param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        for j in range(n_param):
            p_bytes = param_mat[j]
            p_str = "".join([
                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                for c in p_bytes
            ]).strip()

            if p_str:
                if p_str == "PRES" or p_str.startswith("DOWN_IRRADIANCE") or p_str.startswith("UP_RADIANCE"):
                    pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"
                elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY", "NITRATE", "PH_IN_SITU_TOTAL"]:
                    if has_valid_measurements(bgc_file, p_str, iprof):
                        pdm[iprof, j] = b"D" if pdm.dtype.kind in ["S", "U", "O"] else "D"
                        data_mode[iprof] = b"D" if data_mode.dtype.kind in ["S", "U", "O"] else "D"
                    else:
                        pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def get_profile_qc_grade(qc_masked_array, parameter_present=True):
    if not parameter_present:
        return " "

    if hasattr(qc_masked_array, "compressed"):
        unmasked_vals = qc_masked_array.compressed()
    else:
        unmasked_vals = qc_masked_array

    raw_qcs = []
    for q in np.ravel(unmasked_vals):
        if isinstance(q, np.ma.core.MaskedConstant) or q is np.ma.masked:
            continue
        s = q.decode("utf-8", errors="ignore").strip() if isinstance(q, bytes) else str(q).strip()
        if s != "":
            raw_qcs.append(s)

    valid_qcs = [q for q in raw_qcs if q not in ["", " ", "9"]]
    total_points = len(valid_qcs)

    if total_points == 0 or all(q == "4" for q in valid_qcs):
        return "F"

    good_count = sum(1 for q in valid_qcs if q in ["1", "2"])
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


def get_pres_array_for_profile(ds, iprof):
    for p_var in ["PRES_ADJUSTED", "PRES"]:
        if p_var in ds.variables:
            v = ds.variables[p_var]
            arr = v[iprof, :] if v.ndim > 1 else v[:]
            p_vals = arr.filled(99999.0) if hasattr(arr, "filled") else np.array(arr)
            valid = (~np.isnan(p_vals)) & (p_vals != 99999.0)
            if valid.any():
                return np.float32(p_vals)

    if "PRES" in ds.variables:
        v = ds.variables["PRES"]
        arr = v[0, :] if v.ndim > 1 else v[:]
        p_vals = arr.filled(99999.0) if hasattr(arr, "filled") else np.array(arr)
        return np.float32(p_vals)

    return np.full((ds.dimensions["N_LEVELS"].size,), 99999.0, dtype="float32")


def match_csv_row_by_pres(cycle_df, nc_pres):
    if cycle_df.empty or np.isnan(nc_pres) or nc_pres == 99999.0:
        return None

    csv_pres = cycle_df["PRES"].values
    diffs = np.abs(csv_pres - nc_pres)
    min_idx = np.argmin(diffs)

    if nc_pres <= 50.0:
        tol = 3.0
    elif nc_pres <= 200.0:
        tol = 12.0
    elif nc_pres <= 600.0:
        tol = 35.0
    else:
        tol = 80.0

    if diffs[min_idx] <= tol:
        return cycle_df.iloc[min_idx]
    return None


def clean_and_fill_qc_variables(bgc_file, WMOfloatid=None):
    n_prof = bgc_file.dimensions["N_PROF"].size
    skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

    base_params = [
        "DOXY", "CHLA", "CHLA_FLUORESCENCE", "BBP700", "FLUORESCENCE_CHLA",
        "DOWN_IRRADIANCE380", "DOWN_IRRADIANCE443", "DOWN_IRRADIANCE490", "DOWN_IRRADIANCE555"
    ]

    for param in base_params:
        if skip_doxy and "DOXY" in param:
            continue

        adj_param = f"{param}_ADJUSTED"
        adj_qc_param = f"{param}_ADJUSTED_QC"
        raw_qc_param = f"{param}_QC"
        adj_err_param = f"{param}_ADJUSTED_ERROR"

        if param in bgc_file.variables and adj_param in bgc_file.variables:
            raw_var = bgc_file.variables[param]
            adj_var = bgc_file.variables[adj_param]
            raw_qc_var = bgc_file.variables.get(raw_qc_param)
            adj_qc_var = bgc_file.variables.get(adj_qc_param)
            adj_err_var = bgc_file.variables.get(adj_err_param)
            err_default = 10.0 if "DOXY" in param else (0.0005 if "BBP" in param else 0.07)

            for iprof in range(n_prof):
                mode = get_param_data_mode(bgc_file, iprof, param)
                if mode == "R":
                    continue

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

    # Step B: Clean level QC arrays
    for var_name, var in list(bgc_file.variables.items()):
        if not var_name.endswith("_QC") or var_name.startswith("PROFILE_"):
            continue
        if "PH" in var_name or "NITRATE" in var_name or "TEMP_CPU_CHLA" in var_name:
            continue
        if skip_doxy and "DOXY" in var_name:
            continue

        is_adjusted_qc = var_name.endswith("_ADJUSTED_QC")
        base_param = var_name[:-12] if is_adjusted_qc else var_name[:-3]

        for iprof in range(n_prof):
            try:
                mode = get_param_data_mode(bgc_file, iprof, base_param)

                if is_adjusted_qc and mode == "R":
                    n_lev = bgc_file.dimensions["N_LEVELS"].size
                    fill_qc = np.full(shape=(n_lev,), fill_value=b" ", dtype="|S1")
                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = fill_qc
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = b" "
                    else:
                        bgc_file.variables[var_name][:] = fill_qc
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

                    if missing_mask.any():
                        q_char[missing_mask] = b"9"

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

                    invalid_qc_mask = valid_data_mask & (~np.isin(q_char, [b"1", b"2", b"3", b"4"]))
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

    # Step C: Update PROFILE_<PARAM>_QC
    for var_name in list(bgc_file.variables.keys()):
        if not var_name.startswith("PROFILE_") or not var_name.endswith("_QC"):
            continue

        param_base = var_name[8:-3]
        for iprof in range(n_prof):
            if param_base == "PRES":
                target_qc_name = "PRES_QC"
                source_prof = 0
            else:
                mode = get_param_data_mode(bgc_file, iprof, param_base)
                target_qc_name = (
                    f"{param_base}_ADJUSTED_QC"
                    if mode in ["D", "A"] and f"{param_base}_ADJUSTED_QC" in bgc_file.variables
                    else f"{param_base}_QC"
                )
                source_prof = iprof

            param_is_present = has_valid_measurements(bgc_file, param_base, source_prof)

            if target_qc_name in bgc_file.variables:
                target_qc_var = bgc_file.variables[target_qc_name]
                target_qc_array = (
                    target_qc_var[source_prof, :] if target_qc_var.ndim > 1 else target_qc_var[:]
                )
                prof_grade = get_profile_qc_grade(target_qc_array, parameter_present=param_is_present)
                grade_bytes = prof_grade.encode("utf-8")

                if bgc_file.variables[var_name].ndim == 1:
                    bgc_file.variables[var_name][iprof] = grade_bytes
                elif bgc_file.variables[var_name].ndim == 2:
                    bgc_file.variables[var_name][iprof, 0] = grade_bytes


def add_missing_valid_range_attributes(bgc_file):
    for var_name, (vmin, vmax) in ARGO_VALID_RANGES.items():
        if var_name in bgc_file.variables:
            var = bgc_file.variables[var_name]
            if "valid_min" not in var.ncattrs():
                var.setncattr("valid_min", np.float32(vmin))
            if "valid_max" not in var.ncattrs():
                var.setncattr("valid_max", np.float32(vmax))


def remove_forbidden_attributes(bgc_file):
    forbidden_attrs = ["valid_min", "valid_max"]
    for var_name, var in bgc_file.variables.items():
        if var_name not in ARGO_VALID_RANGES:
            for attr in forbidden_attrs:
                if attr in var.ncattrs():
                    var.delncattr(attr)


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


def write_history_metadata(bgc_file, iprof_chla, iprof_doxy, WMOfloatid=None):
    bgc_file.history = datetime.datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ creation")
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

    if WMOfloatid not in DOXY_EXCLUDE_FLOAT_IDS:
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
        bgc_file.variables["DATE_UPDATE"][:] = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))


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
        dark_cols = sorted([col for col in df_bio.columns if col.startswith("MIN_FLUOCHLA_CYCLE")])
        dark_vals = (
            df_bio[dark_cols].iloc[0].dropna().astype(int).tolist()
            if dark_cols
            else []
        )
        dark_str = " ".join(map(str, dark_vals)) if dark_vals else "NA"
        scale_val = df_bio["SCALE_CHLA"].iloc[0] if "SCALE_CHLA" in df_bio.columns else "NA"
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

    SciCalComArray_CHLA = nc.stringtochar(np.array(comment_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalComArray_CHLA_FLU = nc.stringtochar(np.array(comment_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalComArray_BBP = nc.stringtochar(np.array(comment_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalEquArray_CHLA = nc.stringtochar(np.array(equation_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalEquArray_BBP = nc.stringtochar(np.array(equation_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA = nc.stringtochar(np.array(coef_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA_FLU = nc.stringtochar(np.array(coef_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_BBP = nc.stringtochar(np.array(coef_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
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
                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                for c in p_bytes
            ]).strip()

            if param_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "CHLA":
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "BBP700":
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray


def write_chla_BBP_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio, iprof_idx=0):
    n_prof = bgc_file.dimensions["N_PROF"].size
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    fluo_csv_col = next((c for c in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE_ADJUSTED"] if c in df_bio.columns), None)
    fluo_qc_csv_col = next((c for c in ["CHLA_FLUORESCENCE_QC", "FLUORESCENCE_CHLA_QC", "CHLA_FLUORESCENCE_ADJUSTED_QC"] if c in df_bio.columns), None)

    for prof in range(n_prof):
        mode = get_param_data_mode(bgc_file, prof, "CHLA")

        CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_Array[:] = 99999.0
        CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")

        has_chla = "CHLA" in bgc_file.variables
        has_chla_qc = "CHLA_QC" in bgc_file.variables
        raw_chla_var = bgc_file.variables["CHLA"][prof, :] if (has_chla and bgc_file.variables["CHLA"].ndim > 1) else (bgc_file.variables["CHLA"][:] if has_chla else None)
        raw_chla_qc_var = bgc_file.variables["CHLA_QC"][prof, :] if (has_chla_qc and bgc_file.variables["CHLA_QC"].ndim > 1) else (bgc_file.variables["CHLA_QC"][:] if has_chla_qc else None)

        chla_data_arr = np.array(raw_chla_var, copy=True) if raw_chla_var is not None else None
        chla_qc_arr = np.array(raw_chla_qc_var, copy=True) if raw_chla_qc_var is not None else None

        CHLA_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_ERROR_Array[:] = 99999.0

        FLUORESCENCE_CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        FLUORESCENCE_CHLA_Adjusted_Array[:] = 99999.0
        FLUORESCENCE_CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")

        target_fluo_name = "FLUORESCENCE_CHLA" if "FLUORESCENCE_CHLA" in bgc_file.variables else "CHLA_FLUORESCENCE"
        target_fluo_qc_name = f"{target_fluo_name}_QC"
        has_fluo = target_fluo_name in bgc_file.variables
        has_fluo_qc = target_fluo_qc_name in bgc_file.variables

        raw_fluo_var = bgc_file.variables[target_fluo_name][prof, :] if (has_fluo and bgc_file.variables[target_fluo_name].ndim > 1) else (bgc_file.variables[target_fluo_name] if has_fluo else None)
        raw_fluo_qc_var = bgc_file.variables[target_fluo_qc_name][prof, :] if (has_fluo_qc and bgc_file.variables[target_fluo_qc_name].ndim > 1) else (bgc_file.variables[target_fluo_qc_name] if has_fluo_qc else None)

        fluo_data_arr = np.array(raw_fluo_var, copy=True) if raw_fluo_var is not None else None
        fluo_qc_arr = np.array(raw_fluo_qc_var, copy=True) if raw_fluo_qc_var is not None else None

        FLUORESCENCE_CHLA_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[:] = 99999.0

        nc_pres_all = get_pres_array_for_profile(bgc_file, prof)

        if mode in ["D", "A"]:
            for i in range(n_levels):
                nc_pres = nc_pres_all[i]
                is_pres_missing_or_bad = np.isnan(nc_pres) or nc_pres == 99999.0
                is_raw_chla_missing = chla_data_arr is None or pd.isna(chla_data_arr[i]) or chla_data_arr[i] == 99999.0

                matched_row = match_csv_row_by_pres(cycle_df, nc_pres) if not (is_broken or is_pres_missing_or_bad) else None

                raw_chla_qc_char = chla_qc_arr[i].decode("utf-8") if (chla_qc_arr is not None and isinstance(chla_qc_arr[i], bytes)) else str(chla_qc_arr[i]) if chla_qc_arr is not None else "1"
                raw_fluo_qc_char = fluo_qc_arr[i].decode("utf-8") if (fluo_qc_arr is not None and isinstance(fluo_qc_arr[i], bytes)) else str(fluo_qc_arr[i]) if fluo_qc_arr is not None else "1"

                if is_pres_missing_or_bad or is_raw_chla_missing or is_broken:
                    CHLA_Adjusted_Array[i] = 99999.0
                    CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    CHLA_Adjusted_Array.mask[i] = True
                    CHLA_Adjusted_ERROR_Array.mask[i] = True
                    CHLA_AdjustedQC_Array[i] = b"4" if is_broken else b"9"
                elif matched_row is not None:
                    csv_chla_qc = matched_row.get("CHLA_FINAL_QC", "1")
                    chla_qc_str = str(int(csv_chla_qc)) if pd.notna(csv_chla_qc) else "1"
                    chla_final_val = matched_row.get("CHLA_FINAL", np.nan)

                    if raw_chla_qc_char in ["4"] or chla_qc_str in ["4"] or pd.isna(chla_final_val) or chla_final_val == 99999.0:
                        CHLA_Adjusted_Array[i] = 99999.0
                        CHLA_Adjusted_ERROR_Array[i] = 99999.0
                        CHLA_Adjusted_Array.mask[i] = True
                        CHLA_Adjusted_ERROR_Array.mask[i] = True
                        CHLA_AdjustedQC_Array[i] = b"4" if (raw_chla_qc_char == "4" or chla_qc_str == "4") else b"1"
                    else:
                        CHLA_Adjusted_Array[i] = np.float32(chla_final_val)
                        CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                        CHLA_Adjusted_Array.mask[i] = False
                        CHLA_Adjusted_ERROR_Array.mask[i] = False
                        valid_flag = chla_qc_str.encode("utf-8") if chla_qc_str in ["1", "2", "3"] else b"1"
                        CHLA_AdjustedQC_Array[i] = valid_flag
                else:
                    CHLA_Adjusted_Array[i] = np.float32(chla_data_arr[i])
                    CHLA_Adjusted_Array.mask[i] = False
                    CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                    CHLA_Adjusted_ERROR_Array.mask[i] = False
                    fallback_flag = raw_chla_qc_char.encode("utf-8") if raw_chla_qc_char in ["1", "2", "3"] else b"1"
                    CHLA_AdjustedQC_Array[i] = fallback_flag

                if is_pres_missing_or_bad or is_broken:
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"4" if is_broken else b"9"
                elif matched_row is not None and fluo_csv_col is not None and pd.notna(matched_row.get(fluo_csv_col)):
                    fluo_csv_val = matched_row.get(fluo_csv_col)
                    fluo_adj_qc_val = matched_row.get(fluo_qc_csv_col, "1")
                    fluo_adj_qc_str = str(int(fluo_adj_qc_val)) if pd.notna(fluo_adj_qc_val) else "1"

                    if raw_fluo_qc_char in ["4"] or fluo_adj_qc_str in ["4"] or pd.isna(fluo_csv_val) or fluo_csv_val == 99999.0:
                        FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                        FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                        FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"4" if (raw_fluo_qc_char == "4" or fluo_adj_qc_str == "4") else b"1"
                    else:
                        FLUORESCENCE_CHLA_Adjusted_Array[i] = np.float32(fluo_csv_val)
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                        FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = False
                        FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = False
                        valid_fluo_flag = fluo_adj_qc_str.encode("utf-8") if fluo_adj_qc_str in ["1", "2", "3"] else b"1"
                        FLUORESCENCE_CHLA_AdjustedQC_Array[i] = valid_fluo_flag
                elif fluo_data_arr is not None and not (pd.isna(fluo_data_arr[i]) or fluo_data_arr[i] == 99999.0):
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = np.float32(fluo_data_arr[i])
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = False
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = False
                    fallback_flag = raw_fluo_qc_char.encode("utf-8") if raw_fluo_qc_char in ["1", "2", "3"] else b"1"
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = fallback_flag
                else:
                    FLUORESCENCE_CHLA_Adjusted_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = 99999.0
                    FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = True
                    FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = True
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"9"

        bgc_file.variables["CHLA_ADJUSTED"][prof, :] = CHLA_Adjusted_Array
        bgc_file.variables["CHLA_ADJUSTED_QC"][prof, :] = CHLA_AdjustedQC_Array
        bgc_file.variables["CHLA_ADJUSTED_ERROR"][prof, :] = CHLA_Adjusted_ERROR_Array

        for var_prefix in ["FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]:
            adj_var = f"{var_prefix}_ADJUSTED"
            adj_qc_var = f"{var_prefix}_ADJUSTED_QC"
            adj_err_var = f"{var_prefix}_ADJUSTED_ERROR"

            if adj_var in bgc_file.variables:
                bgc_file.variables[adj_var][prof, :] = FLUORESCENCE_CHLA_Adjusted_Array
            if adj_qc_var in bgc_file.variables:
                bgc_file.variables[adj_qc_var][prof, :] = FLUORESCENCE_CHLA_AdjustedQC_Array
            if adj_err_var in bgc_file.variables:
                bgc_file.variables[adj_err_var][prof, :] = FLUORESCENCE_CHLA_Adjusted_ERROR_Array


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
        raw_bbp_var = bgc_file.variables["BBP700"][prof, :] if (has_bbp and bgc_file.variables["BBP700"].ndim > 1) else (bgc_file.variables["BBP700"][:] if has_bbp else None)
        raw_bbp_qc_var = bgc_file.variables["BBP700_QC"][prof, :] if (has_bbp_qc and bgc_file.variables["BBP700_QC"].ndim > 1) else (bgc_file.variables["BBP700"][:] if has_bbp else None)

        bbp_data_arr = np.array(raw_bbp_var, copy=True) if raw_bbp_var is not None else None
        bbp_qc_arr = np.array(raw_bbp_qc_var, copy=True) if raw_bbp_qc_var is not None else None

        BBP700_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        BBP700_Adjusted_Array[:] = 99999.0
        BBP700_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")
        BBP700_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        BBP700_Adjusted_ERROR_Array[:] = 99999.0

        nc_pres_all = get_pres_array_for_profile(bgc_file, prof)

        if mode in ["D", "A"]:
            for i in range(n_levels):
                nc_pres = nc_pres_all[i]
                is_pres_missing_or_bad = np.isnan(nc_pres) or nc_pres == 99999.0
                is_raw_bbp_missing = bbp_data_arr is None or pd.isna(bbp_data_arr[i]) or bbp_data_arr[i] == 99999.0

                matched_row = match_csv_row_by_pres(cycle_df, nc_pres) if not (is_broken or is_pres_missing_or_bad) else None
                raw_qc_char = bbp_qc_arr[i].decode("utf-8") if (bbp_qc_arr is not None and isinstance(bbp_qc_arr[i], bytes)) else str(bbp_qc_arr[i]) if bbp_qc_arr is not None else "1"

                if is_pres_missing_or_bad or is_raw_bbp_missing or is_broken:
                    BBP700_Adjusted_Array[i] = 99999.0
                    BBP700_Adjusted_ERROR_Array[i] = 99999.0
                    BBP700_Adjusted_Array.mask[i] = True
                    BBP700_Adjusted_ERROR_Array.mask[i] = True
                    BBP700_AdjustedQC_Array[i] = b"4" if is_broken else b"9"
                elif matched_row is not None:
                    csv_qc = matched_row.get("BBP700_FINAL_QC", "1")
                    qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "1"
                    raw_bbp = matched_row.get("BBP700_FINAL", np.nan)

                    if raw_qc_char in ["4"] or qc_str in ["4"] or pd.isna(raw_bbp) or raw_bbp == 99999.0:
                        BBP700_Adjusted_Array[i] = 99999.0
                        BBP700_Adjusted_ERROR_Array[i] = 99999.0
                        BBP700_Adjusted_Array.mask[i] = True
                        BBP700_Adjusted_ERROR_Array.mask[i] = True
                        BBP700_AdjustedQC_Array[i] = b"4" if (raw_qc_char == "4" or qc_str == "4") else b"1"
                    else:
                        BBP700_Adjusted_Array[i] = np.float32(raw_bbp)
                        BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                        BBP700_Adjusted_Array.mask[i] = False
                        BBP700_Adjusted_ERROR_Array.mask[i] = False
                        valid_bbp_flag = qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1"
                        BBP700_AdjustedQC_Array[i] = valid_bbp_flag
                else:
                    BBP700_Adjusted_Array[i] = np.float32(bbp_data_arr[i])
                    BBP700_Adjusted_Array.mask[i] = False
                    BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                    BBP700_Adjusted_ERROR_Array.mask[i] = False
                    fallback_flag = raw_qc_char.encode("utf-8") if raw_qc_char in ["1", "2", "3"] else b"1"
                    BBP700_AdjustedQC_Array[i] = fallback_flag

        bgc_file.variables["BBP700_ADJUSTED"][prof, :] = BBP700_Adjusted_Array
        if "BBP700_ADJUSTED_QC" in bgc_file.variables:
            bgc_file.variables["BBP700_ADJUSTED_QC"][prof, :] = BBP700_AdjustedQC_Array
        if "BBP700_ADJUSTED_ERROR" in bgc_file.variables:
            bgc_file.variables["BBP700_ADJUSTED_ERROR"][prof, :] = BBP700_Adjusted_ERROR_Array


def write_DOXY_slope_drift(ds, profile_idx, float_df, target_cycle, WMOfloatid=None):
    if WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS:
        return

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
                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                for c in p_bytes
            ]).strip()

            if param_str == "DOXY":
                if "SCIENTIFIC_CALIB_COEFFICIENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof, 0, j, :] = char_coef
                if "SCIENTIFIC_CALIB_EQUATION" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_EQUATION"][iprof, 0, j, :] = char_equ
                if "SCIENTIFIC_CALIB_COMMENT" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_COMMENT"][iprof, 0, j, :] = char_com
                if "SCIENTIFIC_CALIB_DATE" in ds.variables:
                    ds.variables["SCIENTIFIC_CALIB_DATE"][iprof, 0, j, :] = char_date


def write_DOXY_from_csv(ds, float_df, target_cycle, WMOfloatid=None):
    is_excluded = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

    var_names = ds.variables.keys()
    n_prof = ds.dimensions["N_PROF"].size
    n_levels = ds.dimensions["N_LEVELS"].size

    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()

    for iprof in range(n_prof):
        mode = get_param_data_mode(ds, iprof, "DOXY")

        has_doxy = "DOXY" in ds.variables
        has_doxy_qc = "DOXY_QC" in ds.variables
        raw_doxy_var = ds.variables["DOXY"][iprof, :] if (has_doxy and ds.variables["DOXY"].ndim > 1) else (ds.variables["DOXY"][:] if has_doxy else None)
        raw_doxy_qc_var = ds.variables["DOXY_QC"][iprof, :] if (has_doxy_qc and ds.variables["DOXY_QC"].ndim > 1) else (ds.variables["DOXY"][:] if has_doxy else None)

        doxy_data_arr = np.array(raw_doxy_var, copy=True) if raw_doxy_var is not None else None
        doxy_qc_arr = np.array(raw_doxy_qc_var, copy=True) if raw_doxy_qc_var is not None else None

        if is_excluded:
            if "PROFILE_DOXY_QC" in ds.variables and doxy_qc_arr is not None:
                prof_grade = get_profile_qc_grade(doxy_qc_arr, parameter_present=has_valid_measurements(ds, "DOXY", iprof))
                grade_bytes = prof_grade.encode("utf-8")
                if ds.variables["PROFILE_DOXY_QC"].ndim == 1:
                    ds.variables["PROFILE_DOXY_QC"][iprof] = grade_bytes
                elif ds.variables["PROFILE_DOXY_QC"].ndim == 2:
                    ds.variables["PROFILE_DOXY_QC"][iprof, 0] = grade_bytes
            continue

        if mode == "R":
            continue

        nc_pres_all = get_pres_array_for_profile(ds, iprof)

        DOXY_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        DOXY_Adjusted_Array[:] = 99999.0
        DOXY_Adjusted_Error_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        DOXY_Adjusted_Error_Array[:] = 99999.0
        DOXY_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

        for i in range(n_levels):
            nc_pres = nc_pres_all[i]
            is_pres_missing_or_bad = np.isnan(nc_pres) or nc_pres == 99999.0
            has_valid_raw_data = doxy_data_arr is not None and not pd.isna(doxy_data_arr[i]) and doxy_data_arr[i] != 99999.0
            raw_qc_char = doxy_qc_arr[i].decode("utf-8") if (doxy_qc_arr is not None and isinstance(doxy_qc_arr[i], bytes)) else str(doxy_qc_arr[i]) if doxy_qc_arr is not None else "1"

            if is_pres_missing_or_bad or not has_valid_raw_data:
                DOXY_Adjusted_Array[i] = 99999.0
                DOXY_Adjusted_Array.mask[i] = True
                DOXY_Adjusted_Error_Array[i] = 99999.0
                DOXY_Adjusted_Error_Array.mask[i] = True
                DOXY_AdjustedQC_Array[i] = b"9"
                continue

            matched_row = match_csv_row_by_pres(cycle_df, nc_pres)

            if matched_row is not None:
                csv_qc = matched_row.get("DOXY_FINAL_QC", "1")
                qc_str = str(int(csv_qc)) if pd.notna(csv_qc) else "1"
                raw_doxy_final = matched_row.get("DOXY_FINAL")

                if raw_qc_char in ["4"] or qc_str in ["4"] or pd.isna(raw_doxy_final) or raw_doxy_final == 99999.0:
                    DOXY_Adjusted_Array[i] = 99999.0
                    DOXY_Adjusted_Error_Array[i] = 99999.0
                    DOXY_Adjusted_Array.mask[i] = True
                    DOXY_Adjusted_Error_Array.mask[i] = True
                    DOXY_AdjustedQC_Array[i] = b"4" if (raw_qc_char == "4" or qc_str == "4") else b"1"
                else:
                    DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
                    DOXY_Adjusted_Array.mask[i] = False
                    valid_flag = raw_qc_char.encode("utf-8") if raw_qc_char in ["1", "2", "3"] else (qc_str.encode("utf-8") if qc_str in ["1", "2", "3"] else b"1")
                    DOXY_AdjustedQC_Array[i] = valid_flag

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
                fallback_qc = raw_qc_char.encode("utf-8") if raw_qc_char in ["1", "2", "3"] else b"1"
                DOXY_AdjustedQC_Array[i] = fallback_qc

        if "DOXY_ADJUSTED" in var_names:
            doxy_adj_var = ds.variables["DOXY_ADJUSTED"]
            doxy_adj_qc_var = ds.variables.get("DOXY_ADJUSTED_QC")

            if doxy_adj_var.ndim > 1:
                doxy_adj_var[iprof, :] = DOXY_Adjusted_Array
                if doxy_adj_qc_var is not None:
                    doxy_adj_qc_var[iprof, :] = DOXY_AdjustedQC_Array
            else:
                doxy_adj_var[:] = DOXY_Adjusted_Array
                if doxy_adj_qc_var is not None:
                    doxy_adj_qc_var[:] = DOXY_AdjustedQC_Array

        if "DOXY_ADJUSTED_ERROR" in var_names:
            doxy_adj_err_var = ds.variables["DOXY_ADJUSTED_ERROR"]
            if doxy_adj_err_var.ndim > 1:
                doxy_adj_err_var[iprof, :] = DOXY_Adjusted_Error_Array
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
# SECTION 3: Plotly DMQC Status Plot Generation
# ==============================================================================

data_type_color_map = {"D": "#647ed5", "A": "#c07ab8", "R": "#ff4d6d"}
label_data_type_color_map = {
    "Delayed Mode (D)": "#647ed5",
    "Adjusted Mode (A)": "#c07ab8",
    "Real Time Mode (R)": "#ff4d6d",
}
data_mode_rank_map = {"R": 1, "A": 2, "D": 3}

position_system_color_map = {"GPS": "#af43f0", "ARGOS": "#00bcd4"}
position_qc_color_map = {
    "1": "#1b7837", "2": "#ffffbf", "3": "#fdae61",
    "4": "#d73027", "5": "#404040", "8": "#7fbf7b", "9": "#808080"
}
label_qc_color_map = {
    "1 - Good data": "#1b7837",
    "2 - Probably good": "#ffffbf",
    "3 - Probably bad": "#fdae61",
    "4 - Bad data": "#d73027",
    "5 - Value changed": "#404040",
    "8 - Interpolated": "#7fbf7b",
    "9 - Missing": "#808080",
}

param_qc_color_map = {
    "A": "#1b7837", "B": "#a6dba0", "C": "#ffffbf",
    "D": "#fdae61", "E": "#f46d43", "F": "#d73027", " ": "#ffffff"
}
label_param_qc_color_map = {
    "Grade A (100% Good)": "#1b7837",
    "Grade B (75-99% Good)": "#a6dba0",
    "Grade C (50-75% Good)": "#ffffbf",
    "Grade D (25-50% Good)": "#fdae61",
    "Grade E (0-25% Good)": "#f46d43",
    "Grade F (0% Good)": "#d73027",
}


def get_subsample_prof_idx(temp_list, n):
    if len(temp_list) <= n:
        return temp_list
    idx = np.round(np.linspace(0, len(temp_list) - 1, n)).astype(int)
    return [temp_list[i] for i in idx]


def add_Data_mode_sp(fig, float_map_inst, std_fig_width=1000):
    x_nprof = []
    y_dataType = []
    color_dataType = []
    n_profs = len(float_map_inst["NPROF_list"])

    for n_prof in float_map_inst["NPROF_list"]:
        x_nprof.append(str(n_prof + 1))
        y_dataType.append(" ")
        df_sub = float_map_inst["NPROF_NPARAM_df"]
        modes = df_sub[df_sub["N_PROF"] == n_prof]["PARAMETER_DATA_MODE"].tolist()
        best_mode = min(modes, key=lambda m: data_mode_rank_map.get(m, 99)) if modes else "R"
        color_dataType.append(data_type_color_map.get(best_mode, "#FFFFFF"))

    x_ticks = get_subsample_prof_idx(x_nprof, min(7, n_profs))
    fig.update_xaxes(
        tickmode="array",
        tickvals=x_ticks,
        ticktext=[str(x) for x in x_ticks],
        row="all",
        col=1,
    )

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(n_profs, 1), 1.0)
    fig.add_trace(
        go.Scatter(
            x=x_nprof,
            y=y_dataType,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color_dataType,
            marker_line_color=color_dataType,
            marker_line_width=marker_line_w,
            marker_size=30,
            orientation="h",
            showlegend=False,
        ),
        row=1,
        col=1,
    )


def add_positioning_system_sp(fig, float_map_inst, std_fig_width=1000):
    x1_nprof = []
    y1_positioning_system = []
    y1_position_qc = []
    color1_positioning_system = []
    color1_position_qc = []
    n_profs = len(float_map_inst["NPROF_list"])

    for n_prof in float_map_inst["NPROF_list"]:
        pos_sys = float_map_inst["NPROF_df"]["POSITIONING_SYSTEM"].get(n_prof, "GPS").strip()
        pos_qc = float_map_inst["NPROF_df"]["POSITION_QC"].get(n_prof, "1").strip()
        x1_nprof.append(str(n_prof + 1))
        y1_positioning_system.append(" ")
        y1_position_qc.append("  ")
        color1_positioning_system.append(position_system_color_map.get(pos_sys, "#af43f0"))
        color1_position_qc.append(position_qc_color_map.get(pos_qc, "#808080"))

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(n_profs, 1), 1.0)
    fig.add_trace(
        go.Scatter(
            x=x1_nprof,
            y=y1_position_qc,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color1_position_qc,
            marker_line_color=color1_position_qc,
            marker_line_width=marker_line_w,
            marker_size=15,
            orientation="h",
            showlegend=False,
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x1_nprof,
            y=y1_positioning_system,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color1_positioning_system,
            marker_line_color=color1_positioning_system,
            marker_line_width=marker_line_w,
            marker_size=15,
            orientation="h",
            showlegend=False,
        ),
        row=2,
        col=1,
    )


def add_profile_depth_sp(fig, float_map_inst, std_fig_width=1000):
    x2_nprof = []
    y2_depth = []
    color_depth = "#41AACE"

    for n_prof in float_map_inst["NPROF_list"]:
        x2_nprof.append(str(n_prof + 1))
        df_levels = float_map_inst["NPROF_NLEVEL_df"]
        p_slice = df_levels[df_levels["N_PROF"] == n_prof]["PRES_ADJUSTED"]
        max_d = np.nanmax(p_slice) if len(p_slice) > 0 and not np.all(np.isnan(p_slice)) else 0.0
        y2_depth.append(max_d)

    max_prof_depth = np.max(y2_depth) if y2_depth and not np.all(np.isnan(y2_depth)) else 2000.0
    fig.add_trace(
        go.Bar(
            x=x2_nprof,
            y=y2_depth,
            orientation="v",
            marker=dict(color=color_depth, line=dict(color=color_depth), cornerradius="0%"),
            showlegend=False,
        ),
        row=3,
        col=1,
    )

    fig.update_yaxes(range=[max_prof_depth, 0], visible=False, row=3, col=1)
    return np.rint(max_prof_depth)


def add_param_QC_sp(fig, float_map_inst, PARAM, subplotNum, std_fig_width=1000):
    x3_nprof = []
    y3_p1 = []  # QC Flag track (lower)
    y3_p2 = []  # Data Mode track (upper)
    color3_qc = []
    color3_pdm = []
    hover_texts_qc = []
    hover_texts_pdm = []
    n_profs = len(float_map_inst["NPROF_list"])

    for n_prof in float_map_inst["NPROF_list"]:
        cyc_str = str(n_prof + 1)
        x3_nprof.append(cyc_str)
        y3_p1.append("Grade")
        y3_p2.append("Mode")

        flag = float_map_inst["NPROF_df"].get(f"PROFILE_{PARAM}_QC", {}).get(n_prof, " ")
        color3_qc.append(param_qc_color_map.get(flag, "#FFFFFF"))
        hover_texts_qc.append(f"Cycle {cyc_str}<br>Param: {PARAM}<br>QC Grade: {flag}")

        pdm = float_map_inst["PDM_MAP"].get(PARAM, {}).get(n_prof, "R")
        color3_pdm.append(data_type_color_map.get(pdm, "#FFFFFF"))
        hover_texts_pdm.append(f"Cycle {cyc_str}<br>Param: {PARAM}<br>Data Mode: {pdm}")

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(n_profs, 1), 1.0)

    fig.add_trace(
        go.Scatter(
            x=x3_nprof,
            y=y3_p1,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color3_qc,
            marker_line_color=color3_qc,
            marker_line_width=marker_line_w,
            marker_size=16,
            hoverinfo="text",
            hovertext=hover_texts_qc,
            orientation="h",
            showlegend=False,
        ),
        row=subplotNum,
        col=1,
    )

    fig.add_trace(
        go.Scatter(
            x=x3_nprof,
            y=y3_p2,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color3_pdm,
            marker_line_color=color3_pdm,
            marker_line_width=marker_line_w,
            marker_size=16,
            hoverinfo="text",
            hovertext=hover_texts_pdm,
            orientation="h",
            showlegend=False,
        ),
        row=subplotNum,
        col=1,
    )


def update_layout_DMQC_STATUS(fig, float_map_inst, max_prof_depth, subplotList, num_rows, std_fig_width=1000):
    n_profs = len(float_map_inst["NPROF_list"])
    fig.update_xaxes(range=[-0.05 * n_profs, n_profs * 1.05], dtick=n_profs * 0.33)

    fig.update_layout(
        autosize=False,
        height=max(850, 45 * num_rows + 200),
        width=std_fig_width,
        bargap=0,
        bargroupgap=0,
        font=dict(color=TEXT_COLOR),
        title=dict(
            text=f"WMOID {float_map_inst['wmoID_str']} - BD File Verification & DMQC Status",
            font=dict(color=TEXT_COLOR, size=22),
            x=0.5,
            xanchor="center",
            yanchor="top",
        ),
        annotations=[
            dict(
                text="Data Type (R/A/D)",
                xref="paper", yref="paper",
                x=0, y=0.9,
                xanchor="right", yanchor="top",
                showarrow=False, font=dict(size=11),
            ),
            dict(
                text="Positioning System",
                xref="paper", yref="paper",
                x=0, y=sum(fig.layout["yaxis2"].domain) / 2,
                xanchor="right", yanchor="bottom",
                showarrow=False, font=dict(size=11),
            ),
            dict(
                text="Position QC",
                xref="paper", yref="paper",
                x=0, y=sum(fig.layout["yaxis2"].domain) / 2,
                xanchor="right", yanchor="top",
                showarrow=False, font=dict(size=11),
            ),
            dict(
                text="Profile Depth",
                xref="paper", yref="paper",
                x=0, y=sum(fig.layout["yaxis3"].domain) / 2,
                xanchor="right", yanchor="bottom",
                showarrow=False, font=dict(size=11),
            ),
            dict(
                text=f"Max: ~{max_prof_depth:.0f} dbar",
                xref="paper", yref="paper",
                x=0, y=sum(fig.layout["yaxis3"].domain) / 2,
                xanchor="right", yanchor="top",
                showarrow=False,
                font=dict(size=11, color="#41AACE"),
            ),
        ],
        plot_bgcolor=BG_COLOR,
        paper_bgcolor=BG_COLOR,
        margin=dict(l=175, r=20, t=110, b=80),
    )

    for mode_name, hex_color in label_data_type_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{mode_name}",
                legendgroup="DataMode",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="<b>Data Mode</b>", selector=dict(legendgroup="DataMode"))

    for grade_name, hex_color in label_param_qc_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{grade_name}",
                legendgroup="ProfileQualityFlag",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="<b>QC Profile Grade</b>", selector=dict(legendgroup="ProfileQualityFlag"))

    for qc_name, hex_color in label_qc_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{qc_name}",
                legendgroup="PositioningSystemQC",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="<b>Position QC</b>", selector=dict(legendgroup="PositioningSystemQC"))

    for PARAM, subplotidx in zip(float_map_inst["param_list"], subplotList):
        tempPARAM = "PH" if PARAM == "PH_IN_SITU_TOTAL" else PARAM
        domain_center = sum(fig.layout[f"yaxis{subplotidx}"].domain) / 2

        fig.add_annotation(
            text=f"{tempPARAM} Mode",
            xref="paper", yref="paper",
            x=0, y=domain_center + 0.007,
            xanchor="right", yanchor="bottom",
            showarrow=False, font=dict(size=11),
        )
        fig.add_annotation(
            text=f"Profile {tempPARAM} QC",
            xref="paper", yref="paper",
            x=0, y=domain_center - 0.007,
            xanchor="right", yanchor="top",
            showarrow=False, font=dict(size=11),
        )

    total_slots = 0
    d_mode_count = 0
    for p in float_map_inst["param_list"]:
        for m in float_map_inst["PDM_MAP"].get(p, {}).values():
            total_slots += 1
            if m == "D":
                d_mode_count += 1
    pct_d = (d_mode_count / total_slots * 100) if total_slots > 0 else 0.0

    float_status = float_map_inst.get("activity_status", "Active")
    status_annotation = (
        f'<span style="font-size: 16px; font-weight: bold; color: #16a085;">BD Filler: {pct_d:.0f}% D-Mode</span><br>'
        f'<span style="font-size: 14px; color: #27ae60;">Status: {float_status}</span>'
    )

    fig.add_annotation(
        text=status_annotation,
        xref="paper", yref="paper",
        x=1, y=1,
        xanchor="right", yanchor="top",
        showarrow=False,
        yshift=75, xshift=140,
    )

    curr_date_str = dt.now().strftime("%Y-%m-%d")
    fig.add_annotation(
        text=f"Figure Created: {curr_date_str}",
        xref="paper", yref="paper",
        x=0, y=0,
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11),
        xshift=-150, yshift=-60,
    )

    prof_dates = float_map_inst.get("JULD_list", [])
    if prof_dates:
        first_prof = float_map_inst["NPROF_list"][0]
        last_prof = float_map_inst["NPROF_list"][-1]
        first_d = prof_dates[0].strftime("%Y-%m-%d") if isinstance(prof_dates[0], dt) else str(prof_dates[0])
        last_d = prof_dates[-1].strftime("%Y-%m-%d") if isinstance(prof_dates[-1], dt) else str(prof_dates[-1])

        fig.add_annotation(
            text=f"Cycle {first_prof + 1}:  {first_d}",
            xref="paper", yref="paper",
            x=0, y=0,
            xanchor="left", yanchor="bottom",
            showarrow=False, font=dict(size=11),
            xshift=0, yshift=-50,
        )
        fig.add_annotation(
            text=f"Cycle {last_prof + 1}:  {last_d}",
            xref="paper", yref="paper",
            x=1, y=0,
            xanchor="right", yanchor="bottom",
            showarrow=False, font=dict(size=11),
            xshift=0, yshift=-50,
        )


def make_DMQC_STATUS_fig(float_map_inst, output_filepath, show_grid=True):
    num_rows = 3 + len(float_map_inst["param_list"])
    row_heights = [1.0 / num_rows] * num_rows

    fig = make_subplots(
        rows=num_rows,
        cols=1,
        shared_xaxes=True,
        x_title="Cycle Number",
        shared_yaxes=False,
        vertical_spacing=0.015,
        row_heights=row_heights,
        print_grid=show_grid,
    )

    add_Data_mode_sp(fig, float_map_inst)
    add_positioning_system_sp(fig, float_map_inst)
    max_prof_depth = add_profile_depth_sp(fig, float_map_inst)

    subplotNumList = []
    for i, PARAM in enumerate(float_map_inst["param_list"]):
        subplotNum = i + 4
        try:
            subplotNumList.append(subplotNum)
            add_param_QC_sp(fig, float_map_inst, PARAM, subplotNum)
        except KeyError:
            print(f"Warning: Missing data for {PARAM}. Skipping subplot.")
            continue

    update_layout_DMQC_STATUS(fig, float_map_inst, max_prof_depth, subplotNumList, num_rows)

    logo_path = "/a2/ARGO_OPERATIONS/float_vis/drive/pictures/noaa_25_logo3.png"
    if os.path.exists(logo_path):
        try:
            logo = Image.open(logo_path)
            fig.add_layout_image(
                dict(
                    source=logo,
                    xref="paper", yref="paper",
                    x=0, y=1,
                    opacity=1,
                    xanchor="right", yanchor="bottom",
                    sizex=0.1, sizey=0.1,
                    layer="above",
                )
            )
        except Exception:
            pass

    fig.write_image(output_filepath, format="png", width=1050, height=max(850, 45 * num_rows + 200), scale=2)
    return fig


def generate_dmqc_status_plot(WMOfloatid, processed_bd_files, output_dir):
    """
    Harvests data dynamically across BD files and renders the DMQC status figure.
    """
    if not processed_bd_files:
        return

    plot_out_dir = os.path.join(output_dir, "DMQC_status_plot")
    os.makedirs(plot_out_dir, exist_ok=True)
    out_png = os.path.join(plot_out_dir, f"{WMOfloatid}_DMQC_status.png")

    candidate_params = ["PRES", "TEMP", "PSAL", "DOXY", "BBP700", "CHLA", "PH_IN_SITU_TOTAL", "NITRATE", "CHLA_FLUORESCENCE"]
    present_params = set()

    nprof_list = []
    juld_list = []
    nprof_df_data = {"POSITIONING_SYSTEM": {}, "POSITION_QC": {}}
    pdm_map = {p: {} for p in candidate_params}
    nprof_nparam_records = []
    nprof_nlevel_records = []

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            cycle_num = int(fpath[-6:-3])
            prof_idx = cycle_num - 1
            nprof_list.append(prof_idx)

            c_date = dt.now()
            if "JULD" in ds.variables and ds.variables["JULD"][:].size > 0:
                juld_val = ds.variables["JULD"][0]
                if not np.ma.is_masked(juld_val):
                    ref_date = dt(1950, 1, 1, tzinfo=timezone.utc)
                    c_date = ref_date + datetime.timedelta(days=float(juld_val))
            juld_list.append(c_date)

            pos_sys = "GPS"
            if "POSITIONING_SYSTEM" in ds.variables:
                raw_sys = ds.variables["POSITIONING_SYSTEM"][0]
                pos_sys = raw_sys.decode("utf-8") if isinstance(raw_sys, bytes) else str(raw_sys)
            nprof_df_data["POSITIONING_SYSTEM"][prof_idx] = pos_sys

            pos_qc = "1"
            if "POSITION_QC" in ds.variables:
                raw_qc = ds.variables["POSITION_QC"][0]
                pos_qc = raw_qc.decode("utf-8") if isinstance(raw_qc, bytes) else str(raw_qc)
            nprof_df_data["POSITION_QC"][prof_idx] = pos_qc

            p_pres_idx = detect_parameter_profile(ds, "PRES")
            pres_vals = get_pres_array_for_profile(ds, p_pres_idx)
            for p_v in pres_vals:
                if not np.isnan(p_v) and p_v != 99999.0:
                    nprof_nlevel_records.append({"N_PROF": prof_idx, "PRES_ADJUSTED": p_v})

            for p in candidate_params:
                p_prof = detect_parameter_profile(ds, p)

                var_exists = (
                    p in ds.variables or f"{p}_ADJUSTED" in ds.variables or f"PROFILE_{p}_QC" in ds.variables
                )
                if not var_exists and p == "CHLA_FLUORESCENCE":
                    var_exists = "FLUORESCENCE_CHLA" in ds.variables or "PROFILE_FLUORESCENCE_CHLA_QC" in ds.variables

                if var_exists:
                    present_params.add(p)

                mode_val = get_param_data_mode(ds, p_prof, p)
                pdm_map[p][prof_idx] = mode_val
                nprof_nparam_records.append({"N_PROF": prof_idx, "PARAMETER_DATA_MODE": mode_val})

                target_qc = (
                    f"{p}_ADJUSTED_QC"
                    if mode_val in ["D", "A"] and f"{p}_ADJUSTED_QC" in ds.variables
                    else f"{p}_QC"
                )
                has_data = has_valid_measurements(ds, p, p_prof)

                grade = " "
                if target_qc in ds.variables:
                    qc_slice = ds.variables[target_qc][p_prof, :] if ds.variables[target_qc].ndim > 1 else ds.variables[target_qc][:]
                    grade = get_profile_qc_grade(qc_slice, parameter_present=has_data)
                else:
                    qc_var_names = [f"PROFILE_{p}_QC"]
                    if p == "CHLA_FLUORESCENCE":
                        qc_var_names.append("PROFILE_FLUORESCENCE_CHLA_QC")
                    for qvn in qc_var_names:
                        if qvn in ds.variables:
                            g = ds.variables[qvn][p_prof] if ds.variables[qvn].ndim > 0 else ds.variables[qvn][:]
                            if isinstance(g, np.ndarray) and g.ndim > 0:
                                g = g[0]
                            grade = g.decode("utf-8", errors="ignore") if isinstance(g, bytes) else str(g)
                            if grade.strip() != "":
                                break

                if p == "PRES":
                    if "PRES_QC" in ds.variables:
                        pres_qc_slice = ds.variables["PRES_QC"][p_prof, :] if ds.variables["PRES_QC"].ndim > 1 else ds.variables["PRES_QC"][:]
                        grade = get_profile_qc_grade(pres_qc_slice, parameter_present=True)

                qc_col_key = f"PROFILE_{p}_QC"
                if qc_col_key not in nprof_df_data:
                    nprof_df_data[qc_col_key] = {}
                nprof_df_data[qc_col_key][prof_idx] = grade

            ds.close()
        except Exception:
            traceback.print_exc()

    if not nprof_list:
        return

    ordered_param_list = [p for p in candidate_params if p in present_params]

    float_map_inst = {
        "wmoID_str": str(WMOfloatid),
        "activity_status": "Active",
        "param_list": ordered_param_list,
        "NPROF_list": nprof_list,
        "JULD_list": juld_list,
        "NPROF_df": nprof_df_data,
        "PDM_MAP": pdm_map,
        "NPROF_NPARAM_df": pd.DataFrame(nprof_nparam_records),
        "NPROF_NLEVEL_df": pd.DataFrame(nprof_nlevel_records),
    }

    try:
        make_DMQC_STATUS_fig(float_map_inst, output_filepath=out_png)
        print(f"[{WMOfloatid}] Final Output BD Plotly DMQC Status Plot saved to: {out_png}")
    except Exception:
        traceback.print_exc()


# ==============================================================================
# SECTION 3.5: Section Contour & QC Summary Plots (R_mode & D_mode)
# ==============================================================================

# QC Flag color palette
QC_COLOR_MAP = {
    "0": "#b3b3b3",  # gray70
    "1": "#33A02C",  # Green
    "2": "#B2DF8A",  # Light green
    "3": "#FDBF6F",  # Orange
    "4": "#E31A1C",  # Red
    "5": "#1F78B4",  # Blue
    "8": "#CAB2D6",  # Purple
    "9": "#000000",  # Black
}


def create_interpolated_grid(df, val_col, grid_nx=100, grid_ny=100):
    clean_df = df.dropna(subset=[val_col, "CYCLE_NUMBER", "PRES"]).copy()
    clean_df = clean_df[clean_df[val_col] != 99999.0]

    if clean_df.empty:
        return None

    grouped = clean_df.groupby(["CYCLE_NUMBER", "PRES"])[val_col].mean().reset_index()

    if len(grouped) < 5 or grouped["CYCLE_NUMBER"].nunique() < 2 or grouped["PRES"].nunique() < 2:
        return None

    x = grouped["CYCLE_NUMBER"].values
    y = grouped["PRES"].values
    z = grouped[val_col].values

    if np.nanmax(z) == np.nanmin(z):
        return None

    xi = np.linspace(np.nanmin(x), np.nanmax(x), grid_nx)
    yi = np.linspace(np.nanmin(y), np.nanmax(y), grid_ny)
    xi_grid, yi_grid = np.meshgrid(xi, yi)

    zi_grid = griddata((x, y), z, (xi_grid, yi_grid), method="linear")

    if np.isnan(zi_grid).any():
        zi_nearest = griddata((x, y), z, (xi_grid, yi_grid), method="nearest")
        zi_grid = np.where(np.isnan(zi_grid), zi_nearest, zi_grid)

    return xi_grid, yi_grid, zi_grid


def generate_section_qc_plot(
    df_base, var_name, qc_var_name, title_prefix, file_suffix,
    output_dir, WMOfloatid, col_palette="viridis", unit_str="", is_chla=False
):
    if var_name not in df_base.columns:
        return

    clean_var = df_base[var_name].replace(99999.0, np.nan)
    if clean_var.dropna().empty:
        return

    df_qc = df_base.copy()
    df_qc["VAR_VAL"] = clean_var
    df_qc = df_qc.dropna(subset=["VAR_VAL", "PRES", "CYCLE_NUMBER"]).copy()

    if df_qc.empty:
        return

    if qc_var_name in df_qc.columns:
        df_qc["QC_VAL"] = (
            df_qc[qc_var_name]
            .astype(str)
            .str.extract(r"([0-9])")[0]
            .fillna("1")
        )
    else:
        df_qc["QC_VAL"] = "1"

    interp_res = create_interpolated_grid(df_qc, "VAR_VAL")
    max_pres = np.nanmax(df_base["PRES"].replace(99999.0, np.nan))
    if np.isnan(max_pres) or max_pres <= 0:
        max_pres = 2000.0

    fig, (ax_data, ax_qc) = plt.subplots(
        nrows=2, ncols=1, figsize=(12, 8), dpi=150, sharex=True,
        gridspec_kw={"height_ratios": [1, 1]}
    )
    fig.patch.set_facecolor("white")

    divider_data = make_axes_locatable(ax_data)
    cax_data = divider_data.append_axes("right", size="3%", pad=0.15)

    divider_qc = make_axes_locatable(ax_qc)
    cax_qc = divider_qc.append_axes("right", size="3%", pad=0.15)
    cax_qc.axis("off")

    # --- TOP PANEL: Scientific Contour or Tile Plot ---
    if interp_res is not None:
        xi_grid, yi_grid, zi_grid = interp_res

        if is_chla:
            chla_levels = [-0.1, 0.01, 0.03, 0.05, 0.1, 0.2, 0.35, 0.5, 0.8, 1.2, 2.0, 5.0]
            cmap = plt.get_cmap("turbo").copy()
            norm = mcolors.BoundaryNorm(boundaries=chla_levels, ncolors=cmap.N, extend="both")
            cf = ax_data.contourf(xi_grid, yi_grid, zi_grid, levels=chla_levels, cmap=cmap, norm=norm, extend="both")
            cbar = fig.colorbar(cf, cax=cax_data)
            cbar.set_label(unit_str, fontsize=10)
        else:
            cmap_name = col_palette if col_palette in plt.colormaps() else "viridis"
            cf = ax_data.contourf(xi_grid, yi_grid, zi_grid, levels=15, cmap=cmap_name)
            cbar = fig.colorbar(cf, cax=cax_data)
            cbar.set_label(unit_str, fontsize=10)

        ax_data.set_title(f"Float {WMOfloatid} - {title_prefix} Section Profile", fontsize=12, fontweight="bold")
    else:
        sc = ax_data.scatter(df_qc["CYCLE_NUMBER"], df_qc["PRES"], c=df_qc["VAR_VAL"], cmap=col_palette, s=15, alpha=0.8)
        cbar = fig.colorbar(sc, cax=cax_data)
        cbar.set_label(unit_str, fontsize=10)
        ax_data.set_title(f"Float {WMOfloatid} - {title_prefix} Section Profile", fontsize=12, fontweight="bold")

    ax_data.set_ylim(max_pres, 0)
    ax_data.set_ylabel("Pressure (dbar)", fontsize=11)
    ax_data.grid(True, linestyle="--", alpha=0.4)

    # --- BOTTOM PANEL: Discrete QC Flag Map ---
    unique_qcs = sorted(df_qc["QC_VAL"].unique(), key=lambda q: int(q) if q.isdigit() else 99)
    for q_flag in unique_qcs:
        sub = df_qc[df_qc["QC_VAL"] == q_flag]
        ax_qc.scatter(
            sub["CYCLE_NUMBER"], sub["PRES"],
            c=QC_COLOR_MAP.get(q_flag, "#000000"),
            label=f"Flag {q_flag}",
            s=12, alpha=0.7, edgecolors="none"
        )

    ax_qc.set_ylim(max_pres, 0)
    ax_qc.set_ylabel("Pressure (dbar)", fontsize=11)
    ax_qc.set_xlabel("Cycle Number", fontsize=11)
    ax_qc.set_title(f"{title_prefix} QC Flags Section Map", fontsize=12, fontweight="bold")
    ax_qc.grid(True, linestyle="--", alpha=0.4)

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=QC_COLOR_MAP.get(q, "#000000"),
               markersize=8, label=f"QC {q}")
        for q in unique_qcs
    ]
    ax_qc.legend(
        handles=legend_handles, 
        title="QC Flag", 
        bbox_to_anchor=(0.0, 1.0), 
        loc="upper left", 
        bbox_transform=cax_qc.transAxes, 
        frameon=True, 
        fontsize=9
    )

    today_stamp = dt.now().strftime("%Y-%m-%d")
    plt.suptitle(
        f"BGC-Argo Float {WMOfloatid} - {title_prefix} Section DMQC Summary (Verified BD Output)\nGenerated on: {today_stamp}",
        fontsize=13, fontweight="bold", y=0.98
    )

    out_file = os.path.join(output_dir, f"Float_{WMOfloatid}_{file_suffix}_Section_Summary.png")
    plt.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[{WMOfloatid}] Saved Aligned BD-Verified QC Summary Plot: {out_file}")


def extract_bd_parameter_dataframe(processed_bd_files, param_name, adj_param_name=None):
    """
    Independently extracts levels for a specific parameter paired directly with
    its own pressure array. Decouples CTD, oxygen, and optical sampling depths
    so data between 150m and 2000m is completely preserved.
    """
    records = []

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            cycle_num = int(fpath[-6:-3])
            iprof = detect_parameter_profile(ds, param_name)

            raw_v_name = param_name
            if raw_v_name not in ds.variables and param_name == "CHLA_FLUORESCENCE":
                raw_v_name = "FLUORESCENCE_CHLA"

            adj_v_name = adj_param_name if adj_param_name else f"{param_name}_ADJUSTED"
            if adj_v_name not in ds.variables and "FLUORESCENCE_CHLA" in raw_v_name:
                adj_v_name = "FLUORESCENCE_CHLA_ADJUSTED"

            raw_qc_name = f"{raw_v_name}_QC"
            adj_qc_name = f"{adj_v_name}_QC"

            pres_arr = get_pres_array_for_profile(ds, iprof)
            n_levels = len(pres_arr)

            def get_arr(v_name):
                if v_name not in ds.variables:
                    return None
                v = ds.variables[v_name]
                arr = v[iprof, :] if v.ndim > 1 else v[:]
                if hasattr(arr, "filled"):
                    arr = arr.filled(99999.0 if arr.dtype.kind in ['f', 'i'] else b" ")
                return np.array(arr)

            raw_arr = get_arr(raw_v_name)
            raw_qc_arr = get_arr(raw_qc_name)
            adj_arr = get_arr(adj_v_name)
            adj_qc_arr = get_arr(adj_qc_name)

            for k in range(n_levels):
                p_val = float(pres_arr[k])
                if np.isnan(p_val) or p_val == 99999.0:
                    continue

                def parse_val(arr):
                    if arr is None or k >= len(arr):
                        return np.nan
                    val = arr[k]
                    try:
                        f = float(val)
                        return np.nan if f == 99999.0 else f
                    except (ValueError, TypeError):
                        return np.nan

                def parse_qc(arr):
                    if arr is None or k >= len(arr):
                        return "9"
                    val = arr[k]
                    if isinstance(val, (bytes, np.bytes_)):
                        val = val.decode("utf-8", errors="ignore").strip()
                    else:
                        val = str(val).strip()
                    return val if val in ["0", "1", "2", "3", "4", "5", "8", "9"] else "9"

                records.append({
                    "CYCLE_NUMBER": cycle_num,
                    "PRES": p_val,
                    param_name: parse_val(raw_arr),
                    f"{param_name}_QC": parse_qc(raw_qc_arr),
                    f"{param_name}_ADJUSTED": parse_val(adj_arr),
                    f"{param_name}_ADJUSTED_QC": parse_qc(adj_qc_arr if adj_qc_arr is not None else raw_qc_arr),
                })

            ds.close()
        except Exception:
            traceback.print_exc()

    return pd.DataFrame(records)


def generate_all_qc_summary_plots(processed_bd_files, WMOfloatid, output_root):
    """
    Builds strictly 2 section summary plots per variable:
      1. R_mode: <PARAM> & <PARAM>_QC
      2. D_mode: <PARAM>_ADJUSTED & <PARAM>_ADJUSTED_QC
    Both explicitly labeled in titles and file suffixes.
    """
    if not processed_bd_files:
        return

    qc_summary_dir = os.path.join(output_root, "QC_SUMMARIES")
    os.makedirs(qc_summary_dir, exist_ok=True)

    tasks = [
        ("DOXY", "plasma", r"Oxygen ($\mu$mol/kg)", False),
        ("BBP700", "turbo", r"$b_{bp}(700)$ (m$^{-1}$)", False),
        ("CHLA", "turbo", r"Chlorophyll a (mg/m$^3$)", True),
        ("CHLA_FLUORESCENCE", "mako", "Fluorescence (RU)", False),
    ]

    for param_base, palette, unit, is_chla in tasks:
        try:
            df_param = extract_bd_parameter_dataframe(processed_bd_files, param_base)
            if df_param.empty:
                continue

            # 1. Raw Section Plot labeled explicitly as R_mode
            generate_section_qc_plot(
                df_param,
                param_base,
                f"{param_base}_QC",
                f"{param_base} R_mode",
                f"{param_base}_R_mode",
                qc_summary_dir,
                WMOfloatid,
                col_palette=palette,
                unit_str=unit,
                is_chla=is_chla,
            )

            # 2. Delayed-Mode Section Plot labeled explicitly as D_mode
            generate_section_qc_plot(
                df_param,
                f"{param_base}_ADJUSTED",
                f"{param_base}_ADJUSTED_QC",
                f"{param_base} D_mode",
                f"{param_base}_D_mode",
                qc_summary_dir,
                WMOfloatid,
                col_palette=palette,
                unit_str=unit,
                is_chla=is_chla,
            )
        except Exception:
            traceback.print_exc()


# ==============================================================================
# SECTION 4: Processing Loop
# ==============================================================================

processed_output_paths = []
script_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in locals() else os.getcwd()

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

            # 1. Update PARAMETER_DATA_MODE and DATA_MODE
            write_parameter_data_modes(ds)

            # 2. Process calibrations and CSV adjustments
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_chla)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_bbp)
            write_DOXY_slope_drift(ds, 0, df_bio, idx_profile, WMOfloatid=WMOfloatid)
            write_DOXY_from_csv(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)

            # 3. Carry over missing adjusted levels and clean QC variables
            clean_and_fill_qc_variables(ds, WMOfloatid=WMOfloatid)

            # 4. Metadata and attributes
            write_history_metadata(ds, iprof_chla, iprof_doxy, WMOfloatid=WMOfloatid)
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

    final_processed_files = []
    for w_file in new_bd_files:
        path, name = os.path.split(w_file)

        new_name = name.replace("w_AOML_BD", "BD")
        new_name = new_name.replace("w_AOML_BR", "BD")
        new_name = new_name.replace("w_BD", "BD")
        new_name = new_name.replace("w_BR", "BD")

        new_path = os.path.join(final_out_dir, new_name)
        safe_rename(w_file, new_path)
        processed_output_paths.append(new_path)
        final_processed_files.append(new_path)

    # 5. Generate DMQC Status Plot reading directly from output BD files
    generate_dmqc_status_plot(WMOfloatid, final_processed_files, script_dir)

    # 6. Generate Section Contour & QC Summary Plots directly from output BD files (R_mode & D_mode)
    generate_all_qc_summary_plots(final_processed_files, WMOfloatid, script_dir)