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

VALID_QC_CHARS = {"0", "1", "2", "3", "4", "5", "8", "9"}
MISSING_VAL = -99999.0


def is_missing_value(val):
    if val is None or pd.isna(val):
        return True
    try:
        f = float(val)
        return f == -99999.0 or f == 99999.0
    except (ValueError, TypeError):
        return False


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
        "CHLA_FINAL_QC", "BBP700_FINAL_QC",
        "CHLA_ADJUSTED_QC", "CHLA_FLUORESCENCE_ADJUSTED_QC",
        "FLUORESCENCE_CHLA_ADJUSTED_QC", "BBP700_ADJUSTED_QC",
    ]

    for col in float_cols:
        if col in df.columns:
            df.loc[override_condition, col] = MISSING_VAL

    for col in qc_cols:
        if col in df.columns:
            df.loc[override_condition, col] = 4

    doxy_hardcode_cond = (float_str == "2904010") & (cycle_num == 60)
    if doxy_hardcode_cond.any():
        for doxy_col in ["DOXY_FINAL", "DOXY_ADJUSTED_ERROR"]:
            if doxy_col in df.columns:
                df.loc[
                    doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0) | (df["DOXY"] == MISSING_VAL)),
                    doxy_col,
                ] = MISSING_VAL

        if "DOXY_FINAL_QC" in df.columns:
            df.loc[
                doxy_hardcode_cond & (df["DOXY"].isna() | (df["DOXY"] == 99999.0) | (df["DOXY"] == MISSING_VAL)),
                "DOXY_FINAL_QC",
            ] = 9

    return df


def detect_parameter_profile(ds, param_prefix):
    """
    Robust profile index detection handling padding, byte string decodings,
    and parameter aliases across N_PROF.
    """
    if "STATION_PARAMETERS" not in ds.variables:
        return 0

    n_prof = ds.dimensions["N_PROF"].size
    n_param = ds.dimensions["N_PARAM"].size

    target_prefixes = [param_prefix]
    if param_prefix in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_prefixes = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
    elif param_prefix in ["PH", "PH_IN_SITU_TOTAL"]:
        target_prefixes = ["PH_IN_SITU_TOTAL", "PH"]
    elif param_prefix in ["NITRATE"]:
        target_prefixes = ["NITRATE"]

    # 1. Search for exact prefix match inside STATION_PARAMETERS
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

    # 2. Fallback: inspect actual 2D data variables directly across N_PROF
    for tp in target_prefixes:
        for candidate_var in [tp, f"{tp}_ADJUSTED"]:
            if candidate_var in ds.variables:
                var = ds.variables[candidate_var]
                if var.ndim > 1:
                    for iprof in range(var.shape[0]):
                        sl = var[iprof, :]
                        vals = sl.filled(MISSING_VAL) if hasattr(sl, "filled") else np.array(sl)
                        if (~pd.isna(vals) & (vals != MISSING_VAL) & (vals != 99999.0)).any():
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
    elif param == "NITRATE":
        param_candidates.extend(["NITRATE", "NITRATE_ADJUSTED"])

    for p in param_candidates:
        if p in bgc_file.variables:
            var = bgc_file.variables[p]
            p_slice = var[iprof, :] if var.ndim > 1 else var[:]
            p_vals = p_slice.filled(MISSING_VAL) if hasattr(p_slice, "filled") else np.array(p_slice)
            valid_mask = (~pd.isna(p_vals)) & (p_vals != MISSING_VAL) & (p_vals != 99999.0)
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
            return m.decode("utf-8").strip() if isinstance(m, (bytes, np.bytes_)) else str(m).strip()
        return "R"

    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
    pdm_array = bgc_file.variables["PARAMETER_DATA_MODE"][iprof]

    if isinstance(station_params, np.ma.MaskedArray):
        station_params = station_params.filled(b" ")

    target_names = [param_base]
    if param_base in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_names = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
    elif param_base in ["PH", "PH_IN_SITU_TOTAL"]:
        target_names = ["PH_IN_SITU_TOTAL", "PH"]

    for j, p_bytes in enumerate(station_params):
        p_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
            for c in p_bytes
        ]).strip()

        if p_str in target_names or any(p_str.startswith(tn) for tn in target_names):
            m_val = pdm_array[j]
            return m_val.decode("utf-8").strip() if isinstance(m_val, (bytes, np.bytes_)) else str(m_val).strip()

    dm_var = bgc_file.variables.get("DATA_MODE")
    if dm_var is not None:
        m = dm_var[iprof] if dm_var.ndim > 0 else dm_var[:]
        return m.decode("utf-8").strip() if isinstance(m, (bytes, np.bytes_)) else str(m).strip()
    return "R"


def write_parameter_data_modes(bgc_file):
    """
    Updates PARAMETER_DATA_MODE and DATA_MODE.
    Explicitly ignores PH and NITRATE so their original modes are never modified.
    """
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

            if not p_str or "PH" in p_str or "NITRATE" in p_str:
                continue

            if p_str == "PRES" or p_str.startswith("DOWN_IRRADIANCE") or p_str.startswith("UP_RADIANCE"):
                pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"
            elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY"]:
                if has_valid_measurements(bgc_file, p_str, iprof):
                    pdm[iprof, j] = b"D" if pdm.dtype.kind in ["S", "U", "O"] else "D"
                    data_mode[iprof] = b"D" if data_mode.dtype.kind in ["S", "U", "O"] else "D"
                else:
                    pdm[iprof, j] = b"R" if pdm.dtype.kind in ["S", "U", "O"] else "R"

    bgc_file.variables["PARAMETER_DATA_MODE"][:] = pdm
    bgc_file.variables["DATA_MODE"][:] = data_mode


def get_profile_qc_grade(qc_array, parameter_present=True):
    """
    Implements Argo Users Manual Table 2a (consistent with get_profile_param_qc):
    - Values of 1, 2, 5, 8 are GOOD data.
    - Flag 9 and blanks (' ', '') are NOT USED in computation.
    - All other flags are BAD data.
    - If no valid measurements exist: returns 'F' if parameter exists on float, otherwise ' '.
    """
    if not parameter_present:
        return " "

    raw_qcs = []
    has_any_point = False
    for q in np.ravel(qc_array):
        if isinstance(q, np.ma.core.MaskedConstant) or q is np.ma.masked:
            continue
        s = q.decode("utf-8", errors="ignore").strip() if isinstance(q, (bytes, np.bytes_)) else str(q).strip()
        if s != "":
            has_any_point = True
        if s not in ["", " ", "9"]:
            raw_qcs.append(s)

    total_points = len(raw_qcs)
    if total_points == 0:
        return "F" if has_any_point else " "

    good_qc = {"1", "2", "5", "8"}
    good_count = sum(1 for q in raw_qcs if q in good_qc)
    good_fraction = good_count / total_points

    if good_count == total_points:
        return "A"
    elif good_fraction >= 0.75:
        return "B"
    elif good_fraction >= 0.50:
        return "C"
    elif good_fraction >= 0.25:
        return "D"
    elif good_fraction > 0.0:
        return "E"
    else:
        return "F"


def get_pres_and_qc_for_profile(ds, iprof):
    n_levels = ds.dimensions["N_LEVELS"].size
    pres_vals = np.full((n_levels,), MISSING_VAL, dtype="float32")
    pres_qcs = np.full((n_levels,), b"9", dtype="|S1")

    p_var_name = "PRES_ADJUSTED" if "PRES_ADJUSTED" in ds.variables else ("PRES" if "PRES" in ds.variables else None)
    if p_var_name:
        v = ds.variables[p_var_name]
        arr = v[iprof, :] if v.ndim > 1 else v[:]
        p_raw = arr.filled(MISSING_VAL) if hasattr(arr, "filled") else np.array(arr)
        if (~np.isnan(p_raw) & (p_raw != MISSING_VAL) & (p_raw != 99999.0)).any():
            pres_vals = np.float32(p_raw)
        elif "PRES" in ds.variables:
            arr0 = ds.variables["PRES"][0, :] if ds.variables["PRES"].ndim > 1 else ds.variables["PRES"][:]
            pres_vals = np.float32(arr0.filled(MISSING_VAL) if hasattr(arr0, "filled") else np.array(arr0))

    q_var_name = "PRES_ADJUSTED_QC" if "PRES_ADJUSTED_QC" in ds.variables else ("PRES_QC" if "PRES_QC" in ds.variables else None)
    if q_var_name:
        qv = ds.variables[q_var_name]
        q_slice = qv[iprof, :] if qv.ndim > 1 else qv[:]
        if hasattr(q_slice, "filled"):
            q_slice = q_slice.filled(b"9")
        pres_qcs = np.array(q_slice).astype("|S1")

    return pres_vals, pres_qcs


def match_csv_row_by_pres(cycle_df, nc_pres):
    if cycle_df.empty or is_missing_value(nc_pres):
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


def sanitize_qc_char(val, default="9"):
    if pd.isna(val) or val is None or val == "":
        return default
    s = str(val).strip()
    match = re.search(r"([0-9])", s)
    if match and match.group(1) in VALID_QC_CHARS:
        return match.group(1)
    return default


def parse_raw_qc_char(q_arr, idx):
    if q_arr is None or idx >= len(q_arr):
        return "9"
    val = q_arr[idx]
    if isinstance(val, (bytes, np.bytes_)):
        val = val.decode("utf-8", errors="ignore").strip()
    else:
        val = str(val).strip()
    return val if val in VALID_QC_CHARS else "9"


def standardize_raw_missing_data(bgc_file):
    """
    Standardizes ONLY the raw variables (CHLA, DOXY, CHLA_FLUORESCENCE/FLUORESCENCE_CHLA, BBP700)
    to fill missing values with -99999.0 and assign a QC flag of '9'.
    Explicitly excludes and does NOT touch pH or Nitrate.
    """
    n_prof = bgc_file.dimensions["N_PROF"].size
    raw_params = ["CHLA", "DOXY", "BBP700", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]

    for param in raw_params:
        if param in bgc_file.variables:
            var = bgc_file.variables[param]
            qc_name = f"{param}_QC"
            qc_var = bgc_file.variables.get(qc_name)

            for iprof in range(n_prof):
                arr = var[iprof, :] if var.ndim > 1 else var[:]
                vals = arr.filled(MISSING_VAL) if hasattr(arr, "filled") else np.array(arr, copy=True)
                missing_mask = pd.isna(vals) | (vals == 99999.0) | (vals == MISSING_VAL)

                if missing_mask.any():
                    vals[missing_mask] = MISSING_VAL
                    if var.ndim > 1:
                        bgc_file.variables[param][iprof, :] = vals
                    else:
                        bgc_file.variables[param][:] = vals

                    if qc_var is not None:
                        q_slice = qc_var[iprof, :] if qc_var.ndim > 1 else qc_var[:]
                        q_vals = np.array(q_slice, copy=True)
                        if hasattr(q_vals, "filled"):
                            q_vals = q_vals.filled(b" ")
                        q_char = q_vals.astype("|S1")
                        q_char[missing_mask] = b"9"

                        if qc_var.ndim > 1:
                            bgc_file.variables[qc_name][iprof, :] = q_char
                        else:
                            bgc_file.variables[qc_name][:] = q_char


def clean_and_fill_qc_variables(bgc_file, WMOfloatid=None):
    """
    Cleans adjusted QC flags:
    - Missing/unadjusted levels receive QC '9' and value -99999.0.
    - Broken sensors remain flagged '4'.
    - Updates PROFILE_<PARAM>_QC across its proper profile index.
    - Excludes pH and Nitrate entirely.
    """
    n_prof = bgc_file.dimensions["N_PROF"].size
    skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

    for var_name, var in list(bgc_file.variables.items()):
        if not var_name.endswith("_ADJUSTED_QC") or var_name.startswith("PROFILE_"):
            continue
        if "PH" in var_name or "NITRATE" in var_name:
            continue
        if skip_doxy and "DOXY" in var_name:
            continue

        base_param = var_name[:-12]
        adj_data_var_name = f"{base_param}_ADJUSTED"

        for iprof in range(n_prof):
            try:
                mode = get_param_data_mode(bgc_file, iprof, base_param)

                if mode == "R":
                    n_lev = bgc_file.dimensions["N_LEVELS"].size
                    fill_qc = np.full(shape=(n_lev,), fill_value=b" ", dtype="|S1")
                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = fill_qc
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = b" "
                    else:
                        bgc_file.variables[var_name][:] = fill_qc
                    continue

                if adj_data_var_name in bgc_file.variables:
                    adj_p_var = bgc_file.variables[adj_data_var_name]
                    p_raw = adj_p_var[iprof, :] if adj_p_var.ndim > 1 else adj_p_var[:]
                    p_vals = p_raw.filled(MISSING_VAL) if hasattr(p_raw, "filled") else np.array(p_raw)

                    q_raw = var[iprof, :] if var.ndim > 1 else var[iprof]
                    q_vals = np.array(q_raw, copy=True)
                    if hasattr(q_vals, "filled"):
                        q_vals = q_vals.filled(b" ")
                    q_char = q_vals.astype("|S1")

                    missing_mask = pd.isna(p_vals) | (p_vals == 99999.0) | (p_vals == MISSING_VAL)
                    if missing_mask.any():
                        p_vals[missing_mask] = MISSING_VAL
                        if adj_p_var.ndim > 1:
                            bgc_file.variables[adj_data_var_name][iprof, :] = p_vals
                        else:
                            bgc_file.variables[adj_data_var_name][:] = p_vals

                        keep_four_mask = missing_mask & (q_char == b"4")
                        q_char[missing_mask] = b"9"
                        q_char[keep_four_mask] = b"4"

                    if var.ndim == 2:
                        bgc_file.variables[var_name][iprof, :] = q_char
                    elif var.ndim == 1:
                        bgc_file.variables[var_name][iprof] = q_char if q_char.ndim == 0 else q_char[0]
                    else:
                        bgc_file.variables[var_name][:] = q_char

            except Exception:
                traceback.print_exc()

    # Recalculate and populate PROFILE_<PARAM>_QC on the exact parameter sub-profile
    for var_name in list(bgc_file.variables.keys()):
        if not var_name.startswith("PROFILE_") or not var_name.endswith("_QC"):
            continue
        if "PH" in var_name or "NITRATE" in var_name:
            continue

        param_base = var_name[8:-3]
        target_prof = detect_parameter_profile(bgc_file, param_base)

        mode = get_param_data_mode(bgc_file, target_prof, param_base)
        target_qc_name = (
            f"{param_base}_ADJUSTED_QC"
            if mode in ["D", "A"] and f"{param_base}_ADJUSTED_QC" in bgc_file.variables
            else f"{param_base}_QC"
        )
        if param_base == "PRES":
            target_qc_name = "PRES_QC"

        param_is_present = has_valid_measurements(bgc_file, param_base, target_prof)

        if target_qc_name in bgc_file.variables:
            target_qc_var = bgc_file.variables[target_qc_name]
            target_qc_array = (
                target_qc_var[target_prof, :] if target_qc_var.ndim > 1 else target_qc_var[:]
            )
            prof_grade = get_profile_qc_grade(target_qc_array, parameter_present=param_is_present)
            grade_bytes = prof_grade.encode("utf-8")

            # Update ALL profile slots so any reader looking at iprof=0 or target_prof sees the grade
            for iprof in range(n_prof):
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

    fluo_csv_col = next((c for c in ["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"] if c in df_bio.columns), None)
    fluo_qc_csv_col = next((c for c in ["CHLA_FLUORESCENCE_ADJUSTED_QC", "FLUORESCENCE_CHLA_ADJUSTED_QC", "CHLA_FLUORESCENCE_QC", "FLUORESCENCE_CHLA_QC"] if c in df_bio.columns), None)

    for prof in range(n_prof):
        mode = get_param_data_mode(bgc_file, prof, "CHLA")

        CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_Array[:] = MISSING_VAL
        CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")
        CHLA_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        CHLA_Adjusted_ERROR_Array[:] = MISSING_VAL

        FLUORESCENCE_CHLA_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        FLUORESCENCE_CHLA_Adjusted_Array[:] = MISSING_VAL
        FLUORESCENCE_CHLA_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")
        FLUORESCENCE_CHLA_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        FLUORESCENCE_CHLA_Adjusted_ERROR_Array[:] = MISSING_VAL

        raw_chla_arr = bgc_file.variables["CHLA"][prof, :] if "CHLA" in bgc_file.variables else None
        raw_chla_qc_arr = bgc_file.variables["CHLA_QC"][prof, :] if "CHLA_QC" in bgc_file.variables else None

        fluo_name = "FLUORESCENCE_CHLA" if "FLUORESCENCE_CHLA" in bgc_file.variables else "CHLA_FLUORESCENCE"
        raw_fluo_arr = bgc_file.variables[fluo_name][prof, :] if fluo_name in bgc_file.variables else None
        raw_fluo_qc_arr = bgc_file.variables[f"{fluo_name}_QC"][prof, :] if f"{fluo_name}_QC" in bgc_file.variables else None

        nc_pres_all, nc_pres_qc_all = get_pres_and_qc_for_profile(bgc_file, prof)

        if mode in ["D", "A"]:
            for i in range(n_levels):
                nc_pres = nc_pres_all[i]
                pres_qc_char = nc_pres_qc_all[i].decode("utf-8") if isinstance(nc_pres_qc_all[i], (bytes, np.bytes_)) else str(nc_pres_qc_all[i])

                is_pres_missing = is_missing_value(nc_pres) or pres_qc_char == "9"
                is_pres_bad = pres_qc_char == "4"

                raw_chla_val = raw_chla_arr[i] if raw_chla_arr is not None else MISSING_VAL
                raw_chla_qc = parse_raw_qc_char(raw_chla_qc_arr, i)

                raw_fluo_val = raw_fluo_arr[i] if raw_fluo_arr is not None else MISSING_VAL
                raw_fluo_qc = parse_raw_qc_char(raw_fluo_qc_arr, i)

                # CHLA_ADJUSTED
                if is_pres_missing:
                    CHLA_AdjustedQC_Array[i] = b"9"
                elif is_pres_bad:
                    CHLA_AdjustedQC_Array[i] = b"4"
                elif is_broken or raw_chla_qc == "4":
                    CHLA_AdjustedQC_Array[i] = b"4"
                elif raw_chla_qc == "9" or is_missing_value(raw_chla_val):
                    CHLA_AdjustedQC_Array[i] = b"9"
                else:
                    matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                    if matched_row is not None and "CHLA_FINAL" in matched_row:
                        chla_final_val = matched_row.get("CHLA_FINAL")
                        csv_chla_qc = sanitize_qc_char(matched_row.get("CHLA_FINAL_QC"), default="1")

                        if is_missing_value(chla_final_val) or csv_chla_qc in ["4", "9"]:
                            CHLA_AdjustedQC_Array[i] = csv_chla_qc.encode("utf-8")
                        else:
                            CHLA_Adjusted_Array[i] = np.float32(chla_final_val)
                            CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                            CHLA_Adjusted_Array.mask[i] = False
                            CHLA_Adjusted_ERROR_Array.mask[i] = False
                            CHLA_AdjustedQC_Array[i] = csv_chla_qc.encode("utf-8")
                    else:
                        CHLA_AdjustedQC_Array[i] = b"9"

                # FLUORESCENCE_CHLA_ADJUSTED
                if is_pres_missing:
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"9"
                elif is_pres_bad:
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"4"
                elif is_broken or raw_fluo_qc == "4":
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"4"
                elif raw_fluo_qc == "9" or is_missing_value(raw_fluo_val):
                    FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"9"
                else:
                    matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                    if matched_row is not None and fluo_csv_col is not None and pd.notna(matched_row.get(fluo_csv_col)):
                        fluo_csv_val = matched_row.get(fluo_csv_col)
                        fluo_qc_char = sanitize_qc_char(matched_row.get(fluo_qc_csv_col), default="1")

                        if is_missing_value(fluo_csv_val) or fluo_qc_char in ["4", "9"]:
                            FLUORESCENCE_CHLA_AdjustedQC_Array[i] = fluo_qc_char.encode("utf-8")
                        else:
                            FLUORESCENCE_CHLA_Adjusted_Array[i] = np.float32(fluo_csv_val)
                            FLUORESCENCE_CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)
                            FLUORESCENCE_CHLA_Adjusted_Array.mask[i] = False
                            FLUORESCENCE_CHLA_Adjusted_ERROR_Array.mask[i] = False
                            FLUORESCENCE_CHLA_AdjustedQC_Array[i] = fluo_qc_char.encode("utf-8")
                    else:
                        FLUORESCENCE_CHLA_AdjustedQC_Array[i] = b"9"

        if "CHLA_ADJUSTED" in bgc_file.variables:
            bgc_file.variables["CHLA_ADJUSTED"][prof, :] = CHLA_Adjusted_Array
        if "CHLA_ADJUSTED_QC" in bgc_file.variables:
            bgc_file.variables["CHLA_ADJUSTED_QC"][prof, :] = CHLA_AdjustedQC_Array
        if "CHLA_ADJUSTED_ERROR" in bgc_file.variables:
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

        BBP700_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        BBP700_Adjusted_Array[:] = MISSING_VAL
        BBP700_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9" if mode in ["D", "A"] else b" ", dtype="|S1")
        BBP700_Adjusted_ERROR_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        BBP700_Adjusted_ERROR_Array[:] = MISSING_VAL

        raw_bbp_arr = bgc_file.variables["BBP700"][prof, :] if "BBP700" in bgc_file.variables else None
        raw_bbp_qc_arr = bgc_file.variables["BBP700_QC"][prof, :] if "BBP700_QC" in bgc_file.variables else None

        nc_pres_all, nc_pres_qc_all = get_pres_and_qc_for_profile(bgc_file, prof)

        if mode in ["D", "A"]:
            for i in range(n_levels):
                nc_pres = nc_pres_all[i]
                pres_qc_char = nc_pres_qc_all[i].decode("utf-8") if isinstance(nc_pres_qc_all[i], (bytes, np.bytes_)) else str(nc_pres_qc_all[i])

                is_pres_missing = is_missing_value(nc_pres) or pres_qc_char == "9"
                is_pres_bad = pres_qc_char == "4"

                raw_bbp_val = raw_bbp_arr[i] if raw_bbp_arr is not None else MISSING_VAL
                raw_bbp_qc = parse_raw_qc_char(raw_bbp_qc_arr, i)

                if is_pres_missing:
                    BBP700_AdjustedQC_Array[i] = b"9"
                elif is_pres_bad:
                    BBP700_AdjustedQC_Array[i] = b"4"
                elif is_broken or raw_bbp_qc == "4":
                    BBP700_AdjustedQC_Array[i] = b"4"
                elif raw_bbp_qc == "9" or is_missing_value(raw_bbp_val):
                    BBP700_AdjustedQC_Array[i] = b"9"
                else:
                    matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                    if matched_row is not None and "BBP700_FINAL" in matched_row:
                        raw_bbp_final = matched_row.get("BBP700_FINAL")
                        csv_qc = sanitize_qc_char(matched_row.get("BBP700_FINAL_QC"), default="1")

                        if is_missing_value(raw_bbp_final) or csv_qc in ["4", "9"]:
                            BBP700_AdjustedQC_Array[i] = csv_qc.encode("utf-8")
                        else:
                            BBP700_Adjusted_Array[i] = np.float32(raw_bbp_final)
                            BBP700_Adjusted_ERROR_Array[i] = np.float32(BBP700_Adjusted_ERROR_est)
                            BBP700_Adjusted_Array.mask[i] = False
                            BBP700_Adjusted_ERROR_Array.mask[i] = False
                            BBP700_AdjustedQC_Array[i] = csv_qc.encode("utf-8")
                    else:
                        BBP700_AdjustedQC_Array[i] = b"9"

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

        if is_excluded:
            if "PROFILE_DOXY_QC" in ds.variables and "DOXY_QC" in ds.variables:
                doxy_qc_slice = ds.variables["DOXY_QC"][iprof, :] if ds.variables["DOXY_QC"].ndim > 1 else ds.variables["DOXY_QC"][:]
                prof_grade = get_profile_qc_grade(doxy_qc_slice, parameter_present=has_valid_measurements(ds, "DOXY", iprof))
                grade_bytes = prof_grade.encode("utf-8")
                if ds.variables["PROFILE_DOXY_QC"].ndim == 1:
                    ds.variables["PROFILE_DOXY_QC"][iprof] = grade_bytes
                elif ds.variables["PROFILE_DOXY_QC"].ndim == 2:
                    ds.variables["PROFILE_DOXY_QC"][iprof, 0] = grade_bytes
            continue

        if mode == "R":
            continue

        raw_doxy_arr = ds.variables["DOXY"][iprof, :] if "DOXY" in ds.variables else None
        raw_doxy_qc_arr = ds.variables["DOXY_QC"][iprof, :] if "DOXY_QC" in ds.variables else None

        nc_pres_all, nc_pres_qc_all = get_pres_and_qc_for_profile(ds, iprof)

        DOXY_Adjusted_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        DOXY_Adjusted_Array[:] = MISSING_VAL
        DOXY_Adjusted_Error_Array = np.ma.masked_all(shape=(n_levels,), dtype="float32")
        DOXY_Adjusted_Error_Array[:] = MISSING_VAL
        DOXY_AdjustedQC_Array = np.full(shape=(n_levels,), fill_value=b"9", dtype="|S1")

        for i in range(n_levels):
            nc_pres = nc_pres_all[i]
            pres_qc_char = nc_pres_qc_all[i].decode("utf-8") if isinstance(nc_pres_qc_all[i], (bytes, np.bytes_)) else str(nc_pres_qc_all[i])

            is_pres_missing = is_missing_value(nc_pres) or pres_qc_char == "9"
            is_pres_bad = pres_qc_char == "4"

            raw_doxy_val = raw_doxy_arr[i] if raw_doxy_arr is not None else MISSING_VAL
            raw_doxy_qc = parse_raw_qc_char(raw_doxy_qc_arr, i)

            if is_pres_missing:
                DOXY_AdjustedQC_Array[i] = b"9"
            elif is_pres_bad:
                DOXY_AdjustedQC_Array[i] = b"4"
            elif raw_doxy_qc == "4":
                DOXY_AdjustedQC_Array[i] = b"4"
            elif raw_doxy_qc == "9" or is_missing_value(raw_doxy_val):
                DOXY_AdjustedQC_Array[i] = b"9"
            else:
                matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                if matched_row is not None and "DOXY_FINAL" in matched_row:
                    raw_doxy_final = matched_row.get("DOXY_FINAL")
                    csv_qc = sanitize_qc_char(matched_row.get("DOXY_FINAL_QC"), default="1")

                    if is_missing_value(raw_doxy_final) or csv_qc in ["4", "9"]:
                        DOXY_AdjustedQC_Array[i] = csv_qc.encode("utf-8")
                    else:
                        DOXY_Adjusted_Array[i] = np.float32(raw_doxy_final)
                        DOXY_Adjusted_Array.mask[i] = False
                        DOXY_AdjustedQC_Array[i] = csv_qc.encode("utf-8")

                        raw_doxy_error = matched_row.get("DOXY_ADJUSTED_ERROR")
                        if pd.notna(raw_doxy_error) and not is_missing_value(raw_doxy_error):
                            DOXY_Adjusted_Error_Array[i] = np.float32(raw_doxy_error)
                            DOXY_Adjusted_Error_Array.mask[i] = False
                        else:
                            DOXY_Adjusted_Error_Array[i] = np.float32(DOXY_Adjusted_ERROR_est_default)
                            DOXY_Adjusted_Error_Array.mask[i] = False
                else:
                    DOXY_AdjustedQC_Array[i] = b"9"

        if "DOXY_ADJUSTED" in var_names:
            doxy_adj_var = ds.variables["DOXY_ADJUSTED"]
            if doxy_adj_var.ndim > 1:
                doxy_adj_var[iprof, :] = DOXY_Adjusted_Array
            else:
                doxy_adj_var[:] = DOXY_Adjusted_Array

        if "DOXY_ADJUSTED_QC" in var_names:
            doxy_adj_qc_var = ds.variables["DOXY_ADJUSTED_QC"]
            if doxy_adj_qc_var.ndim > 1:
                doxy_adj_qc_var[iprof, :] = DOXY_AdjustedQC_Array
            else:
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
    idx = np.round(np.linspace(0, len(temp_list) - 1, n)).astype(int)
    return [temp_list[i] for i in idx]


def add_Data_mode_sp(fig, float_map_inst, std_fig_width=1000):
    x_nprof = []
    y_dataType = []
    color_dataType = []

    for n_prof in float_map_inst["NPROF_list"]:
        x_nprof.append(str(n_prof + 1))
        y_dataType.append(" ")
        df_sub = float_map_inst["NPROF_NPARAM_df"]
        modes = df_sub[df_sub["N_PROF"] == n_prof]["PARAMETER_DATA_MODE"].tolist()
        best_mode = min(modes, key=lambda x: data_mode_rank_map.get(x, 99)) if modes else "R"
        color_dataType.append(data_type_color_map.get(best_mode, "#FFFFFF"))

    x_ticks = get_subsample_prof_idx(x_nprof, min(7, len(x_nprof)))
    fig.update_xaxes(
        tickmode="array",
        tickvals=x_ticks,
        ticktext=[str(x) for x in x_ticks],
        row="all",
        col=1,
    )

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(len(float_map_inst["NPROF_list"]), 1), 1.0)
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

    for n_prof in float_map_inst["NPROF_list"]:
        positioning_system = float_map_inst["NPROF_df"]["POSITIONING_SYSTEM"][n_prof].rstrip()
        position_qc = float_map_inst["NPROF_df"]["POSITION_QC"][n_prof]
        x1_nprof.append(str(n_prof + 1))
        y1_positioning_system.append(" ")
        y1_position_qc.append("  ")

        color1_positioning_system.append(position_system_color_map.get(positioning_system, "#ff00f7"))
        color1_position_qc.append(position_qc_color_map.get(position_qc, "#808080"))

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(len(float_map_inst["NPROF_list"]), 1), 1.0)
    
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
        nprof_nlevel_df = float_map_inst["NPROF_NLEVEL_df"]
        depth_profile = nprof_nlevel_df[nprof_nlevel_df["N_PROF"] == n_prof]["PRES_ADJUSTED"]
        max_d = np.nanmax(depth_profile) if len(depth_profile) > 0 and not np.all(np.isnan(depth_profile)) else 0.0
        y2_depth.append(max_d)

    max_prof_depth = np.max(y2_depth) if y2_depth and not np.all(np.isnan(y2_depth)) else 2000.0

    fig.add_trace(
        go.Bar(
            x=x2_nprof,
            y=y2_depth,
            orientation="v",
            marker=dict(
                color=color_depth,
                line=dict(color=color_depth),
                cornerradius="0%",
            ),
            showlegend=False,
        ),
        row=3,
        col=1,
    )

    fig.update_yaxes(
        range=[max_prof_depth, 0],
        visible=False,
        row=3,
        col=1,
    )

    return np.rint(max_prof_depth)


def add_param_QC_sp(fig, float_map_inst, PARAM, subplotNum, std_fig_width=1000):
    x3_nprof = []
    y3_placehodler1 = []
    y3_placehodler2 = []
    color3_profile_flag_QC = []
    color3_param_data_mode = []

    for n_prof in float_map_inst["NPROF_list"]:
        x3_nprof.append(str(n_prof + 1))
        y3_placehodler1.append(" ")
        y3_placehodler2.append("  ")

        # Resolve variable aliasing across both standard and shortened names
        possible_qc_keys = [f"PROFILE_{PARAM}_QC"]
        if PARAM == "CHLA_FLUORESCENCE":
            possible_qc_keys.append("PROFILE_FLUORESCENCE_CHLA_QC")
        elif PARAM == "FLUORESCENCE_CHLA":
            possible_qc_keys.append("PROFILE_CHLA_FLUORESCENCE_QC")
        elif PARAM == "PH_IN_SITU_TOTAL":
            possible_qc_keys.append("PROFILE_PH_QC")
        elif PARAM == "PH":
            possible_qc_keys.append("PROFILE_PH_IN_SITU_TOTAL_QC")

        profile_flag_QC = " "
        for key in possible_qc_keys:
            if key in float_map_inst["NPROF_df"] and n_prof in float_map_inst["NPROF_df"][key]:
                val = float_map_inst["NPROF_df"][key][n_prof]
                if val and str(val).strip() != "":
                    profile_flag_QC = str(val).strip()
                    break

        profile_flag_QC_color = param_qc_color_map.get(profile_flag_QC, "#FFFFFF")
        color3_profile_flag_QC.append(profile_flag_QC_color)

        possible_pdm_keys = [PARAM]
        if PARAM in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
            possible_pdm_keys.extend(["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"])
        elif PARAM in ["PH", "PH_IN_SITU_TOTAL"]:
            possible_pdm_keys.extend(["PH_IN_SITU_TOTAL", "PH"])

        param_data_mode = "R"
        for key in possible_pdm_keys:
            if key in float_map_inst["PDM_MAP"] and n_prof in float_map_inst["PDM_MAP"][key]:
                m = float_map_inst["PDM_MAP"][key][n_prof]
                if m and str(m).strip() != "":
                    param_data_mode = str(m).strip()
                    break

        param_data_mode_color = data_type_color_map.get(param_data_mode, "#FFFFFF")
        color3_param_data_mode.append(param_data_mode_color)

    marker_line_w = max(0.9 * (std_fig_width - 170) / max(len(float_map_inst["NPROF_list"]), 1), 1.0)

    fig.add_trace(
        go.Scatter(
            x=x3_nprof,
            y=y3_placehodler1,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color3_profile_flag_QC,
            marker_line_color=color3_profile_flag_QC,
            marker_line_width=marker_line_w,
            marker_size=15,
            orientation="h",
            showlegend=False,
        ),
        row=subplotNum,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x3_nprof,
            y=y3_placehodler2,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color3_param_data_mode,
            marker_line_color=color3_param_data_mode,
            marker_line_width=marker_line_w,
            marker_size=15,
            orientation="h",
            showlegend=False,
        ),
        row=subplotNum,
        col=1,
    )


def update_layout_DMQC_STATUS(fig, float_map_inst, max_prof_depth, subplotList, num_rows, std_fig_width=1000):
    n_profs = len(float_map_inst["NPROF_list"])
    fig.update_xaxes(
        range=[-0.05 * n_profs, n_profs * 1.05],
        dtick=max(1, int(n_profs * 0.33)),
    )

    calc_height = max(800, int(60 * num_rows + 150))

    fig.update_layout(
        autosize=False,
        height=calc_height,
        width=std_fig_width,
        bargap=0,
        bargroupgap=0,
        font=dict(color=TEXT_COLOR),
        title=dict(
            text=f"WMOID {float_map_inst['wmoID_str']}",
            font=dict(color=TEXT_COLOR, size=24),
            x=0.5,
            xanchor="center",
            yanchor="top",
        ),
        annotations=[
            dict(
                text="Data Type (R/A/D)",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=sum(fig.layout["yaxis1"].domain) / 2,
                yanchor="middle",
                showarrow=False,
                textangle=0,
                font=dict(size=12),
            ),
            dict(
                text="Positioning System",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=sum(fig.layout["yaxis2"].domain) / 2,
                yanchor="bottom",
                showarrow=False,
                textangle=0,
                font=dict(size=12),
            ),
            dict(
                text="Position QC",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=sum(fig.layout["yaxis2"].domain) / 2,
                yanchor="top",
                showarrow=False,
                textangle=0,
                font=dict(size=12),
            ),
            dict(
                text="Profile Depth",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=sum(fig.layout["yaxis3"].domain) / 2,
                yanchor="bottom",
                showarrow=False,
                textangle=0,
                font=dict(size=12),
            ),
            dict(
                text=f"Max: ~{max_prof_depth:.0f} dbar",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=sum(fig.layout["yaxis3"].domain) / 2,
                yanchor="top",
                showarrow=False,
                textangle=0,
                font=dict(size=12, color="#41AACE"),
            ),
        ],
        plot_bgcolor=BG_COLOR,
        paper_bgcolor=BG_COLOR,
        margin=dict(l=170, r=20, t=100, b=80),
    )

    for data_mode, hex_color in label_data_type_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{data_mode}",
                legendgroup="DataMode",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="Data Mode: ", selector=dict(legendgroup="DataMode"))

    for qc_code, hex_color in label_qc_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{qc_code}",
                legendgroup="PositioningSystemQC",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="Positioning System QC: ", selector=dict(legendgroup="PositioningSystemQC"))

    for qc_code, hex_color in label_param_qc_color_map.items():
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=hex_color),
                name=f"{qc_code}",
                legendgroup="ProfileQualityFlag",
                showlegend=True,
            )
        )
    fig.update_traces(legendgrouptitle_text="Profile PARAM Grade: ", selector=dict(legendgroup="ProfileQualityFlag"))

    new_annotations = []
    for PARAM, subplotidx in zip(float_map_inst["param_list"], subplotList):
        tempPARAM = "PH" if PARAM in ["PH_IN_SITU_TOTAL", "PH"] else ("FLUOR" if PARAM in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"] else PARAM)
        domain_center = sum(fig.layout[f"yaxis{subplotidx}"].domain) / 2

        new_annotations.extend([
            dict(
                text=f"{tempPARAM} Data Mode",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=domain_center,
                yanchor="bottom",
                showarrow=False,
                textangle=0,
                font=dict(size=11),
            ),
            dict(
                text=f"Profile {tempPARAM} QC",
                xref="paper", yref="paper",
                x=0,
                xanchor="right",
                y=domain_center,
                yanchor="top",
                showarrow=False,
                textangle=0,
                font=dict(size=11),
            ),
        ])
    fig.layout.annotations = list(fig.layout.annotations) + new_annotations

    float_status = float_map_inst.get("activity_status", "Active")
    if float_status == "Active":
        status_color = "#147628"
        annotation_text = f"Status: {float_status}"
    elif float_status == "Inactive":
        status_color = "#FA0226"
        annotation_text = f"Status: {float_status}"
    else:
        parts = float_status.split(" ", 1)
        base_status = parts[0]
        overdue_detail = parts[1] if len(parts) > 1 else ""
        status_color = "#f5c133"
        annotation_text = (
            f'<span style="font-size: 20px;">Status: {base_status}</span>'
            f'<span style="font-size: 14px; text-align: right;"><br>{overdue_detail}</span>'
        )

    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=1, y=1,
        xanchor="right", yanchor="top",
        showarrow=False,
        font=dict(size=20, color=status_color),
        yshift=60,
        xshift=140,
    )

    current_date = dt.now().strftime("%Y-%m-%d")
    fig.add_annotation(
        text=f"Figure Created: {current_date}",
        xref="paper", yref="paper",
        x=0, y=0,
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=12),
        xshift=-150, yshift=-60,
    )

    prof_date_map = dict(zip(float_map_inst["NPROF_list"], float_map_inst["JULD_list"]))
    if prof_date_map:
        first_prof = next(iter(prof_date_map))
        first_date = prof_date_map[first_prof].strftime("%Y-%m-%d") if isinstance(prof_date_map[first_prof], dt) else str(prof_date_map[first_prof])
        last_prof = next(reversed(prof_date_map))
        last_date = prof_date_map[last_prof].strftime("%Y-%m-%d") if isinstance(prof_date_map[last_prof], dt) else str(prof_date_map[last_prof])

        fig.add_annotation(
            text=f"Cycle {first_prof + 1}:  {first_date}",
            xref="paper", yref="paper",
            x=0, y=0,
            xanchor="left", yanchor="bottom",
            showarrow=False,
            font=dict(size=12),
            xshift=0, yshift=-50,
        )
        fig.add_annotation(
            text=f"Cycle {last_prof + 1}:  {last_date}",
            xref="paper", yref="paper",
            x=1, y=0,
            xanchor="right", yanchor="bottom",
            showarrow=False,
            font=dict(size=12),
            xshift=0, yshift=-50,
        )


def make_DMQC_STATUS_fig(float_map_inst, output_filepath=None, show_grid=True):
    num_rows = 3 + len(float_map_inst["param_list"])
    row_heights_list = [1 / num_rows] * num_rows

    fig = make_subplots(
        rows=num_rows,
        cols=1,
        shared_xaxes=True,
        x_title="Cycle Number",
        shared_yaxes=False,
        vertical_spacing=0.02,
        row_heights=row_heights_list,
        print_grid=show_grid,
    )

    add_Data_mode_sp(fig, float_map_inst)
    add_positioning_system_sp(fig, float_map_inst)
    max_prof_depth = add_profile_depth_sp(fig, float_map_inst)

    subplotNumList = []
    for i in range(len(float_map_inst["param_list"])):
        PARAM = float_map_inst["param_list"][i]
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

    out_file = output_filepath if output_filepath else f"/a2/ARGO_OPERATIONS/float_vis/drive/pictures/{float_map_inst['wmoID_str']}_DMQC_status.png"
    fig.write_image(
        out_file,
        format="png",
        width=1000,
        height=max(800, int(60 * num_rows + 150)),
        scale=2,
    )
    return fig


def generate_dmqc_status_plot(WMOfloatid, processed_bd_files, output_dir):
    """
    Directly harvests output BD NetCDF files to generate the verified DMQC Status Plot.
    Properly extracts all sub-profiles so PH, NITRATE, DOXY, CHLA, and FLUOR grades populate.
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
            pres_vals, _ = get_pres_and_qc_for_profile(ds, p_pres_idx)
            for p_v in pres_vals:
                if not is_missing_value(p_v):
                    nprof_nlevel_records.append({"N_PROF": prof_idx, "PRES_ADJUSTED": p_v})

            for p in candidate_params:
                p_prof = detect_parameter_profile(ds, p)

                var_exists = (
                    p in ds.variables or f"{p}_ADJUSTED" in ds.variables or f"PROFILE_{p}_QC" in ds.variables
                )
                if not var_exists and p == "CHLA_FLUORESCENCE":
                    var_exists = "FLUORESCENCE_CHLA" in ds.variables or "PROFILE_FLUORESCENCE_CHLA_QC" in ds.variables
                elif not var_exists and p == "PH_IN_SITU_TOTAL":
                    var_exists = "PH" in ds.variables or "PROFILE_PH_QC" in ds.variables

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
                qc_var_names = [f"PROFILE_{p}_QC"]
                if p == "CHLA_FLUORESCENCE":
                    qc_var_names.append("PROFILE_FLUORESCENCE_CHLA_QC")
                elif p == "PH_IN_SITU_TOTAL":
                    qc_var_names.append("PROFILE_PH_QC")

                # Harvest profile QC grade directly across sub-profiles
                for qvn in qc_var_names:
                    if qvn in ds.variables:
                        g_var = ds.variables[qvn]
                        raw_g_slice = g_var[:]
                        for candidate_idx in [p_prof, 0]:
                            if candidate_idx < len(raw_g_slice):
                                cand = raw_g_slice[candidate_idx]
                                if isinstance(cand, np.ndarray) and cand.ndim > 0:
                                    cand = cand[0]
                                c_str = cand.decode("utf-8", errors="ignore") if isinstance(cand, (bytes, np.bytes_)) else str(cand)
                                if c_str.strip() != "":
                                    grade = c_str.strip()
                                    break
                    if grade != " ":
                        break

                # If missing or blank, compute grade dynamically using Table 2a
                if grade == " " and target_qc in ds.variables:
                    qc_slice = ds.variables[target_qc][p_prof, :] if ds.variables[target_qc].ndim > 1 else ds.variables[target_qc][:]
                    grade = get_profile_qc_grade(qc_slice, parameter_present=has_data)

                if p == "PRES":
                    if "PRES_QC" in ds.variables:
                        pres_qc_slice = ds.variables["PRES_QC"][p_prof, :] if ds.variables["PRES_QC"].ndim > 1 else ds.variables["PRES_QC"][:]
                        grade = get_profile_qc_grade(pres_qc_slice, parameter_present=True)

                qc_col_key = f"PROFILE_{p}_QC"
                if qc_col_key not in nprof_df_data:
                    nprof_df_data[qc_col_key] = {}
                nprof_df_data[qc_col_key][prof_idx] = grade

                # Mirror aliases
                if p == "CHLA_FLUORESCENCE":
                    if "PROFILE_FLUORESCENCE_CHLA_QC" not in nprof_df_data:
                        nprof_df_data["PROFILE_FLUORESCENCE_CHLA_QC"] = {}
                    nprof_df_data["PROFILE_FLUORESCENCE_CHLA_QC"][prof_idx] = grade
                elif p == "PH_IN_SITU_TOTAL":
                    if "PROFILE_PH_QC" not in nprof_df_data:
                        nprof_df_data["PROFILE_PH_QC"] = {}
                    nprof_df_data["PROFILE_PH_QC"][prof_idx] = grade

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
    clean_df = clean_df[~clean_df[val_col].isin([MISSING_VAL, 99999.0])]

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
    """
    Renders exactly 2 panels:
      Top: Scientific section profile (valid data only) with missing data overlaid as faint 'x' markers.
      Bottom: Discrete QC flags section map (all active flags including Flag 9).
    """
    if var_name not in df_base.columns:
        return

    df_plot = df_base.copy()
    if qc_var_name in df_plot.columns:
        df_plot["QC_VAL"] = (
            df_plot[qc_var_name]
            .astype(str)
            .str.extract(r"([0-9])")[0]
            .fillna("9")
        )
    else:
        df_plot["QC_VAL"] = "9"

    is_missing_mask = (
        df_plot[var_name].isna()
        | (df_plot[var_name] == MISSING_VAL)
        | (df_plot[var_name] == 99999.0)
        | (df_plot["QC_VAL"] == "9")
    )

    df_valid = df_plot[~is_missing_mask].copy()
    df_missing = df_plot[is_missing_mask].copy()

    df_valid = df_valid.dropna(subset=[var_name, "PRES", "CYCLE_NUMBER"]).copy()
    df_valid = df_valid[(df_valid["PRES"] != MISSING_VAL) & (df_valid["PRES"] != 99999.0)]

    df_missing = df_missing.dropna(subset=["PRES", "CYCLE_NUMBER"]).copy()
    df_missing = df_missing[(df_missing["PRES"] != MISSING_VAL) & (df_missing["PRES"] != 99999.0)]

    if df_valid.empty and df_missing.empty:
        return

    interp_res = create_interpolated_grid(df_valid, var_name) if not df_valid.empty else None
    
    pres_all = df_plot["PRES"].replace([MISSING_VAL, 99999.0], np.nan)
    max_pres = np.nanmax(pres_all) if not pres_all.dropna().empty else 2000.0
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

    # --- TOP PANEL: Scientific Contour or Scatter Plot ---
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
    elif not df_valid.empty:
        sc = ax_data.scatter(df_valid["CYCLE_NUMBER"], df_valid["PRES"], c=df_valid[var_name], cmap=col_palette, s=15, alpha=0.8)
        cbar = fig.colorbar(sc, cax=cax_data)
        cbar.set_label(unit_str, fontsize=10)
        ax_data.set_title(f"Float {WMOfloatid} - {title_prefix} Section Profile", fontsize=12, fontweight="bold")
    else:
        cax_data.axis("off")
        ax_data.set_title(f"Float {WMOfloatid} - {title_prefix} Section Profile (No Valid Data)", fontsize=12, fontweight="bold")

    # OVERLAY MISSING DATA (Faint 'x' directly over section plot)
    if not df_missing.empty:
        ax_data.scatter(
            df_missing["CYCLE_NUMBER"],
            df_missing["PRES"],
            marker="x",
            c="#7f7f7f",
            s=12,
            alpha=0.35,
            linewidths=0.7,
            label="Missing Data (-99999)",
        )
        ax_data.legend(loc="upper right", fontsize=8, frameon=True, facecolor="white", edgecolor="none")

    ax_data.set_ylim(max_pres, 0)
    ax_data.set_ylabel("Pressure (dbar)", fontsize=11)
    ax_data.grid(True, linestyle="--", alpha=0.4)

    # --- BOTTOM PANEL: Discrete QC Flag Map ---
    unique_qcs = sorted(df_plot["QC_VAL"].unique(), key=lambda q: int(q) if q.isdigit() else 99)
    for q_flag in unique_qcs:
        sub = df_plot[df_plot["QC_VAL"] == q_flag]
        sub = sub[(sub["PRES"] != MISSING_VAL) & (sub["PRES"] != 99999.0)]
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
    records = []

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            cycle_num = int(fpath[-6:-3])
            iprof = detect_parameter_profile(ds, param_name)

            raw_v_name = param_name
            if raw_v_name not in ds.variables and param_name == "CHLA_FLUORESCENCE":
                raw_v_name = "FLUORESCENCE_CHLA"
            elif raw_v_name not in ds.variables and param_name == "PH_IN_SITU_TOTAL":
                raw_v_name = "PH"

            adj_v_name = adj_param_name if adj_param_name else f"{param_name}_ADJUSTED"
            if adj_v_name not in ds.variables and "FLUORESCENCE_CHLA" in raw_v_name:
                adj_v_name = "FLUORESCENCE_CHLA_ADJUSTED"
            elif adj_v_name not in ds.variables and raw_v_name == "PH":
                adj_v_name = "PH_ADJUSTED"

            raw_qc_name = f"{raw_v_name}_QC"
            adj_qc_name = f"{adj_v_name}_QC"

            pres_arr, _ = get_pres_and_qc_for_profile(ds, iprof)
            n_levels = len(pres_arr)

            def get_arr(v_name):
                if v_name not in ds.variables:
                    return None
                v = ds.variables[v_name]
                arr = v[iprof, :] if v.ndim > 1 else v[:]
                if hasattr(arr, "filled"):
                    arr = arr.filled(MISSING_VAL if arr.dtype.kind in ['f', 'i'] else b" ")
                return np.array(arr)

            raw_arr = get_arr(raw_v_name)
            raw_qc_arr = get_arr(raw_qc_name)
            adj_arr = get_arr(adj_v_name)
            adj_qc_arr = get_arr(adj_qc_name)

            for k in range(n_levels):
                p_val = float(pres_arr[k])
                if is_missing_value(p_val):
                    continue

                def parse_val(arr):
                    if arr is None or k >= len(arr):
                        return np.nan
                    val = arr[k]
                    try:
                        f = float(val)
                        return np.nan if is_missing_value(f) else f
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
                    return val if val in VALID_QC_CHARS else "9"

                records.append({
                    "CYCLE_NUMBER": cycle_num,
                    "PRES": p_val,
                    param_name: parse_val(raw_arr),
                    f"{param_name}_QC": parse_qc(raw_qc_arr),
                    f"{param_name}_ADJUSTED": parse_val(adj_arr),
                    f"{param_name}_ADJUSTED_QC": parse_qc(adj_qc_arr),
                })

            ds.close()
        except Exception:
            traceback.print_exc()

    return pd.DataFrame(records)


def generate_all_qc_summary_plots(processed_bd_files, WMOfloatid, output_root):
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

            # 1. Update PARAMETER_DATA_MODE and DATA_MODE (pH and Nitrate explicitly skipped)
            write_parameter_data_modes(ds)

            # 2. Standardize missing data in raw variables to -99999.0 and QC 9 (only change to raw data)
            standardize_raw_missing_data(ds)

            # 3. Process calibrations and CSV adjustments
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_chla)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio, iprof_bbp)
            write_DOXY_slope_drift(ds, 0, df_bio, idx_profile, WMOfloatid=WMOfloatid)
            write_DOXY_from_csv(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)

            # 4. Clean adjusted QC variables and update PROFILE_<PARAM>_QC (pH and Nitrate skipped)
            clean_and_fill_qc_variables(ds, WMOfloatid=WMOfloatid)

            # 5. Metadata and attributes
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

    # Step 6: Finalize output files
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

    # Step 7: Verification Plotting (reading directly from output BD files)
    generate_dmqc_status_plot(WMOfloatid, final_processed_files, script_dir)
    generate_all_qc_summary_plots(final_processed_files, WMOfloatid, script_dir)