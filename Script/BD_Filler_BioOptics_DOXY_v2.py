#!/usr/bin/env python3
# coding: utf-8

# ==============================================================================
# SECTION 1: Library Imports & Global Constants
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
    # 4903624 excluded from analysis
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

# Central storage location for input B-files
INPUT_STORAGE_DIR = "/a1/ARGO_DELAY/DMQC_BGC/storage/"

# Output destination base directory & dynamic date-stamped folder (MMDDYYYY)
BASE_OUTPUT_DIR = "/data/a1/ARGO_DELAY/DMQC_BGC/BDfiller/"
TODAY_FOLDER_STR = dt.now().strftime("%m%d%Y")
OUTPUT_PLOTS_CSV_DIR = os.path.join(BASE_OUTPUT_DIR, f"BD_Output_Plots_{TODAY_FOLDER_STR}")
os.makedirs(OUTPUT_PLOTS_CSV_DIR, exist_ok=True)

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
VALID_DATA_MODES = {"R", "A", "D"}
MISSING_VAL = 99999.0

TARGET_R_TO_D_SYNC_PARAMS = {
    "CHLA",
    "CHLA_FLUORESCENCE",
    "FLUORESCENCE_CHLA",
    "DOXY",
    "BBP700",
}


def write_nc_char_1d(var, iprof, char_str):
    """Safely writes a 1D string of length N_LEVELS / N_PARAM into a 2D NC_CHAR variable."""
    n_len = var.shape[1]
    padded = char_str.ljust(n_len)[:n_len]
    if var.dtype.kind == "S":
        var[iprof, :] = nc.stringtochar(np.array(padded, dtype=f"S{n_len}"))
    else:
        var[iprof, :] = list(padded)


def is_missing_value(val):
    if val is None or np.ma.is_masked(val) or pd.isna(val):
        return True
    try:
        f = float(val)
        return np.isnan(f) or np.isinf(f) or f == -99999.0 or f == 99999.0 or f == -999.0
    except (ValueError, TypeError):
        return False


def is_missing_qc(qc_val):
    if qc_val is None or np.ma.is_masked(qc_val):
        return True
    if isinstance(qc_val, (bytes, np.bytes_)):
        s = qc_val.decode("utf-8", errors="ignore").strip()
    else:
        s = str(qc_val).strip()
    return s in ["", " ", "9"] or s not in VALID_QC_CHARS


def is_blank_qc(qc_val):
    """Returns True if the QC character is a blank (' ' or '') or masked."""
    if qc_val is None or np.ma.is_masked(qc_val):
        return True
    if isinstance(qc_val, (bytes, np.bytes_)):
        s = qc_val.decode("utf-8", errors="ignore").strip()
    else:
        s = str(qc_val).strip()
    return s == ""


def get_r_mode_missing_mask(ds, iprof, base_param, n_levels):
    """
    Returns a boolean array of length `n_levels` where True indicates that
    the R-mode variable (<PARAM>) OR its R-mode QC (<PARAM>_QC) at that
    pressure level is missing / fill value.
    """
    missing_mask = np.zeros(n_levels, dtype=bool)

    candidate_names = [base_param]
    if base_param in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        candidate_names = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]

    raw_var_name = next((c for c in candidate_names if c in ds.variables), None)
    raw_qc_name = f"{raw_var_name}_QC" if raw_var_name and f"{raw_var_name}_QC" in ds.variables else None

    if raw_var_name is not None:
        v = ds.variables[raw_var_name]
        arr = v[iprof, :] if v.ndim > 1 else v[:]
        vals = arr.filled(MISSING_VAL) if hasattr(arr, "filled") else np.array(arr, dtype="float64")
        val_miss = (
            np.isnan(vals)
            | np.isinf(vals)
            | (vals == MISSING_VAL)
            | (vals == -99999.0)
            | (vals == 99999.0)
            | (vals == -999.0)
        )
        missing_mask |= val_miss[:n_levels]
    else:
        missing_mask[:] = True

    if raw_qc_name is not None:
        qv = ds.variables[raw_qc_name]
        q_arr = qv[iprof, :] if qv.ndim > 1 else qv[:]
        if hasattr(q_arr, "filled"):
            q_arr = q_arr.filled(b"9")
        for k in range(min(n_levels, len(q_arr))):
            if is_missing_qc(q_arr[k]):
                missing_mask[k] = True

    return missing_mask


def get_profile_max_level_idx(ds, iprof, base_param=None):
    """
    Finds the 0-based index (`max_k`) of the deepest actual level (`nlev - 1`)
    for `iprof` across PRES and all parameters on that profile.
    - Levels `0 .. max_k` are inside the profile depth range (missing values -> '9').
    - Levels `max_k + 1 .. N_LEVELS - 1` are array padding beyond nlev (must be blank ' ').
    """
    n_levels = ds.dimensions["N_LEVELS"].size
    pres_vals, pres_qcs = get_pres_and_qc_for_profile(ds, iprof)

    # 1. Valid PRES values or non-blank PRES_QC
    combined_active = (
        (~np.isnan(pres_vals))
        & (pres_vals != MISSING_VAL)
        & (pres_vals != -99999.0)
        & (pres_vals != 99999.0)
        & (pres_vals != -999.0)
    )
    for k in range(min(n_levels, len(pres_qcs))):
        if not is_blank_qc(pres_qcs[k]):
            combined_active[k] = True

    # 2. Check non-blank R-mode <PARAM>_QC flags across ALL parameters in this iprof
    prof_params = get_station_parameters_for_profile(ds, iprof)
    if base_param and base_param not in prof_params:
        prof_params.append(base_param)

    for p in prof_params:
        candidates = (
            ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
            if p in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
            else [p]
        )
        for cand in candidates:
            qc_name = f"{cand}_QC"
            if qc_name in ds.variables:
                qv = ds.variables[qc_name]
                q_arr = qv[iprof, :] if qv.ndim > 1 else qv[:]
                for k in range(min(n_levels, len(q_arr))):
                    if not is_blank_qc(q_arr[k]):
                        combined_active[k] = True

    active_indices = np.where(combined_active)[0]
    return int(active_indices[-1]) if len(active_indices) > 0 else -1


# ==============================================================================
# SECTION 2: Architecture & NetCDF Helper Functions
# ==============================================================================

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


def is_broken_sensor_float(WMOfloatid, cycle_num):
    float_str = str(WMOfloatid)
    cycle_num = int(cycle_num)
    return (
        ((float_str == "2904010") & (cycle_num >= 49))
        | ((float_str == "2904011") & (cycle_num >= 24))
    )


def detect_parameter_profile(ds, param_prefix):
    prof = get_profile_for_parameter(ds, param_prefix)
    if prof is not None:
        return prof

    # Secondary fallback: inspect actual variable arrays if STATION_PARAMETERS is missing/incomplete
    target_prefixes = [param_prefix]
    if param_prefix in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_prefixes = ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
    elif param_prefix in ["PH", "PH_IN_SITU_TOTAL"]:
        target_prefixes = ["PH_IN_SITU_TOTAL", "PH"]
    elif param_prefix in ["NITRATE"]:
        target_prefixes = ["NITRATE"]

    for tp in target_prefixes:
        for candidate_var in [tp, f"{tp}_ADJUSTED"]:
            if candidate_var in ds.variables:
                var = ds.variables[candidate_var]
                if var.ndim > 1:
                    for iprof in range(var.shape[0]):
                        sl = var[iprof, :]
                        vals = sl.filled(MISSING_VAL) if hasattr(sl, "filled") else np.array(sl)
                        if (~pd.isna(vals) & (vals != MISSING_VAL) & (vals != 99999.0) & (vals != -999.0)).any():
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
            valid_mask = (~pd.isna(p_vals)) & (p_vals != MISSING_VAL) & (p_vals != 99999.0) & (p_vals != -999.0)
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
            mode_str = m.decode("utf-8").strip() if isinstance(m, (bytes, np.bytes_)) else str(m).strip()
            return mode_str if mode_str in VALID_DATA_MODES else "R"
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
            mode_str = m_val.decode("utf-8").strip() if isinstance(m_val, (bytes, np.bytes_)) else str(m_val).strip()
            return mode_str if mode_str in VALID_DATA_MODES else "R"

    # Parameter is not listed in STATION_PARAMETERS[iprof] -> must be treated as R/FillValue on this iprof
    return "R"


def sync_data_mode_from_parameter_data_mode(bgc_file):
    """
    Ensures DATA_MODE[iprof] is strictly consistent with PARAMETER_DATA_MODE[iprof, :]:
    - If any parameter on iprof is 'D' -> DATA_MODE[iprof] = 'D'
    - Else if any parameter on iprof is 'A' -> DATA_MODE[iprof] = 'A'
    - Else -> DATA_MODE[iprof] = 'R'
    """
    if "PARAMETER_DATA_MODE" not in bgc_file.variables or "DATA_MODE" not in bgc_file.variables:
        return

    n_prof = bgc_file.dimensions["N_PROF"].size
    pdm_var = bgc_file.variables["PARAMETER_DATA_MODE"]
    dm_var = bgc_file.variables["DATA_MODE"]
    pdm = pdm_var[:]

    for iprof in range(n_prof):
        row_modes = []
        for m in np.ravel(pdm[iprof]):
            s = m.decode("utf-8", errors="ignore").strip() if isinstance(m, (bytes, np.bytes_)) else str(m).strip()
            if s in VALID_DATA_MODES:
                row_modes.append(s)

        if "D" in row_modes:
            best = "D"
        elif "A" in row_modes:
            best = "A"
        else:
            best = "R"

        if dm_var.dtype.kind == "S":
            dm_var[iprof] = nc.stringtochar(np.array(best, dtype="S1"))
        else:
            dm_var[iprof] = best


def write_parameter_data_modes(bgc_file, WMOfloatid=None, cycle_df=None):
    n_param = bgc_file.dimensions["N_PARAM"].size
    n_prof = bgc_file.dimensions["N_PROF"].size
    skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

    pdm_var = bgc_file.variables["PARAMETER_DATA_MODE"]
    pdm = np.array(pdm_var[:], copy=True)

    csv_col_map = {
        "CHLA": ["CHLA_FINAL"],
        "CHLA_FLUORESCENCE": ["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"],
        "FLUORESCENCE_CHLA": ["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"],
        "BBP700": ["BBP700_FINAL"],
        "DOXY": ["DOXY_FINAL"],
    }

    def has_csv_calibration(p_name):
        if cycle_df is None or cycle_df.empty:
            return False
        cols = csv_col_map.get(p_name, [])
        for c in cols:
            if c in cycle_df.columns:
                vals = cycle_df[c].dropna()
                if not vals.empty and (~vals.isin([-99999.0, 99999.0, -999.0])).any():
                    return True
        return False

    for iprof in range(n_prof):
        param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
        if isinstance(param_mat, np.ma.MaskedArray):
            param_mat = param_mat.filled(b" ")

        row_chars = []
        for j in range(n_param):
            p_bytes = param_mat[j]
            p_str = "".join([
                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                for c in p_bytes
            ]).strip()

            cur_m = pdm[iprof, j]
            cur_char = cur_m.decode("utf-8", errors="ignore") if isinstance(cur_m, (bytes, np.bytes_)) else str(cur_m)
            if not cur_char:
                cur_char = " "

            if not p_str or "PH" in p_str or "NITRATE" in p_str:
                row_chars.append(cur_char[0])
                continue

            if skip_doxy and p_str == "DOXY":
                row_chars.append(cur_char[0])
                continue

            if p_str == "PRES" or p_str.startswith("DOWN_IRRADIANCE") or p_str.startswith("UP_RADIANCE"):
                row_chars.append("R")
            elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY"]:
                # Only switch to 'D' if both the NetCDF profile AND the CSV contain valid data
                if has_valid_measurements(bgc_file, p_str, iprof) and has_csv_calibration(p_str):
                    row_chars.append("D")
                else:
                    row_chars.append("R")
            else:
                row_chars.append(cur_char[0])

        write_nc_char_1d(pdm_var, iprof, "".join(row_chars))

    sync_data_mode_from_parameter_data_mode(bgc_file)


def get_profile_qc_grade(qc_array, parameter_present=True):
    if not parameter_present:
        return " "

    valid_levels = []
    for q in np.ravel(qc_array):
        if isinstance(q, np.ma.core.MaskedConstant) or q is np.ma.masked:
            continue
        s = q.decode("utf-8", errors="ignore").strip() if isinstance(q, (bytes, np.bytes_)) else str(q).strip()
        if s not in ["", " ", "9"]:
            valid_levels.append(s)

    total_levels = len(valid_levels)
    if total_levels == 0:
        return " "

    good_qc = {"1", "2", "5"}
    good_count = sum(1 for q in valid_levels if q in good_qc)
    good_pct = (good_count / total_levels) * 100.0

    if good_count == total_levels:
        return "A"
    elif good_pct >= 75.0:
        return "B"
    elif good_pct >= 50.0:
        return "C"
    elif good_pct >= 25.0:
        return "D"
    elif good_pct > 0.0:
        return "E"
    else:
        return "F"


def get_station_parameters_for_profile(ds, iprof):
    """
    Extracts cleaned parameter name strings for a specific profile (N_PROF index)
    directly from STATION_PARAMETERS.
    """
    if "STATION_PARAMETERS" not in ds.variables:
        return []

    param_mat = ds.variables["STATION_PARAMETERS"][iprof]
    if isinstance(param_mat, np.ma.MaskedArray):
        param_mat = param_mat.filled(b" ")

    param_names = []
    n_param = ds.dimensions["N_PARAM"].size
    for j in range(n_param):
        p_bytes = param_mat[j]
        p_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
            for c in p_bytes
        ]).strip()
        if p_str:
            param_names.append(p_str)
    return param_names


def get_profile_for_parameter(ds, param_base):
    """
    Finds the exact N_PROF index that contains the given parameter according
    to STATION_PARAMETERS. If not present in STATION_PARAMETERS, returns None.
    Also handles parameter aliases (e.g. FLUORESCENCE_CHLA vs CHLA_FLUORESCENCE).
    """
    if "STATION_PARAMETERS" not in ds.variables:
        return None

    target_names = {param_base}
    if param_base in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_names = {"CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"}
    elif param_base in ["PH", "PH_IN_SITU_TOTAL"]:
        target_names = {"PH_IN_SITU_TOTAL", "PH"}
    elif param_base == "NITRATE":
        target_names = {"NITRATE"}

    n_prof = ds.dimensions["N_PROF"].size
    for iprof in range(n_prof):
        prof_params = get_station_parameters_for_profile(ds, iprof)
        for p in prof_params:
            if p in target_names or any(p.startswith(t) for t in target_names):
                return iprof
    return None


def get_pres_and_qc_for_profile(ds, iprof):
    """
    Retrieves PRES strictly for `iprof`.
    Returns pres_vals (float32 array, with missing values set to NaN) and pres_qcs.
    """
    n_levels = ds.dimensions["N_LEVELS"].size
    pres_vals = np.full((n_levels,), np.nan, dtype="float32")
    pres_qcs = np.full((n_levels,), b" ", dtype="|S1")

    # Verify PRES is actually listed for this specific profile in STATION_PARAMETERS
    prof_params = get_station_parameters_for_profile(ds, iprof)
    if "PRES" not in prof_params and not any(p.startswith("PRES") for p in prof_params):
        return pres_vals, pres_qcs

    p_var_name = "PRES_ADJUSTED" if "PRES_ADJUSTED" in ds.variables else ("PRES" if "PRES" in ds.variables else None)
    if p_var_name:
        v = ds.variables[p_var_name]
        arr = v[iprof, :] if v.ndim > 1 else v[:]
        p_raw = arr.filled(np.nan) if hasattr(arr, "filled") else np.array(arr, dtype="float32")
        # Ensure NetCDF fill values are converted to NaN
        mask = (p_raw == MISSING_VAL) | (p_raw == 99999.0) | (p_raw == -99999.0) | (p_raw == -999.0) | np.isinf(p_raw)
        p_raw[mask] = np.nan
        pres_vals = np.float32(p_raw)

    q_var_name = "PRES_ADJUSTED_QC" if "PRES_ADJUSTED_QC" in ds.variables else ("PRES_QC" if "PRES_QC" in ds.variables else None)
    if q_var_name:
        qv = ds.variables[q_var_name]
        q_slice = qv[iprof, :] if qv.ndim > 1 else qv[:]
        if hasattr(q_slice, "filled"):
            q_slice = q_slice.filled(b" ")
        pres_qcs = np.array(q_slice).astype("|S1")

    return pres_vals, pres_qcs


def match_csv_row_by_pres(cycle_df, nc_pres):
    """
    Matches the closest pressure in the CSV for the given NetCDF pressure level.
    Uses adaptive tolerance based on depth.
    """
    if cycle_df.empty or is_missing_value(nc_pres):
        return None

    csv_pres = cycle_df["PRES"].values
    diffs = np.abs(csv_pres - nc_pres)
    min_idx = np.argmin(diffs)

    # Allow generous tolerance so discrete sample depths match their calibration rows
    if nc_pres <= 50.0:
        tol = 10.0
    elif nc_pres <= 200.0:
        tol = 25.0
    elif nc_pres <= 600.0:
        tol = 60.0
    else:
        tol = 150.0

    if diffs[min_idx] <= tol:
        return cycle_df.iloc[min_idx]
    return None


def sanitize_qc_char(val, default="1"):
    if pd.isna(val) or val is None or val == "":
        return default
    s = str(val).strip()
    match = re.search(r"([0-9])", s)
    if match and match.group(1) in VALID_QC_CHARS:
        return match.group(1)
    return default


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

    prof_chla = get_profile_for_parameter(bgc_file, "CHLA")
    prof_fluo = get_profile_for_parameter(bgc_file, "CHLA_FLUORESCENCE")
    prof_bbp = get_profile_for_parameter(bgc_file, "BBP700")

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

            if param_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"] and iprof_idx == prof_fluo:
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA_FLU
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "CHLA" and iprof_idx == prof_chla:
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "BBP700" and iprof_idx == prof_bbp:
                bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_BBP
                bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray


def write_chla_BBP_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio):
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    if cycle_df.empty:
        return

    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)
    fluo_csv_col = next((c for c in ["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"] if c in df_bio.columns), None)
    fluo_qc_csv_col = next((c for c in ["CHLA_FLUORESCENCE_ADJUSTED_QC", "FLUORESCENCE_CHLA_ADJUSTED_QC", "CHLA_FLUORESCENCE_QC", "FLUORESCENCE_CHLA_QC"] if c in df_bio.columns), None)

    prof_chla = get_profile_for_parameter(bgc_file, "CHLA")
    if prof_chla is not None:
        mode = get_param_data_mode(bgc_file, prof_chla, "CHLA")
        if mode in ["D", "A"]:
            nc_pres_all, _ = get_pres_and_qc_for_profile(bgc_file, prof_chla)
            r_miss_mask = get_r_mode_missing_mask(bgc_file, prof_chla, "CHLA", n_levels)
            max_k = get_profile_max_level_idx(bgc_file, prof_chla, "CHLA")

            if max_k >= 0:
                CHLA_Adjusted_Array = np.array(bgc_file.variables["CHLA_ADJUSTED"][prof_chla, :], copy=True) if "CHLA_ADJUSTED" in bgc_file.variables else None
                CHLA_AdjustedQC_Array = np.array(bgc_file.variables["CHLA_ADJUSTED_QC"][prof_chla, :], copy=True) if "CHLA_ADJUSTED_QC" in bgc_file.variables else None
                CHLA_Adjusted_ERROR_Array = np.array(bgc_file.variables["CHLA_ADJUSTED_ERROR"][prof_chla, :], copy=True) if "CHLA_ADJUSTED_ERROR" in bgc_file.variables else None

                for i in range(max_k + 1):
                    nc_pres = nc_pres_all[i]
                    if is_missing_value(nc_pres) or np.isnan(nc_pres) or r_miss_mask[i]:
                        if CHLA_Adjusted_Array is not None:
                            CHLA_Adjusted_Array[i] = np.float32(MISSING_VAL)
                        if CHLA_AdjustedQC_Array is not None:
                            CHLA_AdjustedQC_Array[i] = b"9"
                        if CHLA_Adjusted_ERROR_Array is not None:
                            CHLA_Adjusted_ERROR_Array[i] = np.float32(MISSING_VAL)
                        continue

                    matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                    if matched_row is None:
                        continue

                    if "CHLA_FINAL" in matched_row and pd.notna(matched_row.get("CHLA_FINAL")):
                        chla_val = matched_row.get("CHLA_FINAL")
                        csv_qc = "4" if is_broken else sanitize_qc_char(matched_row.get("CHLA_FINAL_QC"), default="1")
                        if is_missing_value(chla_val) or csv_qc in ["4", "9"]:
                            if CHLA_Adjusted_Array is not None:
                                CHLA_Adjusted_Array[i] = np.float32(MISSING_VAL)
                            if CHLA_AdjustedQC_Array is not None:
                                CHLA_AdjustedQC_Array[i] = (csv_qc if csv_qc in ["4", "9"] else "9").encode("utf-8")
                            if CHLA_Adjusted_ERROR_Array is not None:
                                CHLA_Adjusted_ERROR_Array[i] = np.float32(MISSING_VAL)
                        else:
                            if CHLA_Adjusted_Array is not None:
                                CHLA_Adjusted_Array[i] = np.float32(chla_val)
                            if CHLA_AdjustedQC_Array is not None:
                                CHLA_AdjustedQC_Array[i] = csv_qc.encode("utf-8")
                            if CHLA_Adjusted_ERROR_Array is not None:
                                CHLA_Adjusted_ERROR_Array[i] = np.float32(CHLA_Adjusted_ERROR_est)

                if CHLA_Adjusted_Array is not None:
                    bgc_file.variables["CHLA_ADJUSTED"][prof_chla, :] = CHLA_Adjusted_Array
                if CHLA_AdjustedQC_Array is not None:
                    qc_str = "".join([
                        (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))[:1] or " "
                        for c in CHLA_AdjustedQC_Array
                    ])
                    write_nc_char_1d(bgc_file.variables["CHLA_ADJUSTED_QC"], prof_chla, qc_str)
                if CHLA_Adjusted_ERROR_Array is not None:
                    bgc_file.variables["CHLA_ADJUSTED_ERROR"][prof_chla, :] = CHLA_Adjusted_ERROR_Array

    # Same pattern for FLUORESCENCE_CHLA / CHLA_FLUORESCENCE
    prof_fluo = get_profile_for_parameter(bgc_file, "CHLA_FLUORESCENCE")
    if prof_fluo is not None and fluo_csv_col is not None:
        mode = get_param_data_mode(bgc_file, prof_fluo, "CHLA_FLUORESCENCE")
        if mode in ["D", "A"]:
            nc_pres_all, _ = get_pres_and_qc_for_profile(bgc_file, prof_fluo)
            r_miss_mask = get_r_mode_missing_mask(bgc_file, prof_fluo, "CHLA_FLUORESCENCE", n_levels)
            max_k = get_profile_max_level_idx(bgc_file, prof_fluo, "CHLA_FLUORESCENCE")

            if max_k >= 0:
                for var_prefix in ["FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]:
                    adj_var = f"{var_prefix}_ADJUSTED"
                    adj_qc_var = f"{var_prefix}_ADJUSTED_QC"
                    adj_err_var = f"{var_prefix}_ADJUSTED_ERROR"

                    if adj_var in bgc_file.variables:
                        fluo_adj = np.array(bgc_file.variables[adj_var][prof_fluo, :], copy=True)
                        fluo_qc = np.array(bgc_file.variables[adj_qc_var][prof_fluo, :], copy=True) if adj_qc_var in bgc_file.variables else None
                        fluo_err = np.array(bgc_file.variables[adj_err_var][prof_fluo, :], copy=True) if adj_err_var in bgc_file.variables else None

                        for i in range(max_k + 1):
                            nc_pres = nc_pres_all[i]
                            if is_missing_value(nc_pres) or np.isnan(nc_pres) or r_miss_mask[i]:
                                fluo_adj[i] = np.float32(MISSING_VAL)
                                if fluo_qc is not None:
                                    fluo_qc[i] = b"9"
                                if fluo_err is not None:
                                    fluo_err[i] = np.float32(MISSING_VAL)
                                continue

                            matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
                            if matched_row is None:
                                continue

                            val = matched_row.get(fluo_csv_col)
                            if pd.notna(val):
                                q = "4" if is_broken else sanitize_qc_char(matched_row.get(fluo_qc_csv_col), default="1")
                                if is_missing_value(val) or q in ["4", "9"]:
                                    fluo_adj[i] = np.float32(MISSING_VAL)
                                    if fluo_qc is not None:
                                        fluo_qc[i] = (q if q in ["4", "9"] else "9").encode("utf-8")
                                    if fluo_err is not None:
                                        fluo_err[i] = np.float32(MISSING_VAL)
                                else:
                                    fluo_adj[i] = np.float32(val)
                                    if fluo_qc is not None:
                                        fluo_qc[i] = q.encode("utf-8")
                                    if fluo_err is not None:
                                        fluo_err[i] = np.float32(CHLA_Adjusted_ERROR_est)

                        bgc_file.variables[adj_var][prof_fluo, :] = fluo_adj
                        if fluo_qc is not None:
                            qc_str = "".join([
                                (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))[:1] or " "
                                for c in fluo_qc
                            ])
                            write_nc_char_1d(bgc_file.variables[adj_qc_var], prof_fluo, qc_str)
                        if fluo_err is not None:
                            bgc_file.variables[adj_err_var][prof_fluo, :] = fluo_err


def write_BBP700_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio):
    prof_bbp = get_profile_for_parameter(bgc_file, "BBP700")
    if prof_bbp is None or "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    mode = get_param_data_mode(bgc_file, prof_bbp, "BBP700")
    if mode not in ["D", "A"]:
        return

    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    if cycle_df.empty:
        return

    n_levels = bgc_file.dimensions["N_LEVELS"].size
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)
    nc_pres_all, _ = get_pres_and_qc_for_profile(bgc_file, prof_bbp)
    r_miss_mask = get_r_mode_missing_mask(bgc_file, prof_bbp, "BBP700", n_levels)
    max_k = get_profile_max_level_idx(bgc_file, prof_bbp, "BBP700")

    if max_k < 0:
        return

    bbp_adj = np.array(bgc_file.variables["BBP700_ADJUSTED"][prof_bbp, :], copy=True)
    bbp_qc = np.array(bgc_file.variables["BBP700_ADJUSTED_QC"][prof_bbp, :], copy=True) if "BBP700_ADJUSTED_QC" in bgc_file.variables else None
    bbp_err = np.array(bgc_file.variables["BBP700_ADJUSTED_ERROR"][prof_bbp, :], copy=True) if "BBP700_ADJUSTED_ERROR" in bgc_file.variables else None

    for i in range(max_k + 1):
        nc_pres = nc_pres_all[i]
        if is_missing_value(nc_pres) or r_miss_mask[i]:
            bbp_adj[i] = np.float32(MISSING_VAL)
            if bbp_qc is not None:
                bbp_qc[i] = b"9"
            if bbp_err is not None:
                bbp_err[i] = np.float32(MISSING_VAL)
            continue

        matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
        if matched_row is None:
            continue

        if "BBP700_FINAL" in matched_row and pd.notna(matched_row.get("BBP700_FINAL")):
            val = matched_row.get("BBP700_FINAL")
            q = "4" if is_broken else sanitize_qc_char(matched_row.get("BBP700_FINAL_QC"), default="1")
            if is_missing_value(val) or q in ["4", "9"]:
                bbp_adj[i] = np.float32(MISSING_VAL)
                if bbp_qc is not None:
                    bbp_qc[i] = (q if q in ["4", "9"] else "9").encode("utf-8")
                if bbp_err is not None:
                    bbp_err[i] = np.float32(MISSING_VAL)
            else:
                bbp_adj[i] = np.float32(val)
                if bbp_qc is not None:
                    bbp_qc[i] = q.encode("utf-8")
                if bbp_err is not None:
                    bbp_err[i] = np.float32(BBP700_Adjusted_ERROR_est)

    bgc_file.variables["BBP700_ADJUSTED"][prof_bbp, :] = bbp_adj
    if bbp_qc is not None:
        qc_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))[:1] or " "
            for c in bbp_qc
        ])
        write_nc_char_1d(bgc_file.variables["BBP700_ADJUSTED_QC"], prof_bbp, qc_str)
    if bbp_err is not None:
        bgc_file.variables["BBP700_ADJUSTED_ERROR"][prof_bbp, :] = bbp_err


def write_DOXY_slope_drift(ds, float_df, target_cycle, WMOfloatid=None):
    if WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS:
        return

    prof_doxy = get_profile_for_parameter(ds, "DOXY")
    if prof_doxy is None:
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

    n_param = ds.dimensions["N_PARAM"].size
    station_params_mat = ds.variables["STATION_PARAMETERS"][prof_doxy]
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
                ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"][prof_doxy, 0, j, :] = char_coef
            if "SCIENTIFIC_CALIB_EQUATION" in ds.variables:
                ds.variables["SCIENTIFIC_CALIB_EQUATION"][prof_doxy, 0, j, :] = char_equ
            if "SCIENTIFIC_CALIB_COMMENT" in ds.variables:
                ds.variables["SCIENTIFIC_CALIB_COMMENT"][prof_doxy, 0, j, :] = char_com
            if "SCIENTIFIC_CALIB_DATE" in ds.variables:
                ds.variables["SCIENTIFIC_CALIB_DATE"][prof_doxy, 0, j, :] = char_date


def write_DOXY_from_csv(ds, float_df, target_cycle, WMOfloatid=None):
    if WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS:
        return

    prof_doxy = get_profile_for_parameter(ds, "DOXY")
    if prof_doxy is None:
        return

    mode = get_param_data_mode(ds, prof_doxy, "DOXY")
    if mode not in ["D", "A"]:
        return

    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()
    if cycle_df.empty:
        return

    n_levels = ds.dimensions["N_LEVELS"].size
    nc_pres_all, _ = get_pres_and_qc_for_profile(ds, prof_doxy)
    r_miss_mask = get_r_mode_missing_mask(ds, prof_doxy, "DOXY", n_levels)
    max_k = get_profile_max_level_idx(ds, prof_doxy, "DOXY")

    if max_k < 0:
        return

    doxy_adj = np.array(ds.variables["DOXY_ADJUSTED"][prof_doxy, :], copy=True) if "DOXY_ADJUSTED" in ds.variables else None
    doxy_qc = np.array(ds.variables["DOXY_ADJUSTED_QC"][prof_doxy, :], copy=True) if "DOXY_ADJUSTED_QC" in ds.variables else None
    doxy_err = np.array(ds.variables["DOXY_ADJUSTED_ERROR"][prof_doxy, :], copy=True) if "DOXY_ADJUSTED_ERROR" in ds.variables else None

    for i in range(max_k + 1):
        nc_pres = nc_pres_all[i]
        if is_missing_value(nc_pres) or r_miss_mask[i]:
            if doxy_adj is not None:
                doxy_adj[i] = np.float32(MISSING_VAL)
            if doxy_qc is not None:
                doxy_qc[i] = b"9"
            if doxy_err is not None:
                doxy_err[i] = np.float32(MISSING_VAL)
            continue

        matched_row = match_csv_row_by_pres(cycle_df, nc_pres)
        if matched_row is None:
            continue

        if "DOXY_FINAL" in matched_row and pd.notna(matched_row.get("DOXY_FINAL")):
            val = matched_row.get("DOXY_FINAL")
            q = sanitize_qc_char(matched_row.get("DOXY_FINAL_QC"), default="1")
            if is_missing_value(val) or q in ["4", "9"]:
                if doxy_adj is not None:
                    doxy_adj[i] = np.float32(MISSING_VAL)
                if doxy_qc is not None:
                    doxy_qc[i] = (q if q in ["4", "9"] else "9").encode("utf-8")
                if doxy_err is not None:
                    doxy_err[i] = np.float32(MISSING_VAL)
            else:
                if doxy_adj is not None:
                    doxy_adj[i] = np.float32(val)
                if doxy_qc is not None:
                    doxy_qc[i] = q.encode("utf-8")
                if doxy_err is not None:
                    err = matched_row.get("DOXY_ADJUSTED_ERROR")
                    doxy_err[i] = np.float32(err) if pd.notna(err) and not is_missing_value(err) else np.float32(DOXY_Adjusted_ERROR_est_default)

    if doxy_adj is not None:
        ds.variables["DOXY_ADJUSTED"][prof_doxy, :] = doxy_adj
    if doxy_qc is not None:
        qc_str = "".join([
            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))[:1] or " "
            for c in doxy_qc
        ])
        write_nc_char_1d(ds.variables["DOXY_ADJUSTED_QC"], prof_doxy, qc_str)
    if doxy_err is not None:
        ds.variables["DOXY_ADJUSTED_ERROR"][prof_doxy, :] = doxy_err


def update_profile_qc_grades(bgc_file, WMOfloatid=None):
    n_prof = bgc_file.dimensions["N_PROF"].size
    skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

    for var_name in list(bgc_file.variables.keys()):
        if not var_name.startswith("PROFILE_") or not var_name.endswith("_QC"):
            continue
        if "PH" in var_name or "NITRATE" in var_name:
            continue

        param_base = var_name[8:-3]
        if skip_doxy and param_base == "DOXY":
            continue

        target_prof = detect_parameter_profile(bgc_file, param_base)
        mode = get_param_data_mode(bgc_file, target_prof, param_base)

        adj_qc_name = f"{param_base}_ADJUSTED_QC"
        raw_qc_name = "PRES_QC" if param_base == "PRES" else f"{param_base}_QC"

        has_adj_flags = False
        if mode in ["D", "A"] and adj_qc_name in bgc_file.variables:
            adj_raw = bgc_file.variables[adj_qc_name][target_prof, :] if bgc_file.variables[adj_qc_name].ndim > 1 else bgc_file.variables[adj_qc_name][:]
            flags = [q.decode("utf-8", errors="ignore").strip() if isinstance(q, (bytes, np.bytes_)) else str(q).strip() for q in np.ravel(adj_raw)]
            if any(q in ["1", "2", "3", "4", "5", "8"] for q in flags):
                has_adj_flags = True

        target_qc_name = adj_qc_name if has_adj_flags else raw_qc_name
        param_is_present = has_valid_measurements(bgc_file, param_base, target_prof)

        if target_qc_name in bgc_file.variables:
            target_qc_var = bgc_file.variables[target_qc_name]
            target_qc_array = (
                target_qc_var[target_prof, :] if target_qc_var.ndim > 1 else target_qc_var[:]
            )
            prof_grade = get_profile_qc_grade(target_qc_array, parameter_present=param_is_present)
            grade_char = nc.stringtochar(np.array(prof_grade, dtype="S1"))

            for iprof in range(n_prof):
                if bgc_file.variables[var_name].ndim == 1:
                    bgc_file.variables[var_name][iprof] = grade_char if bgc_file.variables[var_name].dtype.kind == "S" else prof_grade
                elif bgc_file.variables[var_name].ndim == 2:
                    bgc_file.variables[var_name][iprof, 0] = grade_char if bgc_file.variables[var_name].dtype.kind == "S" else prof_grade


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


def standardize_missing_data_and_qc(ds):
    """
    Standardizes data and QC flags for all _ADJUSTED variables using nc.stringtochar
    for NC_CHAR variables:
    - Syncs DATA_MODE[iprof] with PARAMETER_DATA_MODE[iprof, :].
    - If DATA_MODE[iprof] is 'R' OR PARAMETER_DATA_MODE for (iprof, base_param) is NOT
      'D' or 'A': forces <PARAM>_ADJUSTED and <PARAM>_ADJUSTED_ERROR to _FillValue
      (99999.0), and <PARAM>_ADJUSTED_QC to blank (' ') across ALL N_LEVELS.
    - If PARAMETER_DATA_MODE for (iprof, base_param) IS 'D' or 'A':
      - k <= max_k (within profile nlev): missing/blank -> '9' and 99999.0; QC '4'/'9' -> 99999.0
      - k >  max_k (beyond profile nlev): padding -> ' ' and 99999.0
    - RAW variables (PRES, DOXY, CHLA, BBP700, etc.) and their _QC flags are NEVER modified.
    """
    ds.set_auto_mask(False)
    sync_data_mode_from_parameter_data_mode(ds)

    n_prof = ds.dimensions["N_PROF"].size
    n_levels = ds.dimensions["N_LEVELS"].size

    # Precompute valid PRES masks per iprof
    prof_pres_valid = {}
    for iprof in range(n_prof):
        pres_vals, _ = get_pres_and_qc_for_profile(ds, iprof)
        prof_pres_valid[iprof] = (
            (~np.isnan(pres_vals))
            & (pres_vals != MISSING_VAL)
            & (pres_vals != -99999.0)
            & (pres_vals != 99999.0)
            & (pres_vals != -999.0)
        )

    # Precompute R-mode missing masks per (iprof, base_param)
    r_mode_missing_cache = {}
    for iprof in range(n_prof):
        for p in TARGET_R_TO_D_SYNC_PARAMS:
            r_mode_missing_cache[(iprof, p)] = get_r_mode_missing_mask(ds, iprof, p, n_levels)

    for v_name in ds.variables.keys():
        # Strictly only operate on _ADJUSTED variables; never touch raw data, raw QC, or metadata
        if "_ADJUSTED" not in v_name or v_name.startswith("PRES_ADJUSTED"):
            continue

        var = ds.variables[v_name]
        if var.ndim != 2:
            continue

        base_param = (
            v_name.replace("_ADJUSTED_QC", "")
            .replace("_ADJUSTED_ERROR", "")
            .replace("_ADJUSTED", "")
        )

        # -------------------------------------------------------------
        # 1. Level QC flags (<PARAM>_ADJUSTED_QC):
        #    - DATA_MODE == 'R' OR Parameter mode not in ('D', 'A'): ALL levels -> ' '
        #    - Parameter mode 'D'/'A':
        #      - k <= max_k (within nlev): missing/blank -> '9'
        #      - k >  max_k (deeper than nlev): padding -> ' '
        # -------------------------------------------------------------
        if v_name.endswith("_ADJUSTED_QC"):
            for iprof in range(var.shape[0]):
                mode = get_param_data_mode(ds, iprof, base_param)
                prof_dm = ""
                if "DATA_MODE" in ds.variables:
                    dm_raw = ds.variables["DATA_MODE"][iprof]
                    prof_dm = (
                        dm_raw.decode("utf-8", errors="ignore").strip()
                        if isinstance(dm_raw, (bytes, np.bytes_))
                        else str(dm_raw).strip()
                    )

                # If DATA_MODE[iprof] is 'R' OR parameter mode is not 'D'/'A', ALL levels MUST be blank (' ')
                if prof_dm == "R" or mode not in ["D", "A"]:
                    write_nc_char_1d(var, iprof, " " * var.shape[1])
                    continue

                q_slice = np.array(var[iprof, :], copy=True)
                max_k = get_profile_max_level_idx(ds, iprof, base_param)
                pres_mask = prof_pres_valid[iprof]
                r_miss = r_mode_missing_cache.get((iprof, base_param))

                out_chars = []
                for k in range(var.shape[1]):
                    if max_k < 0 or k > max_k:
                        # Deeper than deepest pressure (nlev) of this nprof -> blank padding
                        out_chars.append(" ")
                    elif not pres_mask[k] or (r_miss is not None and r_miss[k]):
                        # Shallower than or at deepest pressure, but missing -> '9'
                        out_chars.append("9")
                    else:
                        c = (
                            q_slice[k].decode("utf-8", errors="ignore").strip()
                            if isinstance(q_slice[k], (bytes, np.bytes_))
                            else str(q_slice[k]).strip()
                        )
                        out_chars.append(c if (c in VALID_QC_CHARS and c != "") else "9")

                write_nc_char_1d(var, iprof, "".join(out_chars))
            continue

        # -------------------------------------------------------------
        # 2. Numeric Adjusted variables: _ADJUSTED and _ADJUSTED_ERROR
        # -------------------------------------------------------------
        if var.dtype.kind in ["f", "i"]:
            raw_fill = getattr(var, "_FillValue", MISSING_VAL)
            fill_val = float(np.ravel(raw_fill)[0]) if isinstance(raw_fill, (np.ndarray, list)) else float(raw_fill)
            adj_qc_var_name = f"{base_param}_ADJUSTED_QC"
            adj_qc_var = ds.variables.get(adj_qc_var_name)

            for iprof in range(var.shape[0]):
                mode = get_param_data_mode(ds, iprof, base_param)
                prof_dm = ""
                if "DATA_MODE" in ds.variables:
                    dm_raw = ds.variables["DATA_MODE"][iprof]
                    prof_dm = (
                        dm_raw.decode("utf-8", errors="ignore").strip()
                        if isinstance(dm_raw, (bytes, np.bytes_))
                        else str(dm_raw).strip()
                    )

                # If DATA_MODE[iprof] is 'R' OR parameter mode is not 'D'/'A', ALL levels MUST be _FillValue
                if prof_dm == "R" or mode not in ["D", "A"]:
                    var[iprof, :] = np.full(var.shape[1], fill_val, dtype=var.dtype)
                    continue

                max_k = get_profile_max_level_idx(ds, iprof, base_param)
                pres_mask = prof_pres_valid[iprof]
                r_miss = r_mode_missing_cache.get((iprof, base_param))
                data_slice = np.array(var[iprof, :], dtype=var.dtype, copy=True)
                qc_slice = np.array(adj_qc_var[iprof, :], copy=True) if adj_qc_var is not None else None
                qc_modified = False

                for k in range(len(data_slice)):
                    qc_char = ""
                    if qc_slice is not None and k < len(qc_slice):
                        qc_char = (
                            qc_slice[k].decode("utf-8", errors="ignore").strip()
                            if isinstance(qc_slice[k], (bytes, np.bytes_))
                            else str(qc_slice[k]).strip()
                        )

                    # Argo rule: QC = '4' (bad) or '9' (missing) MUST have _FillValue in _ADJUSTED & _ADJUSTED_ERROR
                    if (
                        max_k < 0
                        or k > max_k
                        or not pres_mask[k]
                        or (r_miss is not None and r_miss[k])
                        or qc_char in ["4", "9", ""]
                    ):
                        data_slice[k] = fill_val
                    else:
                        val = float(data_slice[k])
                        if np.isnan(val) or np.isinf(val) or val == -99999.0 or val == 99999.0 or val == -999.0:
                            data_slice[k] = fill_val
                            if qc_slice is not None and k <= max_k and qc_char not in ["4", "9"]:
                                qc_slice[k] = b"9"
                                qc_modified = True

                var[iprof, :] = data_slice
                if qc_modified and adj_qc_var is not None:
                    qc_str = "".join([
                        (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))[:1] or " "
                        for c in qc_slice
                    ])
                    write_nc_char_1d(adj_qc_var, iprof, qc_str)


def enforce_r_to_d_fill_in_final_bd_files(final_bd_files, WMOfloatid=None):
    """
    Final post-write pass directly on the finalized BD*.nc files in final_out_dir.
    Ensures DATA_MODE[iprof] matches PARAMETER_DATA_MODE[iprof, :], wipes all
    _ADJUSTED* variables to FillValue when DATA_MODE or PARAMETER_DATA_MODE is 'R',
    and enforces '9' inside nlev / ' ' past nlev when PARAMETER_DATA_MODE is 'D' or 'A'.
    """
    for fpath in final_bd_files:
        try:
            ds = nc.Dataset(fpath, "r+")
            ds.set_auto_mask(False)

            standardize_missing_data_and_qc(ds)
            update_profile_qc_grades(ds, WMOfloatid=WMOfloatid)

            ds.sync()
            ds.close()
        except Exception:
            traceback.print_exc()
            if "ds" in locals() and ds.isopen():
                ds.close()


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


# ==============================================================================
# SECTION 2.5: Missing Data & QC Alignment Auditor
# ==============================================================================

def check_bd_missing_data_mismatches(processed_bd_files, WMOfloatid):
    """
    Audits final BD NetCDF files by PRES level per float and cycle.
    Detects and reports mismatches between:
    - CHLA and CHLA_ADJUSTED, CHLA_QC and CHLA_ADJUSTED_QC
    - CHLA_FLUORESCENCE and CHLA_FLUORESCENCE_ADJUSTED, CHLA_FLUORESCENCE_QC and CHLA_FLUORESCENCE_ADJUSTED_QC
    - BBP700 and BBP700_ADJUSTED, BBP700_QC and BBP700_ADJUSTED_QC
    - DOXY and DOXY_ADJUSTED, DOXY_QC and DOXY_ADJUSTED_QC
    - QC = 4 or 9 where _ADJUSTED is not FillValue
    - Trailing padding levels (k > max_k) where QC must be blank (' ') rather than '9'
    - R-mode profiles where _ADJUSTED and _ADJUSTED_QC must be 100% FillValue
    """
    pairs_to_check = [
        ("CHLA", "CHLA_ADJUSTED", "CHLA_QC", "CHLA_ADJUSTED_QC"),
        ("CHLA_FLUORESCENCE", "CHLA_FLUORESCENCE_ADJUSTED", "CHLA_FLUORESCENCE_QC", "CHLA_FLUORESCENCE_ADJUSTED_QC"),
        ("BBP700", "BBP700_ADJUSTED", "BBP700_QC", "BBP700_ADJUSTED_QC"),
        ("DOXY", "DOXY_ADJUSTED", "DOXY_QC", "DOXY_ADJUSTED_QC"),
    ]

    mismatches = []

    def is_qc_missing(qc_char):
        if qc_char is None or pd.isna(qc_char):
            return True
        c = qc_char.decode("utf-8", errors="ignore").strip() if isinstance(qc_char, (bytes, np.bytes_)) else str(qc_char).strip()
        return c in ["", " ", "9"]

    def decode_qc(qc_val):
        if isinstance(qc_val, (bytes, np.bytes_)):
            return qc_val.decode("utf-8", errors="ignore").strip()
        return str(qc_val).strip()

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            ds.set_auto_mask(False)
            cycle_num = int(fpath[-6:-3])

            for raw_var, adj_var, raw_qc_var, adj_qc_var in pairs_to_check:
                # Handle FLUORESCENCE_CHLA naming alias if present
                if raw_var == "CHLA_FLUORESCENCE" and raw_var not in ds.variables and "FLUORESCENCE_CHLA" in ds.variables:
                    raw_var = "FLUORESCENCE_CHLA"
                    adj_var = "FLUORESCENCE_CHLA_ADJUSTED"
                    raw_qc_var = "FLUORESCENCE_CHLA_QC"
                    adj_qc_var = "FLUORESCENCE_CHLA_ADJUSTED_QC"

                iprof = detect_parameter_profile(ds, raw_var)
                mode = get_param_data_mode(ds, iprof, raw_var)

                def get_slice(v_name):
                    if v_name not in ds.variables:
                        return None
                    v = ds.variables[v_name]
                    arr = v[iprof, :] if v.ndim > 1 else v[:]
                    if hasattr(arr, "filled"):
                        arr = arr.filled(MISSING_VAL if arr.dtype.kind in ["f", "i"] else b" ")
                    return np.array(arr)

                raw_data = get_slice(raw_var)
                adj_data = get_slice(adj_var)
                raw_qc_data = get_slice(raw_qc_var)
                adj_qc_data = get_slice(adj_qc_var)

                if raw_data is None and adj_data is None:
                    continue

                pres_vals, _ = get_pres_and_qc_for_profile(ds, iprof)
                n_levels = len(pres_vals)

                # In R-mode, _ADJUSTED must be FillValue and _ADJUSTED_QC must be blank (' ') everywhere
                if mode not in ["D", "A"]:
                    for k in range(n_levels):
                        a_val = adj_data[k] if adj_data is not None and k < len(adj_data) else MISSING_VAL
                        a_qc_str = decode_qc(adj_qc_data[k]) if adj_qc_data is not None and k < len(adj_qc_data) else ""
                        if not is_missing_value(a_val) or a_qc_str != "":
                            mismatches.append({
                                "WMO": WMOfloatid,
                                "CYCLE": cycle_num,
                                "PRES": float(pres_vals[k]) if not is_missing_value(pres_vals[k]) else np.nan,
                                "LEVEL_IDX": k,
                                "COMPARISON_TYPE": "R_MODE_FILL",
                                "PAIR": f"{adj_var} / {adj_qc_var}",
                                "RAW_VAL": a_val,
                                "ADJ_VAL": a_qc_str,
                                "RAW_IS_MISSING": is_missing_value(a_val),
                                "ADJ_IS_MISSING": a_qc_str == "",
                            })
                    continue

                max_k = get_profile_max_level_idx(ds, iprof, raw_var)

                for k in range(n_levels):
                    r_qc = raw_qc_data[k] if raw_qc_data is not None and k < len(raw_qc_data) else b" "
                    a_qc = adj_qc_data[k] if adj_qc_data is not None and k < len(adj_qc_data) else b" "
                    a_qc_str = decode_qc(a_qc)

                    # Check trailing padding levels beyond nlev (k > max_k): must be blank (' ')
                    if max_k < 0 or k > max_k:
                        if a_qc_str != "":
                            mismatches.append({
                                "WMO": WMOfloatid,
                                "CYCLE": cycle_num,
                                "PRES": np.nan,
                                "LEVEL_IDX": k,
                                "COMPARISON_TYPE": "PADDING_QC",
                                "PAIR": f"{raw_qc_var} vs {adj_qc_var}",
                                "RAW_VAL": decode_qc(r_qc),
                                "ADJ_VAL": a_qc_str,
                                "RAW_IS_MISSING": True,
                                "ADJ_IS_MISSING": False,
                            })
                        continue

                    # Inside profile depth range (k <= max_k): missing levels must have QC == '9'
                    pres = pres_vals[k]
                    if is_missing_value(pres) or np.isnan(pres):
                        if a_qc_str != "9":
                            mismatches.append({
                                "WMO": WMOfloatid,
                                "CYCLE": cycle_num,
                                "PRES": np.nan,
                                "LEVEL_IDX": k,
                                "COMPARISON_TYPE": "INTERIOR_QC",
                                "PAIR": f"{raw_qc_var} vs {adj_qc_var}",
                                "RAW_VAL": decode_qc(r_qc),
                                "ADJ_VAL": a_qc_str,
                                "RAW_IS_MISSING": True,
                                "ADJ_IS_MISSING": False,
                            })
                        continue

                    r_val = raw_data[k] if raw_data is not None and k < len(raw_data) else MISSING_VAL
                    a_val = adj_data[k] if adj_data is not None and k < len(adj_data) else MISSING_VAL

                    r_val_miss = is_missing_value(r_val)
                    a_val_miss = is_missing_value(a_val)
                    r_qc_miss = is_qc_missing(r_qc)
                    a_qc_miss = is_qc_missing(a_qc)

                    # Check Argo rule: QC == '4' or '9' must have missing/FillValue in _ADJUSTED
                    if a_qc_str in ["4", "9"] and not a_val_miss:
                        mismatches.append({
                            "WMO": WMOfloatid,
                            "CYCLE": cycle_num,
                            "PRES": float(pres),
                            "LEVEL_IDX": k,
                            "COMPARISON_TYPE": "QC_4_OR_9_VAL",
                            "PAIR": f"{adj_var} vs {adj_qc_var}",
                            "RAW_VAL": a_qc_str,
                            "ADJ_VAL": a_val,
                            "RAW_IS_MISSING": True,
                            "ADJ_IS_MISSING": a_val_miss,
                        })

                    # Check value missingness mismatch (excluding valid QC '4' bad-flagged levels)
                    if a_qc_str != "4" and r_val_miss != a_val_miss:
                        mismatches.append({
                            "WMO": WMOfloatid,
                            "CYCLE": cycle_num,
                            "PRES": float(pres),
                            "LEVEL_IDX": k,
                            "COMPARISON_TYPE": "VALUE",
                            "PAIR": f"{raw_var} vs {adj_var}",
                            "RAW_VAL": r_val,
                            "ADJ_VAL": a_val,
                            "RAW_IS_MISSING": r_val_miss,
                            "ADJ_IS_MISSING": a_val_miss,
                        })

                    # Check QC flag missingness mismatch
                    if r_qc_miss != a_qc_miss:
                        mismatches.append({
                            "WMO": WMOfloatid,
                            "CYCLE": cycle_num,
                            "PRES": float(pres),
                            "LEVEL_IDX": k,
                            "COMPARISON_TYPE": "QC_FLAG",
                            "PAIR": f"{raw_qc_var} vs {adj_qc_var}",
                            "RAW_VAL": decode_qc(r_qc),
                            "ADJ_VAL": a_qc_str,
                            "RAW_IS_MISSING": r_qc_miss,
                            "ADJ_IS_MISSING": a_qc_miss,
                        })

            ds.close()
        except Exception:
            traceback.print_exc()

    df_mismatches = pd.DataFrame(mismatches)

    print("\n" + "=" * 90)
    print(f"MISSING DATA / QC AUDIT SUMMARY: FLOAT {WMOfloatid}")
    print("=" * 90)

    if df_mismatches.empty:
        print(f"[{WMOfloatid}] All missing value masks and QC flags align consistently across all cycles.")
    else:
        print(f"[{WMOfloatid}] Found {len(df_mismatches)} level-by-level mismatches.\n")

        for _, row in df_mismatches.iterrows():
            pres_str = f"{row['PRES']:7.1f}" if pd.notna(row["PRES"]) else "  FILL "
            print(
                f"  Float {row['WMO']} | Cycle {row['CYCLE']:03d} | PRES: {pres_str} dbar (lvl {row['LEVEL_IDX']:03d}) | "
                f"Type: {row['COMPARISON_TYPE']:<13} | Pair: {row['PAIR']:<28} | "
                f"Raw: {str(row['RAW_VAL']):<10} (Miss: {row['RAW_IS_MISSING']}) vs "
                f"Adj: {str(row['ADJ_VAL']):<10} (Miss: {row['ADJ_IS_MISSING']})"
            )

        print("\n--- Summary Count Table by Cycle and Parameter Pair ---")
        summary_table = (
            df_mismatches.groupby(["WMO", "CYCLE", "COMPARISON_TYPE", "PAIR"])
            .size()
            .reset_index(name="MISMATCH_COUNT")
        )
        print(summary_table.to_string(index=False))

    print("=" * 90 + "\n")
    return df_mismatches


# ==============================================================================
# SECTION 3: Plotly DMQC Status Plot Generation
# ==============================================================================

data_type_color_map = {
    "D": "#647ed5",
    "A": "#c07ab8",
    "R": "#ff4d6d",
    "99999": "#808080",
    " ": "#808080"
}
label_data_type_color_map = {
    "Delayed Mode (D)": "#647ed5",
    "Adjusted Mode (A)": "#c07ab8",
    "Real Time Mode (R)": "#ff4d6d",
    "Missing Mode (99999)": "#808080",
}
data_mode_rank_map = {"99999": 0, "R": 1, "A": 2, "D": 3}

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
    "D": "#fdae61", "E": "#f46d43", "F": "#d73027",
    "9": "#808080", " ": "#e0e0e0"
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

        valid_modes = [m for m in modes if m in ["R", "A", "D"]]
        best_mode = max(valid_modes, key=lambda x: data_mode_rank_map.get(x, 0)) if valid_modes else "99999"
        color_dataType.append(data_type_color_map.get(best_mode, "#808080"))

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
            y=y1_position_qc,
            mode="markers",
            marker_symbol="line-ns",
            marker_color=color1_positioning_system,
            marker_line_color=color1_position_qc,
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

        profile_flag_QC_color = param_qc_color_map.get(profile_flag_QC, "#808080")
        color3_profile_flag_QC.append(profile_flag_QC_color)

        possible_pdm_keys = [PARAM]
        if PARAM in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
            possible_pdm_keys.extend(["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"])
        elif PARAM in ["PH", "PH_IN_SITU_TOTAL"]:
            possible_pdm_keys.extend(["PH_IN_SITU_TOTAL", "PH"])

        param_data_mode = "-99999"
        for key in possible_pdm_keys:
            if key in float_map_inst["PDM_MAP"] and n_prof in float_map_inst["PDM_MAP"][key]:
                m = float_map_inst["PDM_MAP"][key][n_prof]
                if m and str(m).strip() in VALID_DATA_MODES:
                    param_data_mode = str(m).strip()
                    break

        param_data_mode_color = data_type_color_map.get(param_data_mode, "#808080")
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
    if not processed_bd_files:
        return

    os.makedirs(output_dir, exist_ok=True)
    out_png = os.path.join(output_dir, f"{WMOfloatid}_DMQC_status.png")

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

                has_data = has_valid_measurements(ds, p, p_prof)
                grade = " "

                qc_var_names = [f"PROFILE_{p}_QC"]
                if p == "CHLA_FLUORESCENCE":
                    qc_var_names.append("PROFILE_FLUORESCENCE_CHLA_QC")
                elif p == "PH_IN_SITU_TOTAL":
                    qc_var_names.append("PROFILE_PH_QC")

                for qvn in qc_var_names:
                    if qvn in ds.variables:
                        raw_val = ds.variables[qvn][:]
                        if hasattr(raw_val, "filled"):
                            raw_val = raw_val.filled(b" ")

                        target_indices = [p_prof, 0] if p_prof < len(raw_val) else [0]
                        for idx_cand in target_indices:
                            if idx_cand < len(raw_val):
                                item = raw_val[idx_cand]
                                while isinstance(item, (np.ndarray, list)) and len(item) > 0:
                                    item = item[0]
                                item_str = item.decode("utf-8", errors="ignore").strip() if isinstance(item, (bytes, np.bytes_)) else str(item).strip()
                                if item_str in ["A", "B", "C", "D", "E", "F"]:
                                    grade = item_str
                                    break
                    if grade != " ":
                        break

                if grade == " ":
                    adj_qc_name = f"{p}_ADJUSTED_QC"
                    raw_qc_name = "PRES_QC" if p == "PRES" else f"{p}_QC"
                    target_qc = raw_qc_name

                    if mode_val in ["D", "A"] and adj_qc_name in ds.variables:
                        adj_q_slice = ds.variables[adj_qc_name][p_prof, :] if ds.variables[adj_qc_name].ndim > 1 else ds.variables[adj_qc_name][:]
                        flags = [q.decode("utf-8", errors="ignore").strip() if isinstance(q, (bytes, np.bytes_)) else str(q).strip() for q in np.ravel(adj_q_slice)]
                        if any(q in ["1", "2", "3", "4", "5", "8"] for q in flags):
                            target_qc = adj_qc_name

                    if target_qc in ds.variables:
                        qc_slice = ds.variables[target_qc][p_prof, :] if ds.variables[target_qc].ndim > 1 else ds.variables[target_qc][:]
                        grade = get_profile_qc_grade(qc_slice, parameter_present=has_data)

                qc_col_key = f"PROFILE_{p}_QC"
                if qc_col_key not in nprof_df_data:
                    nprof_df_data[qc_col_key] = {}
                nprof_df_data[qc_col_key][prof_idx] = grade

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
    "0": "#b3b3b3",
    "1": "#33A02C",
    "2": "#B2DF8A",
    "3": "#FDBF6F",
    "4": "#E31A1C",
    "5": "#1F78B4",
    "8": "#CAB2D6",
    "9": "#000000",
}


def create_interpolated_grid(df, val_col, grid_nx=100, grid_ny=100):
    clean_df = df.dropna(subset=[val_col, "CYCLE_NUMBER", "PRES"]).copy()
    clean_df = clean_df[~clean_df[val_col].isin([MISSING_VAL, 99999.0, -999.0])]

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
        | (df_plot[var_name] == -999.0)
        | (df_plot["QC_VAL"] == "9")
    )

    df_valid = df_plot[~is_missing_mask].copy()
    df_missing = df_plot[is_missing_mask].copy()

    df_valid = df_valid.dropna(subset=[var_name, "PRES", "CYCLE_NUMBER"]).copy()
    df_valid = df_valid[(df_valid["PRES"] != MISSING_VAL) & (df_valid["PRES"] != 99999.0) & (df_valid["PRES"] != -999.0)]

    df_missing = df_missing.dropna(subset=["PRES", "CYCLE_NUMBER"]).copy()
    df_missing = df_missing[(df_missing["PRES"] != MISSING_VAL) & (df_missing["PRES"] != 99999.0) & (df_missing["PRES"] != -999.0)]

    if df_valid.empty and df_missing.empty:
        return

    interp_res = create_interpolated_grid(df_valid, var_name) if not df_valid.empty else None

    pres_all = df_plot["PRES"].replace([MISSING_VAL, 99999.0, -999.0], np.nan)
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

    unique_qcs = sorted(df_plot["QC_VAL"].unique(), key=lambda q: int(q) if q.isdigit() else 99)
    for q_flag in unique_qcs:
        sub = df_plot[df_plot["QC_VAL"] == q_flag]
        sub = sub[(sub["PRES"] != MISSING_VAL) & (sub["PRES"] != 99999.0) & (sub["PRES"] != -999.0)]
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
            ds.set_auto_mask(False)
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

    os.makedirs(output_root, exist_ok=True)

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

            generate_section_qc_plot(
                df_param,
                param_base,
                f"{param_base}_QC",
                f"{param_base} R_mode",
                f"{param_base}_R_mode",
                output_root,
                WMOfloatid,
                col_palette=palette,
                unit_str=unit,
                is_chla=is_chla,
            )

            generate_section_qc_plot(
                df_param,
                f"{param_base}_ADJUSTED",
                f"{param_base}_ADJUSTED_QC",
                f"{param_base} D_mode",
                f"{param_base}_D_mode",
                output_root,
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
print(f"Master plot and CSV output directory created at: {OUTPUT_PLOTS_CSV_DIR}")

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
    main_float_dir = f"/data/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
    bio_dmqc_csv_path = (
        f"/data/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_BBP_DOXY_{WMOfloatid}.csv"
    )

    output_lut_dir = os.path.join(main_float_dir, "LUT")
    today_str = dt.now().strftime("%Y-%m-%d")
    final_out_dir = os.path.join(output_lut_dir, today_str)

    if not os.path.exists(bio_dmqc_csv_path):
        print(f"[{WMOfloatid}] CSV calibration file not found: {bio_dmqc_csv_path}. Skipping.")
        continue

    os.makedirs(main_float_dir, exist_ok=True)
    os.makedirs(output_lut_dir, exist_ok=True)
    os.makedirs(final_out_dir, exist_ok=True)

    df_bio = pd.read_csv(bio_dmqc_csv_path, low_memory=False)

    # --------------------------------------------------------------------------
    # Search for input B-files directly in the central storage location
    # --------------------------------------------------------------------------
    all_bd_files = sorted(
        glob.glob(os.path.join(INPUT_STORAGE_DIR, f"BD*{WMOfloatid}_*.nc"))
    )
    all_br_files = sorted(
        glob.glob(os.path.join(INPUT_STORAGE_DIR, f"AOML_BR*{WMOfloatid}_*.nc"))
    )
    if not all_br_files:
        all_br_files = sorted(
            glob.glob(os.path.join(INPUT_STORAGE_DIR, f"BR*{WMOfloatid}_*.nc"))
        )

    if not all_bd_files and not all_br_files:
        print(f"[{WMOfloatid}] No input BD or BR files found in central storage {INPUT_STORAGE_DIR}. Skipping.")
        continue

    sorted_b_files = organize_b_files(all_bd_files, all_br_files)
    new_bd_files = []

    for bgc_filename in sorted_b_files:
        idx_profile = int(bgc_filename[-6:-3])

        # Stage working BD file in the float's directory for processing
        w_bgc_filename = create_working_bd_file(
            bgc_filename, dest_dir=main_float_dir
        )

        try:
            ds = nc.Dataset(w_bgc_filename, "r+")
            ds.set_auto_mask(False)

            iprof_chla = get_profile_for_parameter(ds, "CHLA")
            iprof_bbp = get_profile_for_parameter(ds, "BBP700")
            iprof_doxy = get_profile_for_parameter(ds, "DOXY")

            # 1. Update PARAMETER_DATA_MODE and sync DATA_MODE
            cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
            write_parameter_data_modes(ds, WMOfloatid=WMOfloatid, cycle_df=cycle_df)

            # 2. Process calibrations and CSV adjustments
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio)
            write_DOXY_slope_drift(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)
            write_DOXY_from_csv(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)

            # 3. Propagate R-mode fill values ('9' inside nlev, ' ' past nlev; FillValue if R-mode)
            standardize_missing_data_and_qc(ds)

            # 4. Update PROFILE_<PARAM>_QC based on evaluated & standardized levels
            update_profile_qc_grades(ds, WMOfloatid=WMOfloatid)

            # 5. Metadata and attributes
            write_history_metadata(ds, iprof_chla if iprof_chla is not None else 0,
                                   iprof_doxy if iprof_doxy is not None else 0,
                                   WMOfloatid=WMOfloatid)
            add_missing_valid_range_attributes(ds)
            remove_forbidden_attributes(ds)

            active_iprofs = [p for p in [iprof_chla, iprof_bbp] if p is not None]
            if WMOfloatid not in DOXY_EXCLUDE_FLOAT_IDS and iprof_doxy is not None:
                active_iprofs.append(iprof_doxy)

            for i in range(ds.dimensions["N_PROF"].size):
                if i in active_iprofs:
                    ds.variables["DATA_STATE_INDICATOR"][i, :] = nc.stringtochar(
                        np.array("2C  ", dtype="S4")
                    )

            ds.sync()
            ds.close()
            new_bd_files.append(w_bgc_filename)

        except Exception:
            traceback.print_exc()
            if "ds" in locals() and ds.isopen():
                ds.close()

    # Step 5: Finalize output NetCDF files into the float's LUT/<date>/ directory
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

    # Step 5.2: Enforce R-mode FillValue and D/A-mode ('9' inside nlev, ' ' past nlev) in final BD output
    enforce_r_to_d_fill_in_final_bd_files(final_processed_files, WMOfloatid=WMOfloatid)

    # --------------------------------------------------------------------------
    # Step 5.5: Audit Missing Value / QC Alignment & Save CSV to BDfiller folder
    # --------------------------------------------------------------------------
    df_missing_mismatches = check_bd_missing_data_mismatches(final_processed_files, WMOfloatid)

    if not df_missing_mismatches.empty:
        audit_csv = os.path.join(OUTPUT_PLOTS_CSV_DIR, f"{WMOfloatid}_missing_data_mismatches.csv")
        df_missing_mismatches.to_csv(audit_csv, index=False)
        print(f"[{WMOfloatid}] Saved missing mismatch audit table to: {audit_csv}")

    # --------------------------------------------------------------------------
    # Step 6: Verification Plotting (Directing all plots to BDfiller folder)
    # --------------------------------------------------------------------------
    generate_dmqc_status_plot(WMOfloatid, final_processed_files, OUTPUT_PLOTS_CSV_DIR)
    generate_all_qc_summary_plots(final_processed_files, WMOfloatid, OUTPUT_PLOTS_CSV_DIR)