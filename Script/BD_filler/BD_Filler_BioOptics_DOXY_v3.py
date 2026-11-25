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
INPUT_STORAGE_DIR = "/data/a1/ARGO_DELAY/DMQC_BGC/storage/"

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
scientific_calibration_equation_CHLA_FLU = (
    "FLUORESCENCE_CHLA_ADJUSTED = (FLUORESCENCE_CHLA - MEDIAN(PRELIM_DARK_CHLA)) * SCALE_CHLA"
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

CSV_COL_MAP = {
    "CHLA": ["CHLA_FINAL", "CHLA_ADJUSTED"],
    "CHLA_FLUORESCENCE": [
        "CHLA_FLUORESCENCE_FINAL",
        "FLUORESCENCE_CHLA_FINAL",
        "CHLA_FLUORSCENCE_FINAL",
        "CHLA_FLUORESCENCE_ADJUSTED",
        "FLUORESCENCE_CHLA_ADJUSTED",
    ],
    "FLUORESCENCE_CHLA": [
        "CHLA_FLUORESCENCE_FINAL",
        "FLUORESCENCE_CHLA_FINAL",
        "CHLA_FLUORSCENCE_FINAL",
        "CHLA_FLUORESCENCE_ADJUSTED",
        "FLUORESCENCE_CHLA_ADJUSTED",
    ],
    "BBP700": ["BBP700_FINAL", "BBP700_ADJUSTED"],
    "DOXY": ["DOXY_FINAL", "DOXY_ADJUSTED", "DOXY_SLOPE"],
}


def is_missing_value(val):
    if val is None or pd.isna(val):
        return True
    try:
        f = float(val)
        return np.isnan(f) or np.isinf(f) or f == -99999.0 or f == 99999.0
    except (ValueError, TypeError):
        return False


def _decode_qc_char(qval):
    if qval is None or np.ma.is_masked(qval):
        return " "
    if isinstance(qval, (bytes, np.bytes_)):
        return qval.decode("utf-8", errors="ignore")
    return str(qval)


def _get_raw_val_and_qc(ds, raw_candidates, iprof, level_idx):
    """Returns (raw_val, raw_qc_char) for the first matching raw variable in raw_candidates."""
    for vname in raw_candidates:
        if vname in ds.variables:
            v = ds.variables[vname]
            fill_v = getattr(v, "_FillValue", MISSING_VAL)
            val = v[iprof, level_idx] if v.ndim > 1 else v[level_idx]
            if np.ma.is_masked(val) or is_missing_value(val) or val == fill_v:
                val_f = np.nan
            else:
                val_f = float(val)

            qc_char = "9"
            qc_name = f"{vname}_QC"
            if qc_name in ds.variables:
                qv = ds.variables[qc_name]
                qval = qv[iprof, level_idx] if qv.ndim > 1 else qv[level_idx]
                qstr = _decode_qc_char(qval)
                if qstr == " ":
                    qc_char = " "
                elif qstr.strip() in VALID_QC_CHARS:
                    qc_char = qstr.strip()
            return val_f, qc_char
    return np.nan, "9"


def _is_raw_level_missing(ds, raw_candidates, iprof, level_idx):
    """
    Per Argo GDAC format checker rules: a raw measurement is missing if and only if
    the raw value is FillValue/NaN OR <PARAM>_QC is '9' (or blank padding).
    """
    raw_val, raw_qc = _get_raw_val_and_qc(ds, raw_candidates, iprof, level_idx)
    return np.isnan(raw_val) or raw_qc.strip() in ["", "9"]


def _get_param_active_max_level(ds, raw_candidates, iprof):
    """
    Finds the highest level index k where <PARAM>_QC is non-blank (i.e. not end-of-profile padding).
    Returns -1 if all levels are blank padding.
    """
    n_levels = ds.dimensions["N_LEVELS"].size
    for vname in raw_candidates:
        qc_name = f"{vname}_QC"
        if qc_name in ds.variables:
            qv = ds.variables[qc_name]
            qarr = qv[iprof, :] if qv.ndim > 1 else qv[:]
            if hasattr(qarr, "filled"):
                qarr = qarr.filled(b" ")
            non_blank = [idx for idx in range(min(n_levels, len(qarr))) if _decode_qc_char(qarr[idx]) != " "]
            if non_blank:
                return non_blank[-1]
        if vname in ds.variables:
            v = ds.variables[vname]
            arr = v[iprof, :] if v.ndim > 1 else v[:]
            vals = arr.filled(MISSING_VAL) if hasattr(arr, "filled") else np.array(arr)
            valid_idx = np.where(~np.isnan(vals) & (vals != MISSING_VAL) & (vals != 99999.0) & (vals != -99999.0))[0]
            if len(valid_idx) > 0:
                return int(valid_idx[-1])
    return -1


# ==============================================================================
# SECTION 2: Architecture & NetCDF Helper Functions
# ==============================================================================

def safe_rename(from_file, to_file):
    gc.collect()
    try:
        os.replace(from_file, to_file)
    except OSError:
        shutil.copyfile(from_file, to_file)
        os.remove(from_file)


def compact_station_parameters_and_nparam(nc_path):
    """
    Fixes two Argo GDAC format checker warnings on BD files:
    1. 'STATION_PARAMETERS[i,*]: Empty entries in list' -> left-aligns non-empty
       STATION_PARAMETERS[iprof, :] and re-indexes PARAMETER_DATA_MODE and
       SCIENTIFIC_CALIB_* along N_PARAM to match.
    2. 'N_PARAM: Larger than necessary' -> shrinks the N_PARAM dimension to the
       exact maximum number of used parameters across all N_PROF profiles.
    """
    ds = nc.Dataset(nc_path, "r")
    ds.set_auto_mask(False)
    try:
        if "N_PARAM" not in ds.dimensions or "STATION_PARAMETERS" not in ds.variables:
            return

        n_prof = ds.dimensions["N_PROF"].size
        old_n_param = ds.dimensions["N_PARAM"].size

        sp_raw = ds.variables["STATION_PARAMETERS"][:]
        per_prof_indices = []
        has_internal_gaps = False

        for iprof in range(n_prof):
            used_j = []
            seen_blank = False
            for j in range(old_n_param):
                p_str = "".join([
                    (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                    for c in sp_raw[iprof, j]
                ]).strip()
                if p_str:
                    if seen_blank:
                        has_internal_gaps = True
                    used_j.append((j, p_str))
                else:
                    seen_blank = True
            per_prof_indices.append(used_j)

        new_n_param = max((len(u) for u in per_prof_indices), default=1)
        new_n_param = max(1, new_n_param)

        # Nothing to do if N_PARAM is already minimal and no profile has internal blank gaps
        if new_n_param == old_n_param and not has_internal_gaps:
            return

        tmp_compact = f"{nc_path}.compact.tmp"
        out = nc.Dataset(tmp_compact, "w", format=ds.data_model)
        out.set_auto_mask(False)

        try:
            # Copy global attributes
            for attr in ds.ncattrs():
                out.setncattr(attr, ds.getncattr(attr))

            # Copy dimensions, resizing N_PARAM -> new_n_param
            for dname, dim in ds.dimensions.items():
                if dname == "N_PARAM":
                    out.createDimension(dname, new_n_param)
                else:
                    out.createDimension(dname, None if dim.isunlimited() else dim.size)

            # Copy all variables, remapping N_PARAM axis where present
            for vname, var in ds.variables.items():
                fill_val = getattr(var, "_FillValue", None)
                kwargs = {}
                if fill_val is not None:
                    kwargs["fill_value"] = fill_val

                out_var = out.createVariable(vname, var.datatype, var.dimensions, **kwargs)
                for attr in var.ncattrs():
                    if attr != "_FillValue":
                        out_var.setncattr(attr, var.getncattr(attr))

                data = np.array(var[:], copy=True)

                if "N_PARAM" not in var.dimensions:
                    out_var[:] = data
                    continue

                param_axis = var.dimensions.index("N_PARAM")
                prof_axis = var.dimensions.index("N_PROF") if "N_PROF" in var.dimensions else None

                # Build blank/fill-initialized output array with resized N_PARAM axis
                new_shape = list(data.shape)
                new_shape[param_axis] = new_n_param
                if data.dtype.kind == "S":
                    new_data = np.full(new_shape, b" ", dtype=data.dtype)
                else:
                    fval = fill_val if fill_val is not None else MISSING_VAL
                    new_data = np.full(new_shape, fval, dtype=data.dtype)

                if prof_axis is not None:
                    for iprof in range(n_prof):
                        for new_j, (old_j, _) in enumerate(per_prof_indices[iprof]):
                            src_idx = [slice(None)] * data.ndim
                            dst_idx = [slice(None)] * data.ndim
                            src_idx[prof_axis] = iprof
                            src_idx[param_axis] = old_j
                            dst_idx[prof_axis] = iprof
                            dst_idx[param_axis] = new_j
                            new_data[tuple(dst_idx)] = data[tuple(src_idx)]
                else:
                    copy_len = min(old_n_param, new_n_param)
                    src_idx = [slice(None)] * data.ndim
                    dst_idx = [slice(None)] * data.ndim
                    src_idx[param_axis] = slice(0, copy_len)
                    dst_idx[param_axis] = slice(0, copy_len)
                    new_data[tuple(dst_idx)] = data[tuple(src_idx)]

                out_var[:] = new_data

            out.sync()
        finally:
            out.close()
    finally:
        ds.close()

    if os.path.exists(tmp_compact):
        safe_rename(tmp_compact, nc_path)


def is_broken_sensor_float(WMOfloatid, cycle_num):
    float_str = str(WMOfloatid)
    cycle_num = int(cycle_num)
    return (
        ((float_str == "2904010") & (cycle_num >= 49))
        | ((float_str == "2904011") & (cycle_num >= 24))
    )


def get_station_parameters_for_profile(ds, iprof):
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


def get_profiles_for_parameter(ds, param_base):
    """
    Returns ALL profile indices (N_PROF) where param_base is present in STATION_PARAMETERS.
    Crucial for floats with multiple N_PROF axes (e.g., N_PROF=0 primary, N_PROF=1 unpumped/secondary).
    """
    if "STATION_PARAMETERS" not in ds.variables:
        return []

    target_names = {param_base}
    if param_base in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
        target_names = {"CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"}
    elif param_base in ["PH", "PH_IN_SITU_TOTAL"]:
        target_names = {"PH_IN_SITU_TOTAL", "PH"}
    elif param_base == "NITRATE":
        target_names = {"NITRATE"}

    matching_profs = []
    n_prof = ds.dimensions["N_PROF"].size
    for iprof in range(n_prof):
        prof_params = get_station_parameters_for_profile(ds, iprof)
        for p in prof_params:
            if p in target_names or any(p.startswith(t) for t in target_names):
                matching_profs.append(iprof)
                break
    return matching_profs


def get_profile_for_parameter(ds, param_base):
    profs = get_profiles_for_parameter(ds, param_base)
    return profs[0] if profs else None


def detect_parameter_profile(ds, param_prefix):
    prof = get_profile_for_parameter(ds, param_prefix)
    if prof is not None:
        return prof

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
            valid_mask = (~pd.isna(p_vals)) & (p_vals != MISSING_VAL) & (p_vals != 99999.0) & (p_vals != -99999.0)
            if valid_mask.any():
                return True
    return False


def has_csv_calibration(cycle_df, p_name):
    if cycle_df is None or cycle_df.empty:
        return False
    cols = CSV_COL_MAP.get(p_name, [])
    for c in cols:
        if c in cycle_df.columns:
            vals = pd.to_numeric(cycle_df[c], errors="coerce").dropna()
            if not vals.empty and (~vals.isin([-99999.0, 99999.0])).any():
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
            mode_str = _decode_qc_char(m).strip()
            return mode_str if mode_str in VALID_DATA_MODES else "99999"
        return "99999"

    station_params = bgc_file.variables["STATION_PARAMETERS"][iprof]
    pdm_array = bgc_file.variables["PARAMETER_DATA_MODE"][iprof]

    if isinstance(station_params, np.ma.MaskedArray):
        station_params = station_params.filled(b" ")
    if isinstance(pdm_array, np.ma.MaskedArray):
        pdm_array = pdm_array.filled(b" ")

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
            mode_str = _decode_qc_char(pdm_array[j]).strip()
            return mode_str if mode_str in VALID_DATA_MODES else "99999"

    dm_var = bgc_file.variables.get("DATA_MODE")
    if dm_var is not None:
        m = dm_var[iprof] if dm_var.ndim > 0 else dm_var[:]
        mode_str = _decode_qc_char(m).strip()
        return mode_str if mode_str in VALID_DATA_MODES else "99999"
    return "99999"


def write_parameter_data_modes(bgc_file, WMOfloatid=None, cycle_num=None, cycle_df=None):
    """
    Sets PARAMETER_DATA_MODE[iprof, j] = 'D' for DMQC-managed BGC parameters
    (CHLA, CHLA_FLUORESCENCE/FLUORESCENCE_CHLA, BBP700, and DOXY unless excluded)
    and derives DATA_MODE[iprof] strictly from PARAMETER_DATA_MODE[iprof, :].
    """
    bgc_file.set_auto_mask(False)
    try:
        n_param = bgc_file.dimensions["N_PARAM"].size
        n_prof = bgc_file.dimensions["N_PROF"].size
        skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

        pdm_var = bgc_file.variables["PARAMETER_DATA_MODE"]
        dm_var = bgc_file.variables["DATA_MODE"]

        pdm = np.array(pdm_var[:], dtype="|S1", copy=True)
        data_mode = np.array(dm_var[:], dtype="|S1", copy=True)

        mode_rank = {"R": 1, "A": 2, "D": 3}

        for iprof in range(n_prof):
            param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
            prof_modes = []
            for j in range(n_param):
                p_bytes = param_mat[j]
                p_str = "".join([
                    (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                    for c in p_bytes
                ]).strip()

                if not p_str:
                    continue

                existing_m = _decode_qc_char(pdm[iprof, j]).strip()
                if existing_m not in VALID_DATA_MODES:
                    existing_m = "R"

                if "PH" in p_str or "NITRATE" in p_str or (skip_doxy and p_str == "DOXY"):
                    prof_modes.append(existing_m)
                    continue

                if p_str == "PRES" or p_str.startswith("DOWN_IRRADIANCE") or p_str.startswith("UP_RADIANCE"):
                    pdm[iprof, j] = b"R"
                    prof_modes.append("R")
                elif p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY"]:
                    pdm[iprof, j] = b"D"
                    prof_modes.append("D")
                else:
                    prof_modes.append(existing_m)

            if prof_modes:
                best_mode = max(prof_modes, key=lambda m: mode_rank.get(m, 1))
                data_mode[iprof] = best_mode.encode("utf-8")

        pdm_var[:] = pdm
        dm_var[:] = data_mode
        bgc_file.sync()
    finally:
        bgc_file.set_auto_mask(True)


def enforce_empty_adjusted_modes_to_realtime(bgc_file):
    """
    Ensures PARAMETER_DATA_MODE and DATA_MODE align with _ADJUSTED / _ADJUSTED_QC:
    - Keeps 'D'/'A' mode whenever _ADJUSTED has valid values OR _ADJUSTED_QC has valid flags
      ('1', '2', '3', '4', '5', '8' — including broken/uncalibrated cycles flagged '4').
    - Reverts a parameter to 'R' only when all levels are missing ('9' or ' ').
    - Recomputes DATA_MODE[iprof] strictly from PARAMETER_DATA_MODE[iprof, :].
    """
    if "STATION_PARAMETERS" not in bgc_file.variables or "PARAMETER_DATA_MODE" not in bgc_file.variables:
        return

    bgc_file.set_auto_mask(False)
    try:
        n_prof = bgc_file.dimensions["N_PROF"].size
        n_param = bgc_file.dimensions["N_PARAM"].size

        pdm_var = bgc_file.variables["PARAMETER_DATA_MODE"]
        dm_var = bgc_file.variables["DATA_MODE"]

        pdm = np.array(pdm_var[:], dtype="|S1", copy=True)
        data_mode = np.array(dm_var[:], dtype="|S1", copy=True)

        mode_rank = {"R": 1, "A": 2, "D": 3}

        for iprof in range(n_prof):
            param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
            prof_modes = []
            for j in range(n_param):
                p_bytes = param_mat[j]
                p_str = "".join([
                    (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                    for c in p_bytes
                ]).strip()

                if not p_str:
                    continue

                cur_mode = _decode_qc_char(pdm[iprof, j]).strip()
                if cur_mode in ["D", "A"] and p_str in ["CHLA", "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA", "BBP700", "DOXY"]:
                    adj_candidates = [f"{p_str}_ADJUSTED"]
                    if p_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
                        adj_candidates = ["CHLA_FLUORESCENCE_ADJUSTED", "FLUORESCENCE_CHLA_ADJUSTED"]

                    has_any_adj = False
                    for ac in adj_candidates:
                        if ac in bgc_file.variables:
                            vals = np.array(bgc_file.variables[ac][iprof, :], dtype="float32")
                            if (~np.isnan(vals) & (vals != MISSING_VAL) & (vals != 99999.0) & (vals != -99999.0)).any():
                                has_any_adj = True
                                break
                        qc_ac = f"{ac}_QC"
                        if qc_ac in bgc_file.variables:
                            qsl = bgc_file.variables[qc_ac][iprof, :]
                            if any(_decode_qc_char(q).strip() in {"1", "2", "3", "4", "5", "8"} for q in qsl):
                                has_any_adj = True
                                break

                    if not has_any_adj:
                        pdm[iprof, j] = b"R"
                        cur_mode = "R"
                        base_list = (
                            ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
                            if p_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]
                            else [p_str]
                        )
                        for base_c in base_list:
                            for suffix, fill_v in [("_ADJUSTED", MISSING_VAL), ("_ADJUSTED_ERROR", MISSING_VAL)]:
                                vname = f"{base_c}{suffix}"
                                if vname in bgc_file.variables:
                                    v = bgc_file.variables[vname]
                                    fval = getattr(v, "_FillValue", fill_v)
                                    v[iprof, :] = np.full((v.shape[1],), np.float32(fval), dtype="float32")
                            qname = f"{base_c}_ADJUSTED_QC"
                            if qname in bgc_file.variables:
                                qv = bgc_file.variables[qname]
                                qv[iprof, :] = np.full((qv.shape[1],), b" ", dtype="|S1")

                if cur_mode in VALID_DATA_MODES:
                    prof_modes.append(cur_mode)

            if prof_modes:
                best_mode = max(prof_modes, key=lambda m: mode_rank.get(m, 1))
                data_mode[iprof] = best_mode.encode("utf-8")

        pdm_var[:] = pdm
        dm_var[:] = data_mode
        bgc_file.sync()
    finally:
        bgc_file.set_auto_mask(True)


def sync_data_mode_with_parameter_modes(bgc_file, WMOfloatid=None):
    """
    Guarantees strict GDAC consistency between DATA_MODE[iprof] and PARAMETER_DATA_MODE[iprof, :]:
    1. If a BGC file has zero 'D' parameters in PARAMETER_DATA_MODE, promotes the primary BGC
       parameter (CHLA, BBP700, or DOXY) in PARAMETER_DATA_MODE[iprof, j] to 'D' and sets
       _ADJUSTED_QC = '4' wherever raw data exist ('9' where raw is missing, ' ' where blank padding).
    2. Recomputes DATA_MODE[iprof] = max(PARAMETER_DATA_MODE[iprof, :]) across every profile so
       DATA_MODE[i] / PARAMETER_DATA_MODE[i, :] can NEVER be 'D'/'R'.
    """
    if "STATION_PARAMETERS" not in bgc_file.variables or "PARAMETER_DATA_MODE" not in bgc_file.variables:
        return

    bgc_file.set_auto_mask(False)
    try:
        n_prof = bgc_file.dimensions["N_PROF"].size
        n_param = bgc_file.dimensions["N_PARAM"].size
        n_levels = bgc_file.dimensions["N_LEVELS"].size
        skip_doxy = WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS if WMOfloatid is not None else False

        pdm_var = bgc_file.variables["PARAMETER_DATA_MODE"]
        dm_var = bgc_file.variables["DATA_MODE"]

        pdm = np.array(pdm_var[:], dtype="|S1", copy=True)
        data_mode = np.array(dm_var[:], dtype="|S1", copy=True)

        mode_rank = {"R": 1, "A": 2, "D": 3}

        any_d_in_file = False
        for iprof in range(n_prof):
            param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
            for j in range(n_param):
                p_str = "".join([
                    (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                    for c in param_mat[j]
                ]).strip()
                if p_str and _decode_qc_char(pdm[iprof, j]).strip() == "D":
                    any_d_in_file = True
                    break

        if not any_d_in_file:
            preferred_targets = ["CHLA", "BBP700"] if skip_doxy else ["CHLA", "BBP700", "DOXY"]
            promoted = False
            for target_p in preferred_targets:
                for iprof in range(n_prof):
                    param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
                    for j in range(n_param):
                        p_str = "".join([
                            (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                            for c in param_mat[j]
                        ]).strip()
                        if p_str == target_p:
                            pdm[iprof, j] = b"D"
                            adj_v = f"{target_p}_ADJUSTED"
                            adj_q = f"{target_p}_ADJUSTED_QC"
                            adj_e = f"{target_p}_ADJUSTED_ERROR"
                            if adj_v in bgc_file.variables:
                                fill_v = getattr(bgc_file.variables[adj_v], "_FillValue", MISSING_VAL)
                                bgc_file.variables[adj_v][iprof, :] = np.full((n_levels,), np.float32(fill_v), dtype="float32")
                            if adj_e in bgc_file.variables:
                                fill_e = getattr(bgc_file.variables[adj_e], "_FillValue", MISSING_VAL)
                                bgc_file.variables[adj_e][iprof, :] = np.full((n_levels,), np.float32(fill_e), dtype="float32")
                            if adj_q in bgc_file.variables:
                                q_out = np.full((n_levels,), b" ", dtype="|S1")
                                for k in range(n_levels):
                                    raw_val, raw_qc = _get_raw_val_and_qc(bgc_file, [target_p], iprof, k)
                                    if raw_qc == " ":
                                        q_out[k] = b" "
                                    elif np.isnan(raw_val) or raw_qc.strip() in ["", "9"]:
                                        q_out[k] = b"9"
                                    else:
                                        # Raw data present while _ADJUSTED is missing -> must be '4' (never '9')
                                        q_out[k] = b"4"
                                bgc_file.variables[adj_q][iprof, :] = q_out
                            promoted = True
                            break
                    if promoted:
                        break
                if promoted:
                    break

        for iprof in range(n_prof):
            param_mat = bgc_file.variables["STATION_PARAMETERS"][iprof]
            prof_modes = []
            for j in range(n_param):
                p_str = "".join([
                    (c.decode("utf-8", errors="ignore") if isinstance(c, (bytes, np.bytes_)) else str(c))
                    for c in param_mat[j]
                ]).strip()
                if not p_str:
                    continue
                m = _decode_qc_char(pdm[iprof, j]).strip()
                if m in VALID_DATA_MODES:
                    prof_modes.append(m)

            if prof_modes:
                best_mode = max(prof_modes, key=lambda m: mode_rank.get(m, 1))
                data_mode[iprof] = best_mode.encode("utf-8")

        pdm_var[:] = pdm
        dm_var[:] = data_mode
        bgc_file.sync()
    finally:
        bgc_file.set_auto_mask(True)


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


def get_pres_and_qc_for_profile(ds, iprof):
    """
    Retrieves pressure and pressure QC for profile `iprof`.
    Merges PRES_ADJUSTED with raw PRES level-by-level (falling back to PRES wherever
    PRES_ADJUSTED is 99999.0/NaN) so deep BGC levels below the adjusted CTD grid still
    receive their valid pressure values.
    """
    n_levels = ds.dimensions["N_LEVELS"].size
    pres_vals = np.full((n_levels,), np.nan, dtype="float32")
    pres_qcs = np.full((n_levels,), b" ", dtype="|S1")

    def _read_p(var_name):
        if var_name not in ds.variables:
            return None
        v = ds.variables[var_name]
        arr = v[iprof, :] if v.ndim > 1 else v[:]
        p_raw = arr.filled(np.nan) if hasattr(arr, "filled") else np.array(arr, dtype="float32")
        mask = (p_raw == MISSING_VAL) | (p_raw == 99999.0) | (p_raw == -99999.0) | np.isinf(p_raw)
        p_raw[mask] = np.nan
        return np.float32(p_raw)

    p_raw = _read_p("PRES")
    p_adj = _read_p("PRES_ADJUSTED")

    if p_adj is not None and p_raw is not None:
        pres_vals = np.where(np.isnan(p_adj), p_raw, p_adj)
    elif p_adj is not None:
        pres_vals = p_adj
    elif p_raw is not None:
        pres_vals = p_raw

    q_var_name = (
        "PRES_ADJUSTED_QC"
        if (p_adj is not None and not np.all(np.isnan(p_adj)) and "PRES_ADJUSTED_QC" in ds.variables)
        else ("PRES_QC" if "PRES_QC" in ds.variables else None)
    )

    if q_var_name and q_var_name in ds.variables:
        qv = ds.variables[q_var_name]
        q_slice = qv[iprof, :] if qv.ndim > 1 else qv[:]
        if hasattr(q_slice, "filled"):
            q_slice = q_slice.filled(b" ")
        pres_qcs = np.array(q_slice).astype("|S1")

    return pres_vals, pres_qcs


def match_csv_row_by_pres(cycle_df, nc_pres):
    if cycle_df.empty or is_missing_value(nc_pres):
        return None

    csv_pres = pd.to_numeric(cycle_df["PRES"], errors="coerce").values
    valid_mask = ~np.isnan(csv_pres)
    if not valid_mask.any():
        return None

    diffs = np.where(valid_mask, np.abs(csv_pres - nc_pres), np.inf)
    min_idx = int(np.argmin(diffs))

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


def _write_param_from_csv(
    ds,
    iprof,
    raw_candidates,
    adj_var_name,
    cycle_df,
    val_csv_col,
    qc_csv_col,
    err_csv_col=None,
    default_err=None,
    is_broken=False,
    fallback_to_raw=False,
    slope_val=None,
    drift_val=None,
):
    """
    Core writer that populates <PARAM>_ADJUSTED, <PARAM>_ADJUSTED_QC, and <PARAM>_ADJUSTED_ERROR
    from <PARAM>_FINAL and <PARAM>_FINAL_QC in the input CSV across all active levels of <PARAM>.
    Enforces Argo GDAC & BGC QC Manual rules with auto_mask disabled so single-char |S1 flags
    are never silently masked by netCDF4.
    """
    if iprof is None or adj_var_name not in ds.variables:
        return

    ds.set_auto_mask(False)
    try:
        n_levels = ds.dimensions["N_LEVELS"].size
        adj_qc_name = f"{adj_var_name}_QC"
        adj_err_name = f"{adj_var_name}_ERROR"

        adj_var = ds.variables[adj_var_name]
        fill_val = getattr(adj_var, "_FillValue", MISSING_VAL)

        adj_arr = np.full((n_levels,), np.float32(fill_val), dtype="float32")
        qc_arr = np.full((n_levels,), b" ", dtype="|S1")
        err_arr = np.full((n_levels,), np.float32(fill_val), dtype="float32") if adj_err_name in ds.variables else None

        max_k = _get_param_active_max_level(ds, raw_candidates, iprof)
        if max_k < 0:
            adj_var[iprof, :] = adj_arr
            if adj_qc_name in ds.variables:
                ds.variables[adj_qc_name][iprof, :] = qc_arr
            if err_arr is not None:
                ds.variables[adj_err_name][iprof, :] = err_arr
            return

        pres_vals, _ = get_pres_and_qc_for_profile(ds, iprof)

        for i in range(max_k + 1):
            raw_val, raw_qc = _get_raw_val_and_qc(ds, raw_candidates, iprof, i)

            # Preserve blank padding only where raw <PARAM>_QC is blank
            if raw_qc == " ":
                adj_arr[i] = np.float32(fill_val)
                qc_arr[i] = b" "
                if err_arr is not None:
                    err_arr[i] = np.float32(fill_val)
                continue

            # Only set _ADJUSTED_QC = '9' when raw <PARAM> itself is missing!
            if np.isnan(raw_val) or raw_qc.strip() in ["", "9"]:
                adj_arr[i] = np.float32(MISSING_VAL)
                qc_arr[i] = b"9"
                if err_arr is not None:
                    err_arr[i] = np.float32(MISSING_VAL)
                continue

            nc_pres = pres_vals[i]
            # Raw data IS present, so if PRES is missing or sensor is broken, mark bad ('4') — NEVER '9'
            if is_missing_value(nc_pres) or is_broken:
                adj_arr[i] = np.float32(MISSING_VAL)
                qc_arr[i] = b"4"
                if err_arr is not None:
                    err_arr[i] = np.float32(MISSING_VAL)
                continue

            matched_row = match_csv_row_by_pres(cycle_df, nc_pres) if (cycle_df is not None and not cycle_df.empty) else None

            if matched_row is not None and val_csv_col is not None:
                csv_val = matched_row.get(val_csv_col)
                csv_qc = sanitize_qc_char(matched_row.get(qc_csv_col), default="1") if qc_csv_col is not None else "1"

                if csv_qc in ["4", "9"] or pd.isna(csv_val) or is_missing_value(csv_val):
                    adj_arr[i] = np.float32(MISSING_VAL)
                    qc_arr[i] = b"4"
                    if err_arr is not None:
                        err_arr[i] = np.float32(MISSING_VAL)
                else:
                    adj_arr[i] = np.float32(csv_val)
                    qc_arr[i] = csv_qc.encode("utf-8")
                    if err_arr is not None:
                        if err_csv_col and pd.notna(matched_row.get(err_csv_col)) and not is_missing_value(matched_row.get(err_csv_col)):
                            err_arr[i] = np.float32(matched_row.get(err_csv_col))
                        elif default_err is not None:
                            err_arr[i] = np.float32(default_err)
            elif slope_val is not None:
                if raw_qc.strip() == "4":
                    adj_arr[i] = np.float32(MISSING_VAL)
                    qc_arr[i] = b"4"
                    if err_arr is not None:
                        err_arr[i] = np.float32(MISSING_VAL)
                else:
                    d_val = 0.0 if drift_val is None else drift_val
                    adj_arr[i] = np.float32(raw_val * slope_val + d_val)
                    qc_arr[i] = b"1"
                    if err_arr is not None and default_err is not None:
                        err_arr[i] = np.float32(default_err)
            elif fallback_to_raw:
                if raw_qc.strip() == "4":
                    adj_arr[i] = np.float32(MISSING_VAL)
                    qc_arr[i] = b"4"
                    if err_arr is not None:
                        err_arr[i] = np.float32(MISSING_VAL)
                else:
                    adj_arr[i] = np.float32(raw_val)
                    qc_arr[i] = raw_qc.strip().encode("utf-8")
                    if err_arr is not None and default_err is not None:
                        err_arr[i] = np.float32(default_err)
            else:
                # Raw data exists while _ADJUSTED is missing (99999.0) -> BGC manual mandates QC = '4' (never '9')
                adj_arr[i] = np.float32(MISSING_VAL)
                qc_arr[i] = b"4"
                if err_arr is not None:
                    err_arr[i] = np.float32(MISSING_VAL)

        adj_var[iprof, :] = adj_arr
        if adj_qc_name in ds.variables:
            ds.variables[adj_qc_name][iprof, :] = qc_arr
        if err_arr is not None:
            ds.variables[adj_err_name][iprof, :] = err_arr
        ds.sync()
    finally:
        ds.set_auto_mask(True)


def write_scientific_calib_chla(bgc_file, WMOfloatid, idx_profile, df_bio):
    is_broken_chla = is_broken_sensor_float(WMOfloatid, idx_profile)
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

    if not is_broken_chla:
        if not cycle_df.empty and "CHLA_FINAL_QC" in cycle_df.columns:
            qcs = pd.to_numeric(cycle_df["CHLA_FINAL_QC"], errors="coerce").dropna().astype(int)
            if not qcs.empty and (qcs.isin([4, 9])).all():
                is_broken_chla = True

    if is_broken_chla or cycle_df.empty:
        comment_chla = scientific_calibration_comment_BROKEN
        comment_flu = scientific_calibration_comment_BROKEN
        comment_bbp = scientific_calibration_comment_BROKEN
        equation_chla = scientific_calibration_equation_BROKEN
        equation_flu = scientific_calibration_equation_BROKEN
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
        equation_flu = scientific_calibration_equation_CHLA_FLU
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
    SciCalEquArray_CHLA_FLU = nc.stringtochar(np.array(equation_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalEquArray_BBP = nc.stringtochar(np.array(equation_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA = nc.stringtochar(np.array(coef_chla.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_CHLA_FLU = nc.stringtochar(np.array(coef_flu.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalCoeArray_BBP = nc.stringtochar(np.array(coef_bbp.ljust(str256_len)[:str256_len], dtype=f"S{str256_len}"))
    SciCalDateArray = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))

    n_prof_size = bgc_file.dimensions["N_PROF"].size
    n_param_size = bgc_file.dimensions["N_PARAM"].size

    profs_chla = set(get_profiles_for_parameter(bgc_file, "CHLA"))
    profs_fluo = set(get_profiles_for_parameter(bgc_file, "CHLA_FLUORESCENCE"))
    profs_bbp = set(get_profiles_for_parameter(bgc_file, "BBP700"))

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

            if param_str in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"] and iprof_idx in profs_fluo:
                if get_param_data_mode(bgc_file, iprof_idx, "CHLA_FLUORESCENCE") in ["D", "A"]:
                    bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA_FLU
                    bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_CHLA_FLU
                    bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA_FLU
                    bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "CHLA" and iprof_idx in profs_chla:
                if get_param_data_mode(bgc_file, iprof_idx, "CHLA") in ["D", "A"]:
                    bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_CHLA
                    bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_CHLA
                    bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA
                    bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str == "BBP700" and iprof_idx in profs_bbp:
                if get_param_data_mode(bgc_file, iprof_idx, "BBP700") in ["D", "A"]:
                    bgc_file.variables["SCIENTIFIC_CALIB_COMMENT"][iprof_idx, 0, j, :] = SciCalComArray_BBP
                    bgc_file.variables["SCIENTIFIC_CALIB_EQUATION"][iprof_idx, 0, j, :] = SciCalEquArray_BBP
                    bgc_file.variables["SCIENTIFIC_CALIB_COEFFICIENT"][iprof_idx, 0, j, :] = SciCalCoeArray_BBP
                    bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof_idx, 0, j, :] = SciCalDateArray


def write_chla_BBP_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio):
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)

    # 1. CHLA <- CHLA_FINAL & CHLA_FINAL_QC across ALL N_PROF entries containing CHLA
    chla_csv_col = next((c for c in ["CHLA_FINAL", "CHLA_ADJUSTED"] if c in df_bio.columns), None)
    chla_qc_csv_col = next((c for c in ["CHLA_FINAL_QC", "CHLA_ADJUSTED_QC"] if c in df_bio.columns), None)
    for prof_chla in get_profiles_for_parameter(bgc_file, "CHLA"):
        if get_param_data_mode(bgc_file, prof_chla, "CHLA") in ["D", "A"]:
            _write_param_from_csv(
                ds=bgc_file,
                iprof=prof_chla,
                raw_candidates=["CHLA"],
                adj_var_name="CHLA_ADJUSTED",
                cycle_df=cycle_df,
                val_csv_col=chla_csv_col,
                qc_csv_col=chla_qc_csv_col,
                default_err=CHLA_Adjusted_ERROR_est,
                is_broken=is_broken,
                fallback_to_raw=False,
            )

    # 2. CHLA_FLUORESCENCE / FLUORESCENCE_CHLA across ALL N_PROF entries
    fluo_csv_col = next(
        (
            c
            for c in [
                "CHLA_FLUORESCENCE_FINAL",
                "FLUORESCENCE_CHLA_FINAL",
                "CHLA_FLUORSCENCE_FINAL",
                "CHLA_FLUORESCENCE_ADJUSTED",
                "FLUORESCENCE_CHLA_ADJUSTED",
            ]
            if c in df_bio.columns
        ),
        None,
    )
    fluo_qc_csv_col = next(
        (
            c
            for c in [
                "CHLA_FLUORESCENCE_FINAL_QC",
                "FLUORESCENCE_CHLA_FINAL_QC",
                "CHLA_FLUORSCENCE_FINAL_QC",
                "CHLA_FLUORESCENCE_ADJUSTED_QC",
                "FLUORESCENCE_CHLA_ADJUSTED_QC",
            ]
            if c in df_bio.columns
        ),
        None,
    )
    for prof_fluo in get_profiles_for_parameter(bgc_file, "CHLA_FLUORESCENCE"):
        if get_param_data_mode(bgc_file, prof_fluo, "CHLA_FLUORESCENCE") in ["D", "A"]:
            for var_prefix in ["FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]:
                adj_var = f"{var_prefix}_ADJUSTED"
                if adj_var in bgc_file.variables:
                    _write_param_from_csv(
                        ds=bgc_file,
                        iprof=prof_fluo,
                        raw_candidates=[var_prefix, "CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"],
                        adj_var_name=adj_var,
                        cycle_df=cycle_df,
                        val_csv_col=fluo_csv_col,
                        qc_csv_col=fluo_qc_csv_col,
                        default_err=CHLA_Adjusted_ERROR_est,
                        is_broken=is_broken,
                        fallback_to_raw=False,
                    )


def write_BBP700_adjusted(bgc_file, WMOfloatid, idx_profile, df_bio):
    if "BBP700_ADJUSTED" not in bgc_file.variables:
        return

    is_broken = is_broken_sensor_float(WMOfloatid, idx_profile)
    cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]
    bbp_csv_col = next((c for c in ["BBP700_FINAL", "BBP700_ADJUSTED"] if c in df_bio.columns), None)
    bbp_qc_csv_col = next((c for c in ["BBP700_FINAL_QC", "BBP700_ADJUSTED_QC"] if c in df_bio.columns), None)

    for prof_bbp in get_profiles_for_parameter(bgc_file, "BBP700"):
        mode = get_param_data_mode(bgc_file, prof_bbp, "BBP700")
        if mode not in ["D", "A"]:
            continue

        _write_param_from_csv(
            ds=bgc_file,
            iprof=prof_bbp,
            raw_candidates=["BBP700"],
            adj_var_name="BBP700_ADJUSTED",
            cycle_df=cycle_df,
            val_csv_col=bbp_csv_col,
            qc_csv_col=bbp_qc_csv_col,
            default_err=BBP700_Adjusted_ERROR_est,
            is_broken=is_broken,
            fallback_to_raw=False,
        )


def write_DOXY_slope_drift(ds, float_df, target_cycle, WMOfloatid=None):
    if WMOfloatid in DOXY_EXCLUDE_FLOAT_IDS:
        return

    profs_doxy = get_profiles_for_parameter(ds, "DOXY")
    if not profs_doxy:
        return

    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)]
    slope_val = cycle_df["DOXY_SLOPE"].iloc[0] if (not cycle_df.empty and "DOXY_SLOPE" in cycle_df.columns) else np.nan
    drift_val = cycle_df["DOXY_DRIFT"].iloc[0] if (not cycle_df.empty and "DOXY_DRIFT" in cycle_df.columns) else np.nan

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
    for prof_doxy in profs_doxy:
        if get_param_data_mode(ds, prof_doxy, "DOXY") not in ["D", "A"]:
            continue

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

    profs_doxy = get_profiles_for_parameter(ds, "DOXY")
    if not profs_doxy:
        return

    cycle_df = float_df[float_df["CYCLE_NUMBER"] == int(target_cycle)].copy()
    doxy_csv_col = next((c for c in ["DOXY_FINAL", "DOXY_ADJUSTED"] if c in float_df.columns), None)
    doxy_qc_csv_col = next((c for c in ["DOXY_FINAL_QC", "DOXY_ADJUSTED_QC"] if c in float_df.columns), None)
    has_slope = "DOXY_SLOPE" in cycle_df.columns and not cycle_df["DOXY_SLOPE"].dropna().empty

    slope_val = float(cycle_df["DOXY_SLOPE"].dropna().iloc[0]) if has_slope else None
    drift_val = (
        float(cycle_df["DOXY_DRIFT"].dropna().iloc[0])
        if ("DOXY_DRIFT" in cycle_df.columns and not cycle_df["DOXY_DRIFT"].dropna().empty)
        else 0.0
    )

    for prof_doxy in profs_doxy:
        mode = get_param_data_mode(ds, prof_doxy, "DOXY")
        if mode not in ["D", "A"]:
            continue

        _write_param_from_csv(
            ds=ds,
            iprof=prof_doxy,
            raw_candidates=["DOXY"],
            adj_var_name="DOXY_ADJUSTED",
            cycle_df=cycle_df,
            val_csv_col=doxy_csv_col,
            qc_csv_col=doxy_qc_csv_col,
            err_csv_col="DOXY_ADJUSTED_ERROR" if "DOXY_ADJUSTED_ERROR" in cycle_df.columns else None,
            default_err=DOXY_Adjusted_ERROR_est_default,
            is_broken=False,
            fallback_to_raw=False,
            slope_val=slope_val,
            drift_val=drift_val,
        )


def update_profile_qc_grades(bgc_file, WMOfloatid=None):
    """
    Computes PROFILE_<PARAM>_QC independently for EACH iprof in N_PROF using that
    profile's own <PARAM>_ADJUSTED_QC[iprof, :] (if mode in D/A) or <PARAM>_QC[iprof, :] (if mode R).
    """
    bgc_file.set_auto_mask(False)
    try:
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

            matching_profs = set(get_profiles_for_parameter(bgc_file, param_base))
            adj_qc_name = f"{param_base}_ADJUSTED_QC"
            raw_qc_name = "PRES_QC" if param_base == "PRES" else f"{param_base}_QC"

            for iprof in range(n_prof):
                if matching_profs and iprof not in matching_profs:
                    grade_bytes = b" "
                else:
                    mode = get_param_data_mode(bgc_file, iprof, param_base)
                    target_qc_name = adj_qc_name if (mode in ["D", "A"] and adj_qc_name in bgc_file.variables) else raw_qc_name
                    param_is_present = has_valid_measurements(bgc_file, param_base, iprof)

                    if target_qc_name in bgc_file.variables:
                        target_qc_var = bgc_file.variables[target_qc_name]
                        target_qc_array = (
                            target_qc_var[iprof, :] if target_qc_var.ndim > 1 else target_qc_var[:]
                        )
                        prof_grade = get_profile_qc_grade(target_qc_array, parameter_present=param_is_present)
                        grade_bytes = prof_grade.encode("utf-8")
                    else:
                        grade_bytes = b" "

                if bgc_file.variables[var_name].ndim == 1:
                    bgc_file.variables[var_name][iprof] = grade_bytes
                elif bgc_file.variables[var_name].ndim == 2:
                    bgc_file.variables[var_name][iprof, 0] = grade_bytes
    finally:
        bgc_file.set_auto_mask(True)


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
    Enforces BGC-Argo Data Management & GDAC format checker rules across ALL N_PROF profiles:
        <PARAM> / <PARAM>_QC -> <PARAM>_ADJUSTED / <PARAM>_ADJUSTED_QC / <PARAM>_ADJUSTED_ERROR
    Runs with ds.set_auto_mask(False) so single-character |S1 flags ('4', '9', ' ') are never
    masked or lost during NetCDF4 read/write operations.
    - In 'R' mode (or if parameter is absent from STATION_PARAMETERS[iprof]):
      - _ADJUSTED = 99999.0, _ADJUSTED_ERROR = 99999.0, _ADJUSTED_QC = ' ' at all levels.
    - In 'D' or 'A' mode:
      - Blank padding (raw <PARAM>_QC == ' '): _ADJUSTED = 99999.0, _ADJUSTED_ERROR = 99999.0, _ADJUSTED_QC = ' '.
      - Raw <PARAM> missing (99999.0 / NaN / raw QC == '9'):
            _ADJUSTED = 99999.0, _ADJUSTED_ERROR = 99999.0, _ADJUSTED_QC = '9'.
      - Raw <PARAM> NOT missing (valid measurement & raw QC != '9'):
          - If _ADJUSTED is missing (99999.0 / NaN) or _ADJUSTED_QC is '9'/' '/invalid:
              _ADJUSTED = 99999.0, _ADJUSTED_ERROR = 99999.0, _ADJUSTED_QC = '4' (NEVER '9').
    """
    ds.set_auto_mask(False)
    try:
        n_prof = ds.dimensions["N_PROF"].size
        n_levels = ds.dimensions["N_LEVELS"].size
        prof_param_map = {iprof: get_station_parameters_for_profile(ds, iprof) for iprof in range(n_prof)}

        # Pass 1: Standardize _ADJUSTED_QC flags using raw <PARAM> and numeric <PARAM>_ADJUSTED
        for v_name in ds.variables.keys():
            if "_ADJUSTED" not in v_name or v_name.startswith("PRES_ADJUSTED") or not v_name.endswith("_ADJUSTED_QC"):
                continue

            var = ds.variables[v_name]
            if var.ndim != 2:
                continue

            base_param = v_name.replace("_ADJUSTED_QC", "")
            raw_candidates = [base_param]
            if base_param in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
                raw_candidates = [base_param, "FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]

            adj_num_name = f"{base_param}_ADJUSTED"

            for iprof in range(n_prof):
                param_in_prof = (
                    base_param in prof_param_map[iprof]
                    or any(p.startswith(base_param) for p in prof_param_map[iprof])
                )
                if not param_in_prof:
                    var[iprof, :] = np.full((n_levels,), b" " if var.dtype.kind == "S" else " ", dtype=var.dtype)
                    continue

                mode = get_param_data_mode(ds, iprof, base_param)
                if mode not in ["D", "A"]:
                    var[iprof, :] = np.full((n_levels,), b" " if var.dtype.kind == "S" else " ", dtype=var.dtype)
                    continue

                arr = np.array(var[iprof, :], copy=True)

                adj_num_slice = None
                adj_fill_val = MISSING_VAL
                if adj_num_name in ds.variables:
                    adj_v_obj = ds.variables[adj_num_name]
                    adj_fill_val = getattr(adj_v_obj, "_FillValue", MISSING_VAL)
                    adj_num_slice = np.array(adj_v_obj[iprof, :], dtype="float32")

                out_qc = np.full((n_levels,), b" " if var.dtype.kind == "S" else " ", dtype=var.dtype)

                for k in range(n_levels):
                    raw_val, raw_qc = _get_raw_val_and_qc(ds, raw_candidates, iprof, k)
                    raw_missing = np.isnan(raw_val) or raw_qc.strip() in ["", "9"]

                    # Blank padding only when raw_qc is blank AND raw value itself is missing
                    if raw_qc == " " and raw_missing:
                        out_qc[k] = b" " if var.dtype.kind == "S" else " "
                        continue

                    c = _decode_qc_char(arr[k]).strip()
                    adj_val_missing = (
                        (is_missing_value(adj_num_slice[k]) or adj_num_slice[k] == adj_fill_val)
                        if (adj_num_slice is not None and k < len(adj_num_slice))
                        else True
                    )

                    # Only set _ADJUSTED_QC = '9' when raw <PARAM> itself is missing!
                    if raw_missing:
                        out_qc[k] = b"9" if var.dtype.kind == "S" else "9"
                    elif adj_val_missing or c not in VALID_QC_CHARS or c in ["", "9"]:
                        # Raw data is NOT missing, so _ADJUSTED_QC can NEVER be '9' or ' ' -> set to '4'
                        out_qc[k] = b"4" if var.dtype.kind == "S" else "4"
                    else:
                        out_qc[k] = c.encode("utf-8") if var.dtype.kind == "S" else c

                var[iprof, :] = out_qc

        # Pass 2: Standardize numeric _ADJUSTED and _ADJUSTED_ERROR arrays to match finalized _ADJUSTED_QC
        for v_name in ds.variables.keys():
            if "_ADJUSTED" not in v_name or v_name.startswith("PRES_ADJUSTED"):
                continue
            if v_name.endswith("_ADJUSTED_QC"):
                continue

            var = ds.variables[v_name]
            if var.ndim != 2 or var.dtype.kind not in ["f", "i"]:
                continue

            base_param = v_name.replace("_ADJUSTED_ERROR", "").replace("_ADJUSTED", "")
            raw_candidates = [base_param]
            if base_param in ["CHLA_FLUORESCENCE", "FLUORESCENCE_CHLA"]:
                raw_candidates = [base_param, "FLUORESCENCE_CHLA", "CHLA_FLUORESCENCE"]

            fill_val = getattr(var, "_FillValue", MISSING_VAL)
            adj_qc_name = f"{base_param}_ADJUSTED_QC"
            has_adj_qc = adj_qc_name in ds.variables

            for iprof in range(n_prof):
                param_in_prof = (
                    base_param in prof_param_map[iprof]
                    or any(p.startswith(base_param) for p in prof_param_map[iprof])
                )
                if not param_in_prof:
                    var[iprof, :] = np.full((n_levels,), np.float32(fill_val), dtype="float32")
                    continue

                mode = get_param_data_mode(ds, iprof, base_param)
                if mode not in ["D", "A"]:
                    var[iprof, :] = np.full((n_levels,), np.float32(fill_val), dtype="float32")
                    continue

                data_slice = np.array(var[iprof, :], dtype="float32", copy=True)

                qc_slice = None
                if has_adj_qc:
                    qc_slice = np.array(ds.variables[adj_qc_name][iprof, :], copy=True)

                for k in range(n_levels):
                    raw_val, raw_qc = _get_raw_val_and_qc(ds, raw_candidates, iprof, k)
                    raw_missing = np.isnan(raw_val) or raw_qc.strip() in ["", "9"]

                    if raw_qc == " " and raw_missing:
                        data_slice[k] = np.float32(fill_val)
                        continue

                    val = data_slice[k]
                    qc_char = _decode_qc_char(qc_slice[k]).strip() if qc_slice is not None else ""

                    if (
                        raw_missing
                        or qc_char in ["4", "9", ""]
                        or is_missing_value(val)
                        or val == fill_val
                    ):
                        data_slice[k] = np.float32(MISSING_VAL)

                var[iprof, :] = data_slice

        ds.sync()
    finally:
        ds.set_auto_mask(True)


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


def write_history_metadata(bgc_file, profs_chla, profs_doxy, WMOfloatid=None):
    bgc_file.history = datetime.datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ creation")
    bgc_file.setncattr("comment_dmqc_operator", comment_dmqc_operator_chla)

    UTCcurrent = datetime.datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")

    for iprof_c in (profs_chla if profs_chla else [0]):
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
            iprof_c,
        )

    if WMOfloatid not in DOXY_EXCLUDE_FLOAT_IDS:
        for iprof_d in (profs_doxy if profs_doxy else [0]):
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
                iprof_d,
            )

    if "DATE_UPDATE" in bgc_file.variables:
        bgc_file.variables["DATE_UPDATE"][:] = nc.stringtochar(np.array(UTCcurrent, dtype="S14"))


# ==============================================================================
# SECTION 2.5: Mode-Aware Missing Data & Blank Alignment Auditor
# ==============================================================================

def check_bd_missing_data_mismatches(processed_bd_files, WMOfloatid):
    """
    Audits final BD NetCDF files across ALL N_PROF profiles according to Argo User's Manual
    & BGC QC Manual rules for <PARAM> in [DOXY, CHLA, CHLA_FLUORESCENCE, BBP700].
    """
    target_params = ["DOXY", "CHLA", "CHLA_FLUORESCENCE", "BBP700"]
    mismatches = []

    def get_slice(ds, v_name, iprof):
        if v_name not in ds.variables:
            return None
        v = ds.variables[v_name]
        arr = v[iprof, :] if v.ndim > 1 else v[:]
        return np.array(arr)

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            ds.set_auto_mask(False)
            cycle_num = int(fpath[-6:-3])
            n_levels = ds.dimensions["N_LEVELS"].size

            for param_base in target_params:
                matching_profs = get_profiles_for_parameter(ds, param_base)
                if not matching_profs:
                    continue

                raw_v = param_base
                if raw_v not in ds.variables and param_base == "CHLA_FLUORESCENCE" and "FLUORESCENCE_CHLA" in ds.variables:
                    raw_v = "FLUORESCENCE_CHLA"

                adj_v = f"{param_base}_ADJUSTED"
                if adj_v not in ds.variables and raw_v == "FLUORESCENCE_CHLA" and "FLUORESCENCE_CHLA_ADJUSTED" in ds.variables:
                    adj_v = "FLUORESCENCE_CHLA_ADJUSTED"

                raw_qc_v = f"{raw_v}_QC"
                adj_qc_v = f"{adj_v}_QC"

                for iprof in matching_profs:
                    mode = get_param_data_mode(ds, iprof, param_base)

                    pres_vals, _ = get_pres_and_qc_for_profile(ds, iprof)
                    raw_data = get_slice(ds, raw_v, iprof)
                    adj_data = get_slice(ds, adj_v, iprof)
                    raw_qc_data = get_slice(ds, raw_qc_v, iprof)
                    adj_qc_data = get_slice(ds, adj_qc_v, iprof)

                    if raw_data is None and adj_data is None:
                        continue

                    for k in range(n_levels):
                        p_val = pres_vals[k]
                        p_miss = is_missing_value(p_val)

                        r_val = raw_data[k] if raw_data is not None and k < len(raw_data) else MISSING_VAL
                        a_val = adj_data[k] if adj_data is not None and k < len(adj_data) else MISSING_VAL
                        r_qc = _decode_qc_char(raw_qc_data[k]) if raw_qc_data is not None and k < len(raw_qc_data) else " "
                        a_qc = _decode_qc_char(adj_qc_data[k]) if adj_qc_data is not None and k < len(adj_qc_data) else " "

                        r_val_miss = is_missing_value(r_val) or (r_qc.strip() == "9")
                        a_val_miss = is_missing_value(a_val)
                        r_qc_blank = (r_qc == " ")
                        a_qc_blank = (a_qc == " ")
                        a_qc_miss = (a_qc.strip() in ["", "9"])

                        if mode == "R":
                            if not a_val_miss or not a_qc_blank:
                                mismatches.append({
                                    "WMO": WMOfloatid,
                                    "CYCLE": cycle_num,
                                    "N_PROF": iprof,
                                    "LEVEL_IDX": k,
                                    "PRES": float(p_val) if not p_miss else np.nan,
                                    "PARAM": param_base,
                                    "MODE": mode,
                                    "ISSUE": "R_MODE_NON_BLANK_ADJUSTED",
                                    "RAW_QC": repr(r_qc),
                                    "ADJ_QC": repr(a_qc),
                                })
                            continue

                        if r_qc_blank != a_qc_blank:
                            mismatches.append({
                                "WMO": WMOfloatid,
                                "CYCLE": cycle_num,
                                "N_PROF": iprof,
                                "LEVEL_IDX": k,
                                "PRES": float(p_val) if not p_miss else np.nan,
                                "PARAM": param_base,
                                "MODE": mode,
                                "ISSUE": "BLANK_PADDING_MISMATCH",
                                "RAW_QC": repr(r_qc),
                                "ADJ_QC": repr(a_qc),
                            })
                        elif not r_qc_blank and not a_qc_blank:
                            if r_val_miss and (not a_val_miss or a_qc.strip() != "9"):
                                mismatches.append({
                                    "WMO": WMOfloatid,
                                    "CYCLE": cycle_num,
                                    "N_PROF": iprof,
                                    "LEVEL_IDX": k,
                                    "PRES": float(p_val) if not p_miss else np.nan,
                                    "PARAM": param_base,
                                    "MODE": mode,
                                    "ISSUE": "RAW_MISSING_ADJ_NOT_9",
                                    "RAW_QC": repr(r_qc),
                                    "ADJ_QC": repr(a_qc),
                                })
                            elif not r_val_miss and a_qc_miss:
                                mismatches.append({
                                    "WMO": WMOfloatid,
                                    "CYCLE": cycle_num,
                                    "N_PROF": iprof,
                                    "LEVEL_IDX": k,
                                    "PRES": float(p_val) if not p_miss else np.nan,
                                    "PARAM": param_base,
                                    "MODE": mode,
                                    "ISSUE": "RAW_PRESENT_ADJ_QC_9",
                                    "RAW_QC": repr(r_qc),
                                    "ADJ_QC": repr(a_qc),
                                })
                            elif not r_val_miss and a_val_miss and a_qc.strip() != "4":
                                mismatches.append({
                                    "WMO": WMOfloatid,
                                    "CYCLE": cycle_num,
                                    "N_PROF": iprof,
                                    "LEVEL_IDX": k,
                                    "PRES": float(p_val) if not p_miss else np.nan,
                                    "PARAM": param_base,
                                    "MODE": mode,
                                    "ISSUE": "RAW_PRESENT_ADJ_MISSING_NOT_4",
                                    "RAW_QC": repr(r_qc),
                                    "ADJ_QC": repr(a_qc),
                                })

            ds.close()
        except Exception:
            traceback.print_exc()

    df_mismatches = pd.DataFrame(mismatches)

    if not df_mismatches.empty:
        print(f"\n[{WMOfloatid}] Found {len(df_mismatches)} missing/blank alignment mismatches:")
        for _, row in df_mismatches.iterrows():
            print(
                f"  Float {row['WMO']} | Cycle {row['CYCLE']:03d} | N_PROF={row['N_PROF']} | "
                f"lvl={row['LEVEL_IDX']:03d} | PRES={row['PRES']} | {row['PARAM']} (mode={row['MODE']}) | "
                f"{row['ISSUE']} | {row['PARAM']}_QC={row['RAW_QC']} vs {row['PARAM']}_ADJUSTED_QC={row['ADJ_QC']}"
            )

    return df_mismatches


# ==============================================================================
# SECTION 3: Plotly DMQC Status Plot Generation
# ==============================================================================

data_type_color_map = {
    "D": "#647ed5",
    "A": "#c07ab8",
    "R": "#ff4d6d",
    "99999": "#808080",
    " ": "#808080",
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
    "4": "#d73027", "5": "#404040", "8": "#7fbf7b", "9": "#808080",
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
    "9": "#808080", " ": "#e0e0e0",
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
        print(f"[{WMOfloatid}] Saved DMQC status plot to {out_png}")
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
            label="Missing Data (99999)",
        )
        ax_data.legend(loc="upper right", fontsize=8, frameon=True, facecolor="white", edgecolor="none")

    ax_data.set_ylim(max_pres, 0)
    ax_data.set_ylabel("Pressure (dbar)", fontsize=11)
    ax_data.grid(True, linestyle="--", alpha=0.4)

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
        fontsize=9,
    )

    today_stamp = dt.now().strftime("%Y-%m-%d")
    plt.suptitle(
        f"BGC-Argo Float {WMOfloatid} - {title_prefix} Section DMQC Summary (Verified BD Output)\nGenerated on: {today_stamp}",
        fontsize=13, fontweight="bold", y=0.98,
    )

    out_file = os.path.join(output_dir, f"Float_{WMOfloatid}_{file_suffix}_Section_Summary.png")
    plt.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close()


def extract_bd_parameter_dataframe(processed_bd_files, param_name, adj_param_name=None):
    records = []

    for fpath in processed_bd_files:
        try:
            ds = nc.Dataset(fpath, "r")
            ds.set_auto_mask(False)
            cycle_num = int(fpath[-6:-3])

            matching_profs = get_profiles_for_parameter(ds, param_name)
            if not matching_profs:
                matching_profs = [detect_parameter_profile(ds, param_name)]

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

            for iprof in matching_profs:
                pres_arr, _ = get_pres_and_qc_for_profile(ds, iprof)
                n_levels = len(pres_arr)

                def get_arr(v_name):
                    if v_name not in ds.variables:
                        return None
                    v = ds.variables[v_name]
                    arr = v[iprof, :] if v.ndim > 1 else v[:]
                    if hasattr(arr, "filled"):
                        arr = arr.filled(MISSING_VAL if arr.dtype.kind in ["f", "i"] else b" ")
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
                        "N_PROF": iprof,
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
            print(f"[{WMOfloatid}] Saved {param_base} R_mode & D_mode section plots to {output_root}")
        except Exception:
            traceback.print_exc()


# ==============================================================================
# SECTION 4: Processing Loop
# ==============================================================================

processed_output_paths = []
print(f"Master plots & CSV output directory created at: {OUTPUT_PLOTS_CSV_DIR}")

for idx_f, WMOfloatid in enumerate(WMO_FLOAT_IDS, start=1):
    main_float_dir = f"/data/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
    bio_dmqc_csv_path = (
        f"/data/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_BBP_DOXY_{WMOfloatid}.csv"
    )

    output_lut_dir = os.path.join(main_float_dir, "LUT")
    today_str = dt.now().strftime("%Y-%m-%d")
    final_out_dir = os.path.join(output_lut_dir, today_str)

    if not os.path.exists(bio_dmqc_csv_path):
        print(f"[{WMOfloatid}] Skipping: CSV not found at {bio_dmqc_csv_path}")
        continue

    os.makedirs(main_float_dir, exist_ok=True)
    os.makedirs(output_lut_dir, exist_ok=True)
    os.makedirs(final_out_dir, exist_ok=True)

    df_bio = pd.read_csv(bio_dmqc_csv_path, low_memory=False)
    df_bio["CYCLE_NUMBER"] = pd.to_numeric(df_bio["CYCLE_NUMBER"], errors="coerce")

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
        print(f"[{WMOfloatid}] Skipping: No BD or BR files found in {INPUT_STORAGE_DIR}")
        continue

    sorted_b_files = organize_b_files(all_bd_files, all_br_files)
    new_bd_files = []

    for bgc_filename in sorted_b_files:
        idx_profile = int(bgc_filename[-6:-3])

        w_bgc_filename = create_working_bd_file(
            bgc_filename, dest_dir=main_float_dir
        )

        try:
            ds = nc.Dataset(w_bgc_filename, "r+")

            profs_chla = get_profiles_for_parameter(ds, "CHLA")
            profs_doxy = get_profiles_for_parameter(ds, "DOXY")

            cycle_df = df_bio.loc[df_bio["CYCLE_NUMBER"] == idx_profile]

            # 1. Set PARAMETER_DATA_MODE and DATA_MODE to 'D' for all target BGC parameters across all N_PROF
            write_parameter_data_modes(
                ds,
                WMOfloatid=WMOfloatid,
                cycle_num=idx_profile,
                cycle_df=cycle_df,
            )

            # 2. Populate <PARAM>_ADJUSTED and <PARAM>_ADJUSTED_QC across ALL matching N_PROF entries
            write_scientific_calib_chla(ds, WMOfloatid, idx_profile, df_bio)
            write_chla_BBP_adjusted(ds, WMOfloatid, idx_profile, df_bio)
            write_BBP700_adjusted(ds, WMOfloatid, idx_profile, df_bio)
            write_DOXY_slope_drift(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)
            write_DOXY_from_csv(ds, df_bio, idx_profile, WMOfloatid=WMOfloatid)

            # 3. Revert unpopulated 'D' parameters back to 'R'
            enforce_empty_adjusted_modes_to_realtime(ds)

            # 4. Enforce strict per-parameter alignment between <PARAM>/<PARAM>_QC and <PARAM>_ADJUSTED/<PARAM>_ADJUSTED_QC
            standardize_missing_data_and_qc(ds)

            # 5. Guarantee strict GDAC consistency between DATA_MODE and PARAMETER_DATA_MODE (never 'D'/'R')
            sync_data_mode_with_parameter_modes(ds, WMOfloatid=WMOfloatid)

            # 6. Update PROFILE_<PARAM>_QC independently for each N_PROF
            update_profile_qc_grades(ds, WMOfloatid=WMOfloatid)

            # 7. Metadata and attributes
            dm_after = ds.variables["DATA_MODE"][:]
            dm_after_arr = dm_after.filled(b" ") if hasattr(dm_after, "filled") else dm_after
            modes_after = [_decode_qc_char(m).strip() for m in np.ravel(dm_after_arr)]
            prof_is_d = [m == "D" for m in modes_after]

            write_history_metadata(
                ds,
                profs_chla,
                profs_doxy,
                WMOfloatid=WMOfloatid,
            )
            for i, is_d in enumerate(prof_is_d):
                if is_d:
                    ds.variables["DATA_STATE_INDICATOR"][i, :] = nc.stringtochar(
                        np.array("2C  ", dtype="S4")
                    )

            add_missing_valid_range_attributes(ds)
            remove_forbidden_attributes(ds)

            ds.sync()
            ds.close()

            # Compact STATION_PARAMETERS (remove internal empty slots) & shrink oversized N_PARAM
            compact_station_parameters_and_nparam(w_bgc_filename)
            new_bd_files.append(w_bgc_filename)

        except Exception:
            traceback.print_exc()
            if "ds" in locals() and ds.isopen():
                ds.close()

    # Step 8: Always finalize output NetCDF files as BD*.nc
    final_processed_files = []
    for w_file in new_bd_files:
        path, name = os.path.split(w_file)
        new_name = re.sub(r"^w_(AOML_)?B[DR]", "BD", name)
        new_path = os.path.join(final_out_dir, new_name)
        safe_rename(w_file, new_path)
        processed_output_paths.append(new_path)
        final_processed_files.append(new_path)

    print(f"[{WMOfloatid}] Wrote {len(final_processed_files)} BD files to {final_out_dir}")

    # --------------------------------------------------------------------------
    # Step 8.5: Audit Missing Value / Blank Padding Alignment (DOXY, CHLA, CHLA_FLUORESCENCE, BBP700)
    # --------------------------------------------------------------------------
    df_missing_mismatches = check_bd_missing_data_mismatches(final_processed_files, WMOfloatid)

    if not df_missing_mismatches.empty:
        audit_csv = os.path.join(OUTPUT_PLOTS_CSV_DIR, f"{WMOfloatid}_missing_data_mismatches.csv")
        df_missing_mismatches.to_csv(audit_csv, index=False)

    # --------------------------------------------------------------------------
    # Step 9: Verification Plotting (DMQC Status + R_mode / D_mode Section Summaries)
    # --------------------------------------------------------------------------
    generate_dmqc_status_plot(WMOfloatid, final_processed_files, OUTPUT_PLOTS_CSV_DIR)
    generate_all_qc_summary_plots(final_processed_files, WMOfloatid, OUTPUT_PLOTS_CSV_DIR)