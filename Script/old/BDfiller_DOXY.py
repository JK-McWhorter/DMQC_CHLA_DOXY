#!/usr/bin/env python
# coding: utf-8

# In[3]:


# import of libraries
import bgcArgoDMQC
import datetime
import glob
import gsw
import os
import netCDF4 
import shutil
import numpy as np
import pandas as pd
from bgcArgoDMQC.io import netcdf as CGnetcdf


# In[4]:


# Global variables & settings
# library settings
np.set_printoptions(threshold=np.inf);


# In[5]:


# select institution and float type
#inst_float = 'pmel_apex'
inst_float = 'aoml_apex'
#inst_float = 'aoml_apex'
floatid = -999999 # flag
# change this if the profile is not "overall good":
profile_doxy_qc = 'A'


# In[6]:


# institution and float specific settings
# PMEL, APEX floats (testing only)
if inst_float == 'pmel_apex':
    floatid = 5906293
    # directory settings
    float_dir = f"/raid/frenzel/DMQC/ARGO_PROCESSING/DATA/ARGO_REPO/{floatid}/QC/"
    #format: PRIMARY | https://orcid.org/0000-0001-5766-1668 | Tanya Maurer, MBARI
    comment_dmqc_operator = "PRIMARY | https://orcid.org/0000-0002-3481-0867 | Hartmut Frenzel, PMEL" 

    # must be 256 characters or less:
    scientific_calibration_comment = 'Polynomial calibration coefficients were used. G determined by comparison of O2sat surface measurement to WOA18; DOXY_ADJUSTED_ERROR is set to 2 mbar. See Takeshita et al., 2013, doi:10.1002/jgrc.20399'

    history_institution = "PM"  # PM for PMEL
    history_reference = "WOA18" # FIXME not really true since they have in-air calibration, but this is for testing only!

    doxy_adj_err = 2 # in mbar, for floats with in-air calibration
    
    # FIXME - is this always true?
    iprof = 1 # use the second entry in the N_PROF dimension


# In[7]:


# institution and float specific settings
# PMEL, Navis floats
if inst_float == 'pmel_navis':
    floatid = 4903499
    #floatid = 4903500
    float_dir = f"/raid/frenzel/DMQC/ARGO_PROCESSING/DATA/ARGO_REPO/{floatid}/QC/"
    #format: PRIMARY | https://orcid.org/0000-0001-5766-1668 | Tanya Maurer, MBARI
    comment_dmqc_operator = "PRIMARY | https://orcid.org/0000-0002-3481-0867 | Hartmut Frenzel, PMEL"
    # must be 256 characters or less:
    scientific_calibration_comment = 'Polynomial calibration coefficients were used. G determined by comparison of O2sat surface measurement to WOA18; DOXY_ADJUSTED_ERROR is set to 5 mbar. See Takeshita et al., 2013, doi:10.1002/jgrc.20399'

    history_institution = "PM"  # PM for PMEL
    history_reference = "WOA18"
    
    doxy_adj_err = 5 # in mbar, for floats without in-air calibration
    
    # FIXME - is this always true?
    iprof = 0 # use the first entry in the N_PROF dimension


# In[8]:


# institution and float specific settings
# AOML, APEX floats
if inst_float == 'aoml_apex':
    floatid = 4903622
    
    float_dir = f"/Users/madison.soden/BGCARGO_work/data/{floatid}/"
    #DnetCDF_path = /Users/madison.soden/BGCARGO_work/data/reference/netcdf/5905988/BD5905988_012.nc"
    #RnetCDF_path1 = "/Users/madison.soden/BGCARGO_work/data/4903622/BR4903622_001.nc"
    #RnetCDF_path2 = "/Users/madison.soden/BGCARGO_work/data/4903622/BR4903622_002.nc"
    #format: PRIMARY | https://orcid.org/0000-0001-5766-1668 | Tanya Maurer, MBARI
    comment_dmqc_operator = "PRIMARY | https://orcid.org/0009-0001-9960-129X | Madison Soden, AOML" 
    # must be 256 characters or less
    scientific_calibration_comment = 'Polynomial calibration coefficients were used. G determined by surface measurement comparison to NCEP'
    
    history_institution = "AO" # AO for AOML
    history_reference = "WOA2023"
    
    doxy_adj_err = 2 # in mbar, for floats with in-air calibration
    
    # FIXME - is this always true?
    iprof = 1 # use the second entry in the N_PROF dimension


# In[9]:


# institution and float specific settings
# AOML, APEX floats
if inst_float == 'aoml_navis':
    floatid = 7901009
    
    float_dir = f"/Users/madison.soden/BGCARGO_work/data/{floatid}/"
    #DnetCDF_path = /Users/madison.soden/BGCARGO_work/data/reference/netcdf/5905988/BD5905988_012.nc"
    #RnetCDF_path1 = "/Users/madison.soden/BGCARGO_work/data/4903622/BR4903622_001.nc"
    #RnetCDF_path2 = "/Users/madison.soden/BGCARGO_work/data/4903622/BR4903622_002.nc"
    #format: PRIMARY | https://orcid.org/0000-0001-5766-1668 | Tanya Maurer, MBARI
    comment_dmqc_operator = "PRIMARY | https://orcid.org/0009-0001-9960-129X | Madison Soden, AOML" 
    # must be 256 characters or less
    scientific_calibration_comment = 'Polynomial calibration coefficients were used. G determined by surface measurement comparison to NCEP'
    
    history_institution = "AO" # AO for AOML
    history_reference = "WOA2023"
    
    doxy_adj_err = 5 # in mbar, for floats with in-air calibration
    
    # FIXME - is this always true?
    iprof = 0 # use the second entry in the N_PROF dimension


# In[10]:


if floatid < 0:
    raise ValueError('You must select one of the pre-defined institution/float cases')


# In[11]:


# other settings that should apply to all institutions and floats
odv_filename = float_dir + f"ODV{floatid}QC.TXT"

data_state_indicator = ['2','C','','']
    # (Ref: Argo Users Manual pgs 85-86) Effectively most cases will change this value from 2B -> 2C for Dmode processing

parameter_data_mode = "D" 
    # either "R" for real time, or "D" for delayed mode
    
scientific_calibration_equation = 'DOXY_ADJUSTED = DOXY*G; G = G_INIT + G_DRIFT*(JULD_PROF - JULD_INIT)/365' # MBARI version place holder
    # must be 256 characters or less


# In[12]:


def get_gain_odvqc(filename, variable='DOXY'):
    '''This function retrieves the initial gain and drift values
    for the specified variable.
    Returns gain and drift.'''
    # FIXME written and tested only for DOXY for now!
    if variable == 'DOXY':
        odv_var = 'Oxygen'
    else:
        raise ValueError(f'not coded for variable {variable}')
    with open(filename, 'r') as f_in:
        lines = f_in.readlines()
    # find the segment with QC correction values
    for count, line in enumerate(lines):
        if line.startswith('//QUALITY CONTROLLED DATA CORRECTIONS:'):
            break
    # now find the line(s) for the specified variable
    # FIXME currently only written for exactly one line with gain factor
    for line in lines[count+2:]:
        if line.startswith('//' + odv_var):
            values = line.split()
            break
    if int(values[1]) != 1 or abs(float(values[3])) > 1e-10:
        raise ValueError(f'unexpected values: {values}')
    return float(values[2]), float(values[4]) # gain, drift


# In[13]:


def get_juld(filename):
    '''Retrieve and return the JULD value from the given file.
    NOTE: It is always taken from the value with index "iprof".'''
    b_file = netCDF4.Dataset(filename, 'r')
    juld = b_file.variables['JULD'][iprof]
    b_file.close()
    return juld


# In[14]:


def get_iprof_phys(pres_raw, pres_bgc):
    '''Determine which profile from the physical PRES values matches
    those from the bgc file. Returns the (zero-offset) value if found,
    throws an IOError if not found.'''
    iprof_phys = -1 # "not found" flag
    for col in range(pres_raw.shape[0]):
        # raw PRES values from phys and bgc files should be identical;
        # we'll allow some numerical error
        if max(abs(pres_raw[col,:] - pres_bgc[iprof,:])) < 0.1:
            iprof_phys = col
    if iprof_phys < 0:
        raise IOError('Matching PRES values not found')
    else:
        print(f'Using profile {iprof_phys+1} of {phys_filename} to determine density.')
        return iprof_phys


# In[44]:


def create_working_bd_file(filename):
    '''Create a working copy of the B file with a leading 'w_' in the name.
    Note that a possibly existing working copy will be overwritten without
    warning.
    Return the name (with full path) of the working copy.'''
    path, name = os.path.split(filename)
    bd_name = name.replace('BR', 'BD') # no effect if file is named BD*nc already
    w_filename = f'{path}/w_{bd_name}'
    shutil.copyfile(filename, w_filename)
    return w_filename


# In[73]:


def get_phys_filename(bgc_filename):
    '''Determine the name of the physical profile file that corresponds
    to the given bgc profile file and return it.
    Throws an IOError if no corresponding file is found.'''
    path, name = os.path.split(bgc_filename)
    # it doesn't matter here if bgc file is named BR* or BD*
    phys_filename = float_dir + '/D/' + name.replace(name[0:2], 'D') # try D file first
    if not os.path.exists(phys_filename):
        phys_filename = float_dir + '/R/' + name.replace(name[0:2], 'R') # try R file if necessary
        print(f'Using phys R file:{phys_filename}')
        if not os.path.exists(phys_filename):
            raise IOError(f'No corresponding phys file found for {bgc_filename}')
    return phys_filename        


# In[74]:


def get_phys_raw_pres(phys_filename):
    '''Extract the raw PRES values from the given physical profile. 
    Return the full array.
    Pre: File must exist; no check is performed.'''
    phys_file = netCDF4.Dataset(phys_filename, 'r') # read-only access is enough
    pres_raw = phys_file.variables["PRES"][:,:] # for comparison with PRES from BD file
    phys_file.close()
    return pres_raw


# In[75]:


def get_dens(phys_filename, verbose=False):
    '''Extract pTS (adjusted if available) from the physical profile file
    with the given name. Compute and return density. Also return salinity and temperature.
    Pre: File must exist; no check is performed.
    FIXME: So far, only existence of adjusted pTS is checked, not the validity
    of the values. (If at least one of them is all NaNs, computed density
    will be all NaNs as well.)'''
    phys_file = netCDF4.Dataset(phys_filename, 'r') # read-only access is enough    
    num_adj = 0
    try:
        temp_var = phys_file.variables["TEMP_ADJUSTED"]
        num_adj += 1
    except KeyError:
        print('TEMP_ADJUSTED not found, using TEMP instead to determine density')
        temp_var = phys_file.variables["TEMP"]
    try:
        psal_var = phys_file.variables["PSAL_ADJUSTED"]
        num_adj += 1
    except KeyError:
        print('PSAL_ADJUSTED not found, using PSAL instead to determine density')
        psal_var = phys_file.variables["PSAL"]
    try:
        pres_var = phys_file.variables["PRES_ADJUSTED"]
        num_adj += 1
    except KeyError:
        print('PRES_ADJUSTED not found, using PRES instead to determine density')
        temp_var = phys_file.variables["PRES"]
    if num_adj == 3:
        if verbose:
            print('Using ADJUSTED values of p,T,S to determine density')
    elif num_adj == 0:
        print('Using RAW values of p,T,S to determine density')
    else:
        print(f'Found {num_adj} variables - case not coded!')
        raise IOError('Unexpected file contents')
    temp = temp_var[:,:]
    psal = psal_var[:,:]
    pres = pres_var[:,:]
    phys_file.close()
    return gsw.rho_t_exact(psal, temp, pres), psal, temp


# In[76]:


def set_sci_calib_coefficient(gain, juld, juld_init):
    '''Fills the given values into the scientific calibration coefficient,
    which must be a string at most 256 characters long.'''
    return f'G_INIT = {gain:.4f}; G_DRIFT = {drift:.4f}; JULD_PROF = {juld:.4f}; JULD_INIT = {juld_init:.4f}'


# In[77]:


# Using code from Christopher Gordon's library,
# but with an added third argument iprof:
def update_history(nc, dct, iprof):
        '''
       Update HISTORY_<PARAM> values in an Argo netCDF file for the specified profile
        '''

        hix = nc.dimensions['N_HISTORY'].size
        for name, value in dct.items():
            nc[name][hix,iprof,:] = CGnetcdf.string_to_array(value, nc.dimensions[nc[name].dimensions[-1]])


# In[78]:


def write_history(bgc_file):
    '''Write global attributes and HISTORY* variables to the BD file with the given handle.
    Also writes out DATE_CREATION and DATE_UPDATE variables.'''
    # history_institution and history_reference must be set as global variables
    
    # Setting Global Attributes for BD file 
    # Argo Users Manual V3.41.1 §2.2.1

    # format: 'YYYY-MM-DDTHH:MM:SSZ creation' EG: '2023-03-13T19:25:43Z creation'
    bgc_file.history = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ creation")
    bgc_file.setncattr('comment_dmqc_operator', comment_dmqc_operator)
        
    history_step = "ARSQ" 
        # specific to Dmode. See ARGO User Manual section 3.12 Reference table 12: history steps codes 
    history_software= "SAGE" 
        # specific to Dmode. This code is institution dependent.
    history_software_release = "2024"
        # specific to Dmode. This code is institution dependent. 
        # Probably not even necessary to use for SAGE. I ended up using the year of the most recent update for SAGE.
        # SAGE does not appeare to have release numbers, saving update year as HISTORY_SOFTWARE_RELEASE 

    # Reference data that is being fed into SAGE. This code is institution dependent.
    history_action = "IP"
        # See ARGO User Manual section 3.7 Reference table 7: history action codes 
    history_parameter = "DOXY" # specific to Dmode DOXY

    # Setting History Variables  
    # Using code from Christopher Gordon's library,
    # but with an added third argument iprof:
    '''
    Ref https://github.com/ArgoCanada/bgcArgoDMQC/blob/main/bgcArgoDMQC/io/netcdf.py
    Christopher Gordon
    Fisheries and Oceans Canada
    chris.gordon@dfo-mpo.gc.ca

    def update_history(nc, dct, iprof):
        ''
       Update HISTORY_<PARAM> values in an Argo netCDF file for the specified profile
        ''

        hix = nc.dimensions['N_HISTORY'].size
        for name, value in dct.items():
            nc[name][hix,iprof,:] = string_to_array(value, nc.dimensions[nc[name].dimensions[-1]])
    '''

    update_history(bgc_file, {"HISTORY_INSTITUTION": history_institution, # PM
                            "HISTORY_STEP": history_step, # specific to Dmode
                            "HISTORY_SOFTWARE": history_software, # specific to Dmode
                            "HISTORY_SOFTWARE_RELEASE": history_software_release, # specific to Dmode
                            "HISTORY_REFERENCE": history_reference, # Reference data that is being fed into SAGE, Institution dependent
                            "HISTORY_DATE": datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S"),
                            "HISTORY_ACTION": history_action,
                            "HISTORY_PARAMETER": history_parameter}, iprof) # specific to Dmode DOXY

    # Setting DATE_CREATION and DATE_UPDATE variables; format: YYYYMMDDHHMISS (UTC)
    UTCcurrent= datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")
    bgc_file.variables["DATE_CREATION"][:] = CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])
    bgc_file.variables["DATE_UPDATE"][:] = CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])


# In[79]:


def write_parameter_data_mode(bgc_file):
    '''Modify the PARAMETER_DATA_MODE variable in the BD file with the given handle.
    Set it to D for DOXY.
    Note: Uses iprof (assigned above).'''
    ParameterList = bgc_file.variables["STATION_PARAMETERS"][iprof].data.astype(str)
    PDMarray = bgc_file.variables["PARAMETER_DATA_MODE"][iprof,:]

    for j in range(bgc_file.dimensions["N_PARAM"].size): 
        PARAMstr = ''.join(ParameterList[j]);
        if PARAMstr[:4] == 'DOXY': 
            PDMarray[j] = parameter_data_mode
      
    bgc_file.variables["PARAMETER_DATA_MODE"][iprof,:] = PDMarray
    # adjust the overall DATA_MODE as well (it may be 'D' already)
    bgc_file.variables["DATA_MODE"][iprof] = 'D'


# In[80]:


def write_scientific_calib(bgc_file, scientific_calibration_coefficient):
    '''Write the SCIENTIFIC_CALIB_* (DATE, COMMENT, EQUATION, and COEFFICIENT)
    variables in the BD file with the given handle.
    scientific_calibration_comment and scientific_calibration_equation must
    be assigned globally.
    Note: Uses iprof (assigned globally).'''
    # SCIENTIFIC_CALIB_DATE
    UTCcurrent = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")
    #MS bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof,:] = CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])

    # SCIENTIFIC_CALIB_COMMENT
    SciCalComArray = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    SciCalComArray[:] = ''
    SciCalComArray.mask = True
    SciCalComArray[:len(scientific_calibration_comment)] = list(scientific_calibration_comment)

    # SCIENTIFIC_CALIB_EQUATION
    SciCalEquArray =  np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    SciCalEquArray[:] = ''
    SciCalEquArray.mask = True
    SciCalEquArray[:len(scientific_calibration_equation)] = list(scientific_calibration_equation)

    # SCIENTIFIC_CALIB_COEFFICIENT
    SciCalCoeArray = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    SciCalCoeArray[:] = ''
    SciCalCoeArray.mask = True
    SciCalCoeArray[:len(scientific_calibration_coefficient)] = list(scientific_calibration_coefficient)

    ParameterList= bgc_file.variables["STATION_PARAMETERS"][iprof].data.astype(str)
    for j in range(bgc_file.dimensions["N_PARAM"].size): 
        PARAMstr= ''.join(ParameterList[j]);
        if PARAMstr[:4] == 'DOXY': 
            bgc_file.variables['SCIENTIFIC_CALIB_COMMENT'][iprof,0,j,] = SciCalComArray
            bgc_file.variables['SCIENTIFIC_CALIB_EQUATION'][iprof, 0, j,] = SciCalEquArray
            bgc_file.variables['SCIENTIFIC_CALIB_COEFFICIENT'][iprof, 0, j,] = SciCalCoeArray
            bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof,0,j,:] = \
                CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])


# In[81]:


def write_doxy_adjusted(bgc_file, gain):
    '''Read DOXY values, compute DOXY_ADJUSTED = DOXY * gain, 
    and write them to the BD file with the given file handle. 
    Also write DOXY_ADJUSTED_QC values.
    Return the masked array DOXY_Adjusted_Array.'''    
    doxy = bgc_file.variables['DOXY'][iprof,:]
    doxy_adjusted = doxy * gain;
    doxy_adjusted_qc = np.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), dtype='|S1')
    doxy_adjusted_qc[:] = ' ' # missing value
    doxy_adjusted_qc[~doxy_adjusted.mask] = '1'
    
    bgc_file.variables["DOXY_ADJUSTED"][iprof,:] = doxy_adjusted
    bgc_file.variables["DOXY_ADJUSTED_QC"][iprof,:] = doxy_adjusted_qc
    return doxy_adjusted


# In[82]:


def write_doxy_adjusted_error(bgc_file, doxy_adj_err_mbar, psal, temp, pres, dens, doxy_adj):
    '''Given the (scalar) value of doxy_adj_err_mbar (in mbar), compute it for all
    depth levels in umol/kg and write it to the given bgc_file.
    Return value is in umol/l.'''
    idx_good = np.where(~np.ma.getmask(psal) & ~np.ma.getmask(doxy_adj))
    # initialize the error with its full size
    doxy_adj_error = np.ma.masked_array(psal.shape[0]*[0.0], mask=np.ma.getmask(psal) | np.ma.getmask(doxy_adj), fill_value=99999.0)
    error = bgcArgoDMQC.calc_fixed_doxy_adjusted_error(psal[idx_good], temp[idx_good], pres[idx_good], 
                                                       fix_err=doxy_adj_err_mbar)
    doxy_adj_error[idx_good] = error
    # convert to umol/kg
    doxy_adj_error_umol_kg = doxy_adj_error * 1e3 / dens
    bgc_file.variables["DOXY_ADJUSTED_ERROR"][iprof,:] = doxy_adj_error_umol_kg


# In[83]:


def get_profile(filename):
    '''Extract the profile index from filename, return it as an int.'''
    profile = filename[-6:-3]
    return int(profile)


# In[84]:


def organize_b_files(bd_files, br_files):
    '''Sort the B*nc files by profile. If both BR and BD files exist for a given profile,
    consider only the BD file. Return the sorted list.'''
    ptr_bd = 0
    ptr_br = 0
    max_br = get_profile(br_files[-1])
    if not bd_files: 
        max_prof= max_br
    else: 
        max_bd = get_profile(bd_files[-1])
        max_prof = max(max_bd, max_br)
    sorted_b_files = []
    # 001 is the first actual profile file
    for idx in range(1, max_prof+1):
        file_found = False
        for ptr in range(ptr_bd, len(bd_files)):
            if get_profile(bd_files[ptr]) == idx:
                sorted_b_files.append(bd_files[ptr])
                ptr_bd = ptr+1
                file_found = True
                break
        if not file_found:
            for ptr in range(ptr_br, len(br_files)):
                if get_profile(br_files[ptr]) == idx:
                    sorted_b_files.append(br_files[ptr])
                    ptr_br = ptr+1
                    break
        # it is possible that neither BD nor BR file exist for an index
    return sorted_b_files


# In[85]:


# extract gain and drift values from ODV*QC.TXT file
gain, drift = get_gain_odvqc(odv_filename) # same values for all BD output files
print(f'Using gain = {gain}, drift = {drift}')


# In[94]:


# loop over all BR files - they must be present in "float_dir" already
# FIXME: if both BR and BD file exist for the same profile, both
# will be processed (only the BD file should be used!)
all_bd_files   = glob.glob(f'{float_dir}/BD/BD*{floatid}_*.nc')
all_bd_files.sort()
all_br_files   = glob.glob(f'{float_dir}/BR/BR*{floatid}_*.nc')
all_br_files.sort()
all_d_files    = glob.glob(f'{float_dir}/D/BR*{floatid}_*.nc')
all_d_files.sort()
sorted_b_files = organize_b_files(all_bd_files, all_br_files)

print(f'{len(sorted_b_files)} relevant B files found in {float_dir}')

juld_init = get_juld(sorted_b_files[0])
    
new_bd_files = list()                                 

for bgc_filename in sorted_b_files:
    print(f'Processing {bgc_filename}')
    idx_profile = int(bgc_filename[-6:-3])
    w_bgc_filename = create_working_bd_file(bgc_filename)    
    # 'a' (append) mode automatically reads previous netcdf format
    bgc_file = netCDF4.Dataset(w_bgc_filename, 'a')
    write_history(bgc_file)
    write_parameter_data_mode(bgc_file)
    juld = get_juld(bgc_filename)
    scientific_calibration_coefficient = set_sci_calib_coefficient(gain, juld, juld_init)
    write_scientific_calib(bgc_file, scientific_calibration_coefficient)
    
    # setting DATA_STATE_INDICATOR variable
    for i in range(bgc_file.dimensions["N_PROF"].size):
        if i != iprof: # FIXME is this always correct?
            continue
        bgc_file.variables["DATA_STATE_INDICATOR"][i] = np.ma.array(data_state_indicator, mask=[False, False, True, True], 
                                                                    dtype='|S1')
    
    # compute DOXY_ADJUSTED from DOXY and gain value and write it to BD file    
    doxy_adjusted = write_doxy_adjusted(bgc_file, gain)
    
    # extract variables from corresponding physical profile file
    phys_filename = get_phys_filename(bgc_filename)
    pres_phys_raw = get_phys_raw_pres(phys_filename)
    dens_phys, psal, temp = get_dens(phys_filename)

    # need to determine which pTS profile index to use
    pres_bgc = bgc_file.variables['PRES'][:,:]
    iprof_phys = get_iprof_phys(pres_phys_raw, pres_bgc)
    
    # calculate DOXY_ADJUSTED_ERROR in umol/kg    
    write_doxy_adjusted_error(bgc_file, doxy_adj_err, psal[iprof_phys,:], temp[iprof_phys,:], 
                              pres_bgc[iprof,:], dens_phys[iprof_phys,:], doxy_adjusted)
    
    # assign a new overall profile QC flag
    bgc_file.variables['PROFILE_DOXY_QC'][iprof] = profile_doxy_qc
    
    bgc_file.close()
    new_bd_files.append(w_bgc_filename)


# In[95]:


# change w_BDfileName to BDfileName after everything is done
for file in new_bd_files:
    path, name = os.path.split(file)
    new_name = name.replace('w_BD', 'BD')
    new_path = f'{float_dir}/BD/{new_name}'
    print(f'Renaming {file} to {new_path}')
    os.rename(file, new_path)
    


# In[ ]:


Users/madison.soden/BGCARGO_work/data/4903622/BD_output()


# In[ ]:





# In[ ]:




