#!/usr/bin/env python
# coding: utf-8

### UPDATES
#    
#    11-Sep-26 BN: Updated to add flags, coefficients, comments, equations
#    10-OCT-25 MS: Program creation date?
#
#
#



### To run me make sure you have the BR files in the input directory /a1/ARGO_DELAY/DMQC_BGC/data/{float you are running}
#
### upload the csv file to /a1/ARGO_DELAY/DMQC_BGC/data/csv with this format CHLA_LUT_7901009.csv
#
### UPDATE THIS LINE WMOfloatid = 6999992 #make dynamic
#
### DEPENDING ON THE FLOAT TYPE YOU MIGHT NEED TO UPDATE THIS TOO iprof = 1
#
### THIS COMMAND ACTUALLY RUNS THE PROGRAM python3.9 BDfiller_CHLA_LUT.py 




# In[19]:


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


# In[20]:


#25 - 0.07 -error estimates VIIRS vs raw CHLA 
#22 - 0.18 -error estimates VIIRS vs raw CHLA 
#09 - 0.42 -error estimates VIIRS vs raw CHLA 


# In[21]:


#WMOfloatid = 7901009
#float_dir = f"/Users/madison.soden/BGCARGO_work/data/{WMOfloatid}"
#CHLA_Adjusted_ERROR_est = 0.42
#BBP_Adjusted_ERROR = 9999.0
#iprof = 0


# In[22]:


#WMOfloatid = 4903622
#CHLA_Adjusted_ERROR_est = 0.18
#BBP_Adjusted_ERROR = 9999.0
#iprof = 1


# In[23]:


#WMOfloatid = 4903625
#CHLA_Adjusted_ERROR_est = 0.07
#BBP_Adjusted_ERROR = 9999.0
#iprof = 1


# In[45]:


WMOfloatid = 6999992 #make dynamic
float_dir = f"/a1/ARGO_DELAY/DMQC_BGC/data/{WMOfloatid}/"
profile_chla_qc = 'A'
profile_bbp_qc =  'A'


# other settings that should apply to all institutions and floats
#odv_filename = float_dir + f"ODV{WMOfloatid}QC.TXT" #unnecessary maybe? commenting out for now
bio_dmqc_csv_path = f"/a1/ARGO_DELAY/DMQC_BGC/data/csv/CHLA_DOXY_{WMOfloatid}.csv"
#bio_dmqc_csv_path = f"/a1/ARGO_DELAY/DMQC_BGC/data/csv/OUTPUT_{WMOfloatid}_dmqc_chla_bbp.csv"
#'/Users/madison.soden/Downloads/OUTPUT_4903622_dmqc_chla_bbp.csv'
comment_dmqc_operator = "PRIMARY | https://orcid.org/0009-0006-7862-6267 | Brandon Navarro, NOAA/AOML;" 

history_parameter = "CHLA"
history_institution = "AO" # AO for AOML
history_reference = "MULT" #datasets? HPLC? 
history_software= "RBIO" 
history_software_release = "2025"

#scientific_calibration_comment_CHLA = "Dark correction, NPQ with ML validation (Xing 2012, MLD criterion only), no slope adjustment"
scientific_calibration_comment_CHLA = "CHLA adjustment (Slope: specified in http://dx.doi.org/10.13155/35385 and computed with MLD_LIMIT = 0.03, and following recommendations of Sauzede et al., 2025 (https://doi.org/10.17882/105732), Quenching: Xing et al., 2018, Terrats et al., 2020)"
scientific_calibration_comment_CHLA_FLU = "CHLA_FLUORESCENCE  (specified in http://dx.doi.org/10.13155/35385 and computed with MLD_LIMIT = 0.03)"
scientific_calibration_equation_CHLA = "CHLA_ADJUSTED = CHLA_NPQ for PRES in [0, ZMaxFluo ], CHLA_ADJUSTED = ((FLUORESCENCE_CHLA-MEDIAN(PRELIM_DARK_CHLA)*SCALE_CHLA)/PHYSIO_RATIO"
scientific_calibration_coefficient_CHLA = "PHYSIO_RATIO=1.0"
#scientific_calibration_coefficient_CHLA ="PRELIM_DARK_CHLA = [57 59 60 60 60], SCALE_CHLA = 0.0073, PHYSIO_RATIO = 2.1"
 
CHLA_Adjusted_ERROR_est = 0.07

#scientific_calibration_comment_BBP700 = "TBD"
#scientific_calibration_equation_BBP700 = "TBD" 
#scientific_calibration_coefficient_BBP700 = "TBD" 
#BBP_Adjusted_ERROR = 9999.0 #file checker rejects if you use fill value. It rejects if you are missing the variable, as well. 

iprof = 0
data_state_indicator = ['2','C','','']
parameter_data_mode = "D" 


# In[46]:


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


# In[47]:


def get_juld(filename):
    '''Retrieve and return the JULD value from the given file.
    NOTE: It is always taken from the value with index "iprof".'''
    b_file = netCDF4.Dataset(filename, 'r')
    juld = b_file.variables['JULD'][iprof]
    b_file.close()
    return juld


# In[48]:


def get_profile(filename):
    '''Extract the profile index from filename, return it as an int.'''
    profile = filename[-6:-3]
    return int(profile)


# In[49]:


def organize_b_files(bd_files, br_files):
    '''Sort the B*nc files by profile. If both BR and BD files exist for a given profile,
    consider only the BD file. Return the sorted list.'''
    ptr_bd = 0
    ptr_br = 0
    

    if not br_files: 
        max_br = 0
    else: 
        max_br = get_profile(br_files[-1])
    
    if not bd_files: 
        max_prof= max_br
    else: 
        max_bd = get_profile(bd_files[-1])
        max_prof = max(max_bd, max_br)
        print(max_bd, max_br)
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


# In[ ]:


def update_history(nc, dct, iprof):
    # Using code from Christopher Gordon's library,
    # but with an added third argument iprof:
        '''
       Update HISTORY_<PARAM> values in an Argo netCDF file for the specified profile
        '''

        hix = nc.dimensions['N_HISTORY'].size
        for name, value in dct.items():
            nc[name][hix,iprof,:] = CGnetcdf.string_to_array(value, nc.dimensions[nc[name].dimensions[-1]])
            #print(name)


# In[ ]:


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

        # specific to Dmode. This code is institution dependent.

        # specific to Dmode. This code is institution dependent. 
        # Probably not even necessary to use for SAGE. I ended up using the year of the most recent update for SAGE.
        # SAGE does not appeare to have release numbers, saving update year as HISTORY_SOFTWARE_RELEASE 

    # Reference data that is being fed into SAGE. This code is institution dependent.
    history_action = "IP"
        # See ARGO User Manual section 3.7 Reference table 7: history action codes 

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

    update_history(bgc_file, {    "HISTORY_INSTITUTION": history_institution, # PM
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


# In[52]:


def write_parameter_data_mode(bgc_file):
    '''Modify the PARAMETER_DATA_MODE variable in the BD file with the given handle.
    Set it to D for DOXY.
    Note: Uses iprof (assigned above).'''
    ParameterList = bgc_file.variables["STATION_PARAMETERS"][iprof].data.astype(str)
    PDMarray = bgc_file.variables["PARAMETER_DATA_MODE"][iprof,:]
    
    for j in range(bgc_file.dimensions["N_PARAM"].size): 
        PARAMstr = ''.join(ParameterList[j]);
        if PARAMstr[:4] == "CHLA": 
            #print(PARAMstr)
            PDMarray[j] = parameter_data_mode
        #elif PARAMstr[:6] == "BBP700": 
        #    #print(PARAMstr)
        #    PDMarray[j] = parameter_data_mode
            
    bgc_file.variables["PARAMETER_DATA_MODE"][iprof,:] = PDMarray
    # adjust the overall DATA_MODE as well (it may be 'D' already)
    bgc_file.variables["DATA_MODE"][iprof] = 'D'


# In[53]:

def write_scientific_calib(bgc_file, idx_profile):
    '''Write the SCIENTIFIC_CALIB_* variables for CHLA in the NetCDF file
    using the cycle-specific physio_ratio from the CSV.'''
    
    # 1. Fetch physio_ratio from the CSV for this cycle
    df_bio = pd.read_csv(bio_dmqc_csv_path)
     
    dark_cols = [f'MIN_FLUOCHLA_CYCLE{i}' for i in range(1, 6)]
    if all(col in df_bio.columns for col in dark_cols):
        dark_vals = df_bio[dark_cols].iloc[0].astype(int).tolist()
        dark_str = " ".join(map(str, dark_vals))
    else:
        dark_str = "NA"
    
    if 'SCALE_CHLA' in df_bio.columns:
        scale_val = df_bio['SCALE_CHLA'].iloc[0]
    else:
        scale_val = "NA"
    
    cycle_df = df_bio.loc[df_bio['CYCLE_NUMBER'] == idx_profile]

    if not cycle_df.empty and 'physio_ratio' in cycle_df.columns:
        physio_val = cycle_df['PHYSIO_RATIO'].iloc[0]
        #calib_coefficient = f"physio_ratio = {physio_val}"
    else:
        #calib_coefficient = scientific_calibration_coefficient_CHLA
        physio_val = "1"

    calib_coefficient = f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}, PHYSIO_RATIO = {physio_val}"
    calib_coefficient_flu =f"PRELIM_DARK_CHLA = [{dark_str}], SCALE_CHLA = {scale_val}"
    
    UTCcurrent = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")

    # 2. Format character arrays matching STRING256 and DATE_TIME dimensions
    str256_len = bgc_file.dimensions["STRING256"].size
    date_len = bgc_file.dimensions["DATE_TIME"].size

    # SCIENTIFIC_CALIB_COMMENT
    SciCalComArray_CHLA = np.ma.empty(shape=(str256_len), dtype='|S1')
    SciCalComArray_CHLA[:] = ''
    SciCalComArray_CHLA.mask = True
    SciCalComArray_CHLA[:len(scientific_calibration_comment_CHLA)] = list(scientific_calibration_comment_CHLA)
   
    # SCIENTIFIC_CALIB_COMMENT
    SciCalComArray_CHLA_FLU = np.ma.empty(shape=(str256_len), dtype='|S1')
    SciCalComArray_CHLA_FLU[:] = ''
    SciCalComArray_CHLA_FLU.mask = True
    SciCalComArray_CHLA_FLU[:len(scientific_calibration_comment_CHLA_FLU)] = list(scientific_calibration_comment_CHLA_FLU)


    # SCIENTIFIC_CALIB_EQUATION
    SciCalEquArray_CHLA = np.ma.empty(shape=(str256_len), dtype='|S1')
    SciCalEquArray_CHLA[:] = ''
    SciCalEquArray_CHLA.mask = True
    SciCalEquArray_CHLA[:len(scientific_calibration_equation_CHLA)] = list(scientific_calibration_equation_CHLA)

    # SCIENTIFIC_CALIB_COEFFICIENT
    SciCalCoeArray_CHLA = np.ma.empty(shape=(str256_len), dtype='|S1')
    SciCalCoeArray_CHLA[:] = ''
    SciCalCoeArray_CHLA.mask = True
    SciCalCoeArray_CHLA[:len(calib_coefficient)] = list(calib_coefficient)


    # SCIENTIFIC_CALIB_COEFFICIENT
    SciCalCoeArray_CHLA_FLU = np.ma.empty(shape=(str256_len), dtype='|S1')
    SciCalCoeArray_CHLA_FLU[:] = ''
    SciCalCoeArray_CHLA_FLU.mask = True
    SciCalCoeArray_CHLA_FLU[:len(calib_coefficient_flu)] = list(calib_coefficient_flu)


    # SCIENTIFIC_CALIB_DATE
    SciCalDateArray = CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])

    # 3. Search across ALL profiles (iprof_idx) to locate CHLA
    n_prof_size = bgc_file.dimensions["N_PROF"].size
    n_param_size = bgc_file.dimensions["N_PARAM"].size

    chla_found = False
    for iprof_idx in range(n_prof_size):
        param_list = bgc_file.variables["STATION_PARAMETERS"][iprof_idx].data.astype(str)
        
        for j in range(n_param_size):
            param_str = ''.join(param_list[j]).strip()
            #print(param_str) 
            if param_str.startswith('CHLA_F'):
                # Write directly to the matching NetCDF indices [iprof_idx, 0, j, :]
                #print('did I work')
                bgc_file.variables['SCIENTIFIC_CALIB_COMMENT'][iprof_idx, 0, j, :] = SciCalComArray_CHLA_FLU 
                #bgc_file.variables['SCIENTIFIC_CALIB_EQUATION'][iprof_idx, 0, j, :] = SciCalEquArray_CHLA_FLU
                bgc_file.variables['SCIENTIFIC_CALIB_COEFFICIENT'][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA_FLU
                bgc_file.variables['SCIENTIFIC_CALIB_DATE'][iprof_idx, 0, j, :] = SciCalDateArray
            elif param_str.startswith('CHLA'):
                # Write directly to the matching NetCDF indices [iprof_idx, 0, j, :]
                bgc_file.variables['SCIENTIFIC_CALIB_COMMENT'][iprof_idx, 0, j, :] = SciCalComArray_CHLA 
                bgc_file.variables['SCIENTIFIC_CALIB_EQUATION'][iprof_idx, 0, j, :] = SciCalEquArray_CHLA 
                bgc_file.variables['SCIENTIFIC_CALIB_COEFFICIENT'][iprof_idx, 0, j, :] = SciCalCoeArray_CHLA 
                bgc_file.variables['SCIENTIFIC_CALIB_DATE'][iprof_idx, 0, j, :] = SciCalDateArray
                
                #print(f"Updated SCIENTIFIC_CALIB_COEFFICIENT for CHLA at profile index {iprof_idx}, param slot {j}: '{calib_coefficient}'")
                chla_found = True
                #break
        #if chla_found:
        #    break

    if not chla_found:
        print(f"Warning: CHLA parameter was not found in file for profile cycle {idx_profile}.")

#    # SCIENTIFIC_CALIB_COMMENT
#    SciCalComArray_CHLA = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
#    SciCalComArray_CHLA[:] = ''
#    SciCalComArray_CHLA.mask = True
#    SciCalComArray_CHLA[:len(calib_comment)] = list(calib_comment)
##
#    # SCIENTIFIC_CALIB_EQUATION
#    SciCalEquArray_CHLA  =  np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
#    SciCalEquArray_CHLA[:] = ''
#    SciCalEquArray_CHLA.mask = True
#    SciCalEquArray_CHLA[:len(scientific_calibration_equation_CHLA )] = list(scientific_calibration_equation_CHLA )
#
#    # SCIENTIFIC_CALIB_COEFFICIENT
#    SciCalCoeArray_CHLA  = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
#    SciCalCoeArray_CHLA[:] = ''
#    SciCalCoeArray_CHLA.mask = True
#    SciCalCoeArray_CHLA[:len(scientific_calibration_coefficient_CHLA )] = list(scientific_calibration_coefficient_CHLA )
#
#    
##BBP700
    # SCIENTIFIC_CALIB_COMMENT
    #SciCalComArray_BBP700 = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    #SciCalComArray_BBP700 [:] = ''
    #SciCalComArray_BBP700 .mask = True
    #SciCalComArray_BBP700 [:len(scientific_calibration_comment_BBP700 )] = list(scientific_calibration_comment_BBP700)

    # SCIENTIFIC_CALIB_EQUATION
    #SciCalEquArray_BBP700  =  np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    #SciCalEquArray_BBP700 [:] = ''
    #SciCalEquArray_BBP700 .mask = True
    #SciCalEquArray_BBP700 [:len(scientific_calibration_equation_BBP700 )] = list(scientific_calibration_equation_BBP700 )

    # SCIENTIFIC_CALIB_COEFFICIENT
    #SciCalCoeArray_BBP700  = np.ma.empty(shape = (bgc_file.dimensions["STRING256"].size), dtype='|S1')
    #SciCalCoeArray_BBP700 [:] = ''
    #SciCalCoeArray_BBP700 .mask = True
    #SciCalCoeArray_BBP700 [:len(scientific_calibration_coefficient_BBP700 )] = list(scientific_calibration_coefficient_BBP700 )

    

#    ParameterList= bgc_file.variables["STATION_PARAMETERS"][iprof].data.astype(str)
#    for j in range(bgc_file.dimensions["N_PARAM"].size#): 
#        PARAMstr= ''.join(ParameterList[j]);
#        #print(f"PARAMSTRING: {PARAMstr}")
#        if PARAMstr[:4] == 'CHLA': 
#            bgc_file.variables['SCIENTIFIC_CALIB_COMMENT'][iprof,0,j, :] = SciCalComArray_CHLA 
#            bgc_file.variables['SCIENTIFIC_CALIB_EQUATION'][iprof, 0, j, :] = SciCalEquArray_CHLA 
#            bgc_file.variables['SCIENTIFIC_CALIB_COEFFICIENT'][iprof, 0, j, :] = SciCalCoeArray_CHLA 
#            bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof,0,j,:] = \
#                CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])
#        #elif PARAMstr[:6] == 'BBP700':
#        #    bgc_file.variables['SCIENTIFIC_CALIB_COMMENT'][iprof,0,j,] = SciCalComArray_BBP700
#        #    bgc_file.variables['SCIENTIFIC_CALIB_EQUATION'][iprof, 0, j,] = SciCalEquArray_BBP700
#        #    bgc_file.variables['SCIENTIFIC_CALIB_COEFFICIENT'][iprof, 0, j,] = SciCalCoeArray_BBP700
#        #    bgc_file.variables["SCIENTIFIC_CALIB_DATE"][iprof,0,j,:] = \
#        #        CGnetcdf.string_to_array(UTCcurrent, bgc_file.dimensions["DATE_TIME"])


# In[54]:


def write_chla_BBP_adjusted(bgc_file, idx_profile): 
    # setting PH_ADJUSTED and CHLA_ADJUSTED_QC 
    # BBP700_ADJUSTED and BBP700_ADJUSTED_QC
    df_bio = pd.read_csv( bio_dmqc_csv_path)
    #df_bio.columns = df_bio.columns.str.strip().str.lower()
    print(f"idx_profile: {idx_profile}, Type: {type(idx_profile)}")
    df_bio = df_bio.loc[df_bio['CYCLE_NUMBER'] == idx_profile]
    #df_bio.drop(columns=["Unnamed: 0", "float_num", "LONGITUDE", "LATITUDE"], inplace=True)
    print(df_bio)
    
    #initializing masked arrays
    CHLA_Adjusted_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    CHLA_Adjusted_Array[:] = 99999.0
    CHLA_Adjusted_Array.mask = True
    
    CHLA_AdjustedQC_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), dtype='|S1')
    CHLA_AdjustedQC_Array[:] = ''
    CHLA_AdjustedQC_Array.mask = True

    CHLA_Adjusted_ERROR_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    CHLA_Adjusted_ERROR_Array[:] = 99999.0
    CHLA_Adjusted_ERROR_Array.mask = True
    
    CHLA_FLUORESCENCE_Adjusted_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    CHLA_FLUORESCENCE_Adjusted_Array[:] = 99999.0
    CHLA_FLUORESCENCE_Adjusted_Array.mask = True
    
    CHLA_FLUORESCENCE_AdjustedQC_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), dtype='|S1')
    CHLA_FLUORESCENCE_AdjustedQC_Array[:] = ''
    CHLA_FLUORESCENCE_AdjustedQC_Array.mask = True

    CHLA_FLUORESCENCE_Adjusted_ERROR_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    CHLA_FLUORESCENCE_Adjusted_ERROR_Array[:] = 99999.0
    CHLA_FLUORESCENCE_Adjusted_ERROR_Array.mask = True
    
    CHLA_Adjusted_physio_ratio = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), dtype='float32')
    CHLA_Adjusted_physio_ratio[:] = 99999.0
    CHLA_Adjusted_physio_ratio.mask = True



    #BBP700_Adjusted_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    #BBP700_Adjusted_Array[:] = 99999.0
    #BBP700_Adjusted_Array.mask = True
    
    #BBP700_AdjustedQC_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), dtype='|S1')
    #BBP700_AdjustedQC_Array[:] = ''
    #BBP700_AdjustedQC_Array.mask = True

    #BBP700_Adjusted_ERROR_Array = np.ma.empty(shape = (bgc_file.dimensions["N_LEVELS"].size), fill_value=99999.0,dtype='float32')
    #BBP700_Adjusted_ERROR_Array[:] = 99999.0
    #BBP700_Adjusted_ERROR_Array.mask = True

    assigned_nc_pres_vals = set()
    bad_match_nc_pres_vals = set()
    
    for row in range(len(df_bio)):
        for i in range(bgc_file.dimensions["N_LEVELS"].size):
            csv_pres = np.float32(df_bio.iloc[row]['PRES'])
            nc_pres = np.float32(bgc_file.variables['PRES'][iprof,i]) 
            #print(f"v_#{row}  CSV_Pres = _{csv_pres}_ and BD_Pres = _{nc_pres}_")

            if nc_pres in assigned_nc_pres_vals: 
                #print(f"v_#{row}  Skipped {nc_pres}")
                continue
                
            if csv_pres == nc_pres:
                CHLA_Adjusted_Array[i]         = np.float32(df_bio.iloc[row]['CHLA_FINAL'])
                #CHLA_Adjusted_Array[i]         = np.float32(df_bio.iloc[row]['CHLA_LUT'])
                #print(type(df_bio.iloc[row]['CHLA_FINAL_QC']))
                CHLA_AdjustedQC_Array[i]       = df_bio.iloc[row]['CHLA_FINAL_QC']
                #CHLA_AdjustedQC_Array[i]       = df_bio.iloc[row]['CHLA_FINAL_QC'] #df_ODV.iloc[row,OxygenQC_i]
                CHLA_Adjusted_ERROR_Array[i]   = np.float32(CHLA_Adjusted_ERROR_est)
                CHLA_Adjusted_physio_ratio[i]  = np.float32(df_bio.iloc[row]['PHYSIO_RATIO'])
                #CHLA_Adjusted_scale_chla[i]    = np.float32(df_bio.iloc[row]['scale_chla'])
                #print(CHLA_Adjusted_physio_ratio[i])
                CHLA_FLUORESCENCE_Adjusted_Array[i]         = np.float32(df_bio.iloc[row]['CHLA_FLUORESCENCE'])
                #CHLA_FLUORESCENCE_AdjustedQC_Array[i]       = df_bio.iloc[row]['CHLA_FINAL_QC']
                CHLA_FLUORESCENCE_Adjusted_ERROR_Array[i]   = np.float32(CHLA_Adjusted_ERROR_est)
                #BBP700_Adjusted_Array[i]       = np.float32(df_bio.iloc[row]['BBP700'])
                #BBP700_AdjustedQC_Array[i]     = '1'
                #BBP700_Adjusted_ERROR_Array[i] = 9999.0
                #print(f"v_#{row}  CSV_Pres = _{csv_pres}_ and BD_Pres = _{nc_pres}_")
                #print(f"v_#{row}  CHLA_ADJ = {CHLA_Adjusted_Array[i]} set")
                #print(f"v_#{row}  BBP_ADJ  = {BBP700_Adjusted_Array[i]} set")
                assigned_nc_pres_vals.add(nc_pres)
                break
            elif nc_pres not in bad_match_nc_pres_vals: 
                #print(f"v_#{row} bad match CSV_Pres = _{csv_pres}_ and BD_Pres = _{nc_pres}_")
                bad_match_nc_pres_vals.add(nc_pres)
                
    bgc_file.variables["CHLA_ADJUSTED"][iprof]         = CHLA_Adjusted_Array
    bgc_file.variables["CHLA_ADJUSTED_QC"][iprof]      = CHLA_AdjustedQC_Array
    bgc_file.variables["CHLA_ADJUSTED_ERROR"][iprof]   = CHLA_Adjusted_ERROR_Array
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED"][iprof]         = CHLA_FLUORESCENCE_Adjusted_Array
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_QC"][iprof]      = CHLA_FLUORESCENCE_AdjustedQC_Array
    bgc_file.variables["CHLA_FLUORESCENCE_ADJUSTED_ERROR"][iprof]   = CHLA_FLUORESCENCE_Adjusted_ERROR_Array
    #bgc_file.variables["BBP700_ADJUSTED"][iprof]       = BBP700_Adjusted_Array
    #bgc_file.variables["BBP700_ADJUSTED_QC"][iprof]    = BBP700_AdjustedQC_Array
    #bgc_file.variables["BBP700_ADJUSTED_ERROR"][iprof] = BBP700_Adjusted_ERROR_Array
    print(f"WMO{WMOfloatid} BAD PRES MATCHS: \n    {bad_match_nc_pres_vals}")
    


# In[55]:


#loop over all BR files - they must be present in "float_dir" already
# FIXME: if both BR and BD file exist for the same profile, both
# will be processed (only the BD file should be used!)
all_bd_files   = glob.glob(f'{float_dir}/BD*{WMOfloatid}_*.nc')
all_bd_files.sort()
all_br_files   = glob.glob(f'{float_dir}/AOML_BR*{WMOfloatid}_*.nc')
all_br_files.sort()
all_d_files    = glob.glob(f'{float_dir}/D*{WMOfloatid}_*.nc')
all_d_files.sort()
sorted_b_files = organize_b_files(all_bd_files, all_br_files)
print( f'{len(sorted_b_files)} relevant B files found in {float_dir}')


# In[56]:


juld_init = get_juld(sorted_b_files[0])    
new_bd_files = list()                                 

for bgc_filename in sorted_b_files:
    print(f'Processing {bgc_filename}')
    idx_profile = int(bgc_filename[-6:-3])
    w_bgc_filename = create_working_bd_file(bgc_filename)    
    
    bgc_file = netCDF4.Dataset(w_bgc_filename, 'a') # 'a' (append) mode automatically reads previous netcdf format
    
    #write_history(bgc_file)
    write_parameter_data_mode(bgc_file)
    
    juld = get_juld(bgc_filename)
    write_scientific_calib(bgc_file, idx_profile)
    #REMOVE SCIENTIFIC CALIB COMMENTS or just WriteNULL? 

    
    #setting DATA_STATE_INDICATOR variable
    for i in range(bgc_file.dimensions["N_PROF"].size):
        if i != iprof: # FIXME is this always correct?
            continue
        bgc_file.variables["DATA_STATE_INDICATOR"][i] = np.ma.array(data_state_indicator, mask=[False, False, True, True], 
                                                                    dtype='|S1')
    
    #populate CHLA_ADJUSTED and BBP700_ADJUSTED  
    write_chla_BBP_adjusted(bgc_file, idx_profile)
    
    # extract variables from corresponding physical profile file
    #phys_filename = get_phys_filename(bgc_filename)
    #pres_phys_raw = get_phys_raw_pres(phys_filename)
    #dens_phys, psal, temp = get_dens(phys_filename)

    # need to determine which pTS profile index to use
    #pres_bgc = bgc_file.variables['PRES'][:,:]
    #iprof_phys = get_iprof_phys(pres_phys_raw, pres_bgc)
    
    # calculate DOXY_ADJUSTED_ERROR in umol/kg    
    #write_doxy_adjusted_error(bgc_file, doxy_adj_err, psal[iprof_phys,:], temp[iprof_phys,:], 
                             # pres_bgc[iprof,:], dens_phys[iprof_phys,:], doxy_adjusted)
    
    # assign a new overall profile QC flag
    bgc_file.variables['PROFILE_CHLA_QC'][iprof] = profile_chla_qc
    #bgc_file.variables['PROFILE_BBP700_QC'][iprof] = profile_bbp_qc
    
    bgc_file.close()
    new_bd_files.append(w_bgc_filename)


# In[57]:


# change w_BDfileName to BDfileName after everything is done
for file in new_bd_files:
    path, name = os.path.split(file)
    new_name = name.replace('w_AOML_BD', 'BD')
    new_path = f'{float_dir}/LUT/{new_name}'
    print(f'Renaming {file} to {new_path}')
    os.rename(file, new_path)
    
#for file in new_bd_files:
#    path, name = os.path.split(file)
#    new_names = name.replace('w_AOML_BD', 'BD')
#    new_paths = f'{float_dir}/LUT/{new_names}'
#    print(f'Renaming {file} to {new_paths}')
#    os.rename(file, new_paths)
