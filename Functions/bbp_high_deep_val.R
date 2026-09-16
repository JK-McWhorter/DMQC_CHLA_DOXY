###############################################################################################
###############################################################################################
######  Function for High BBP value at depth test ###########
###############################################################################################
###############################################################################################


# Transcripted by M.Cornec (2023) from BBP doc :
# Giorgio DALL'OLMO et al (2017). Real-time quality control of optical 
# backscattering data from BGC-Argo floats
# https://doi.org/10.12688/openreseurope.15047.1

# "To flag profiles with anomalously high BBP values at depth. 
# High values at deeper depths could indicate a variety of problems, 
# including biofouling,incorrect calibration coefficients, sensor malfunctioning.
# Note that high deep BBP values could also be valid data, 
# for example in the case of sediment-resuspension events."

### INPUTS :
# bbp : vector of the bbp measurements
# dep_bbp : vector of the depths of bbp measurements

### REQUIREMENT:
# Need RunningFilter function

### OUTPUTS :
# bbp_qc : vector of values of QC for the whole profile
# 1 : no high deep values
# 3 : high deep values --> need further DM QC

bbp_high_deep_val<-function(bbp,dep_bbp) {
  
  # Initialize QC : default = good data 
  bbp_qc<-rep(1,length(bbp))
  
  # Calculate the median filtered bbp (over 5 points)
  bbp_filt<-RunningFilter(2,bbp,ends.fill=T, Method="Median")
  
  # define min. treshold value : 5*10-4m-1 = half of the value
  #typical for surface BBP in the oligotrophic ocean (Dall’Olmo et al., 2012)
  tresh<-0.0005
  
  # define the depth below which to track high values
  depth_id<-700
  
  # attribute QC 3 to the whole profile if more that 5 values above treshold
  # below the deep depth
  if(length(bbp_filt[which(dep_bbp>=depth_id &
                           bbp_filt>tresh)])>5){
    bbp_qc<-rep(3,length(bbp))
  }
  
  return(bbp_qc)
  
}