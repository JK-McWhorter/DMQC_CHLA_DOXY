###############################################################################################
###############################################################################################
######  Function for noisy BBP profile test ###########
###############################################################################################
###############################################################################################


# Transcripted by M.Cornec (2023) from BBP doc :
# Giorgio DALL'OLMO et al (2017). Real-time quality control of optical 
# backscattering data from BGC-Argo floats
# https://doi.org/10.12688/openreseurope.15047.1

# "To flag profiles that are affected by noisy data. 
# This noise could indicate sensor malfunctioning, spikes caused by organisms 
# attracted to the light emittedby the BBP sensor (Haëntjens et al., 2020), 
# or other anomalous conditions."

### INPUTS :
# bbp : vector of the bbp measurements
# dep_bbp : vector of the depths of bbp measurements

### REQUIREMENT:
# Need RunningFilter function

### OUTPUTS :
# bbp_qc : values of QC vector for the whole profile
# 1 : no noisy profile
# 3 : noisy profile 

bbp_noisy<-function(bbp,dep_bbp) {
  # Initialize QC : default = good data
  bbp_qc<-rep(1,length(bbp))
  
  # Calculate the median filtered bbp (over 5 points)
  bbp_filt<-RunningFilter(2,bbp,ends.fill=T, Method="Median")
  
  # Calculate the absolute residuals between raw and filtered data
  res<-abs(bbp-bbp_filt)
  
  # subset the residuals to data below 100 m (above can be noisy for other 
  # reasons)
  res_sub<-res[which(dep_bbp>100)]
  
  # test if more than 10% of the values are above the treshold of 0.0005 m-1
  if(length(which(res_sub>0.0005))>(0.1*length(res_sub))){
    # if so, profile QC = 3, need further DMQC
    bbp_qc<-rep(3,length(bbp))
  }
  
  return(bbp_qc)
  
}