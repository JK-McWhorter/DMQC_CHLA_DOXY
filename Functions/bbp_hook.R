###############################################################################################
###############################################################################################
######  Function for BBP Parking Hook test ###########
###############################################################################################
###############################################################################################


# Transcripted by M.Cornec (2023) from BBP doc :
# Giorgio DALL'OLMO et al (2017). Real-time quality control of optical 
# backscattering data from BGC-Argo floats
# https://doi.org/10.12688/openreseurope.15047.1

# "When the float is drifting with the currents while at its parking pressure 
# (typically 1000 dbar), particles may be depositing on the float and BBP sensor. 
# These accumulated particles are likely released back into the water when the 
# float descends to its maximum pressure (typically 2000 dbar), before starting 
# the ascending profile during which data are collected. 
# However, if the float does not descend to 2000 dbar before starting the BBP
# measurements, but immediately starts ascending towards the surface and 
# measuring, then the accumulated particles might be measured by the BBP sensor 
# as they are released back into the water. 
# This is the likely cause of an increase in BBP at the start of the profile, 
# when the parking pressure is close to the maximum pressure." 

### INPUTS :
# bbp : vector of the bbp measurements
# dep_bbp : vector of the depths of bbp measurements
# dep_park : value of parking depth (PARK_PRES, extracted from the mission 
# configuration )

### REQUIREMENT:
# Need RunningFilter function

### OUTPUTS :
# bbp_qc : vector of values of QC
# 1 : no hook values
# 4 : bad values due to hook effect

bbp_hook<-function(bbp,dep_bbp,dep_park) {
  
  # Initialize QC : default = good data 
  bbp_qc<-rep(1,length(bbp))
  
  # define max depth
  max_depth<-max(dep_bbp,na.rm=T)
  
  # Test if the test can be applied : next non-NA measurement above last 
  # depth must be closer than 20m
  if(min(abs(dep_bbp[which(is.na(bbp)==F)]-max_depth),na.rm=T)<=20){
    
    # test if the profile starts from the parking depth
    if(abs(max_depth-dep_park)<100){
      
      # definition of a pressure range to determine the baseline
      baseline_top<-max_depth - 50
      baseline_bot<-max_depth - 20
      
      # calculate the baseline; 0.0002 m-1 corresponds to the deviation
      baseline<-median(
        bbp[which(dep_bbp>=baseline_top & dep_bbp<baseline_bot)],na.rm=T
        )+0.0002
      
      # Apply QC 4 to values below the baseline layer that are above baseline 
      # value
      bbp_qc[which(dep_bbp>baseline_bot & bbp>baseline)]<-4
    }
    
  }
  
  return(bbp_qc)
  
}