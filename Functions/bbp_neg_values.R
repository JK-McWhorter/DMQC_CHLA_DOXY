###############################################################################################
###############################################################################################
######  Function for negative BBP values test ###########
###############################################################################################
###############################################################################################


# Transcripted by M.Cornec (2023) from BBP doc :
# Giorgio DALL'OLMO et al (2017). Real-time quality control of optical 
# backscattering data from BGC-Argo floats
# https://doi.org/10.12688/openreseurope.15047.1

# "To flag data points or profiles with negative BBP values due to a variety 
# of reasons including: sensor drift or malfunctioning, inaccurate calibration
# coefficients, or BBP sensor exposed to air." 

### INPUTS :
# bbp : vector of the bbp measurements
# dep_bbp : vector of the depths of bbp measurements

### OUTPUTS :
# bbp_qc : vector of QC values for the whole profile
# 1 : no negative values
# 3 : less than 10% of the profiles values below 5 m are negative
# 4 : more than 10% of the profiles values below 5 m are negative and/or in the
# surface


bbp_neg_val<-function(bbp,dep_bbp) {
  # Initialize QC : default = good data
  bbp_qc<-rep(1,length(bbp))
  
  # test if negative values in the profile
  if (length(which(bbp<0))>0){
    # if only negative values in the 5 first meters : values flagged as bad
    if(max(dep_bbp[which(bbp<0)],na.rm=T)<=5){
      bbp_qc[which(bbp<0)]<-4
    } else
    # if more that 10% of values below 5m are negative : all profile is flagged 
    # as bad 
    if(length(which(dep_bbp[which(bbp<0 & dep_bbp>5)] >5)) >= 
       (0.1*length(bbp[which(dep_bbp>5)]))){
      bbp_qc<-rep(4,length(bbp))
    } else
    # if less than 10% of values below 5m are negative : all profile is flagged 
    # as probably bad : need further DM QC
    if(length(which(dep_bbp[which(bbp<0 & dep_bbp>5)] >5)) < 
       (0.1*length(bbp[which(dep_bbp>5)]))){
      bbp_qc<-rep(3,length(bbp))
    }
  }
  
  return(bbp_qc)
  
}