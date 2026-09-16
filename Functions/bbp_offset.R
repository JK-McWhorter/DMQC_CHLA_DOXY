###############################################################################################
###############################################################################################
######  Function for BBP Offset test ###########
###############################################################################################
###############################################################################################


# Written by M.Cornec (2018) adapted from BBP doc :
# Catherine SCHMECHTIG, Emmanuel BOSS, Nathan BRIGGS, Hervé CLAUSTRE, Giorgio
# DALL'OLMO, Antoine POTEAU (2017). BGC Argo quality control manual for particles
# backscattering. https://doi.org/10.13155/60262

### INPUTS :
# bbp : vector of the bbp measurements

### REQUIREMENT:
# Need RunningFilter function

### OUTPUTS :
# bbp_qc : QC adjusted values after despike correction

bbp_offset<-function(bbp) {
  # Initialize QC 
  bbp_qc<-rep(1,length(bbp))
  
  # Calculate the median filtered bbp
  bbp_filt<-RunningFilter(2,bbp,ends.fill=T, Method="Median")
  
  # Calculate min value of the filtered bbp
  min_bbp<-min(bbp_filt,na.rm=T)
  
  # attribute QC 2 to the values -20*bbpmin if the min value is out of range
  if(min_bbp< -2.5*10^-5){
    bbq_qc[which(bbp < -20*min_bbp)]<-3
    bbq_qc[which(bbp > -20*min_bbp)]<-2
  }
  
  return(bbp_qc)
  
}