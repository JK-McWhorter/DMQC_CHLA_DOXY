###############################################################################################
###############################################################################################
######  Function for bbp range test ###########
###############################################################################################
###############################################################################################


# Written by M.Cornec (2018) adapted from BBP doc :
#Catherine SCHMECHTIG, Emmanuel BOSS, Nathan BRIGGS, Hervé CLAUSTRE, Giorgio
#DALL'OLMO, Antoine POTEAU (2017). BGC Argo quality control manual for particles
#backscattering
#https://doi.org/10.13155/60262

### INPUTS :
# bbp : vector of the bbp measurements

### OUTPUTS :
# bbp_range: vector of the bbp measurements
# bbp_qc : QC adjusted values after range test

Range_test_bbp<-function(bbp) {
  # Initialize QC and bbp_range default
  bbp_range<-bbp
  bbp_qc<-rep(1,length(bbp))
  
  #Put out of range values to NA and corresponding qc=3 (probably bad data)
  bbp_range[which(bbp > 0.1 | bbp < -2.5*10^-5)]<-NA
  bbp_qc[which(bbp > 0.1 | bbp < -2.5*10^-5)]<-3
  
  return(list("bbp_range" = bbp_range, "bbp_qc" = bbp_qc))
  
}