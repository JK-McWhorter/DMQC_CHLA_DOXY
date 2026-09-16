###############################################################################################
###############################################################################################
######  Function for Negative Spikes Removal ###########
###############################################################################################
###############################################################################################


# Written by M.Cornec (2018) adapted from CHLA doc :
# Schmechtig Catherine, Claustre Herve, Poteau Antoine, D'Ortenzio Fabrizio (2014). Bio-Argo quality
# control manual for the Chlorophyll-A concentration. http://doi.org/10.13155/35385
# Catherine SCHMECHTIG, Emmanuel BOSS, Nathan BRIGGS, Hervé CLAUSTRE, Giorgio
# DALL'OLMO, Antoine POTEAU (2017). BGC Argo quality control manual for particles
# backscattering. https://doi.org/10.13155/60262

### INPUTS :
# param : vector of the param measurements

### REQUIREMENT:
# Need RunningFilter function

### OUTPUTS :
# param_negsp: vector of the param measurements
# param_qc : QC adjusted values after despike correction

Neg_spike<-function(param) {
  # Initialize QC and param_negsp default
  param_negsp<-param
  param_qc<-rep(1,length(param))
  
  # Calculate the median filtered param
  param_filt<-RunningFilter(2,param,ends.fill=T, Method="Median")
 
  # Calculate the residual
  param_res<-param-param_filt
  
  # calculate the 2*10th percentile of the residuals
  crit<-2*quantile(param_res,probs=0.1,na.rm=T)
  
  #Put negative spike params to NA and corresponding qc=4 (bad data)
  param_negsp[which(param_res<crit)]<-NA
  param_qc[which(param_res<crit)]<-4

  return(list("param_negsp" = param_negsp, "param_qc" = param_qc))
  
}