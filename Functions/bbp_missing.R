###############################################################################################
###############################################################################################
######  Function for Missing data for BBP values test ###########
###############################################################################################
###############################################################################################


# Transcripted by M.Cornec (2023) from BBP doc :
# Giorgio DALL'OLMO et al (2017). Real-time quality control of optical 
# backscattering data from BGC-Argo floats
# https://doi.org/10.12688/openreseurope.15047.1

# "To detect and flag profiles that have a large fraction of missing data. 
# Missing data could indicate shallow profiles (caused by a specific float 
# mission and/or bathymetry) or incomplete profiles due to a malfunctioning
# sensor."


### INPUTS :
# bbp : vector of the bbp measurements
# dep_bbp : vector of the depths of bbp measurements

### OUTPUTS :
# bbp_qc : vector of values of QC for the whole profile
# 1 : no mussing data
# 3 : some missing data
# 4 : only data in one layer
# 9 : no data

bbp_missing<-function(bbp,dep_bbp) {
  
  # initiate default bbp QC (1 = data)
  bbp_qc<-rep(1,length(bbp))
  
  # determine bins of depth layers
  bins<-c(50 ,156, 261, 367, 472, 578, 683, 789, 894, 1000)
  
  # initiate empty summary dtfr
  data_bin<-NULL
  
  #initiate first bin lim
  bin_lim<-0
  
  # iterate test over each bin
  for(i in bins){
    
    # intiate bin value default : "yes" there is data
    bin<-"yes"
    
    # if only NA values or no value in the bin --> bin = "no" data
    if(all(is.na(bbp[which(dep_bbp < i &
                           dep_bbp >= bin_lim)]))){
      bin<-"no"
    }
    
    # compile the information of the current bin to the summary
    data_bin<-c(data_bin,bin)
    
    # change the limit depth for the deeper bin
    bin_lim<-i
  }
  
  # Attribute QC to the profile
  
  # Check if some bins contain no data
  if(length(which(data_bin=="no"))!=0){
    # more than one bin has dat:QC 3 : need further DM
    if(length(which(data_bin=="yes"))>1) {
      bbp_qc<-rep(3,length(bbp))
    } 
    # only one bin has data:QC 4 : sensor malfunction
    if(length(which(data_bin=="yes"))==1) {
      bbp_qc<-rep(4,length(bbp))
    } 
    # no data :QC 9 : sensor malfunction
    if(length(which(data_bin=="yes"))==0) {
      bbp_qc<-rep(9,length(bbp))
    } 
  }

  return(bbp_qc)
  
}