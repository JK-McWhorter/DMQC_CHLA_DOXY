###############################################################################################
######  Function for NPQ correction, Xing 2012 (FIXED) #######################################
###############################################################################################

NPQ_MLD <- function (chla,dep_chla,MLD) {
  
  #04/01/19
  MLD<-0.9*MLD
  
  #vertical resolution of the Fchl profile (above 250m)
  RESO<-NA
  RESO<-mean(abs(diff(dep_chla[which(dep_chla<=250)])),na.rm=T) 
  
  # SAFE FALLBACK: If RESO is NA, NaN, or 0, set it to a default of 1
  if (is.na(RESO) || is.nan(RESO) || RESO <= 0) {
    RESO <- 1 
  }
  
  #Declare variable
  chla_npq<-chla
  
  if (length(chla[!is.na(chla)])>5) { # at least 5 no-NA values needed in Fchl profile
    
    #declare empty variables
    chla_interp_res<-NA
    chla_filt_res<-NA
    
    #create interpolated depth profile (function of vert RESO)
    dep_reg_res<-NA
    dep_reg_res <-seq(0,500,RESO)
    chla_interp_res = approx(dep_chla,chla,dep_reg_res)$y
    
    # Filter to remove spike (according to vertical resolution)
    if (RESO < 3) {
      chla_filt_res = RunningFilter(2,chla_interp_res, Method="Median")
    }
    
    if (RESO < 1) {
      chla_filt_res = RunningFilter(5,chla_interp_res, Method="Median")
    }
    
    if (RESO > 3) {
      chla_filt_res = RunningFilter(1,chla_interp_res, Method="Median")
    }
    
    # Identify max Fchla depth
    max_Fchla_depth<-NA
    max_Fchla_depth<-median(dep_reg_res[which(chla_filt_res==max(chla_filt_res,na.rm=T))],na.rm=T)
    
    if (is.na(MLD)==F) { 
      
      chla_max_npq<-NA
      chla_max_npq<-max(chla_filt_res[which(dep_reg_res<=MLD)], na.rm=T)
      
      #identify the depth related to the max Fchla
      chla_max_npq_depth<-NA
      chla_max_npq_depth<-dep_reg_res[which(chla_filt_res==chla_max_npq & dep_reg_res<=MLD)]
      chla_max_npq_depth<-max(chla_max_npq_depth,na.rm=T)
      
      # correct the profile with the max Fchla
      if (chla_max_npq!=-Inf) {
        chla_npq[which(dep_chla <= chla_max_npq_depth )]<-chla_max_npq
      }
    }
    
  }
  return(chla_npq)
}


###############################################################################################
######  Function for NPQ correction, Xing 2018 (FIXED) #######################################
###############################################################################################

NPQ_cor_X12_XB18 <- function (chl,dep_chl,dep_light,light,MLD) {
  
  #04/01/19
  MLD<-0.9*MLD
  
  #vertical resolution of the Fchl profile (above 250m)
  RESO<-NA
  RESO<-mean(abs(diff(dep_chl[which(dep_chl<=250)])),na.rm=T) 
  
  # SAFE FALLBACK: If RESO is NA, NaN, or 0, set it to a default of 1
  if (is.na(RESO) || is.nan(RESO) || RESO <= 0) {
    RESO <- 1 
  }
  
  #Declare variable
  chl_npq<-chl
  
  if (length(chl[!is.na(chl)])>5) { # at least 5 no-NA values needed in Fchl profile
    
    par_15_depth<-NA
    par_15_depth<-max(dep_light[which(light >= 15)],na.rm = T) # identify max depth where PAR is above 15E
    
    #declare empty variables
    chl_interp_res<-NA
    chl_filt_res<-NA
    
    #create interpolated depth profile (fonction of vert RESO)
    dep_reg_res<-NA
    dep_reg_res <-seq(0,500,RESO)
    chl_interp_res = approx(dep_chl,chl,dep_reg_res)$y
    
    # Filter to remove spike (according to vertical resolution)
    if (RESO < 3) {
      chl_filt_res = RunningFilter(5,chl_interp_res, Method="Median")
    }
    
    if (RESO > 3) {
      chl_filt_res = RunningFilter(1,chl_interp_res, Method="Median")
    }
    
    # Identify max Fchl depth
    max_Fchl_depth<-NA
    max_Fchl_depth<-median(dep_reg_res[which(chl_filt_res==max(chl_filt_res,na.rm=T))],na.rm=T)
    
    if (par_15_depth!=Inf & par_15_depth!=-Inf & is.na(MLD)==F) { # if PAR_15 is not beyond/above bound depths
      
      if (min(MLD,par_15_depth,na.rm=T)==MLD &
          max_Fchl_depth > MLD) {  #XB18 (case for stratified shallow mixing wrs w/o DCM)
        
        r=0.092
        iPARmid= 261
        e=2.2
        
        # interpolation of the PAR profile at the Fchl depth
        light_interp<-approx(dep_light,light,dep_chl)$y
        
        chl_npq<-chl*NA
        
        #04/02  Modifiy xb to put values where light is na
        chl_npq<-chl/(r+(1-r)/(1+(light_interp/iPARmid)^e))
        
        chl_npq[which(dep_chl<10)]<-
          chl_npq[which(dep_chl==min(dep_chl[which(dep_chl >= 10)], na.rm=T))]
        
        chl_npq[which(is.na(chl_npq) & dep_chl >=10)]<-
          chl[which(is.na(chl_npq) & dep_chl >=10)]
        
      } else {  #X12 (case for w/o DCM, and DCM with deep mixing)
        
        #identify minimum depth between MLD and PAR 15 depth
        min_depth_npq<-NA
        min_depth_npq<-min(MLD,par_15_depth,na.rm=T)
        
        # identify the value of max Fchla bewteen the surface and the min depth 
        chl_max_npq<-NA
        chl_max_npq<-max(chl_filt_res[which(dep_reg_res<=min_depth_npq)], na.rm=T)
        
        #identify the depth related to the max Fchla
        chl_max_npq_depth<-NA
        chl_max_npq_depth<-dep_reg_res[which(chl_filt_res==chl_max_npq & dep_reg_res<=min_depth_npq)]
        chl_max_npq_depth<-max(chl_max_npq_depth,na.rm=T)
        
        # correct the profile with the max Fchla
        if (chl_max_npq!=-Inf) {
          chl_npq[which(dep_chl <= chl_max_npq_depth )]<-chl_max_npq
        }
      }
    }
  }
  return(chl_npq)
}