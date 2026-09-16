############## MLD CALCULATION
MLD_calc <- function(sigma,dep_sigma) { 
  MLD<-NA
  if(length(sigma[!is.na(sigma)==TRUE])>=2) {
    sigmaSurface<-NA
    sigmaSurface <- approx(dep_sigma,sigma,10, ties=min)$y # identification de la densit? ? 10m ("surface", pour variations journali?res)
    if (is.na(sigmaSurface)==FALSE) {
      MLD <- max(dep_sigma[sigma <= (sigmaSurface + 0.03)],na.rm = T)
    }
  }
  return(MLD)
}