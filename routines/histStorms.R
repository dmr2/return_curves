histStorms <- function( fil, threshold, group ){
  
  # From observations of daily max flood heights (relative to MHHW),
  # determine the frequency of the historical flood events.
  
  dat <- read.csv(fil,header=T,sep=",")
  gt <- function(x,y) length(which(x>y))/round(length(unique(dat$Year)))
  ct <- NA
  for(iii in 1:length(z) ){
    ct[iii] <- gt(dat$Height,z[iii])
  }
  
  subct <- which(diff(ct)<0)
  subct <- intersect(subct,which(z>threshold))
  
  df <- data.frame(z=z[subct], freq=ct[subct], group=group)
  
  return ( df )
}