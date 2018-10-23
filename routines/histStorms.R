histStorms <- function( fil, threshold ){
  
  # From observations of daily max flood heights (relative to MHHW),
  # determine the frequency of the historical flood events.
  z <- seq(0,10,.01) # some flood heights (meters above tide gauge MHHW)
  
  dat <- read.csv(fil,header=F,sep="\t")
  gt <- function(x,y) length(which(x>y))/(length(x)/365.25)
  ct <- NA
  for(iii in 1:length(z) ){
    ct[iii] <- gt(dat$V1,z[iii])
  }
  
  subct <- which(diff(ct)<0)
  subct <- intersect(subct,which(z>threshold))
  
  df <- data.frame(z=z[subct], freq=ct[subct], group="Observed")
  
  return ( df )
}