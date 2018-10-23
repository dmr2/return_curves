#!/opt/local/bin/R

# Inputs:
#   - Monte Carlo SL samples for each tide gauge (e.g., from LocalizeSL and ProjectSL programs)
#   - GPD paramters
#   - Text file with heights of observed flood events
#
# Outputs:
#   - Flood return frequency curves for each tide gauge
#      * Historical and Projected Future flood return curves
#

rm(list=ls(all=TRUE))
setwd("/Users/dmr/Dropbox/IPCC\ Sea\ Level/allowance")

source("routines/GPDsample.R")
source("routines/GPDNExceed.R")
source("routines/openSLR.R")
source("routines/GPDENExceed.R")
source("routines/histStorms.R")
source("routines/plotGPDNexceed.R")

# Climate scenario
scenarios = c("1.5C","2.0C","2.5C")
slab = c("0p8degree","1p3degree","1p8degree")

# Location of local sea level rise Monte Carlo samples
dir = "/Users/dmr/Dropbox/IPCC\ Sea\ Level/slr_samples"

# Specify the years of interest
targ_years = c(2050,2100)

scenarios = c("1.5C","2.0C","2.5C")
slab = c("1p5","2p0","2p5")


# Open GPD parameters and historical flood data
dat <- read.csv("GPDfits_uhawaii_projectLSLfmt.tsv",header=T,sep="\t")
sites <- dat$Site # tide gauge sites

for(j in 1:length(sites)){
#for(j in which(sites=="Rabaul"):which(sites=="Rabaul")){ # Only process "The Battery, NYC"

  site <- dat$Site[j] # tide gauge site name
  scale <- dat$Scale50[j] # median scale parameter
  shape <- dat$Shape50[j] # median shape parameter
  UHid <- dat$UHAWAII_ID[j]
  threshold <- dat$Q99[j] # GPD threshold
  lambda <- dat$Lambda[j] # mean Poisson arrival rate of threshold
  shapescaleV <- dat$Vscaleshape[j] # covariance of GPD scale and shape parameters
  shapeV <- dat$Vshape[j] # variance of shape parameter
  scaleV <- dat$Vscale[j] # variance of scale parameter
  gauge <- dat$PSMSL_ID[j] # Tide gauge ID
  basin <- dat$Basin[j]
  
  # add something that gets start and end year for tide gauge data
  
 # Account for GPD parameter uncertainty by making draws from a
 # bivariate normal distribution using Latin hypercube sampling
  GPD <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)

  z <- seq(0,10,.01) # some flood heights (meters above tide gauge MHHW)
  
  # Expected historical flood height curve (No SLR) (GPD uncertainty)
  qqq <- matrix(NaN, nrow=length(GPD[,2]), ncol=length(z))
  for(iii in 1:length(GPD[,2]) ){
    qqq[iii,] <- GPDNExceed(z-threshold,lambda,-threshold,GPD[iii,2],GPD[iii,1])
  }
  
  # GPD confidence intervals
  gpdCI <- apply(qqq,2,quantile,probs=c(0.1666,0.5,0.8333),na.rm=T)

  name <- c(rep("q16",length(z)),rep("q50",length(z)),rep("q83",length(z)))
  gpdcurve <- data.frame(freq=c(gpdCI[1,],gpdCI[2,],gpdCI[3,]),name=name,
    height=rep(z,3)) # 50th & 17/83 confidence intervals
    
  
  Ne_hist <- apply(qqq,2,mean,na.rm=T) 
  
  # Future Sea Level
  curve50SLR <- array(NaN,dim=c(length(targ_years),length(scenarios),length(z))) # 50th percentile SLR
  curveEff <- array(NaN,dim=c(length(targ_years),length(scenarios),length(z))) # Effective SLR
  
  for( s in 1:length(scenarios) ){
    
    # Get sea level rise Monte Carlo Samples
    fil <- paste( dir,"/LSLproj_MC_",gauge,"_",slab[s],"degree.tsv",sep="")
    if (!file.exists(fil)){
      fil <- paste( dir,"/LSLprojMC_",gauge,"_",slab[s],"degree.tsv",sep="")
    }
    SLR <- getRSLsamps( fil )
    
    SLRMC <- SLR$samples
    years <- SLR$years 
    
    for( t in 1:length(targ_years) ){
      
      indx = which(years==targ_years[t])
    
    # Future flood curves that include GPD parameter and SLR uncertainty
      zminusSLR <- mapply(function(x) x - SLRMC[,indx], z)
    
      ans <- GPDENExceed(zminusSLR-threshold,lambda,-threshold,GPD[,2],GPD[,1],1000)
      curveEff[t,s,] <- apply(matrix(exp(ans),nrow=10000,ncol=1001),2,mean)
      
    # Shift historic flood return curves
    # Median SLR
      slr_q50 <- quantile(SLRMC[,indx],probs=c(.5))
      qqq <- matrix(NaN, nrow=length(GPD[,2]), ncol=length(z))
      for(iii in 1:length(GPD[,2]) ){
        qqq[iii,] <- GPDNExceed(z-slr_q50-threshold,lambda,-threshold,GPD[iii,2],GPD[iii,1])
      }
      curve50SLR[t,s,]  <- apply(qqq,2,mean,na.rm=T)
    
    } # each target year
  } # each scenario
  
  
  # Plot flood curves
  # Get historical flood event data 
  fil <- paste("/Users/dmr/Dropbox/IPCC\ Sea\ Level/GPDfit/obs_q99exceed_",basin,"_gauge_",UHid,".csv",sep="")
  obs <- histStorms(fil)
  
  histcurve <- data.frame(height=z,N=Ne_hist,loc=site,name)

  # Put data into data frame for plotting...
          
  for( t in 1:length(targ_years) ){
    
    freq <- c(Ne_hist,gpdCI[1,],gpdCI[2,],gpdCI[3,],
              curve50SLR[t,1,],curve50SLR[t,2,],curve50SLR[t,3,],
              curveEff[t,1,],curveEff[t,2,],curveEff[t,3,])   
                      
    freq[is.na(freq)] <- 10e-10
  
    name <- c(rep("Historic Flood Return Curve",length(z)),
              rep("q16",length(z)),rep("q50",length(z)),rep("q83",length(z)),
              rep("N+SL50 1.5°C",length(z)),rep("N+SL50 2.0°C",length(z)),
              rep("N+SL50 2.5°C",length(z)),
              rep("Ne 1.5°C",length(z)),
              rep("Ne 2.0°C",length(z)),
              rep("Ne 2.5°C",length(z)))

    df <- data.frame(height=z,freq=freq,name=name,loc=site,
                     basin=basin,year=targ_years[t],id=UHid)
                   
    plotGPDNExceed(df, obs)
    insert_minor

  } # each target year
}# each tide gauge
