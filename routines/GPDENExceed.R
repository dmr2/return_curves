GPDENExceed <- function( z, lambda, MHHW, scale, shape, maxinterpN){
  
  # Function "GPDENExceed"
  #
  # Calculate log of the number of exceedances of z from a
  # Poisson-Generalized Pareto Distribution with Poisson mean
  # lambda and the specified shape and scale factor. Assumes
  # the exceedance threshold has already been removed from z.
  #
  # This function ensures that values returned stay within the
  # support of the GPD.
  #
  # For values of z below zero, the function will treat as though
  # z = 0 unless MHHW is specified. If MHHW is specified
  # value, exceedances below zero will be assumed to fall on a
  # Gumbel distribution between lambda exceedances at zero and a
  # specified value at z = MHHW(1) < 0. If MHHW(2) exists, it is
  # the value of exceedances at z = MHHW(1); otherwise, will default
  # to 365.25/2.
  
  # Reference:
  # M.K. Buchanan, R.E. Kopp, M. Oppenheimer, and C. Tebaldi. (2016).
  # Allowances for evolving coastal flood risk under uncertain local
  # sea-level rise. Climatic Change.
  
  # History
  # 8/11/2017 (DMR): Function created
  #
  
  eps <- .Machine$double.eps # Machine epsilon

  if(length(shape) == 1){ # Don't sample GPD uncertainty 
    
    result <- GPDNExceed(z, lambda, MHHW, scale, shape)
    
  }else{
    
    minz <- .99999*abs(scale/shape) # (1000)
    
    # interpolate to handle z of different sizes
    testz <- unique(sort(c(0, as.vector(z)))) 
    if( length(testz) > maxinterpN ){
      testz <- seq(testz[1]-eps,testz[length(testz)]+eps,length=maxinterpN)
    }else if( length(testz) == 1 ){
      testz <- c(testz-eps, testz+eps)
    }
  
    testz0 <- testz
    a <- mapply(FUN=function(x)(shape>=0)*x,testz) 
    b <- matrix(mapply(FUN=function(x,y)min(x,y),mapply(rep,times=1000,testz),minz),ncol=1000) 
    testz <- a + b*(shape<0) 
    
    a2 <- (shape+eps)/scale 
    a3 <- matrix(pmax(0,testz),nrow=1000,ncol=1000)
    a <- 1 + (a3*a2) 
    b <- -1/(shape+eps)
    logNref <- log(lambda*mapply(FUN=function(x,y) (abs(x)^y)*sign(x),a,b)) 
    logENref = log(apply(exp(matrix(logNref,nrow=1000,ncol=1000)),2,mean))

    f <- approxfun(testz0,logENref)
    result <- f(pmax(0,zminusSLR))
    
    # Put flood heights below 99.9th threshold and greater than MHHW on a Gumbel distribution
    if( length(MHHW) >= 1 ){
      z <- pmax(MHHW[1],zminusSLR)
      if( length( MHHW ) >= 2 ){
        exceedMHHW <- MHHW[2]
      }else{
        exceedMHHW = 365.25/2
      }
      
      result = (z<0)*(log(lambda)+(log(exceedMHHW)-log(lambda))*z/MHHW[1])+(z>=0)*result
      return( result )
    }
  }
 
  
 }