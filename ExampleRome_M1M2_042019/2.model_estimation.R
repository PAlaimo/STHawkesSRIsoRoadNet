
##########################
### MODEL ESTIMATION #####
##########################
# Dr. Marco Mingione and Dr. Pierfrancesco Alaimo Di Loro, 20/11/2023

# This code estimates a nonparametric spatio-temporal periodic Hawkes process for car accidents on a road network.
# Example: Rome, municipalities I and II, April 2019

# Loading Packages --------------------------------------------------------
require(tidyverse)
require(magrittr)
require(lubridate)
require(MASS)
require(mvtnorm)
require(raster)
require(sf)
require(Matrix)

Rcpp::sourceCpp("auxiliary.cpp")


# Auxiliary ---------------------------------------------------------------

NegLogLikehood <- function(x, bENoMu, bANoMu, tLNoA, trig.Int.NoA, inside, Xmat, Xmat_all, secInside, phiin, rhoin){
  mulik <- exp(x[1]) # Contraining positivity
  betalik <- as.matrix(x[2:(ncol(Xmat)+1)])
  
  # Weighted logs on the events
  slbE <- sum(phiin[inside]*log(mulik * bENoMu))
  sltE <- sum(rhoin[secInside]*log(exp(Xmat[secInside,]%*%betalik)*tLNoA))
  sloglambdaEv <- slbE + sltE
  
  # Integral on on the space-time
  tA <- c(exp(Xmat_all%*%betalik)*trig.Int.NoA)
  lambdaAll <-  mulik* bANoMu + sum(tA)
  
  # Output
  return(-sloglambdaEv + lambdaAll)
}


# Load data and pre-processing --------------------------------------------

load("WS/Real_PrepToEst_42019.RData")
load("WS/Real_EstPre_42019.RData")

# Initialization ----------------------------------------------------------
# Setting initial mu and A
mu0 <- 1
beta0 <- c(-2, rep(0,ncol(Xmat)-1))

# Setting initial backgroung
# Space base
value.background <- rep(1, nrow(grid.background))

# Trend
value.trend <- base.trend*0+1
fun.trend <- approxfun(base.trend, value.trend, yleft=0, yright=0)

# Daily
value.daily <- base.daily*0+1
fun.daily <- approxfun(base.daily, value.daily, yleft=0, yright=0)

# weekly periodicity
value.weekly <- base.weekly*0+1
fun.weekly <- approxfun(base.weekly, value.weekly, yleft=0, yright=0)


# Setting initial excitation
# Temporal
value.tempRes <- dnorm(base.tempRes, 0, 0.05)
value.tempRes <- value.tempRes / sum(value.tempRes*step.tempRes) # Make integrate to 1
fun.tempRes <- approxfun(base.tempRes, value.tempRes, yleft = 0, yright=0)

# Spatial
value.spatRes <- dmvnorm(cbind(base.spatRes, rep(0, length(base.spatRes))), mean = c(0,0), sigma = diag(.01, 2))
value.spatRes <- value.spatRes/sum(2*pi*base.spatRes*value.spatRes*step.spatRes)
fun.spatRes <- approxfun(base.spatRes, value.spatRes, yleft=0, yright=0)

# EXPECTATION STEP --------------------------------------------------------

# Weights computation -----------------------------------------------------

fun_back <- value.background[dataTogrid]
back.Events.NoMu <- (fun_back*fun.trend(dat$t)*fun.daily(dat$timeday)*fun.weekly(dat$timeweek))# background importance of each event
# of the trend weight
trig.Lags.NoA <- (fun.tempRes(lags$timeLags)*fun.spatRes(lags$dxy)) # triggering importance of each pair

trig.Events <- tapply(c(exp(Xmat%*%beta0)*trig.Lags.NoA, 0*(1:nrow(dat))), c(lags$idSecond, 1:nrow(dat)), sum) # triggering importance on each event
lambda.events <- (mu0*back.Events.NoMu+trig.Events) # intensity on each event

# Weights
weight.trend <- c((fun_back*fun.trend(dat$t))/lambda.events)
weight.daily <- c((fun_back*fun.daily(dat$timeday))/lambda.events)
weight.weekly <- c((fun_back*fun.weekly(dat$timeweek))/lambda.events)
phi <- c(mu0*back.Events.NoMu/lambda.events)
rho <- c(c(exp(Xmat%*%beta0)*trig.Lags.NoA)/lambda.events[lags$idSecond]) # Putting the "triggered" lambda at the denominator

# Smoothing trend and background-------------------------------------------------------

# Trend smoothing
smoothed.trend <- Rfast::colsums(kernelWithEdgeCorr.trend*weight.trend)
buf_trend_smooth <- base.trend>=0 & base.trend<=ceiling(max(dat$t[!dat$InBuf]))
value.trend <- smoothed.trend/mean(smoothed.trend[buf_trend_smooth]) # Standardization
fun.trend <- approxfun(base.trend, value.trend, yleft=0, yright=0)

rm(smoothed.trend)

# Daily smoothing
smoothed.daily <- Rfast::colsums(kernel.dailyMat*weight.daily[!dat$InBuf])
value.daily <- smoothed.daily/mean(smoothed.daily) # Standardization
fun.daily <- approxfun(base.daily, value.daily, yleft=0, yright=0)

rm(smoothed.daily)

# Weekly smoothing
smoothed.weekly <- Rfast::colsums(kernel.weeklyMat*weight.weekly[!dat$InBuf])
value.weekly <- smoothed.weekly/mean(smoothed.weekly) # Standardization
fun.weekly <- approxfun(base.weekly, value.weekly, yleft=0, yright=0)

rm(smoothed.weekly)

# Background
gridbaseid <- cbind(match(grid.background[,1], unique(base.x)), 
                    match(grid.background[,2], unique(base.y)))-1
###### Smoothed background ######
if(sparse == 1){
  tictoc::tic("Sparse")
  smoothed.background <- smooth_background_sparse(kernelWithEdgeCorr.backgroundx, kernelWithEdgeCorr.backgroundy, phi, 
                                                  gridbaseid, nrow(gridbaseid))
  tictoc::toc()
}else{
  tictoc::tic("Non Sparse")
  smoothed.background <- smooth_background(kernelWithEdgeCorr.backgroundx, kernelWithEdgeCorr.backgroundy, phi, 
                                           gridbaseid, nrow(gridbaseid))
  tictoc::toc()
}

# Standardization only over the polygon
value.background <- smoothed.background/mean(smoothed.background[grid.Inpoly], na.rm = T)
rm(smoothed.background)

# Smoothing the excitation---------------------------------------------------------------

# Temporal excitation smoothing
smoothed.tempRes <- Rfast::colsums(kernelWithEdgeCorr.tempRes*rho[!lags$InBuf])/repTempTrig
value.tempRes <- smoothed.tempRes/sum(smoothed.tempRes*step.tempRes) # Normalizing
fun.tempRes <- approxfun(base.tempRes, value.tempRes, yleft=0, yright=0)

rm(smoothed.tempRes)

# Spatial excitation smoothing
smoothed.spatRes <- Rfast::colsums(kernelWithEdgeCorr.spatRes*rho[!lags$InBuf])/repSpatTrig
smoothed.spatRes[1] <- smoothed.spatRes[2]
# Making integrate to 1
value.spatRes <- smoothed.spatRes/sum(2*pi*base.spatRes*smoothed.spatRes*step.spatRes)
fun.spatRes <- approxfun(base.spatRes, value.spatRes, yleft=0, yright=0)
rm(smoothed.spatRes)

# Updating weights
fun_back <- value.background[dataTogrid]
# Numerator of the trend weight
back.Events.NoMu <- (fun_back*fun.trend(dat$t)*fun.daily(dat$timeday)*fun.weekly(dat$timeweek))
# Numerator of the pair triggering weights. 
trig.Lags.NoA <- (fun.tempRes(lags$timeLags)*fun.spatRes(lags$dxy)) 
# Second part of lambda
trig.Events <- tapply(c(exp(Xmat%*%beta0)*trig.Lags.NoA, 0*(1:nrow(dat))), c(lags$idSecond, 1:nrow(dat)), sum) 
lambda.events <- (mu0*back.Events.NoMu+trig.Events) # Lambda
phi <- c(mu0*back.Events.NoMu/lambda.events)
rho <- c(c(exp(Xmat%*%beta0)*trig.Lags.NoA)/lambda.events[lags$idSecond]) # Putting the "triggered" lambda at the denominator

# MAXIMIZATION STEP -------------------------------------------------------

# Lambda on the events: first term of the loglikelihood -------------------

# Background on the events
fun_back <- value.background[dataTogrid][!dat$InBuf]
back.Events.NoMu <- (fun_back*fun.trend(dat$t[!dat$InBuf])*fun.daily(dat$timeday[!dat$InBuf])*fun.weekly(dat$timeweek[!dat$InBuf]))

# Computation of the triggering on each pair
trig.Lags.NoA <- (fun.tempRes(lags$timeLags[!lags$InBuf])*fun.spatRes(lags$dxy[!lags$InBuf]))

# Lambda integrated: second term of log-likelihood-------------------
# Integration of the background

tempbgseq <- seq(0, max(base.trend[buf_trend_smooth]), by=min(step.daily, step.weekly, step.trend))
back.All.NoMu <- sum(value.background[grid.Inpoly]*step.x*step.y, na.rm = T)*
  sum(fun.trend(tempbgseq)*fun.daily(tempbgseq-floor(tempbgseq))*fun.weekly(tempbgseq%%7)*step.daily)

# Set up where to store integrals of excitation
totalEffects <- data.frame(eventID=1:nrow(dat))

# Temporal triggering
# Building an approximation to the integral of temporal response from 0 to each t in [0,TT]
tempk <- cumsum(fun.tempRes(base.tempRes)*step.tempRes)
fun.intTempRes <- approxfun(base.tempRes, tempk, yleft=0, yright=1)
rm(tempk)
# Computing the integral from 0 to T-ti (equivalent to ti to T, given that we start in 0)
totalEffects$intTempRes <- fun.intTempRes(ceiling(max(base.trend[buf_trend_smooth]))-dat$t)-
  fun.intTempRes(floor(min(base.trend[buf_trend_smooth]))-dat$t)

# Integration of the spatial triggering
spatmids <- c(fun.spatRes(base.spatRes[1]), fun.spatRes(base.spatRes[-length(base.spatRes)])+diff(fun.spatRes(base.spatRes))/2)
# Integration of spatial response for each event
totalEffects$intSpatRes <- intPointsSpaTrig%*%spatmids
# Combining temporal and spatial triggering
totalEffects$intRes <- totalEffects$intTempRes*totalEffects$intSpatRes

# Likelihood optimization --------------------------------------------------------------

# Initial mu on log-scale
lmu0 <- log(mu0)

######## Likelihood Max ########

# Start from previous value
res.optim1 <- optim(par=c(lmu0, beta0), NegLogLikehood, 
                    bENoMu = back.Events.NoMu, 
                    bANoMu = back.All.NoMu, 
                    tLNoA = trig.Lags.NoA,
                    trig.Int.NoA = totalEffects$intRes,
                    Xmat = Xmat, Xmat_all = Xmat_all,
                    inside = !dat$InBuf, secInside = !lags$InBuf, 
                    phiin = phi, rhoin = rho,
                    control=list(trace=0), hessian = T, method = "BFGS")
# Random start
res.optim2 <- optim(par=c(log(0.5), rnorm(ncol(Xmat), beta0, 2)), NegLogLikehood, 
                    bENoMu = back.Events.NoMu, 
                    bANoMu = back.All.NoMu, 
                    tLNoA = trig.Lags.NoA,
                    trig.Int.NoA = totalEffects$intRes,
                    Xmat=Xmat, Xmat_all = Xmat_all,
                    inside = !dat$InBuf, secInside = !lags$InBuf, 
                    phiin = phi, rhoin = rho,
                    control=list(trace=0))
if(res.optim1$value <= res.optim2$value){
  res.optim <- res.optim1
}else{
  res.optim <- res.optim2
}
mustar <- c()
betastar <- list()
mustar[1] <- exp(res.optim$par[1]) # The squared values were the ones to produce the best
betastar[[1]] <- res.optim$par[2:(ncol(Xmat)+1)]  # likelihood

# Looping -----------------------------------------------------------------

M <- 50
print(paste("Loops=", M))
ll <- rep(-Inf, M+1)
ll[1] <- res.optim$value

for(i in 2:M){
  print(paste("Start loop",i))
  # Weights computation -----------------------------------------------------
  
  fun_back <- value.background[dataTogrid]
  # Numerator of the trend weight
  back.Events.NoMu <- (fun_back*fun.trend(dat$t)*fun.daily(dat$timeday)*fun.weekly(dat$timeweek))
  # Numerator of the pair triggering weights
  trig.Lags.NoA <- (fun.tempRes(lags$timeLags)*fun.spatRes(lags$dxy)) 
  # Second part of lambda
  trig.Events <- tapply(c(exp(Xmat%*%betastar[[i-1]])*trig.Lags.NoA, 0*(1:nrow(dat))), c(lags$idSecond, 1:nrow(dat)), sum) 
  lambda.events <- (mustar[i-1]*back.Events.NoMu+trig.Events) # Lambda
  
  weight.trend <- c((fun_back*fun.trend(dat$t))/lambda.events)
  weight.daily <- c((fun_back*fun.daily(dat$timeday))/lambda.events)
  weight.weekly <- c((fun_back*fun.weekly(dat$timeweek))/lambda.events)
  phi <- c(mustar[i-1]*back.Events.NoMu/lambda.events)
  # Putting the "triggered" lambda at the denominator
  rho <- c(c(exp(Xmat%*%betastar[[i-1]])*trig.Lags.NoA)/lambda.events[lags$idSecond]) 
  
  
  # Smoothing ---------------------------------------------------------------
  
  # Smoothing trend and background-------------------------------------------------------
  
  # Trend smoothing
  smoothed.trend <- Rfast::colsums(kernelWithEdgeCorr.trend*weight.trend)
  value.trend <- smoothed.trend/mean(smoothed.trend[buf_trend_smooth]) # Standardization
  fun.trend <- approxfun(base.trend, value.trend, yleft=0, yright=0)
  
  rm(smoothed.trend)
  
  # Daily smoothing
  smoothed.daily <- Rfast::colsums(kernel.dailyMat*weight.daily[!dat$InBuf])
  value.daily <- smoothed.daily/mean(smoothed.daily) # Standardization
  fun.daily <- approxfun(base.daily, value.daily, yleft=0, yright=0)
  
  rm(smoothed.daily)
  
  # Weekly smoothing
  smoothed.weekly <- Rfast::colsums(kernel.weeklyMat*weight.weekly[!dat$InBuf])
  value.weekly <- smoothed.weekly/mean(smoothed.weekly) # Standardization
  fun.weekly <- approxfun(base.weekly, value.weekly, yleft=0, yright=0)
  
  rm(smoothed.weekly)
  
  # Background
  if(sparse == 1){
    tictoc::tic("Sparse")
    smoothed.background <- smooth_background_sparse(kernelWithEdgeCorr.backgroundx, kernelWithEdgeCorr.backgroundy, phi, 
                                                    gridbaseid, nrow(gridbaseid))
    tictoc::toc()
  }else{
    tictoc::tic("Non Sparse")
    smoothed.background <- smooth_background(kernelWithEdgeCorr.backgroundx, kernelWithEdgeCorr.backgroundy, phi, 
                                             gridbaseid, nrow(gridbaseid))
    tictoc::toc()
  }
  # Standardization only over the polygon
  value.background <- smoothed.background/mean(smoothed.background[grid.Inpoly], na.rm = T)
  rm(smoothed.background)
  
  # Temporal excitation smoothing
  smoothed.tempRes <- Rfast::colsums(kernelWithEdgeCorr.tempRes*rho[!lags$InBuf])/repTempTrig
  value.tempRes <- smoothed.tempRes/sum(smoothed.tempRes*step.tempRes)
  fun.tempRes <- approxfun(base.tempRes, value.tempRes, yleft=0, yright=0)
  
  rm(smoothed.tempRes)
  
  # Spatial excitation smoothings
  smoothed.spatRes <- Rfast::colsums(kernelWithEdgeCorr.spatRes*rho[!lags$InBuf])/repSpatTrig
  smoothed.spatRes[1] <- smoothed.spatRes[2]
  
  # Making integrate to 1
  value.spatRes <- smoothed.spatRes/sum(2*pi*base.spatRes*smoothed.spatRes*step.spatRes)
  fun.spatRes <- approxfun(base.spatRes, value.spatRes, yleft=0, yright=0)
  
  rm(smoothed.spatRes)
  
  # Updating weights
  fun_back <- value.background[dataTogrid]
  # Numerator of the trend weight
  back.Events.NoMu <- (fun_back*fun.trend(dat$t)*fun.daily(dat$timeday)*fun.weekly(dat$timeweek))
  # Numerator of the pair triggering weights
  trig.Lags.NoA <- (fun.tempRes(lags$timeLags)*fun.spatRes(lags$dxy)) 
  # Second part of lambda
  trig.Events <- tapply(c(exp(Xmat%*%betastar[[i-1]])*trig.Lags.NoA, 0*(1:nrow(dat))), c(lags$idSecond, 1:nrow(dat)), sum) 
  lambda.events <- (mustar[i-1]*back.Events.NoMu+trig.Events) # Lambda
  phi <- c(mustar[i-1]*back.Events.NoMu/lambda.events)
  rho <- c(c(exp(Xmat%*%betastar[[i-1]])*trig.Lags.NoA)/lambda.events[lags$idSecond])
  
  # Lambda on the events: first term of the loglikelihood -------------------
  
  # Background on the events
  fun_back <- value.background[dataTogrid][!dat$InBuf]
  back.Events.NoMu <- (fun_back*fun.trend(dat$t[!dat$InBuf])*fun.daily(dat$timeday[!dat$InBuf])*fun.weekly(dat$timeweek[!dat$InBuf]))
  
  # Computation of the triggering on each pair
  trig.Lags.NoA <- (fun.tempRes(lags$timeLags[!lags$InBuf])*fun.spatRes(lags$dxy[!lags$InBuf]))
  
  # Lambda integrated: second term of log-likelihood-------------------
  
  # Integration of the background
  back.All.NoMu <- sum(value.background[grid.Inpoly]*step.x*step.y, na.rm = T)*
    sum(fun.trend(tempbgseq)*fun.daily(tempbgseq-floor(tempbgseq))*fun.weekly(tempbgseq%%7)*step.daily)
  

  # # Integration of the triggering
  # Temporal triggering
  # Building an approximation to the integral of temporal response from 0 to each t in [0,TT]
  tempk <- cumsum(fun.tempRes(base.tempRes)*step.tempRes)
  fun.intTempRes <- approxfun(base.tempRes, tempk, yleft=0, yright=1)
  rm(tempk)
  # Computing the integral from 0 to T-ti (equivalent to ti to T, given that we start in 0)
  totalEffects$intTempRes <- fun.intTempRes(ceiling(max(base.trend[buf_trend_smooth]))-dat$t)-
    fun.intTempRes(floor(min(base.trend[buf_trend_smooth]))-dat$t)
  # Computing the integral from 0 to T-ti (equivalent to ti to T, given that we start in 0)
  # Integration of spatial response for each event
  spatmids <- c(fun.spatRes(base.spatRes[1]), fun.spatRes(base.spatRes[-length(base.spatRes)])+diff(fun.spatRes(base.spatRes))/2)
  totalEffects$intSpatRes <- intPointsSpaTrig%*%spatmids
  # Combining temporal and spatial response
  totalEffects$intRes <- totalEffects$intTempRes*totalEffects$intSpatRes
  
  gc()
  
  # Setting initial
  lmu0 <- log(mustar[i-1])
  beta0prec <- betastar[[i-1]]
  
  # Optimization
  res.optim <- optim(par=c(lmu0, beta0prec), NegLogLikehood, 
                     bENoMu = back.Events.NoMu, 
                     bANoMu = back.All.NoMu, 
                     tLNoA = trig.Lags.NoA,
                     trig.Int.NoA = totalEffects$intRes,
                     Xmat=Xmat, Xmat_all = Xmat_all,
                     inside = !dat$InBuf, secInside = !lags$InBuf, 
                     phiin = phi, rhoin = rho,
                     control=list(trace=0), hessian = T, method = "BFGS")
  
  mustar[i] <- exp(res.optim$par[1]) # The squared values were the ones to produce the best
  betastar[[i]] <- res.optim$par[2:(ncol(Xmat)+1)]  # likelihood
  logLik <- res.optim$value    # Storing the maximum ll value (minimum of negll)
  
  print(paste("mu=", round(mustar[i], 4), "beta=(", paste0(round(betastar[[i]], 2), collapse=", "), ") at Loop", i))
  print(paste("NegLoglik=", logLik))
  ll[i] <- logLik
  if(abs(ll[i] - ll[i-1]) < 1e-04) break;
}

gc()

save(dat, dataTogrid, ll, lags, mustar, betastar, 
     phi, rho, lambda.events, res.optim,
     back.Events.NoMu, back.All.NoMu, trig.Lags.NoA, totalEffects,
     weight.trend, fun.trend, 
     weight.daily, fun.daily, 
     weight.weekly, fun.weekly, 
     base.x, base.y, value.background, 
     fun.tempRes, fun.spatRes, NegLogLikehood,
     file = "WS/model_results_42019.RData")
