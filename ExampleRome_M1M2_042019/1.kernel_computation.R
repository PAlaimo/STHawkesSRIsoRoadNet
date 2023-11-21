
##########################
### KERNEL COMPUTATION ###
##########################
# Dr. Marco Mingione and Dr. Pierfrancesco Alaimo Di Loro, 20/11/2023

# This code creates all necessary ingredients for the nonparametric estimation of the spatio-temporal periodic Hawkes process for car accidents on a road network. In particular, it is used to compute the kernels.
# Example: Rome, municipalities I and II, April 2019

# Loading Packages --------------------------------------------------------
require(tidyverse)
require(magrittr)
require(MASS)
require(sf)
require(mvtnorm)
require(Matrix)

# Auxiliary ---------------------------------------------------------------

Rcpp::sourceCpp("auxiliary.cpp")

# This function compute the integral over the polygon defined by the city boundary
polyCub_mia <- function(pts, polyID, mean, sigma, stepl){
  foo <- dmvnorm(pts, mean = mean, sigma = diag(sigma, length(mean)))*stepl^2
  return(sum(foo[polyID]))
}

# Loading data ------------------------------------------------------------

load("WS/Real_PrepToEstPre_42019.RData")

# Kernels temporal background-----------------------------------------------------------------

print("Start kernels bg")

# Trend kernels
bw.trend <- 7 # We want a pretty smooth function over the time
kernelWithEdgeCorr.trend <- kernelWithEdgeCorr_trend(base.trend, dat$t, bw.trend, start.trend, end.trend)

# Daily kernels
bw.daily <- 0.05
kernel.dailyMat <- kernel_periodic(base.daily, dat$timeday[!dat$InBuf], bw.daily, 1)

# Weekly kernels
bw.weekly <- 1
kernel.weeklyMat <- kernel_periodic(base.weekly, dat$timeweek[!dat$InBuf], bw.weekly, 7)

# Kernels spatial background ----------------------------------------------
# The following function returns the dataset with indiviudual bandwidth
# and edge correction for the background smoothing of each observation
# Defining the minimum bandwidth
dat %<>% mutate(bw=0.1)
a2 <- bwCalc_cpp(dat$x[!dat$InBuf], dat$y[!dat$InBuf], 5)
dat$bwProp <- 0.1
dat$bwProp[!dat$InBuf] <- a2
dat %<>% mutate(bw = case_when(bwProp>bw ~ bwProp, TRUE ~ bw)) %>% dplyr::select(-bwProp)
rm(a2)

# Setting up parallel computation
mc.cores <- parallel::detectCores(logical = T)
n_cores <- min(mc.cores, nrow(dat)) - 2

##### ONLY LINUX #####
##### IT DOES NOT WORK in parallel for Windows/Mac users. Just change parallel::mclapply with lapply #####
dat$int <- parallel::mclapply(1:nrow(dat), function(i){ # can be time consuming
  # print(i)
  polyCub_mia(pts = grid.background, 
              polyID = rep(T, nrow(grid.background)), mean = as.numeric(dat[i,c("x", "y")]), 
              sigma = dat$bw[i], stepl = step.x)
}, mc.cores = n_cores ) %>% unlist()
##### ONLY LINUX #####

# ##### UNCOMMENT THIS FOR WINDOWS #####
# ##### WINDOWS #####
# require(doParallel)
# cl <- makeCluster(n_cores-2)
# registerDoParallel(cl)
# clusterExport(cl = cl, c("grid.background", "dat", "step.x", "polyCub_mia"),
#               envir = .GlobalEnv)
# clusterEvalQ(cl, library("tidyverse"))
# clusterEvalQ(cl, library("magrittr"))
# clusterEvalQ(cl, library("mvtnorm"))
# dat$int <- parLapply(cl=cl, 1:nrow(dat), function(i){ # can be time consuming
#   # print(i)
#   polyCub_mia(pts = grid.background, 
#               polyID = rep(T, nrow(grid.background)), mean = as.numeric(dat[i,c("x", "y")]), 
#               sigma = dat$bw[i], stepl = step.x)
# }) %>% unlist()
# stopCluster(cl)
# ##### WINDOWS #####


######### Spatial background kernels with edge correction ######### 
# Background kernels with edge correction
kernelWithEdgeCorr.backgroundx <- kernelWithEdgeCorr_backgroundx(unique(base.x), dat$x, 
                                                                 dat$bw, dat$int)
kernelWithEdgeCorr.backgroundy <- kernelWithEdgeCorr_backgroundy(unique(base.y), dat$y, 
                                                                 dat$bw, dat$int)
gc()

sparse <- 0
if(sparse == 1){
  kernelWithEdgeCorr.backgroundx[kernelWithEdgeCorr.backgroundx<0.1e-6] <- 0
  kernelWithEdgeCorr.backgroundy[kernelWithEdgeCorr.backgroundy<0.1e-6] <- 0
  kernelWithEdgeCorr.backgroundx <- as(kernelWithEdgeCorr.backgroundx, "sparseMatrix")
  kernelWithEdgeCorr.backgroundy <- as(kernelWithEdgeCorr.backgroundy, "sparseMatrix")
}

# Kernels temporal triggering-----------------------------------------------------------------
# Temporal excitation kernels
bw.tempRes <- 0.03 # We want the bandwidth to be larger than the step (which is 0.01)
kernelWithEdgeCorr.tempRes <- kerneledgeCorr_tempRes(base.tempRes, lags$timeLags[!lags$InBuf], 
                                                     bw.tempRes, end.tempRes)

# Kernels spatial triggering-----------------------------------------------------------------
# Spatial excitation kernels
bw.spatRes <- 0.05 # We want the bandwidth to be larger than the step (which is 0.05)
kernelWithEdgeCorr.spatRes <- kerneledgeCorr_spatResISOrep(base.spatRes, lags$dxy[!lags$InBuf], 
                                                           bw.spatRes, end.spatRes)

# Saving ------------------------------------------------------------------

save(dat, lags, grid.background, sparse,
     step.x, step.y, base.x, base.y,
     step.trend, base.trend, start.trend, end.trend,
     step.daily, base.daily, step.weekly, base.weekly,
     step.tempRes, base.tempRes, start.tempRes, end.tempRes,
     step.spatRes, base.spatRes, start.spatRes, end.spatRes,
     kernelWithEdgeCorr.trend, kernel.dailyMat, kernel.weeklyMat,
     kernelWithEdgeCorr.backgroundx, kernelWithEdgeCorr.backgroundy,
     kernelWithEdgeCorr.tempRes, kernelWithEdgeCorr.spatRes,
     bw.spatRes, bw.tempRes, bw.trend, bw.daily, bw.weekly,
     file = "WS/Real_EstPre_42019.RData")
