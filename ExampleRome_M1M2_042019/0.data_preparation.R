#########################
### DATA PREPARATION  ###
#########################
# Dr. Marco Mingione and Dr. Pierfrancesco Alaimo Di Loro, 20/11/2023

# This code creates all necessary ingredients for the nonparametric estimation of the spatio-temporal periodic Hawkes process for car accidents on a road network.
# Example: Rome, municipalities I and II, April 2019

# Loading Packages --------------------------------------------------------
require(tidyverse)
require(magrittr)
require(lubridate)
require(mvtnorm)
require(raster)
require(sf)
require(stars)

# Auxiliary ---------------------------------------------------------------

Rcpp::sourceCpp("auxiliary.cpp")

# Define spatial grid -----------------------------------------------------

# Load street network
street_net <- st_read("Data/network_m1m2.shp")
crs_streetnet <- st_crs(street_net)

# Load shapefile municipalities
shp_munic <- st_read("Data/municipi_I_II.shp")
shp_munic_buf <- st_read("Data/municipi_I_II_buffer.shp")

# Boundaries
pmap <- shp_munic %>% st_cast("LINESTRING") %>% suppressWarnings()
pmap_all <- shp_munic_buf %>% st_cast("LINESTRING") %>% suppressWarnings()

# Define spatial resolution
step.x <- 0.01
step.y <- 0.01
bbox_buf <- st_bbox(street_net)

# Network as a raster (include NA)
streetnet_rast <- st_rasterize(street_net,
                               template = st_as_stars(bbox_buf, dx = step.x, dy=step.y, values = NA_real_))

# Keep as points only the cells containing a road (no NA)
roadcells <- streetnet_rast %>% st_xy2sfc(as_points = T) %>% st_geometry() %>% st_set_crs(value = crs_streetnet)
roadcellsp <- streetnet_rast %>% st_xy2sfc(as_points = F) %>% st_geometry() %>% st_set_crs(value = crs_streetnet)
streetnet_roma_in <- st_intersection(street_net, shp_munic)
streetnet_rast_in <- st_rasterize(streetnet_roma_in,
                                  template = st_as_stars(bbox_buf, dx = step.x, dy=step.y, values = NA_real_))
roadcellspIn <- streetnet_rast_in %>% st_xy2sfc(as_points = F) %>% st_geometry() %>% st_set_crs(value = crs_streetnet)
rm(streetnet_rast_in)
gc()

# Build up the grid
grid.background <- roadcells %>% st_coordinates() %>% as_tibble() %>% as.matrix()
base.x <- unique(grid.background[,1,drop=T])
base.y <- unique(grid.background[,2,drop=T])

# Grid points inside the GRA boundaries
grid.Inpoly <- st_intersects(roadcells, shp_munic, sparse = F, model = "open")

# Loading Data ------------------------------------------------------------

dat <- read_csv("Data/roadaccidents_m1m2_42019.csv")
# Rename everything and scale times to be in days
dat <- dat %>%
  mutate(tNumModScale = (as.numeric(DataOraIncidente) - as.numeric(ymd_hms("2019-04-01 00:00:00"))),
         t = tNumModScale/(3600*24)+runif(n(), -0.01, 0.01),
         xutm = xutm + runif(n(), -0.02, 0.02),
         yutm = yutm + runif(n(), -0.02, 0.02)) %>% rename(x = xutm, y = yutm) %>% arrange(t)

# Define temporal background grid ----------------------------------------------------------

# Trend base
step.trend <- 0.02 # ogni giorno
start.trend <- floor(min(dat$t))
end.trend <- ceiling(max(dat$t))
base.trend <- seq(start.trend, end.trend, by=step.trend)

# Daily periodicity
dat$timeday <- dat$t - floor(dat$t) # Creating the "time into day" variable
step.daily <- 0.02 # ogni mezz'ora
base.daily <- seq(0, 1, step.daily)

# weekly periodicity
dat$timeweek <- (dat$t)%%7 # Creating the "time into week" variable
step.weekly <- 0.02 # ogni mezza giornata
base.weekly <- seq(0, 7, step.weekly)


# Define spatial background grid ----------------------------------------------------------

# Make data spatial data
foo <- dat[, c("x", "y")] %>% st_as_sf(coords = c("x", "y"), crs = crs_streetnet)

# Set-up parallel computation
mc.cores <- parallel::detectCores(logical = T)
n_cores <- min(mc.cores, nrow(dat)) - 2

##### ONLY LINUX #####
##### IT DOES NOT WORK in parallel for Windows/Mac users. Just change parallel::mclapply with  lapply #####
# Snap points onto network
crds_proj <- parallel::mclapply(seq(nrow(foo)), \(i){
  (st_nearest_points(foo[i,], roadcellsp[st_nearest_feature(foo[i,], roadcellsp)]) %>% st_cast("POINT"))[2]
}, mc.cores = n_cores) %>% do.call(what = c, args = .)
##### ONLY LINUX #####

# ##### UNCOMMENT THIS FOR WINDOWS #####
# ##### WINDOWS #####
# require(doParallel)
# cl <- makeCluster(n_cores-2)
# registerDoParallel(cl)
# clusterExport(cl = cl, c("roadcellsp", "foo"), envir = .GlobalEnv)
# clusterEvalQ(cl, library("tidyverse"))
# clusterEvalQ(cl, library("magrittr"))
# clusterEvalQ(cl, library("sf"))
# # Snap points onto network
# crds_proj <- parLapply(cl = cl, seq(nrow(foo)), \(i){
#   (st_nearest_points(foo[i,], roadcellsp[st_nearest_feature(foo[i,], roadcellsp)]) %>% st_cast("POINT"))[2]
# }) %>% do.call(what = c, args = .)
# stopCluster(cl)
# ##### WINDOWS #####

# Update coordinates after snapping onto network
dat$x <- st_coordinates(crds_proj)[,1,drop=T]
dat$y <- st_coordinates(crds_proj)[,2,drop=T]
foo <- dat[, c("x", "y")] %>% st_as_sf(coords = c("x", "y"), crs = crs_streetnet)

# Space base: Associate background cell to each data-point
dataTogrid <- st_nearest_feature(foo, roadcells)

# Excitation functions -----------------------------------------------------

# Time excitation base
step.tempRes <- 0.01
start.tempRes <- 0
end.tempRes <- .3
base.tempRes <- seq(start.tempRes, end.tempRes, step.tempRes)+0.6e-7 # Add a small value to have 0 in 0

# Space excitation base
step.spatRes <- 0.01
start.spatRes <- 0
end.spatRes <- .5
base.spatRes <- seq(start.spatRes, end.spatRes, by=step.spatRes)+0.6e-7

# Set up excitation lags -------------------------------------------------------

# Lag matrix
lags <- (data.table::CJ(1:nrow(dat), 1:nrow(dat)) %>% as.matrix())[,2:1]
lags <- lags %>% as.data.frame %>% 
  left_join(dat %>% mutate(ID = 1:nrow(.)) %>% dplyr::select(ID, InBuf), by = c("V1" = "ID")) %>% 
  left_join(dat %>% mutate(ID = 1:nrow(.)) %>% dplyr::select(ID, NaturaIncidenteNew), by = c("V2" = "ID"))
# Dropping pairs in which triggered happened before triggering
lags <- lags[lags[,2]>lags[,1],] 
colnames(lags)[1:2] <- c("idFirst", "idSecond")
lags$timeLags <- dat$t[lags$idSecond]-dat$t[lags$idFirst]
lags <- lags[lags$timeLags<end.tempRes & lags$timeLags>0,] # Dropping pairs too distant in time
lags$dxy <- sqrt((dat$x[lags$idSecond]-dat$x[lags$idFirst])^2+(dat$y[lags$idSecond]-dat$y[lags$idFirst])^2)
lags <- lags[lags$dxy<end.spatRes & lags$dxy>0,] # Dropping pairs too distant over x

gc()

#----------------- Repetition correction ------------------
# Temporal repetition correction
repTempTrig <- sapply(base.tempRes, function(i) sum((i+dat$t)<end.trend))

# Spatial repetition correction
##### ONLY LINUX #####
##### IT DOES NOT WORK in parallel for Windows/Mac users. Just change parallel::mclapply with  lapply #####
repSpatTrig <- parallel::mclapply(base.spatRes, function(i)
{
  temp <- st_buffer(foo, i) %>% st_cast("LINESTRING")
  sum(st_length(st_intersection(roadcellsp, temp)))
}, mc.cores = n_cores) %>% unlist()
##### ONLY LINUX #####


# ##### UNCOMMENT THIS FOR WINDOWS #####
# ##### WINDOWS #####
# cl <- makeCluster(n_cores-2)
# registerDoParallel(cl)
# clusterExport(cl = cl, c("roadcellsp", "foo"),
#               envir = .GlobalEnv)
# clusterEvalQ(cl, library("tidyverse"))
# clusterEvalQ(cl, library("magrittr"))
# clusterEvalQ(cl, library("sf"))
# repSpatTrig <- parLapply(cl=cl, base.spatRes, function(i)
# {
#   temp <- st_buffer(foo, i) %>% st_cast("LINESTRING")
#   sum(st_length(st_intersection(roadcellsp, temp)))
# }) %>% unlist()
# stopCluster(cl)
# ##### WINDOWS #####

# Calculate spatial excitation integration points of all events
dfac <- factor(base.spatRes)
oop <- tibble(x=0, y=0) %>% st_as_sf(coords=c("x", "y"), crs=crs_streetnet)
spatCrowns <- map_dfr(1:length(base.spatRes), \(j) 
                      st_difference(st_buffer(oop, base.spatRes[j]), st_buffer(oop, c(0, base.spatRes)[j]))) %>% 
  mutate(d=dfac)

##### ONLY LINUX #####
##### IT DOES NOT WORK in parallel for Windows/Mac users. Just change parallel::mclapply with lapply #####
intPointsSpaTrig <- parallel::mclapply(1:nrow(foo), function(i){ # Time consuming...
  crownsnow <- spatCrowns
  st_geometry(crownsnow) <- (spatCrowns$geometry + foo[i,]$geometry)
  st_crs(crownsnow) <- crs_streetnet
  temp <- st_intersection(crownsnow, roadcellspIn)
  temp %>% mutate(A=as.numeric(st_area(.))) %>% as_tibble() %>% group_by(d, .drop=F) %>% summarise(A=sum(A)) %$% A
  
}, mc.cores = n_cores) %>% do.call(what = rbind, args = .)
##### ONLY LINUX #####

# ##### UNCOMMENT THIS FOR WINDOWS #####
# ##### WINDOWS #####
# cl <- makeCluster(n_cores-2)
# registerDoParallel(cl)
# clusterExport(cl = cl, c("dfac", "roadcellspIn", "foo", "base.spatRes",
#                          "spatCrowns", "crs_streetnet"),
#               envir = .GlobalEnv)
# clusterEvalQ(cl, library("tidyverse"))
# clusterEvalQ(cl, library("magrittr"))
# clusterEvalQ(cl, library("sf"))
# intPointsSpaTrig <- parLapply(cl=cl, 1:nrow(foo), function(i){ # Time consuming...
#   crownsnow <- spatCrowns
#   st_geometry(crownsnow) <- (spatCrowns$geometry + foo[i,]$geometry)
#   st_crs(crownsnow) <- crs_streetnet
#   temp <- st_intersection(crownsnow, roadcellspIn)
#   temp %>% mutate(A=as.numeric(st_area(.))) %>% as_tibble() %>% group_by(d, .drop=F) %>%
#     summarise(A=sum(A)) %$% A
# }) %>% do.call(what = rbind, args = .)
# stopCluster(cl)
# ##### WINDOWS #####

rm(foo)
rm(roadcellsp)
rm(roadcellspIn)
gc()

# Set up covariates -------------------------------------------------------

# Design matrices
Xmat_all <- model.matrix(~., data = dat %>% dplyr::select(NaturaIncidenteNew))
Xmat <- model.matrix(~., data = lags %>% dplyr::select(NaturaIncidenteNew))

# Saving ------------------------------------------------------------------

save(dat, lags, grid.background, shp_munic_buf,
     crs_streetnet, bbox_buf,
     step.x, step.y, base.x, base.y, 
     step.trend, base.trend, start.trend, end.trend,
     step.daily, base.daily, step.weekly, base.weekly, 
     step.tempRes, base.tempRes, start.tempRes, end.tempRes, 
     step.spatRes, base.spatRes, start.spatRes, end.spatRes, 
     file = "WS/Real_PrepToEstPre_42019.RData")

save(Xmat, Xmat_all, grid.Inpoly, pmap, pmap_all, 
     roadcells, dataTogrid, crs_streetnet, 
     repTempTrig, repSpatTrig, intPointsSpaTrig, 
     file = "WS/Real_PrepToEst_42019.RData")

