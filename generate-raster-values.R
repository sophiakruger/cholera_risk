library(dplyr)
library(raster)
library(sf)
library(tools)

source("RasterManipulation.R")

## Read in coordinate data (long, lat) from serosurvey data (Azman et al.)
setwd("../CholeraRiskProject")
## Load in serosurvey coordinates on version of R 3.6 or above
serosurvey_dat <- readRDS("../CholeraRiskProject/serosurvey-with-rf-preds.rds") 
sero_coords <- cbind(serosurvey_dat$lon, serosurvey_dat$lat)  
## Load in serosurvey coordinates on version of R 3.5.1 (due to error: 'cannot unserialize ALTVEC object of class 'wrap_logical' from package 'base'; returning length zero vector')
serosurvey_dat <- read.table("../CholeraRiskProject/serosurvey-w-preds-table")
sero_coords <- cbind(serosurvey_dat$lon, serosurvey_dat$lat)

## Coerce sero_coords into SpatialPoints object with same CRS as env rasters
coords<-SpatialPoints(sero_coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) #ESPG code 4326 is for the WGS84 geodetic CRS

## Bring in 'env_rasters', 'pt_rasters', and 'ptm_rasters' lists from RasterManipulation.R that store desired rasters, join to one list
env_vars <- list(env_rasters[[1]],env_rasters[[2]],env_rasters[[3]],env_rasters[[4]], env_rasters[[5]], env_rasters[[6]], env_rasters[[7]], pt_rasters[[1]],pt_rasters[[2]],pt_rasters[[3]], ptm_rasters[[1]], ptm_rasters[[2]], ptm_rasters[[3]])

## For each raster in list, extract the grid cell value and grid cell number for which each coordinate point belongs
rv<-list()
rc<-list()
k<-1
for(i in env_vars) 
{
  cell_vals <- raster::extract(i,coords, method='simple', cell_numbers=TRUE) #returns large vector of 2930 values (correlated to 2930 coords)
  cell_num <- raster::cellFromXY(i, coords) #cell numbers attributed to each value of the raster cell that contains coordinates
  rv[[k]]<- cell_vals
  rc[[k]]<- cell_num
  k<-k+1
}
## Bind all grid cell values to sero_coords, adding one column of cell numbers after lon/lat; convert landcover values to categories
rvalues<-cbind(sero_coords, rc[[1]], as.factor(rv[[1]]), rv[[2]], rv[[3]], rv[[4]], rv[[5]], rv[[6]], rv[[7]], rv[[8]], rv[[9]], rv[[10]], rv[[11]], rv[[12]], rv[[13]])
colnames(rvalues) <- c("lon","lat","cell_num", "landcover", "elev", "dist_to_water", "drp", "tss", "pop_dens", "dist_to_coast", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m")
rm(rc, rv)
gc()

## Bind pred100:pred365 cols of serosurvey_dat to rvalues and add col names
env_data<-cbind(rvalues, serosurvey_dat$pred_100, serosurvey_dat$pred_200, serosurvey_dat$pred_365)
colnames(env_data)<-c("lon","lat","cell_num", "landcover", "elev", "dist_to_water", "drp", "tss", "pop_dens", "dist_to_coast", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m", "pred100", "pred200", "pred365")
rm(rvalues, serosurvey_dat)
gc()

## Coerce env_data into a dataframe, to use dplyr and tibble functions
env_data<-as.data.frame(env_data)

## Add Presence col to df; add 1 to rows where there is >= one predicted case of cholera across 365 days, else 0
env_data$Presence <- 1
#env_data$Presence[env_data$pred100>0 & env_data$pred200>0 & env_data$pred365>0] <-1
env_data$Presence[env_data$pred100==0 & env_data$pred200==0 & env_data$pred365==0] <-0 

saveRDS(env_data, 'env_data.rds')

## Create a matrix of grid cell values for covariates associated to each grid cell containing a survey point
## Add a Count column to total persons surveyed and Prevalence column of the freq of positive tests per cell
i<-1
gcells<-unique(env_data$cell_num)
grid_data<-matrix(nrow=length(gcells), ncol=16)
for(j in gcells)
{
  grid_data[i,1]<-j
  k<-env_data$landcover[(env_data$cell_num==j)]
  grid_data[i,2]<-k[1]
  m<-env_data$elev[(env_data$cell_num==j)]
  grid_data[i,3]<-m[1]
  n<-env_data$dist_to_water[(env_data$cell_num==j)]
  grid_data[i,4]<-n[1]
  #new drp and tss rasters #8/4/21
  p<-env_data$drp[(env_data$cell_num==j)]
  grid_data[i,5]<-p[1]
  q<-env_data$tss[(env_data$cell_num==j)]
  grid_data[i,6]<-q[1]
  #new pop_dens and dist_to_coast rasters #8/9/21
  r<-env_data$pop_dens[(env_data$cell_num==j)]
  grid_data[i,7]<-r[1]
  s<-env_data$dist_to_coast[(env_data$cell_num==j)]
  grid_data[i,8]<-s[1]
  #temp and ppt rasters
  u<-env_data$avg_ppt[(env_data$cell_num==j)]
  grid_data[i,9]<-u[1]
  v<-env_data$avg_tmax[(env_data$cell_num==j)]
  grid_data[i,10]<-v[1]
  x<-env_data$avg_tmin[(env_data$cell_num==j)]
  grid_data[i,11]<-x[1]
  u2<-env_data$avg_ppt_m[(env_data$cell_num==j)]
  grid_data[i,12]<-u2[1]
  v2<-env_data$avg_tmax_m[(env_data$cell_num==j)]
  grid_data[i,13]<-v2[1]
  x2<-env_data$avg_tmin_m[(env_data$cell_num==j)]
  grid_data[i,14]<-x2[1]
  
  #creating count and cases cols
  y<-env_data$cell_num==j # a vector of TRUE/FALSE where cell_num = value of j (a unique value of cell_num)
  grid_data[i,15]<-sum(y) #count = sum of all the elements in vector n that are TRUE = # of total people tested since each row aligns to one coordinate 
  grid_data[i,16]<-sum(env_data$Presence[(env_data$cell_num==j)]) #cases = sum the values of Presence attributed to unique cell
  i<-i+1 #increase to next row number in grid_data
}
colnames(grid_data) <- c("cell_num", "landcover", "elev", "dist_to_water", "drp", "tss", "pop_dens","dist_to_coast", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m", "count", "cases")

## Add a prevalence (frequency) column to grid_data
grid_data<-as.data.frame(grid_data)
grid_data$prevalence <- grid_data$cases/grid_data$count
rm(i, j, k, m, n, p, q, r, s, u, v, x, u2, v2, x2, y, cell_num, cell_vals, gcells, i)

