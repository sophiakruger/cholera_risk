setwd("../CholeraRiskProject/cholera_risk_mapping")

library(arcgisbinding)
arc.check_product()
library(raster)

## Create SpatialPolygonDataFrame for the boundary of Bangladesh provided by GADM
BGD <- getData("GADM", country="BGD", level=0) #level=0 returns country boundary, level=1 returns district boundaries, extent= 88.01057, 92.67366, 20.74111, 26.6344  (xmin, xmax, ymin, ymax)

## Read in Environmental Rasters from ArcGIS Project folder
landcover_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/Landcover.tif"))     
elev_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/Elevation250m_2.tif")) #tif from ArcGIS modelbuilder manipulations, should be stretched values of dataset now, not colormap values
dist_to_water_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/Dist_to_Water.tif")) #tif from ArcGIS modelbuilder manipulations, should be stretched values of dataset now, not colormap values
drp_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/DRP.tif"))
tss_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/TSS.tif"))
pop_dens_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/pop_density_GPWv4.tif"))
dist_to_coast_arc <- arc.raster(arc.open('../CholeraRiskProject/cholera_risk_mapping/Dist_to_coast.tif'))

raster_list <- list(landcover_arc, elev_arc, dist_to_water_arc, drp_arc, tss_arc, pop_dens_arc, dist_to_coast_arc)

## Convert each raster_list item to RasterLayer Object, then Resample to 250m resolution
elev_r <- crop(as.raster(elev_arc), BGD, snap='near') #250m resolution
env_rasters <- list()
x<-1
for(i in raster_list)
{
  raster_r <- as.raster(i)
  crs(raster_r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  raster_rs <- resample(raster_r, elev_r, method="ngb") 
  raster_c <- crop(raster_rs, BGD, snap='near') #crop to resolve memory issues
  env_rasters[[x]] <- raster_c
  x<-x+1
}
names(env_rasters)<- c("landcover", "elev", "dist_to_water", "drp", "tss", "pop_dens", "dist_to_coast")
rm(raster_r, raster_rs, raster_c, landcover_arc, elev_arc, dist_to_water_arc, drp_arc, tss_arc, pop_dens_arc, dist_to_coast_arc)
gc()

## Read TIF files for precip and temp created in 'nc_to_tiff.R' script

# Create RasterStacks for layers of months Oct 2015-Jan 2016
#source("nc_to_tiff.R") 

ppt<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/precip_acc_2015.tif"), bands=c(10,11,12))
ppt2<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/precip_acc_2016.tif"), bands=1)
ppt_r <- as.raster(ppt) 
ppt2_r<- as.raster(ppt2)
ppt_s <- stack(ppt_r, ppt2_r)
rm(ppt, ppt2, ppt_r, ppt2_r)
gc()

tmax<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_max_2015.tif"), bands=c(10,11,12))
tmax2<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_max_2016.tif"), bands=1)
tmax_r <- as.raster(tmax)
tmax2_r<- as.raster(tmax2)
tmax_s <- stack(tmax_r, tmax2_r)
rm(tmax, tmax2, tmax_r, tmax2_r)
gc()

tmin<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_min_2015.tif"), bands=c(10,11,12))
tmin2<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_min_2016.tif"), bands=1)
tmin_r <- as.raster(tmin) 
tmin2_r<- as.raster(tmin2)
tmin_s <- stack(tmin_r, tmin2_r)
rm(tmin, tmin2, tmin_r, tmin2_r)
gc()

## For the ppt and temp raster stacks, project to WGS84, crop to BGD, and resample to resolution of pop_density layer, then create raster of mean values of stack rasters 
stacks <- list(ppt_s, tmax_s, tmin_s)
remove(ppt_s, tmax_s, tmin_s)
gc()
pt_rasters<-list()
y<-1
for (k in stacks) 
{
  #pt_prj <- projectRaster(k, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #need if projections are not WGS84
  crs(k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  pt_rs<- resample(k, elev_r, method="ngb")
  print("test2")
  pt_c <- crop(pt_rs, BGD, snap='near')
  rm(pt_rs)
  gc()
  print("test3")
  avg <- mean(pt_c) #mean of the re-projected, resampled, and cropped 4 layers
  print("test4")
  pt_rasters[[y]] <- avg
  y<-y+1
}
names(pt_rasters) <- c("avg_ppt","avg_tmax","avg_tmin")
rm(pt_c, avg)
gc()

## Create RasterStacks for layers of "Monsoon Months": June-September (techincally mid-october but october was counted in other avgs)
#source("nc_to_tiff.R") 

pptm<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/precip_acc_2015.tif"), bands=c(6, 7, 8, 9))
pptm_r <- as.raster(pptm) 
pptm_s <- stack(pptm_r)
rm(pptm, pptm_r)
gc()

tmaxm<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_max_2015.tif"), bands=c(6, 7, 8, 9))
tmaxm_r <- as.raster(tmaxm)
tmaxm_s <- stack(tmaxm_r)
rm(tmaxm, tmaxm_r)
gc()

tminm<- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/temp_min_2015.tif"), bands=c(6, 7, 8, 9))
tminm_r <- as.raster(tminm) 
tminm_s <- stack(tminm_r)
rm(tminm, tminm_r)
gc()

## For the monsoon ppt and temp raster stacks, project to WGS84, crop to BGD, and resample to resolution of pop_density layer, then create raster of mean values of stack rasters 
stacks_m <- list(pptm_s, tmaxm_s, tminm_s)
remove(pptm_s, tmaxm_s, tminm_s)
gc()
ptm_rasters<-list()
y<-1
for (k in stacks_m) 
{
  #pt_prj <- projectRaster(k, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #need if projections are not WGS84
  crs(k) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  pt_rs<- resample(k, elev_r, method="ngb")
  print("test2")
  pt_c <- crop(pt_rs, BGD, snap='near')
  rm(pt_rs)
  gc()
  print("test3")
  avg <- mean(pt_c) #mean of the re-projected, resampled, and cropped 4 layers
  print("test4")
  ptm_rasters[[y]] <- avg
  y<-y+1
}
names(ptm_rasters) <- c("avg_ppt_m","avg_tmax_m","avg_tmin_m")
rm(pt_c, avg)
gc()

