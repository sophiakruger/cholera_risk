## Manipulating Future (Ftr) Climate and Environmental Variable Data ##

library(arcgisbinding)
arc.check_product()
library(raster)

source("RasterManipulation.R")
setwd("../CholeraRiskProject/cholera_risk_mapping")

##### Imputing missing values of CoastalDEM 2050 projected layer from 2015 Elevation Raster Layer #####

## Convert future CoastalDEM/Elevation Data to RasterLayer, apply CRS, resample to 250m resolution, and crop to BGD
elev_ftr_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/elev2050.tif")) #this has the uppper missing chunk as values 0
elev_ftr <- as.raster(elev_ftr_arc)
crs(elev_ftr) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
elev_ftr <- resample(elev_ftr, elev_r, method="ngb") #resample to 250 res of elev current raster from RasterManipulation.R
elev_ftr <- crop(elev_ftr, BGD, snap='near')

## Mask out the extent values
elev_cr.m <- mask(env_rasters$elev, BGD, inverse=FALSE, updatevalue=NA, updateNA=FALSE)
elev_ftr.m <- mask(elev_ftr, BGD, inverse=FALSE, updatevalue=NA, updateNA=FALSE)

## Extract raster cell values and cell numbers from future/current data
elev_fr <- extract(elev_ftr.m, BGD, cellnumbers=TRUE)
elev_cr <- extract(elev_cr.m, BGD, cellnumbers=TRUE) 

## Create matrix of extracted cell values for current, future
elev_vals <- cbind(as.vector(elev_cr[[1]][,2]), as.vector(elev_fr[[1]][,2]))

## Regress current elevation to known area to future elevation's available data
elev_mod <- lm(elev_vals[,2] ~ elev_vals[,1] + I(elev_vals[,1]^2) + I(elev_vals[,1]^3) +I(elev_vals[,1]^4)) #R^2 = 0.952; y <- 8.019107e-01*x + 1.374737e-03*x^2 + -2.853147e-06*x^3 + 1.763818e-09*x^4

## Restrict cr & ftr rasters to extent of missing "chunk" of future elev data
ext_chunk <- extent(88.01069, 89.1, 25.9, 26.6337) 
elev_cr.ch <- crop(elev_cr.m, ext_chunk)
elev_ftr.ch <- crop(elev_ftr.m, ext_chunk)
NAvalue(elev_ftr.ch) <- 0

## Extract values for extent of missing chunk
elev_ftr.vals <- as.data.frame(extract(elev_ftr.ch, BGD, cellnumbers=TRUE))
elev_cr.vals <- as.data.frame(extract(elev_cr.ch, BGD, cellnumbers=TRUE))

## Find missing Values (NAs)in future raster's missing chunk
no_data <- which(is.na(elev_ftr.vals[,2])) #row numbers of elev_ftr.ch

## Find matching cell numbers and values of cr data for rows of ftr data's missing values -- this is not all cells from the chunk, just NA cells
no_data_cells <- as.matrix(elev_cr.vals[no_data,])

## Apply the elev_mod linear function to the cr chunk's values, imputing future chunk's values
no_data_vals <- as.array(no_data_cells[,2]) #same #values of no_data
imputing_missing_vals <- function(x) {
  y <- 8.019107e-01*x + 1.374737e-03*x^2 + -2.853147e-06*x^3 + 1.763818e-09*x^4
  return(y)
}
imputed_vals <- apply(no_data_vals, 1, imputing_missing_vals) #60,963 cells that had missing values imputed (for elev_ftr.ch raster)

## Fill the NA values of ftr data's chunk with the imputed values
filled_data_ch <- as.matrix(elev_ftr.vals[no_data,]) #grab the cell numbers of elev_ftr.vals (chunk values) NA cells
filled_data_ch[,2] <- imputed_vals #60,963 rows, each representing one cell of the chunk that was missing data

## Create New Raster for the missing "chunk" of elev_ftr.m with values of filled_data.ch (imputed values)
points <- xyFromCell(elev_ftr.ch, filled_data_ch[,1]) #to get coordinates from the cell  numbers, with proj and cell size determined by elev_ftr.ch raster
chunk_pts <- SpatialPoints(points, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) #create spatial points
filled_chunk.r <- rasterize(chunk_pts, field=filled_data_ch[,2], elev_ftr.ch, update=TRUE) #set imputed values to points of cell numbers that were NA, then rasterize points to cells

## Create final elevation raster for 2050, mosaicking the future elevation layer and missing chunk
elev2050 <- mosaic(filled_chunk.r, elev_ftr.m, fun=max) #pick max value between the two rasters so that the filled chunk is not 0, rather that largest value (which in comparison will be the value in filled_chunk.r for that cell)

## Add elevation 2050 raster to ArcGIS Project
write.loc <- file.path('cholera_risk_mapping.gdb', 'elev2050_raster')
arc.write(write.loc, elev2050)

rm(elev_ftr_arc, elev_ftr, elev_cr.m, elev_fr, elev_cr, elev_mod, ext_chunk, elev_cr.ch, elev_ftr.ch, elev_ftr.vals, elev_cr.vals, no_data_cells, no_data_vals, imputing_missing_vals, imputed_vals, filled_data_ch, points, chunk_data_points, chunk_pts, chunk_sp)
gc()

##### Manipulating Other Future Raster Data #####

## Create SpatialPolygonDataFrame for the boundary of Bangladesh provided by GADM
BGD <- getData("GADM", country="BGD", level=0) #level=0 returns country boundary, level=1 returns district boundaries, extent= 88.01057, 92.67366, 20.74111, 26.6344  (xmin, xmax, ymin, ymax)

## Pop, landcover, dist to water rasters manipulated manually in arcGIS, then imported as tif into R to save memory
lc2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/lc2050.tif")) #300m resolution originally
pop2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/pop2050.tif"))  #originally 1km cells, so any raw population size will be same as raw pop density (per km) in 2015 data
d2w2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/Dist_to_water2050.tif")) #raster values are raw values, not stretch pixels
d2w2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/testD2W_2050.tif")) #gets d2w raster in 0-255 pixel value range as 2015 raster

## Temp and PPT rasters include 12 bands, each with the monthly average values for the period 2041-2060, 4.5km res original
# Average of Jan-Oct months for this period, to match avg for 2015
ppt2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_prec_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_prec_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(1, 10, 11, 12))
tmx2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_tmax_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_tmax_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(1, 10, 11, 12))
tmn2050_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_tmin_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_tmin_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(1, 10, 11, 12))
# Average of June-Sept (monsoon) months for this period, to match avg for 2015
ppt2050m_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_prec_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_prec_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(6, 7, 8, 9))
tmx2050m_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_tmax_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_tmax_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(6, 7, 8, 9))
tmn2050m_arc <- arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/wc2.1_2.5m_tmin_BCC-CSM2-MR_ssp245_2021-2040/share/spatial03/worldclim/cmip6/7_fut/2.5m/BCC-CSM2-MR/ssp245/wc2.1_2.5m_tmin_BCC-CSM2-MR_ssp245_2021-2040.tif"), bands=c(6, 7, 8, 9))
#Importing elev2050 raster from ArcGIS, after it was created in R and exported to ArcGIS
elev2050 <- as.raster(arc.raster(arc.open("../CholeraRiskProject/cholera_risk_mapping/elev2050_raster.tif")))

ftr.raster_list <- list(lc2050_arc, pop2050_arc, d2w2050_arc) #ftr = future

## Convert each Arc Raster to RasterLayer Object, then Resample to 250m resolution of 2015 elevation layer, to standardize cell size for all current/future rasters
ftr_rasters <- list()
x<-1
for(i in ftr.raster_list)
{
  raster_r <- as.raster(i)
  crs(raster_r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  raster_rs <- resample(raster_r, elev_r, method="ngb") 
  raster_c <- crop(raster_rs, BGD, snap='near') #crop to resolve memory issues
  ftr_rasters[[x]] <- raster_c
  x<-x+1
}
## Add Elevation Raster, ensure it is projected, resampled, and cropped accurately 
extent(elev2050) <- extent(elev_r) #ensure extent of elev is similar to that in for loop above: elev_r
crs(elev2050) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
elev2050 <- resample(elev2050, elev_r, method="ngb")
elev2050 <- crop(elev2050, BGD, snap='near')
ftr_rasters[[4]] <- elev2050

rm(lc2050_arc, pop2050_arc, d2w2050_arc, raster_r, raster_rs, raster_c)
gc()

## Manipulating Oct-Jan & Monsoon month data from Temp/PPT rasters 
## Convert Arc Rasters to RasterLayers
ppt_rs <- as.raster(ppt2050_arc)
tmx_rs <- as.raster(tmx2050_arc)
tmn_rs <- as.raster(tmn2050_arc)
pptm_rs <- as.raster(ppt2050m_arc)
tmxm_rs <- as.raster(tmx2050m_arc)
tmnm_rs <- as.raster(tmn2050m_arc)

stacks <- list(ppt_rs, tmx_rs, tmn_rs, pptm_rs, tmxm_rs, tmnm_rs)
rm(ppt2050_arc, tmx2050_arc, tmn2050_arc, ppt_rs, tmx_rs, tmn_rs, pptm_rs, tmxm_rs, tmnm_rs)
gc()
x<-5
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
  ftr_rasters[[x]] <- avg
  x<-x+1
}
names(ftr_rasters) <- c("landcover","pop_dens", "dist_to_water", "elev", "avg_ppt","avg_tmax","avg_tmin", "avg_ppt_m","avg_tmax_m","avg_tmin_m")
rm(pt_c, avg)
gc()

#ftr_rasters is final list of climate rasters to use as predictors in predicting each SDM 
