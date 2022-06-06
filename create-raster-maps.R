## Create Raster Maps from Prediction RDS files from running RF model in MSI

library(raster)
library(arcgisbinding)
arc.check_product()

setwd("../CholeraRiskProject/rfpredictions/rfpredictions")

covariates <- readRDS("../CholeraRiskProject/cholera_risk_mapping/covariates.rds")
ftr_rasters <- readRDS("../CholeraRiskProject/cholera_risk_mapping/ftr_rasters.rds")
cellnumbers <- cellsFromExtent(covariates$elev, extent=extent(covariates$elev))
cell_coords <- cbind(cellnumbers, xyFromCell(covariates$elev, cell=cellnumbers))

## 2015 Map Creation
mean_cells <- matrix(data=NA, ncol=3, nrow=5959968) 
lower_cells <- matrix(data=NA, ncol=3, nrow=5959968)
upper_cells <- matrix(data=NA, ncol=3, nrow=5959968)

for(i in 1:5)
{
  cr.probs <- matrix(nrow=1000, ncol=length((1+(i-1)*1e6):(i*1e6)))
  mean_cellvals <- matrix(ncol=3, nrow=ncol(cr.probs)) 
  lower_cellvals <- matrix(ncol=3, nrow=ncol(cr.probs))
  upper_cellvals <- matrix(ncol=3, nrow=ncol(cr.probs))
  # Create truncated matrix of cells (cols) with values from RF model's 1000 iterations
  for(j in 1:1000)
  {
    setwd("../CholeraRiskProject/rfpredictions/rfpredictions")
    cr.probs[j,] <- readRDS(paste('Replicate',j,'/2015preds.rds',sep=''))[,(1+(i-1)*1e6):(i*1e6)]
    #setwd('..')
  }
  # Find mean and percentiles for each cell, x (col of cr.probs)
  for (x in 1:ncol(cr.probs))
  {
    mean_cellvals[x,3] <- mean(cr.probs[,x]) #cr probs goes 1-1 million, mean_cellvals needs to go to specific rows (e.g., 1000001-2000000)
    mean_cellvals[x,1:2] <- cell_coords[x,2:3]
    lower_cellvals[x,3] <- quantile(cr.probs[,x], 0.025) 
    lower_cellvals[x,1:2] <- cell_coords[x,2:3]
    upper_cellvals[x,3] <- quantile(cr.probs[,x], 0.975)
    upper_cellvals[x,1:2] <- cell_coords[x,2:3]
  }
  mean_cells[(1+(i-1)*1e6):(i*1e6),] <- mean_cellvals
  lower_cells[(1+(i-1)*1e6):(i*1e6),] <- lower_cellvals
  upper_cells[(1+(i-1)*1e6):(i*1e6),] <- upper_cellvals
}

# Complete the cellvals matrices with the remaining cells (<1e6)
cr.probs <- matrix(nrow=1000, ncol=length(c(5000001:5959968)))
mean_cellvals <- matrix(data=NA, ncol=3, nrow=ncol(cr.probs)) 
lower_cellvals <- matrix(data=NA, ncol=3, nrow=ncol(cr.probs))
upper_cellvals <- matrix(data=NA, ncol=3, nrow=ncol(cr.probs))

# Create truncated matrix of cells (cols) with values from RF model's 1000 iterations
for(j in 1:1000)
{
  setwd("../CholeraRiskProject/rfpredictions/rfpredictions")
  cr.probs[j,] <- readRDS(paste('Replicate',j,'/2015preds.rds',sep=''))[,c(5000001:5959968)]
  #setwd('..')
}

# Find mean and percentiles for each cell, x (col of fr.probs)
for(x in 1:length(5000001:5959968))
{
  mean_cellvals[x,3] <- mean(cr.probs[,x]) 
  mean_cellvals[x, 1:2] <- cell_coords[x,2:3]
  lower_cellvals[x,3] <- quantile(cr.probs[,x], 0.025) 
  lower_cellvals[x,1:2] <- cell_coords[x,2:3]
  upper_cellvals[x,3] <- quantile(cr.probs[,x], 0.975)
  upper_cellvals[x,1:2] <- cell_coords[x,2:3]
}
mean_cells[5000001:5959968,] <- mean_cellvals
lower_cells[5000001:5959968,] <- lower_cellvals
upper_cells[5000001:5959968,] <- upper_cellvals

saveRDS(mean_cells, 'mean_cells_2015.rds')
saveRDS(lower_cells, 'lower_cells_2015.rds')
saveRDS(upper_cells, 'upper_cells_2015.rds')

# Create SpatialPoints object from coords, then rasterize to same spatial extent as ftr_rasters$elev with given cell values
mean_pts <- SpatialPoints(mean_cells[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
mean_OP_r <- rasterize(mean_pts, covariates$elev, field=mean_cells[,3])
lower_extreme_pts <- SpatialPoints(lower_cells[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
lower_extreme_OP_r <- rasterize(lower_extreme_pts, covariates$elev, field=lower_cells)
upper_extreme_pts <- SpatialPoints(upper_cells[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
upper_extreme_OP_r <- rasterize(upper_extreme_pts, covariates$elev, field=upper_cells[,3])

# Set NAs <- 1 so color ramp is scaled 0 to 1
mean_OP_r[is.na(mean_OP_r)]<- 1
lower_extreme_OP_r[is.na(lower_extreme_OP_r)] <- 1
upper_extreme_OP_r[is.na(upper_extreme_OP_r)] <- 1

# Mask rasters to extent of Bangladesh, BGD, removing background values
BGD <- getData("GADM", country="BGD", level=0) 
cr_map_mean <- mask(mean_OP_r, BGD, inverse=FALSE)
cr_map_2.5 <- mask(lower_extreme_OP_r, BGD, inverse=FALSE)
cr_map_97.5 <- mask(upper_extreme_OP_r, BGD, inverse=FALSE)

# Ratify rasters to create Raster Attribute Table (RAT)
mean2015_RAT <- levels(ratify(cr_map_mean, count=TRUE))[[1]]
q2.5_2015_RAT <- levels(ratify(cr_map_2.5, count=TRUE))[[1]]
q97.5_2015_RAT <- levels(ratify(cr_map_97.5, count=TRUE))[[1]]

# Write rasters as TIFF files to directory
setwd("../CholeraRiskProject/cholera_risk_mapping")
writeRaster(cr_map_mean, filename="mean_occurrence_p.tif", format="GTiff", overwrite=TRUE)
writeRaster(cr_map_2.5, filename="lower_occurrence_p.tif", format="GTiff", overwrite=TRUE)
writeRaster(cr_map_97.5, filename="upper_occurrence_p.tif", format="GTiff", overwrite=TRUE)

# Write rasters to ArcGIS directly 
#write.loc <- file.path('cholera_risk_mapping.gdb', 'cr_map_mean')
#arc.write(write.loc, cr_map_mean)
#write.loc <- file.path('cholera_risk_mapping.gdb', 'cr_map_2.5')
#arc.write(write.loc, cr_map_2.5)
#write.loc <- file.path('cholera_risk_mapping.gdb', 'cr_map_97.5')
#arc.write(write.loc, cr_map_97.5)


## 2050 Map Creation
mean_cells2 <- matrix(data=NA, ncol=3, nrow=5959968)  
lower_cells2 <- matrix(data=NA, ncol=3, nrow=5959968)
upper_cells2 <- matrix(data=NA, ncol=3, nrow=5959968)

for(i in 1:5)
{
  fr.probs <- matrix(nrow=1000, ncol=length((1+(i-1)*1e6):(i*1e6)))
  mean_cellvals2 <- matrix(ncol=3, nrow=ncol(fr.probs)) 
  lower_cellvals2 <- matrix(ncol=3, nrow=ncol(fr.probs))
  upper_cellvals2 <- matrix(ncol=3, nrow=ncol(fr.probs))
  # Create truncated matrix of cells (cols) with values from RF model's 1000 iterations
  for(j in 1:1000)
  {
    setwd("../CholeraRiskProject/rfpredictions/rfpredictions")
    fr.probs[j,] <- readRDS(paste('Replicate',j,'/2050preds.rds',sep=''))[,(1+(i-1)*1e6):(i*1e6)]
    #setwd('..')
  }
  # Find mean and percentiles for each cell, y (col of cr.probs)
  y <- 1
  for (x in (1+(i-1)*1e6):(i*1e6))
  {
    mean_cellvals2[y,3] <- mean(fr.probs[,y]) 
    mean_cellvals2[y,1:2] <- cell_coords[x,2:3]
    lower_cellvals2[y,3] <- quantile(fr.probs[,y], 0.025) 
    lower_cellvals2[y,1:2] <- cell_coords[x,2:3]
    upper_cellvals2[y,3] <- quantile(fr.probs[,y], 0.975)
    upper_cellvals2[y,1:2] <- cell_coords[x,2:3]
    y <- y+1
  }
  mean_cells2[(1+(i-1)*1e6):(i*1e6),] <- mean_cellvals2
  lower_cells2[(1+(i-1)*1e6):(i*1e6),] <- lower_cellvals2
  upper_cells2[(1+(i-1)*1e6):(i*1e6),] <- upper_cellvals2
}

# Complete the cellvals matrices with the remaining cells (<1e6)
fr.probs <- matrix(nrow=1000, ncol=length(c(5000001:5959968)))
mean_cellvals2 <- matrix(data=NA, ncol=3, nrow=ncol(fr.probs)) 
lower_cellvals2 <- matrix(data=NA, ncol=3, nrow=ncol(fr.probs))
upper_cellvals2 <- matrix(data=NA, ncol=3, nrow=ncol(fr.probs))

# Create truncated matrix of cells (cols) with values from RF model's 1000 iterations
for(j in 1:1000)
{
  setwd("../CholeraRiskProject/rfpredictions/rfpredictions")
  fr.probs[j,] <- readRDS(paste('Replicate',j,'/2050preds.rds',sep=''))[,c(5000001:5959968)]
  #setwd('..')
}

# Find mean and percentiles for each cell, x (col of fr.probs)
y <- 1
for(x in 5000001:5959968)
{
  mean_cellvals2[y,3] <- mean(fr.probs[,y]) 
  mean_cellvals2[y, 1:2] <- cell_coords[x,2:3]
  lower_cellvals2[y,3] <- quantile(fr.probs[,y], 0.025) 
  lower_cellvals2[y,1:2] <- cell_coords[x,2:3]
  upper_cellvals2[y,3] <- quantile(fr.probs[,y], 0.975)
  upper_cellvals2[y,1:2] <- cell_coords[x,2:3]
  y <- y+1
}
mean_cells2[5000001:5959968,] <- mean_cellvals2
lower_cells2[5000001:5959968,] <- lower_cellvals2
upper_cells2[5000001:5959968,] <- upper_cellvals2

saveRDS(mean_cells2, 'mean_cells_2050.rds')
saveRDS(lower_cells2, 'lower_cells_2050.rds')
saveRDS(upper_cells2, 'upper_cells_2050.rds')

# Create SpatialPoints object from coords, then rasterize to same spatial extent as ftr_rasters$elev with given cell values
mean_pts <- SpatialPoints(mean_cells2[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
mean_OP_r <- rasterize(mean_pts, ftr_rasters$elev, field=mean_cells2[,3])
lower_extreme_pts <- SpatialPoints(lower_cells2[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
lower_extreme_OP_r <- rasterize(lower_extreme_pts, ftr_rasters$elev, field=lower_cells2)
upper_extreme_pts <- SpatialPoints(upper_cells2[,1:2], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
upper_extreme_OP_r <- rasterize(upper_extreme_pts, ftr_rasters$elev, field=upper_cells2[,3])

# Set NAs <- 1 so color ramp is scaled 0 to 1
mean_OP_r[is.na(mean_OP_r)]<- 1
lower_extreme_OP_r[is.na(lower_extreme_OP_r)] <- 1
upper_extreme_OP_r[is.na(upper_extreme_OP_r)] <- 1

# Mask rasters to extent of Bangladesh, BGD, removing background values
BGD <- getData("GADM", country="BGD", level=0)
fr_map_mean <- mask(mean_OP_r, BGD, inverse=FALSE)
fr_map_2.5 <- mask(lower_extreme_OP_r, BGD, inverse=FALSE)
fr_map_97.5 <- mask(upper_extreme_OP_r, BGD, inverse=FALSE)

# Ratify rasters to create Raster Attribute Table (RAT)
mean2050_RAT <- levels(ratify(fr_map_mean, count=TRUE))[[1]]
q2.5_2050_RAT <- levels(ratify(fr_map_2.5, count=TRUE))[[1]]
q97.5_2050_RAT <- levels(ratify(fr_map_97.5, count=TRUE))[[1]]

setwd("../CholeraRiskProject/cholera_risk_mapping")
writeRaster(fr_map_mean, filename="mean_occurrence_p_2050.tif", format="GTiff", overwrite=TRUE)
writeRaster(fr_map_2.5, filename="lower_occurrence_p_2050.tif", format="GTiff", overwrite=TRUE)
writeRaster(fr_map_97.5, filename="upper_occurrence_p_2050.tif", format="GTiff", overwrite=TRUE)

# Write rasters to ArcGIS directly
#write.loc <- file.path('cholera_risk_mapping.gdb', 'fr_map_mean')
#arc.write(write.loc, fr_map_mean)
#write.loc <- file.path('cholera_risk_mapping.gdb', 'fr_map_2.5')
#arc.write(write.loc, fr_map_2.5)
#write.loc <- file.path('cholera_risk_mapping.gdb', 'fr_map_97.5')
#arc.write(write.loc, fr_map_97.5)

## Plotting Risk Diff Rasters w New Color Ramp
BuWhRd <- colorRampPalette(c("Blue", "White", "Red"))
BuWhRd_cont <- BuWhRd(100)

risk_diff_mean <- fr_map_mean - cr_map_mean
risk_diff_2.5 <- fr_map_2.5 - cr_map_2.5
risk_diff_97.5 <- fr_map_97.5 - cr_map_97.5

raster::plot(risk_diff_mean, main="Change in Infection Risk from 2015 to 2050\n Average Risk", col=BuWhRd_cont)
raster::plot(risk_diff_2.5, main="Change in Infection Risk from 2015 to 2050\n Quantile 2.5 Risk", col=BuWhRd_cont)
raster::plot(risk_diff_97.5, main="Change in Infection Risk from 2015 to 2050\n Quantile 97.5 Risk", col=BuWhRd_cont)

saveRDS(risk_diff_mean, "../CholeraRiskProject/risk_diff_mean.rds")
saveRDS(risk_diff_2.5, "../CholeraRiskProject/risk_diff_2.5.rds")
saveRDS(risk_diff_97.5, "../CholeraRiskProject/risk_diff_97.5.rds")

setwd("../CholeraRiskProject/cholera_risk_mapping")
writeRaster(risk_diff_mean, filename="risk_diff_mean.tif", format="GTiff", overwrite=TRUE)
writeRaster(risk_diff_2.5, filename="risk_diff_2.5.tif", format="GTiff", overwrite=TRUE)
writeRaster(risk_diff_97.5, filename="risk_diff_97.5.tif", format="GTiff", overwrite=TRUE)


## Plot rasters highlighting OP values >= 0.50
#To export rasters to ArcGIS for cartography
tmpfilter <- cr_map_mean >= 0.50
tmpfilter2 <- fr_map_mean >= 0.50
highrisk_2015 <- mask(cr_map_mean, tmpfilter, maskvalue=0)
highrisk_2050 <- mask(fr_map_mean, tmpfilter2, maskvalue=0)
write.loc <- file.path('cholera_risk_mapping.gdb', 'highrisk_2015')
arc.write(write.loc, highrisk_2015)
write.loc <- file.path('cholera_risk_mapping.gdb', 'highrisk_2050')
arc.write(write.loc, highrisk_2050)
#To plot with breaks highlighting OP values >=0.5 as red
districts <- getData("GADM", country="BGD", level=1)
breaks <- c(0.5, 1)
plot(cr_map_mean, breaks, col=c("white", "red"), xlab= "Longitude", ylab= "Latitude", main="Cell with an Occurrence Probability >= 0.50\n Average, 2015")
lines(districts)
plot(fr_map_mean, breaks, col=c("white", "red"), xlab= "Longitude", ylab= "Latitude", main="Cell with an Occurrence Probability >= 0.50\n Average, 2050")
lines(districts)


## Add segments & letters to identify district names on risk diff & 0.50> figures

# Use the following functions to select coordinates w cursor on the plot
loc1 <- function()
{
  inputs <- locator(1)
  cat(paste("text(",round(inputs$x[1],4),",",round(inputs$y[1],4),"\n",sep=""))
}
seg2 <- function()
{
  inputs <- locator(2)
  cat(paste("segments(",round(inputs$x[1],4),",",round(inputs$y[1],4),",",round(inputs$x[2],4),",",round(inputs$y[2],4),")\n",sep=""))
}

raster::plot(risk_diff_mean, xlab= "Longitude", ylab= "Latitude", main="Change in Infection Risk from 2015 to 2050\n Average Risk", col=BuWhRd_cont)
segments(90.40799,23.707016, 90.7458,25.7378)
text(90.7864,25.8494, "Dhaka, C")
text(90.40799,23.707016, "*") #Marks Dhaka as capital
segments(88.9999,25.8246,88.0389,26.2956)
text(87.9171,26.3947, "A")
segments(88.8916,24.6719,88.1472,25.0685)
text(88.0389,25.1429, "B")
segments(91.7609,24.7091,92.2888,25.5395)
text(92.3971,25.6635, "D")
segments(89.1893,22.9367,88.4179,22.6764)
text(88.3231,22.6144, "E")
segments(90.3263,22.4285,90.2992,21.4122)
text(90.2992,21.2635, "F")
segments(91.6121,22.8623,91.788,23.7299)
text(91.8828,23.8291, "G")

raster::plot(risk_diff_2.5, xlab= "Longitude", ylab= "Latitude", main="Change in Infection Risk from 2015 to 2050\n Quantile 2.5 Risk", col=BuWhRd_cont)
segments(90.40799,23.707016, 90.7458,25.7378)
text(90.7864,25.8494, "Dhaka, C")
text(90.40799,23.707016, "*") #Marks Dhaka as capital
segments(88.9999,25.8246,88.0389,26.2956)
text(87.9171,26.3947, "A")
segments(88.8916,24.6719,88.1472,25.0685)
text(88.0389,25.1429, "B")
segments(91.7609,24.7091,92.2888,25.5395)
text(92.3971,25.6635, "D")
segments(89.1893,22.9367,88.4179,22.6764)
text(88.3231,22.6144, "E")
segments(90.3263,22.4285,90.2992,21.4122)
text(90.2992,21.2635, "F")
segments(91.6121,22.8623,91.788,23.7299)
text(91.8828,23.8291, "G")

raster::plot(risk_diff_97.5, xlab= "Longitude", ylab= "Latitude", main="Change in Infection Risk from 2015 to 2050\n Quantile 97.5 Risk", col=BuWhRd_cont)
segments(90.40799,23.707016, 90.7458,25.7378)
text(90.7864,25.8494, "Dhaka, C")
text(90.40799,23.707016, "*") #Marks Dhaka as capital
segments(88.9999,25.8246,88.0389,26.2956)
text(87.9171,26.3947, "A")
segments(88.8916,24.6719,88.1472,25.0685)
text(88.0389,25.1429, "B")
segments(91.7609,24.7091,92.2888,25.5395)
text(92.3971,25.6635, "D")
segments(89.1893,22.9367,88.4179,22.6764)
text(88.3231,22.6144, "E")
segments(90.3263,22.4285,90.2992,21.4122)
text(90.2992,21.2635, "F")
segments(91.6121,22.8623,91.788,23.7299)
text(91.8828,23.8291, "G")

plot(cr_map_mean, breaks, col=c("white", "red"), xlab= "Longitude", ylab= "Latitude", main="Cell with an Occurrence Probability >= 0.50\n Average, 2015")
lines(districts)
segments(90.40799,23.707016, 90.7458,25.7378)
text(90.7864,25.8494, "Dhaka, C")
text(90.40799,23.707016, "*") #Marks Dhaka as capital
segments(88.9999,25.8246,88.0389,26.2956)
text(87.9171,26.3947, "A")
segments(88.8916,24.6719,88.1472,25.0685)
text(88.0389,25.1429, "B")
segments(91.7609,24.7091,92.2888,25.5395)
text(92.3971,25.6635, "D")
segments(89.1893,22.9367,88.4179,22.6764)
text(88.3231,22.6144, "E")
segments(90.3263,22.4285,90.2992,21.4122)
text(90.2992,21.2635, "F")
segments(91.6121,22.8623,91.788,23.7299)
text(91.8828,23.8291, "G")

plot(fr_map_mean, breaks, col=c("white", "red"), xlab= "Longitude", ylab= "Latitude", main="Cell with an Occurrence Probability >= 0.50\n Average, 2050")
lines(districts)
segments(90.40799,23.707016, 90.7458,25.7378)
text(90.7864,25.8494, "Dhaka, C")
text(90.40799,23.707016, "*") #Marks Dhaka as capital
segments(88.9999,25.8246,88.0389,26.2956)
text(87.9171,26.3947, "A")
segments(88.8916,24.6719,88.1472,25.0685)
text(88.0389,25.1429, "B")
segments(91.7609,24.7091,92.2888,25.5395)
text(92.3971,25.6635, "D")
segments(89.1893,22.9367,88.4179,22.6764)
text(88.3231,22.6144, "E")
segments(90.3263,22.4285,90.2992,21.4122)
text(90.2992,21.2635, "F")
segments(91.6121,22.8623,91.788,23.7299)
text(91.8828,23.8291, "G")
