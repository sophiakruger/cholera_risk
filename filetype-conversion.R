## Convert NC spatial files to Tiff files for use in ArcGIS ##

install.packages("ncdf4")
library(raster)
library(ncdf4)

## Converting .nc files to geoTiff files for Precipitation and Temperature TERRA Data
setwd("../CholeraRiskProject")

# Review nc files
precip_acc15 <- nc_open('TerraClimate_ppt_2015.nc')
temp_max15 <- nc_open('TerraClimate_tmax_2015.nc')
temp_min15 <- nc_open('TerraClimate_tmin_2015.nc')

precip_acc16 <- nc_open('TerraClimate_ppt_2016.nc')
temp_max16 <- nc_open('TerraClimate_tmax_2016.nc')
temp_min16 <- nc_open('TerraClimate_tmin_2016.nc')

# Create RasterBrick for desired variables with extent to BGD
ppt15 <- brick('TerraClimate_ppt_2015.nc', varname="ppt", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)
tmax15 <- brick('TerraClimate_tmax_2015.nc', varname="tmax", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)
tmin15 <- brick('TerraClimate_tmin_2015.nc', varname="tmin", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)
ppt16 <- brick('TerraClimate_ppt_2016.nc', varname="ppt", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)
tmax16 <- brick('TerraClimate_tmax_2016.nc', varname="tmax", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)
tmin16 <- brick('TerraClimate_tmin_2016.nc', varname="tmin", xmn=88.01057, xmx=92.67366, ymn=20.74111, ymx=26.6344)

# Create Tif files from RasterBricks, save to ArcGIS Project Folder
setwd("../CholeraRiskProject/cholera_risk_mapping")
writeRaster(x=ppt15, filename='precip_acc_2015.tif', format='GTiff', overwrite=TRUE)
writeRaster(x=tmax15, filename='temp_max_2015.tif', format='GTiff', overwrite=TRUE)
writeRaster(x=tmin15, filename='temp_min_2015.tif', format='GTiff', overwrite=TRUE)
writeRaster(x=ppt16, filename='precip_acc_2016.tif', format='GTiff', overwrite=TRUE)
writeRaster(x=tmax16, filename='temp_max_2016.tif', format='GTiff', overwrite=TRUE)
writeRaster(x=tmin16, filename='temp_min_2016.tif', format='GTiff', overwrite=TRUE)





