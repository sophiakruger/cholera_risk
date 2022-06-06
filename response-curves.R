## Creating Response Curve Plot for 1000 iterations of the best-fitting RF model (elev + d2w)
library(raster)

# Retrieve 2015 covariates and the map of occurrence probability values for 2015 predicted by the RF model
covariates <- readRDS("../CholeraRiskProject/cholera_risk_mapping/covariates.rds")
cr_map_mean_preds <- readRDS("../CholeraRiskProject/cr_map_mean.rds")

# Ensure that the cooordinates for each cell # of elev, d2w, and 2015 preds are the same, ensuring that the values we extract from each cell # will be for the same grid cell 
xyFromCell(covariates$elev, cell=1) == xyFromCell(covariates$dist_to_water, cell=1)
xyFromCell(covariates$elev, cell=1) == xyFromCell(cr_map_mean_preds, cell=1)

xyFromCell(covariates$elev, cell=2979984) == xyFromCell(covariates$dist_to_water, cell=2979984)
xyFromCell(covariates$elev, cell=2979984) == xyFromCell(cr_map_mean_preds, cell=2979984)

xyFromCell(covariates$elev, cell=5959968) == xyFromCell(covariates$dist_to_water, cell=5959968)
xyFromCell(covariates$elev, cell=5959968) == xyFromCell(cr_map_mean_preds, cell=5959968)

# Create vectors containing elev & dist_to_water covariates' values for 2015 and the occurrence probabilities (prediction values) for 2015
elev_2015 <- as.vector(covariates$elev@data@values)
d2w_2015 <- as.vector(covariates$dist_to_water@data@values)
preds_2015 <- as.vector(cr_map_mean_preds@data@values)

# Plot response curves w/ each covariate's values along x-axis, predictions along y-axis
plot(x=elev_2015, y=preds_2015, xlab= "Elevation (m)", ylab="Occurence Probability", main="Predicted Cholera Occurrence Probability\n Given Elevation, 2015")
plot(x=(d2w_2015*250), y=preds_2015, xlab= "Distance to Water (m)", ylab="Occurence Probability", main="Predicted Cholera Occurrence Probability\n Given Distance to Water, 2015") #d2w orginal values are pixels, multiple by 250 to get meters

