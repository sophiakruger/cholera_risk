## To Generate 1000 Iterations of the Best-fitting Model on a Supercomputer ##

library(randomForest)
library(dismo) 
library(pROC)
library(raster)

#env_rasters <- readRDS("env_rasters.rds")
env_rasters <- readRDS("/../rfpredictions/env_rasters.rds")
pt_rasters <- readRDS("/../rfpredictions/pt_rasters.rds")
ptm_rasters <- readRDS("/../rfpredictions/ptm_rasters.rds")
ftr_rasters <- readRDS("/../rfpredictions/ftr_rasters.rds")
env_data <- readRDS("/../rfpredictions/env_data.rds")
covariates <- readRDS("/../rfpredictions/covariates.rds")

serosurvey_mat <- readRDS("/../rfpredictions/serosurvey-dat.rds") #cols 17-24 of serosurvey_mat are categorical 
covariates_rf <- as.data.frame(cbind(env_data[,1:2], serosurvey_mat[,16:24], env_data[,4:16], env_data[,20])) #each row = one individual surveyed, full sero matrix combined with covariates of interest
colnames(covariates_rf) <- c("lon", "lat", "age", "sex", "elec", "ownerHH", "land", "educationSumm", "incomeSumm", "travel", "urban", "landcover","elev", "dist_to_water", "drp", "tss", "pop_dens", "dist_to_coast", "avg_ppt", "avg_tmax", "avg_tmin","avg_ppt_m", "avg_tmax_m", "avg_tmin_m", "presence")

## Characterize necessary variables as factors/categorical
cat.vars <- c("sex","elec", "ownerHH", "land", "educationSumm", "incomeSumm", "travel", "urban", "landcover")
covariates_rf[cat.vars]<-lapply(covariates_rf[cat.vars], factor)

## Create presence and absence matrices of environmental conditions
pres_vals <- covariates_rf[covariates_rf$presence==1, ]
abs_vals <- covariates_rf[covariates_rf$presence==0, ]

## Withold samples for training data and testing data of presence and absence data
set.seed(2)

fld <- kfold(pres_vals, k=5)
pres_tr <- pres_vals[fld!=1, ]
pres_ts <- pres_vals[fld==1, ]

fld <- kfold(abs_vals, k=5)
abs_tr <- abs_vals[fld!=1, ]
abs_ts <- abs_vals[fld==1, ]

## Partition the presence/absence environmental matrices into training and testing datasets 
pv_tr <- pres_vals[row.names(pres_tr),] #select rows of pres vals that are the rows taken for training
pv_ts <- pres_vals[row.names(pres_ts),]
av_tr <- abs_vals[row.names(abs_tr),]
av_ts <- abs_vals[row.names(abs_ts),]

## Vector response vector, where presence training points are 1 and absence training points are 0
pa_tr <- c(rep(1, nrow(pv_tr)), rep(0, nrow(av_tr)))

## Create data frame of presence and absence training points' environmental values 
envtrain <- as.data.frame(rbind(pv_tr, av_tr))
envtrain <- data.frame(cbind(pa_tr, envtrain))
envtrain <- envtrain[,-c(2,3,26)] #remove lon, lat, and presence columns
testset <- rbind(pv_ts[,3:25], av_ts[,3:25])
testset <- na.omit(testset)
envtrain <- na.omit(envtrain)

## Stack spatial covariates
#covariates <- stack(env_rasters[[1]],env_rasters[[2]],env_rasters[[3]], env_rasters[[4]], env_rasters[[5]], env_rasters[[6]], env_rasters[[7]], pt_rasters[[1]],pt_rasters[[2]],pt_rasters[[3]], ptm_rasters[[1]], ptm_rasters[[2]], ptm_rasters[[3]])
names(covariates) <- c("landcover", "elev", "dist_to_water", "drp", "tss", "pop_dens", "dist_to_coast", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m")  

## Create matrices for cell predictions
cellnumbers <- cellsFromExtent(covariates$elev, extent=extent(covariates$elev))
cell_coords <- cbind(cellnumbers, xyFromCell(covariates$elev, cell=cellnumbers))

set.seed(3)

cr.probs <- matrix(data=NA, nrow=1, ncol=length(cellnumbers))
fr.probs <- matrix(data=NA, nrow=1, ncol=length(cellnumbers))
train <- cbind(envtrain[,c(12:13)])
test <- cbind(testset[,c(11:12)])
rf_model <- randomForest(x=train, y=factor(envtrain[,1]), keep.forest=TRUE)
rf_test_pred <- predict(rf_model, newdata=test, type="prob")
rf_pred_cr <- raster::predict(stack(covariates$elev, covariates$dist_to_water), rf_model, type="prob", index=2, na.rm=TRUE)
rs <- stack(na.omit(ftr_rasters$elev), na.omit(ftr_rasters$dist_to_water))
names(rs) <- c("elev", "dist_to_water")
rf_pred_fr <- raster::predict(object=rs, model=rf_model, type="prob", index=2, na.rm=TRUE)
print(i)
for(x in 1:10)
{
  cr.probs[1,x] <- extract(rf_pred_cr, cellnumbers[x])
  fr.probs[1,x] <- extract(rf_pred_fr, cellnumbers[x])
}
saveRDS(cr.probs, paste("2015preds",i,".rds",sep=''))
saveRDS(fr.probs, paste("2050preds",i,".rds",sep=''))




