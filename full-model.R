##### Random Forest Method to Generate 'Full Model' #####
library(randomForest)
library(dismo) 
library(pROC) 

setwd("../CholeraRiskProject")
source('GenerateRasterValues.R')
source('Ftr_RasterManipulation.R')

## Create data frame of lon, lat, environmental predictors, and presence values
serosurvey_mat <- readRDS('serosurvey-dat.rds') #cols 17-24 of serosurvey_mat are categorical 
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
envtrain <- na.omit(envtrain)

## Create test data to use as "new data" in predicting the fitted RF model
testset <- rbind(pv_ts[,3:25], av_ts[,3:25])
testset <- na.omit(testset)

## Runnning the Full Model ##

## Creating Confidence Intervals for the Variable Importance and AUC of Full Model (1000 runs)
varImport <- matrix(nrow=1000, ncol=10) 
colnames(varImport) <- c("landcover", "elev", "dist_to_water", "pop_dens", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m")
aucs <- c()

for(i in 1:1000)
{
  rf_env <- randomForest(x=envtrain[,11:20], y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred <- predict(rf_env, newdata=testset[,10:19], type="prob") #predict to testing data, to create ROC curve out of matrix of class probabilities
  rf_env_pred_ <- predict(rf_env, newdata=testset[,10:19], type="response") #same process, just identifying predicted classes (0s and 1s)
  test_preds_roc <- roc(testset[,20], rf_env_pred[,2]) #create ROC curve
  aucs[i] <- auc(test_preds_roc) - 0.50
  for(x in 1:10)
  {
    varImport[i,x] <- rf_env$importance[x,]
  }
}

auc_confint <- matrix(data = NA, nrow=1, ncol=3)
var_import <- matrix(nrow=10, ncol=4)
auc_confint[1,1] <- "Full Model"
auc_confint[1,2] <- quantile(aucs, 0.025)
auc_confint[1,3] <- quantile(aucs, 0.975)
saveRDS(auc_confint, "FullModel_AUCConfInt.rds")

var_import[,1] <- c("landcover", "elev", "dist_to_water", "pop_dens", "avg_ppt", "avg_tmax", "avg_tmin", "avg_ppt_m", "avg_tmax_m", "avg_tmin_m")
colnames(var_import) <- c("variable", "2.5%", "97.5%", "mean")
for(y in 1:10)
{
  var_import[y,2] <- quantile(varImport[,y], 0.025)
  var_import[y,3] <- quantile(varImport[,y], 0.975)
  var_import[y,4] <- mean(varImport[,y])
  
}
saveRDS(var_import, "FullModel_VarImport.rds")

## Create table of conf. ints of full model's variables' importance
means <- matrix(ncol=2, nrow=10000)
names <- colnames(varImport)
for(i in 1:10)
{
  means[((1+(i-1)*1000):(i*1000)), 1] <- as.numeric(varImport[,i])
  means[((1+(i-1)*1000):(i*1000)), 2] <- as.character(names[i])
}
colnames(means) <- c("variable_importance", "variable")
means <- as.data.frame(means)


