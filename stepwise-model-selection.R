##### Stepwise Model Selection #####

## Creating confidence intervals for partial models to determine if variable offers significant gain in AUC to justify its addition to the model
auc_comp <- matrix(ncol=3, nrow=40, data=NA)

#elev only model
aucs <- c()
aucs1<- c()
for(i in 1:1000)
{
  envtrain31 <- cbind(envtrain[,c(12)]) 
  testset31 <- cbind(testset[,c(11)])
  rf_env31 <- randomForest(x=envtrain31, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred31 <- predict(rf_env31, newdata=testset31, type="prob")
  rf_env_pred31_ <- predict(rf_env31, newdata=testset31, type="response") 
  test_preds_roc31 <- roc(testset[,23], rf_env_pred31[,2]) 
  aucs1[i] <- auc(test_preds_roc31)
  aucs[i] <- auc(test_preds_roc31) - 0.50
}
mean_elev_auc <- mean(aucs1)
auc_comp[1,1] <- "elev"
auc_comp[1,2] <- quantile(aucs, 0.025)
auc_comp[1,3] <- quantile(aucs, 0.975)

## elev+d2w
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{
  #set.seed(21)
  envtrain18 <- cbind(envtrain[,c(12:13)]) 
  testset18 <- cbind(testset[,c(11:12)])
  rf_env18 <- randomForest(x=envtrain18, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred18 <- predict(rf_env18, newdata=testset18, type="prob")
  rf_env_pred18_ <- predict(rf_env18, newdata=testset18, type="response") 
  test_preds_roc18 <- roc(testset[,23], rf_env_pred18[,2]) 
  aucs1[i] <- auc(test_preds_roc18)
  aucs[i] <- auc(test_preds_roc18) - mean_elev_auc #get the increase from elev model's mean AUC #AUC=0.6135
}
mean_edw_auc <- mean(aucs1)
auc_comp[2,1] <- "elev+d2w"
auc_comp[2,2] <- quantile(aucs, 0.025)
auc_comp[2,3] <- quantile(aucs, 0.975)

## d2w only
aucs <- c()
aucs1<- c()
for(i in 1:1000)
{
  #set.seed(31)
  envtrain31 <- cbind(envtrain[,c(13)]) 
  testset31 <- cbind(testset[,c(12)])
  rf_env31 <- randomForest(x=envtrain31, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred31 <- predict(rf_env31, newdata=testset31, type="prob")
  rf_env_pred31_ <- predict(rf_env31, newdata=testset31, type="response") 
  test_preds_roc31 <- roc(testset[,23], rf_env_pred31[,2]) 
  aucs1[i] <- auc(test_preds_roc31)
  aucs[i] <- auc(test_preds_roc31) - 0.50
}
mean_d2w_auc <- mean(aucs1)
auc_comp[3,1] <- "d2w"
auc_comp[3,2] <- quantile(aucs, 0.025)
auc_comp[3,3] <- quantile(aucs, 0.975)

## elev+d2W+pop_dens
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{
  #set.seed(6)
  envtrain3 <- cbind(envtrain[,c(12:13, 16)]) 
  testset3 <- cbind(testset[,c(11:12, 15)])
  rf_env3 <- randomForest(x=envtrain3, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred3 <- predict(rf_env3, newdata=testset3, type="prob")
  rf_env_pred3_ <- predict(rf_env3, newdata=testset3, type="response")  
  test_preds_roc3 <- roc(testset[,23], rf_env_pred3[,2])
  aucs1[i] <- auc(test_preds_roc3) 
  aucs[i] <- auc(test_preds_roc3) - mean_edw_auc
}
mean_edwp_auc <- mean(aucs1)
auc_comp[4,1] <- "elev+d2W + popdens"
auc_comp[4,2] <- quantile(aucs, 0.025)
auc_comp[4,3] <- quantile(aucs, 0.975)

## elev+d2W+landcover
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain32 <- cbind(envtrain[,c(11, 12:13)]) 
  testset32 <- cbind(testset[,c(10, 11:13)])
  rf_env32 <- randomForest(x=envtrain32, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred32 <- predict(rf_env32, newdata=testset32, type="prob")
  rf_env_pred32_ <- predict(rf_env32, newdata=testset32, type="response")  
  test_preds_roc32 <- roc(testset[,23], rf_env_pred32[,2])
  aucs1[i] <- auc(test_preds_roc32) 
  aucs[i] <- auc(test_preds_roc32) - mean_edw_auc
}
mean_edwl_auc <- mean(aucs1)
auc_comp[5,1] <- "elev+d2w+landcover"
auc_comp[5,2] <- quantile(aucs, 0.025)
auc_comp[5,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_tmin_m
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain33 <- cbind(envtrain[,c(12:13, 23)]) 
  testset33 <- cbind(testset[,c(11:12, 22)])
  rf_env33 <- randomForest(x=envtrain33, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred33 <- predict(rf_env33, newdata=testset33, type="prob")
  rf_env_pred33_ <- predict(rf_env33, newdata=testset33, type="response")  
  test_preds_roc33 <- roc(testset[,23], rf_env_pred33[,2])
  aucs1[i] <- auc(test_preds_roc33) 
  aucs[i] <- auc(test_preds_roc33) - mean_edw_auc
}
mean_edwtmnm_auc <- mean(aucs1)
auc_comp[6,1] <- "elev+d2w+avg_tmin_m"
auc_comp[6,2] <- quantile(aucs, 0.025)
auc_comp[6,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_ppt
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain34 <- cbind(envtrain[,c(12:13, 18)]) 
  testset34 <- cbind(testset[,c(11:12, 17)])
  rf_env34 <- randomForest(x=envtrain34, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred34 <- predict(rf_env34, newdata=testset34, type="prob")
  rf_env_pred34_ <- predict(rf_env34, newdata=testset34, type="response")  
  test_preds_roc34 <- roc(testset[,23], rf_env_pred34[,2])
  aucs1[i] <- auc(test_preds_roc34) 
  aucs[i] <- auc(test_preds_roc34) - mean_edw_auc
}
mean_edwpt_auc <- mean(aucs1)
auc_comp[7,1] <- "elev+d2w+avg_ppt"
auc_comp[7,2] <- quantile(aucs, 0.025)
auc_comp[7,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_tmin
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain35 <- cbind(envtrain[,c(12:13, 20)]) 
  testset35 <- cbind(testset[,c(11:12, 19)])
  rf_env35 <- randomForest(x=envtrain35, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred35 <- predict(rf_env35, newdata=testset35, type="prob")
  rf_env_pred35_ <- predict(rf_env35, newdata=testset35, type="response")  
  test_preds_roc35 <- roc(testset[,23], rf_env_pred35[,2])
  aucs1[i] <- auc(test_preds_roc35) 
  aucs[i] <- auc(test_preds_roc35) - mean_edw_auc
}
mean_edwtmn_auc <- mean(aucs1)
auc_comp[8,1] <- "elev+d2w+avg_tmin"
auc_comp[8,2] <- quantile(aucs, 0.025)
auc_comp[8,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_tmax_m
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain36 <- cbind(envtrain[,c(12:13, 22)]) 
  testset36 <- cbind(testset[,c(11:12, 21)])
  rf_env36 <- randomForest(x=envtrain36, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred36 <- predict(rf_env36, newdata=testset36, type="prob")
  rf_env_pred36_ <- predict(rf_env36, newdata=testset36, type="response")  
  test_preds_roc36 <- roc(testset[,23], rf_env_pred36[,2])
  aucs1[i] <- auc(test_preds_roc36) 
  aucs[i] <- auc(test_preds_roc36) - mean_edw_auc
}
mean_edwtmxm_auc <- mean(aucs1)
auc_comp[9,1] <- "elev+d2w+avg_tmax_m"
auc_comp[9,2] <- quantile(aucs, 0.025)
auc_comp[9,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_ppt_m
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain37 <- cbind(envtrain[,c(12:13, 21)]) 
  testset37 <- cbind(testset[,c(11:12, 20)])
  rf_env37 <- randomForest(x=envtrain37, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred37 <- predict(rf_env37, newdata=testset37, type="prob")
  rf_env_pred37_ <- predict(rf_env37, newdata=testset37, type="response")  
  test_preds_roc37 <- roc(testset[,23], rf_env_pred37[,2])
  aucs1[i] <- auc(test_preds_roc37) 
  aucs[i] <- auc(test_preds_roc37) - mean_edw_auc
}
mean_edwpptm_auc <- mean(aucs1)
auc_comp[10,1] <- "elev+d2w+avg_ppt_m"
auc_comp[10,2] <- quantile(aucs, 0.025)
auc_comp[10,3] <- quantile(aucs, 0.975)

## elev+d2w+avg_tmax
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain38 <- cbind(envtrain[,c(12:13, 19)]) 
  testset38 <- cbind(testset[,c(11:12, 18)])
  rf_env38 <- randomForest(x=envtrain38, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred38 <- predict(rf_env38, newdata=testset38, type="prob")
  rf_env_pred38_ <- predict(rf_env38, newdata=testset38, type="response")  
  test_preds_roc38 <- roc(testset[,23], rf_env_pred38[,2])
  aucs1[i] <- auc(test_preds_roc38) 
  aucs[i] <- auc(test_preds_roc38) - mean_edw_auc
}
mean_edwtmx_auc <- mean(aucs1)
auc_comp[11,1] <- "elev+d2w+avg_tmax"
auc_comp[11,2] <- quantile(aucs, 0.025)
auc_comp[11,3] <- quantile(aucs, 0.975)

## elev+d2w+tss
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain39 <- cbind(envtrain[,c(12:13, 15)]) 
  testset39 <- cbind(testset[,c(11:12, 14)])
  rf_env39 <- randomForest(x=envtrain39, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred39 <- predict(rf_env39, newdata=testset39, type="prob")
  rf_env_pred39_ <- predict(rf_env39, newdata=testset39, type="response")  
  test_preds_roc39 <- roc(testset[,23], rf_env_pred39[,2])
  aucs1[i] <- auc(test_preds_roc39) 
  aucs[i] <- auc(test_preds_roc39) - mean_edw_auc
}
mean_edwtss_auc <- mean(aucs1)
auc_comp[12,1] <- "elev+d2w+tss"
auc_comp[12,2] <- quantile(aucs, 0.025)
auc_comp[12,3] <- quantile(aucs, 0.975)

auc_comp <- rbind(auc_comp, 13)

## elev+d2w+drp
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain40 <- cbind(envtrain[,c(12:13, 16)]) 
  testset40 <- cbind(testset[,c(11:12, 15)])
  rf_env40 <- randomForest(x=envtrain40, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred40 <- predict(rf_env40, newdata=testset40, type="prob")
  rf_env_pred40_ <- predict(rf_env40, newdata=testset40, type="response")  
  test_preds_roc40 <- roc(testset[,23], rf_env_pred40[,2])
  aucs1[i] <- auc(test_preds_roc40) 
  aucs[i] <- auc(test_preds_roc40) - mean_edw_auc
}
mean_edwdrp_auc <- mean(aucs1)
auc_comp[13,1] <- "elev+d2w+drp"
auc_comp[13,2] <- quantile(aucs, 0.025)
auc_comp[13,3] <- quantile(aucs, 0.975)

## d2w + pop_dens
aucs <- c()
aucs1 <- c()
for(i in 1:1000)
{ 
  envtrain41 <- cbind(envtrain[,c(13, 16)]) 
  testset41 <- cbind(testset[,c(12, 15 )])
  rf_env41 <- randomForest(x=envtrain41, y=factor(envtrain[,1]), keep.forest=TRUE)
  rf_env_pred41 <- predict(rf_env41, newdata=testset41, type="prob")
  rf_env_pred41_ <- predict(rf_env41, newdata=testset41, type="response")  
  test_preds_roc41 <- roc(testset[,23], rf_env_pred41[,2])
  aucs1[i] <- auc(test_preds_roc41) 
  aucs[i] <- auc(test_preds_roc41) - mean_elev_auc
}
mean_d2wp_auc <- mean(aucs1)
auc_comp[13,1] <- "d2w+popdens"
auc_comp[13,2] <- quantile(aucs, 0.025)
auc_comp[13,3] <- quantile(aucs, 0.975)

## AUC Comp: 'auc_comp' shows us elev+d2w model is best fitting
saveRDS(auc_comp, "../CholeraRiskProject/Final_Figures/AUCConfIntComparison.rds")
saveRDS(rf_env18, "../CholeraRiskProject/Final_Figures/elev+d2w_rfmodel.rds")
