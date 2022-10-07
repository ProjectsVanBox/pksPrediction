##### load random forest library
library(randomForest)


#### load data 
#### data from EGA archive, processed as described in processing scripts
load("trainingDataSNV_organoids.Rdata")

#### train first model 
set.seed(210209)

output.forest <- randomForest(as.factor(VAL) ~ ., na.action=na.roughfix, data=train.data.SNV, mtry = 3, ntree = 1000, localImp = T)

#### predict the patient data with the first random forest model 
train2 = train.data.SNV.patients
predictionTrain2 = predict(output.forest, train2, type = "prob")

### check false negatives 
predictedScore = predictionTrain2[,2]

falsePositives = which(predictedScore > 0.8 & train2$VAL == "PKS")
truePositives = which(predictedScore > 0.5 & train2$VAL == "noPKS")

trainNew = train2[c(falsePositives, truePositives),]

forest2 <- randomForest(as.factor(VAL) ~ ., na.action=na.roughfix, data=trainNew, mtry = 3, ntree = 1000, localImp = T)
