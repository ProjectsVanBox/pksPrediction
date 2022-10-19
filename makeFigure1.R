setwd("~/surfdrive/Projects/pksClassifier/Figures/paper/Figure1")

library(randomForest)
library(MutationalPatterns)
####1a 
####mutational pattern 
existing = get_known_signatures(muttype = "snv", incl_poss_artifacts = F)

####1b 
####performance both random forests
#### load models and test data 
load("trainData.Rdata")
load("forest1.Rdata")
load("forest2.Rdata")

#### 
pred1 = predict(output.forest, train.data, type = "prob")[,2]
pred2 = predict(forest2, train.data, type = "prob")[,2]
predTog = pred1*pred2

#### pool samples 
newNames = NULL

for(i in 1:nrow(train.data.high)){
  newNames[i] = strsplit(rownames(train.data.high)[i], "[.]")[[1]][1]
}

#### score per patient 
patientScore1 = NULL

for(i in 1:length(unique(newNames))){
  patient = unique(newNames)[i]
  patientScore1[i] = length(which(pred1[which(newNames == patient)] > 0.5))/length(which(newNames == patient))
}

names(patientScore1) = unique(newNames)


patientScore2 = NULL

for(i in 1:length(unique(newNames))){
  patient = unique(newNames)[i]
  patientScore2[i] = length(which(pred2[which(newNames == patient)] > 0.5))/length(which(newNames == patient))
}

names(patientScore2) = unique(newNames)


patientScore = NULL

for(i in 1:length(unique(newNames))){
  patient = unique(newNames)[i]
  patientScore3[i] = length(which(predTog[which(newNames == patient)] > 0.5))/length(which(newNames == patient))
}

names(patientScore3) = unique(newNames)

#### plot all scores together 
allScores = rbind(patientScore1, patientScore3)
pdf("barplot_blood_bothRF.pdf", width = 4.5)
par(mar = c(8,5,4,2))
barplot(allScores[,sort(patientScore1, index = T)$ix], beside = T, las =3, col = c(rep("azure3", 12), rep("aquamarine3",6)), ylab = "Fraction PKS positive", cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", c("Blood", "PKS+ tumor"), col = c("azure3", "aquamarine3"), lty = 0, pch = 15, bty = "n", cex  =1.5)
dev.off()

####1c
####feature importance both random forests 
FI = rbind(output.forest$importance[,3], forest2$importance[,3])
pdf("featureImportance_bothModels.pdf", width = 8)
par(mar = c(8,5,4,2))
barplot(FI, beside = T, las = 3, ylab = "Mean decrease in accuracy", col = c("gray81", "gray30"), cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", c("RF model 1", "RF model 2"), col = c("gray81", "gray30"), pch = 15, lty = 0, bty = "n", cex = 1.5)
dev.off()


#####1d
#####sequence motif TP
#### extended sequence context --> organoids and FP/TP from patients 
library(motifStack)
TP = train.data[which(train.data$VAL == "PKS" & predTog > 0.5),6:26]
context = matrix(NA,4,21)
for(i in 1:21){
  context[1,i] = length(which(TP[,i] == "A"))
  context[2,i] = length(which(TP[,i] == "C"))
  context[3,i] = length(which(TP[,i] == "G"))
  context[4,i] = length(which(TP[,i] == "T"))
  
}

rownames(context) <- c("A","C","G","T")
colnames(context) = colnames(TP)
motifTP <- new("pcm", mat=as.matrix(context), name="True Positives")


FP = train.data.high[which(train.data.high$VAL == "noPKS"  & pred1 > 0.5),6:26]
contextFP = matrix(NA,4,21)
for(i in 1:21){
  contextFP[1,i] = length(which(FP[,i] == "A"))
  contextFP[2,i] = length(which(FP[,i] == "C"))
  contextFP[3,i] = length(which(FP[,i] == "G"))
  contextFP[4,i] = length(which(FP[,i] == "T"))
  
}

rownames(contextFP) <- c("A","C","G","T")
colnames(contextFP) = colnames(TP)
motifFP <- new("pcm", mat=as.matrix(contextFP), name="False Positives")
plot(motifFP)

#### true motif 
#### load organoid data 
load("trainingData_allOrganoids.Rdata")
 
expP = train.data.SNV[which(train.data.SNV$VAL == "PKS"),6:26]

newNames = NULL

for(i in 1:nrow(expP)){
  newNames[i] = strsplit(rownames(expP)[i], "[.]")[[1]][1]
}


contextTrue = matrix(NA,4,21)
for(i in 1:21){
  contextTrue[1,i] = length(which(expP[,i] == "A"))
  contextTrue[2,i] = length(which(expP[,i] == "C"))
  contextTrue[3,i] = length(which(expP[,i] == "G"))
  contextTrue[4,i] = length(which(expP[,i] == "T"))
  
}

rownames(contextTrue) <- c("A","C","G","T")
colnames(contextTrue) = colnames(expP)
motifTrue <- new("pcm", mat=as.matrix(contextTrue), name="Organoid pattern")
plot(motifTrue)

pdf("MotifTogether.pdf", width = 12)

plot(motifTrue, cex.axis = 1.5, cex.lab = 1.5)
plot(motifTP, cex.axis = 1.5, cex.lab = 1.5)
plot(motifFP, cex.axis = 1.5, cex.lab = 1.5)
dev.off()


