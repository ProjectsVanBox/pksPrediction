### load data 
library(randomForest)

#### load random forest models 
load("forest1.Rdata")
load("forest2.Rdata")

### load data 
load("testDataScreeningCohort.Rdata")
clinInfo = read.csv("41586_2019_1672_MOESM4_ESM.csv")
rownames(clinInfo) = clinInfo$patient_label_in_files

#### predict pks status 

predict1 = predict(output.forest, test.data, type = "prob")[,2]
predict2 = predict(forest2, test.data, type = "prob")[,2]

finalProbNoM = predictNoM1*predictNoM2 #### probability of pks+



pksClass = rep("noPKS", nrow(test.data))

pksClass[which(finalProbNoM > 0.5)] = "PKS"


test.data$pks = pksClass
test.data$prob = finalProb
test.data = test.data.LS[which(!is.na(finalProb)),]

#### 
newNames = NULL
for(i in 1:nrow(test.data.LS)){
  newNames[i] = strsplit(rownames(test.data.LS)[i], "[.]")[[1]][1]
}

samples = unique(newNames)

ratios = NULL


for(i in 1:length(samples)){
  ratios[i] = table(test.data$pks[which(newNames== samples[i])])[2]/length(which(newNames == samples[i]))
}

names(ratios) = samples

pksStatus = rep("noPKS", length(samples))
pksStatus[which(ratios > 0.1)] = "PKS"

#### plot age & cancer status 
pdf("ageVsRatio.pdf")
plot(clinInfo$ratio[which(clinInfo$bowel_cancer_diagnosis == "No")],clinInfo$age[which(clinInfo$bowel_cancer_diagnosis == "No")], ylim =c(0,80), col = "azure3", pch =4, lwd =  5, ylab = "Age", xlab = "Fraction PKS+", cex.axis = 1.5, cex.lab = 1.5, xlim = c(0,0.55))
points(clinInfo$ratio[which(clinInfo$bowel_cancer_diagnosis == "Yes")],clinInfo$age[which(clinInfo$bowel_cancer_diagnosis == "Yes")], pch = 3, col = "aquamarine3", lwd = 5)
legend("bottomleft", c("No colorectal cancer", "Colorectal cancer"), col = c("azure3","aquamarine3"), pch = c(4,3), bty = "n", lwd = 5, lty = 0, cex = 1.5)
dev.off()

### check screening cohort
clinInfo = clinInfo[samples,]
clinInfo$ratio = ratios
screeningCohort = samples[which(clinInfo$cohort == "bowel_cancer_screening_programme_cohort")]

clinInfoScreening = clinInfo[screeningCohort,]

pdf("Figures/paper/Figure3/AgeAtDiagnosis.pdf", width = 4.5)
boxplot(clinInfoScreening$age[which(clinInfoScreening$ratio > 0.1 & clinInfoScreening$bowel_cancer_diagnosis == "Yes")], clinInfoScreening$age[which(clinInfoScreening$ratio <= 0.1 & clinInfoScreening$bowel_cancer_diagnosis == "Yes")], names = c("PKS+", "PKS-"),ylab = "Age at diagnosis", cex.lab = 1.5, cex.axis = 1.5, cex.names = 1.5)
dev.off()


#### check driver genes 
driverGenes = c("APC", "TP53", "KRAS", "BRAF", "PIK3CA", "SMAD4", "FBXW7", "TCF7L2","FAT4", "AMER1", "LRP1B", "SOX9", "ATM")
noDriver = NULL
for(i in 1:length(ratios)){
  patient = test.data[which(!is.na(match(newNames,names(ratios)[i]))),]
  noDriver[i] = length(which(!is.na(match(driverGenes, patient$GENE))))
}

names(noDriver) = names(ratios)

noDriverScreening = noDriver[screeningCohort]

clinInfo$ratio = ratiosLeeSix[rownames(clinInfo)]
pdf("noDriverVsRatio.pdf")
plot(clinInfo$ratio[which(clinInfo$bowel_cancer_diagnosis == "No")],noDriver[which(clinInfo$bowel_cancer_diagnosis == "No")], ylim =c(0,9), col = "azure3", pch =4, lwd =  5, ylab = "Number of drivers", xlab = "Fraction PKS+", cex.axis = 1.5, cex.lab = 1.5, xlim = c(0,0.55))
points(clinInfo$ratio[which(clinInfo$bowel_cancer_diagnosis == "Yes")],noDriver[which(clinInfo$bowel_cancer_diagnosis == "Yes")], pch = 3, col = "aquamarine3", lwd = 5)
legend("topright", c("No colorectal cancer", "Colorectal cancer"), col = c("azure3","aquamarine3"), pch = c(4,3), bty = "n", lwd = 5, lty = 0, cex = 1.5)
dev.off()