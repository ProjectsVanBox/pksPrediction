library(randomForest)

#### load scores
load("crcMuts.Rdata") ### data frame as obtained from parseBedFiles and random forest posterior probabilities
colnames(crcMuts) = c("CHROM", "START", "END", "GENE", "REPLISEQ", "pred1", "pred2")
crcMuts$finalProb = as.numeric(crcMuts$pred1)*as.numeric(crcMuts$pred2)



newNames = NULL

for(i in 1:nrow(crcMuts)){
  newNames[i] = strsplit(rownames(crcMuts)[i], "[.]")[[1]][1]
}

#### score per sample 
score = NULL

for(i in 1:length(unique(newNames))){
  patient = unique(newNames)[i]
  score[i] = length(which(crcMuts$finalProb[which(newNames == patient)] > 0.5))/length(which(newNames == patient))
}

names(score) = unique(newNames)


#### calculate driver genes 


posPatients = crcMuts[which(!is.na(match(sampleNames,names(score)[which(score > 0.1)]))),]
negPatients = crcMuts[which(!is.na(match(sampleNames,names(scoreColon)[which(score <= 0.1)]))),]

driverGenes = c("APC", "TP53", "KRAS", "BRAF", "PIK3CA", "SMAD4", "FBXW7", "TCF7L2","FAT4", "ATM")

diffPerGene = matrix(NA, 2,length(driverGenes))
set.seed(220718)

pValuesDriver = rep(NA, length(driverGenes))
pValuesDriver_neg = rep(NA, length(driverGenes))
sampledPos = list()
sampledNeg = list()
probPos = list()
probNeg = list()
differencesPos = list()
differencesNeg = list()

#### get non-driver genes to sample
genesToSample = intersect(names(which(table(posPatients$GENE) > 0 & table(posPatients$GENE) < 4)),names(which(table(negPatients$GENE) >0  & table(negPatients$GENE) < 20)) )

randomGenes = sample(genesToSample, 1000, replace = F)
meanPos = rep(NA, length(randomGenes))
meanNeg = rep(NA, length(randomGenes))

for(s in 1:length(randomGenes)){
  
  meanPos[s] = mean(posPatients$finalProb[which(posPatients$GENE == randomGenes[s])], na.rm = T)
  meanNeg[s] = mean(negPatients$finalProb[which(negPatients$GENE == randomGenes[s])], na.rm = T)
  
}

diffRandom = meanPos - meanNeg
for(g in 1:length(driverGenes)){
    diffPos = rep(NA,length(which(posPatients$GENE == driverGenes[g])))
    diffNeg = rep(NA,length(which(negPatients$GENE == driverGenes[g])))
    
    
    
    meanPosUse = mean(meanPos, na.rm = T)
    meanNegUse = mean(meanNeg, na.rm = T)
    
    probPos[[g]] = posPatients$finalProb[which(posPatients$GENE == driverGenes[g])]
    probNeg[[g]] = negPatients$finalProb[which(negPatients$GENE == driverGenes[g])]
    
    for(i in 1:length(which(posPatients$GENE == driverGenes[g]))){
      diffPos[i] = posPatients$finalProb[which(posPatients$GENE == driverGenes[g])[i]] - meanPosUse
    }
    
    for(i in 1:length(which(negPatients$GENE == driverGenes[g]))){
      diffNeg[i] = negPatients$finalProb[which(negPatients$GENE == driverGenes[g])[i]] - meanNegUse
    }
    
  
    differencesPos[[g]] = diffPos
    differencesNeg[[g]] = diffNeg
    
    
    diffPerGene[1,g] = mean(differencesPos[[g]])
    diffPerGene[2,g] = mean(differencesNeg[[g]])
    
    
}
  
diffProb = NULL

for(i in 1:length(driverGenes)){
  diffProb[i] = mean(probPos[[i]]) - mean(probNeg[[i]])
}

#### 

pValues_t = rep(NA, length(driverGenes))
pValues_neg = rep(NA, length(driverGenes))


for(i in 1:length(driverGenes)){
    pValues_t[i] = t.test(posPatients$finalProb[which(posPatients$GENE == driverGenes[i])], meanPos)$p.value
    pValues_neg[i] = t.test(negPatients$finalProb[which(negPatients$GENE == driverGenes[i])], meanNeg)$p.value
    
}

names(pValues_t) = driverGenes

diffPerGeneT = matrix(NA, 2, length(driverGenes))

for(i in 1:length(driverGenes)){
  diffPerGeneT[1,i] = mean(diffPerGene[1,i], na.rm = T)
  diffPerGeneT[2,i] = mean(diffPerGene[2,i], na.rm = T)
}

diffPerGeneT_sd = matrix(NA, 2, length(driverGenes))

for(i in 1:length(driverGenes)){
  diffPerGeneT_sd[1,i] = sd(diffPerGene[1,i], na.rm = T)
  diffPerGeneT_sd[2,i] = sd(diffPerGene[2,i], na.rm = T)
}

#### plot 
#### boxplot 
library(ggplot2)

genes = NULL
label = NULL
prob = NULL

driverGenes_ordered = driverGenes[sort(diffProb, decreasing = T, index = T)$ix]
for(i in 1:length(driverGenes_ordered)){
  genes = c(genes, rep(driverGenes_ordered[i], length(which(posPatients$GENE == driverGenes_ordered[i]))), rep(driverGenes_ordered[i], length(which(negPatients$GENE == driverGenes_ordered[i]))))
  label = c(label, rep("PKS+", length(which(posPatients$GENE == driverGenes_ordered[i]))), rep("PKS-",length(which(negPatients$GENE == driverGenes_ordered[i]))))
  prob = c(prob, posPatients$finalProb[which(posPatients$GENE == driverGenes_ordered[i])], negPatients$finalProb[which(negPatients$GENE == driverGenes_ordered[i])])
  
}

genes = c(genes[1:1224], rep("Random", 2000), genes[1225:2032])
label = c(label[1:1224], rep("PKS+", 1000), rep("PKS-", 1000), label[1225:2032])
prob = c(prob[1:1224], meanPos, meanNeg, prob[1225:2032])

level_order = factor(c(driverGenes_ordered[1:6], "Random", driverGenes_ordered[7:10]))
 dataToPlot = data.frame(genes, label, prob)

 pdf("Boxplot_genes_prob.pdf", width = 12)
 ggplot(dataToPlot, aes(x=factor(genes, level = level_order), y=prob, fill=label)) + ylab("Posterior probability") + xlab(NA)+
   geom_boxplot() +  scale_fill_manual(values = c("aquamarine3", "azure3")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16, angle = 90))
 dev.off()
 
 #### plot APC and p53 probabilities 
 ### APC 
 
 pdf("probDistAPC.pdf", width = 6, height = 4)
 hist(crcMuts$finalProb[which(crcMuts$GENE == "APC")], main = NA, breaks = 100, xlim = c(0,0.9), xlab = "Posterior probability", cex.lab = 1.5, cex.axis = 1.5)
 dev.off()
 
 pdf("probDistP53.pdf", width = 6, height = 4)
 hist(crcMuts$finalProb[which(crcMuts$GENE == "TP53")], ylim = c(0,25), main = NA, breaks = 20, xlim = c(0,0.9), xlab = "Posterior probability", cex.lab = 1.5, cex.axis = 1.5)
 dev.off()
 
 
 #### mutational signature analysis 
#### mutation matrices as produced by mutational patterns function mut_matrix()
load("MutMatAPC.Rdata")
load("MutMatBRAF.Rdata")
load("MutMatSMAD4.Rdata")
load("MutMatFBXW7.Rdata")

library(MutationalPatterns)

existing = get_known_signatures(muttype = "snv", incl_poss_artifacts = F, genome = "GRCh37")

allMuts = cbind((allMutsFBXW7[,1] +allMutsSMAD4[,1] + allMutsAPC[,1] +allMutsBRAF[,1]),(allMutsFBXW7[,2] + allMutsSMAD4[,2] +allMutsAPC[,2] + allMutsBRAF[,2]))
colnames(allMuts) = c("PKS+", "PKS-")

checkContr = fit_to_signatures(allMuts, existing)

relativePosAPC = checkContr$contribution[,1]/sum(checkContr$contribution[,1]) 
relativeNegAPC = checkContr$contribution[,2]/sum(checkContr$contribution[,2]) 

diffAPC = relativePosAPC - relativeNegAPC
names(diffAPC) = rownames(checkContr$contribution)

#### plot profile 
pdf("96profile_drivers.pdf", width = 6, height = 4)
plot_96_profile(allMuts, ymax = 0.1)
dev.off()

#### plot bootstrapped contribution
sigsToUse = unique(c(names(which(relativeNegAPC > 0.02)),names(which(relativePosAPC > 0.02)) ))
existingUse = existing[,sigsToUse]

checkContr_sel = fit_to_signatures(allMuts, existingUse)

checkContr_bootstrapped = fit_to_signatures_bootstrapped(allMuts, existingUse)

pdf("Bootstrapped_contribution.pdf", width = 6, height = 4)
plot_bootstrapped_contribution(checkContr_bootstrapped, mode = "relative", plot_type = "jitter")
dev.off()

relativePosAPC_sel = checkContr_sel$contribution[,1]/sum(checkContr_sel$contribution[,1]) 
relativeNegAPC_sel = checkContr_sel$contribution[,2]/sum(checkContr_sel$contribution[,2]) 

diffAPC = relativePosAPC_sel - relativeNegAPC_sel
names(diffAPC) = rownames(checkContr_sel$contribution)

sigsToPlot = names(diffAPC)[which(abs(diffAPC) > 0.05)]

pdf("barplotSigDifferences.pdf", width =6, height = 4)
barplot(diffAPC[sigsToPlot[c(1,4,5,7,8,2,3,6)]], ylim = c(-0.10, 0.15),las = 3, ylab = "Difference relative contribution", cex.axis = 1.3, cex.lab = 1.3, cex.names = 1.3)
dev.off()
