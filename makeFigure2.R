### score == fraction of mutations that is predicted to be PKS-positive
### tumorTypes == annotation per sample of diagnosed type of cancer

numberPositive = NULL

for(i in 1:length(unique(tumorTypes))){
  numberPositive[i] =length(which(score[which(tumorTypes == unique(tumorTypes)[i])] >= 0.1))/length(which(tumorTypes == unique(tumorTypes)[i]))
}

names(numberPositive) = unique(tumorTypes)

pdf("HMF_percpositive.pdf", width = 11)
par(mar = c(9,6,4,2))
barplot(sort(numberPositive, decreasing = T), las = 3, ylab = "PKS positive", cex.axis = 1.5, cex.lab = 1.5) 
dev.off()


pdf("correlationBurdenPks.pdf")
plot(ratios,as.numeric(table(newNames)),  ylim = c(0,45000), ylab = "T>N mutations", xlab = "Fraction pks-positive", cex.lab = 1.5, cex.axis = 1.3)
abline(v = 0.1, lty = 2, lwd = 3, col = "red")
dev.off()



##### correlate with SBS88 contribution
library(MutationalPatterns)
mutMatColon = read.table("HMF_mut_mat_colon_PASS.txt") ### matrix with 96-mutation types of colon cancers in HMF dataset

#### make sure matrix is in same order as the pks score
mutMatColon = mutMatColon[,names(score)]

existing = get_known_signatures(muttype = "snv", incl_poss_artifacts = F, genome = "GRCh37")
checkContr = fit_to_signatures(mutMatColon, existing)
relContr88 = checkContr$contribution[54,]/apply(checkContr$contribution,2,sum)

pdf("Correlation_score_SBS88_another.pdf")
plot(scoreColon, relContr88, xlab = "Fraction pks-positive", ylab = "Relative contribution SBS88", lwd = 2,cex.lab = 1.5, cex.axis = 1.5)
abline(v = 0.1, col = "red", lty = 2, lwd = 3)
dev.off()

### check in breast cancer 
mutMatHMF = read.table("/Users/joskeubels/Desktop/HMF_mut_mat.txt")
mutMatBreast = mutMatHMF[,sampleName[which(tumorTypes == "Breast")]]
scoreBreast = score[which(tumorTypes == "Breast")]

checkContrBreast = fit_to_signatures(mutMatBreast, existing)

relContr88 = checkContrBreast$contribution[54,]/apply(checkContrBreast$contribution,2,sum)

pdf("Correlation_score_SBS88_Breast.pdf")
plot(scoreBreast, relContr88, xlab = "Fraction pks-positive", ylab = "Relative contribution SBS88", lwd = 2,cex.lab = 1.5, cex.axis = 1.5)
dev.off()

### read HMF metdata
metaHMF = read.delim("metadata.tsv", header = T, fill = TRUE, sep = "\t")
rownames(metaHMF)= metaHMF$sampleId
metaColon = metaHMF[names(score),]


