### parse train data 

files_snv_train <- list.files(path="path/to/files",pattern="merged.bed",full.names=T, recursive = T) ### output BED files from 2_parseFeatures.sh

train.data <- list()
val <- c()
samples <- c()

for ( PKS_file in files_snv_train) {
  PKS.data <- read.table(PKS_file)
  sample <- gsub(".+/(.+)_snvs.+merged.bed","\\1",PKS_file)
  feature <- gsub(".+sorted.(.+).merged.bed","\\1",PKS_file)
  if ( feature == "TSB" ) { 
    colnames(PKS.data) <- c("CHROM","START","END","NAME","SCORE","STRAND","GENE",feature)
  } else {
    colnames(PKS.data) <- c("CHROM","START","END","NAME","SCORE","STRAND",feature)  
  }
  PKS.data <- PKS.data[PKS.data$CHROM != 'M',]
  
  if ( ! sample %in% names(train.data) ) {
    train.data[[sample]] <- PKS.data
    train.data[[sample]]$SAMPLE <- sample
    train.data[[sample]]$VAL <- "PKS"
  } else{
    train.data[[sample]] <- merge(train.data[[sample]], PKS.data, by=c("CHROM","START","END","NAME","SCORE","STRAND"))
  }
}

##### combine all samples in train.data objects
train.data <- do.call("rbind",train.data)
train.data <- train.data[!duplicated(train.data[,"START"]),]
train.data$TSB2 <- "Unknown"
train.data[train.data$STRAND == "+" & train.data$TSB == "+",]$TSB2 <- "Untranscribed"
train.data[train.data$STRAND == "+" & train.data$TSB == "-",]$TSB2 <- "Transcribed"
train.data[train.data$STRAND == "-" & train.data$TSB == "+",]$TSB2 <- "Transcribed"
train.data[train.data$STRAND == "-" & train.data$TSB == "-",]$TSB2 <- "Untranscribed"

df.context <- data.frame(do.call('rbind', strsplit(as.character(train.data$NAME),'',fixed=TRUE)))
colnames(df.context) <- c(paste("POSm",c(10:1),sep=""), "POS0",  paste("POSp",c(1:10),sep=""))
train.data <- cbind(train.data, df.context)

train.data.sub <- train.data[,c("GENEBODY","VAL","REPLISEQ","SIMPLEREPEAT","TSB2",paste("POSm",c(10:1),sep=""),"POS0",paste("POSp",c(1:10),sep=""))]

train.data.sub$GENEBODY <- as.numeric(as.vector(train.data.sub$GENEBODY))
train.data.sub$REPLISEQ <- as.numeric(as.vector(train.data.sub$REPLISEQ))
train.data.sub$SIMPLEREPEAT <- as.numeric(as.vector(train.data.sub$SIMPLEREPEAT))

#### prepare as input data for random forest 
train.data.sub$VAL <- as.factor(train.data.sub$VAL)
train.data.sub$TSB2 <- as.factor(train.data.sub$TSB2)
train.data.sub$POSm10 <- as.factor(train.data.sub$POSm10)
train.data.sub$POSm9 <- as.factor(train.data.sub$POSm9)
train.data.sub$POSm8 <- as.factor(train.data.sub$POSm8)
train.data.sub$POSm7 <- as.factor(train.data.sub$POSm7)
train.data.sub$POSm6 <- as.factor(train.data.sub$POSm6)
train.data.sub$POSm5 <- as.factor(train.data.sub$POSm5)
train.data.sub$POSm4 <- as.factor(train.data.sub$POSm4)
train.data.sub$POSm3 <- as.factor(train.data.sub$POSm3)
train.data.sub$POSm2 <- as.factor(train.data.sub$POSm2)
train.data.sub$POSm1 <- as.factor(train.data.sub$POSm1)
train.data.sub$POS0 <- as.factor(train.data.sub$POS0)
train.data.sub$POSp10 <- as.factor(train.data.sub$POSp10)
train.data.sub$POSp9 <- as.factor(train.data.sub$POSp9)
train.data.sub$POSp8 <- as.factor(train.data.sub$POSp8)
train.data.sub$POSp7 <- as.factor(train.data.sub$POSp7)
train.data.sub$POSp6 <- as.factor(train.data.sub$POSp6)
train.data.sub$POSp5 <- as.factor(train.data.sub$POSp5)
train.data.sub$POSp4 <- as.factor(train.data.sub$POSp4)
train.data.sub$POSp3 <- as.factor(train.data.sub$POSp3)
train.data.sub$POSp2 <- as.factor(train.data.sub$POSp2)
train.data.sub$POSp1 <- as.factor(train.data.sub$POSp1)

save(train.data.sub, file = "trainingDataSNV.Rdata")
