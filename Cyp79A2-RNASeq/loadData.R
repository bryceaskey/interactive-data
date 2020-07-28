# Script to load FPKM and gene description data from .csv files, merge them, calculate FPKM mean
# and standard error for each set of samples, and save the prepared data as a .RData file.

# Install and load packages, and set working directory
setwd("C:/Users/Bryce/Documents/interactive-data/Cyp79A2-RNASeq/")

# Load FPKM and gene description data from .csv files, and merge into a single object
fpkm <- read.csv("data/fpkm_multimapping.csv")
colnames(fpkm)[1] <- "locus"
araport11 <- read.csv("data/Araport11.csv", header=FALSE)
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")
rawData <- merge(araport11, fpkm, by="locus", all.y=TRUE)

# Mutate/transmute functions from dplyr don't work with mean/sd functions
# Mean and sd need to be calculated with for loop
WT_FPKM <- vector(mode="numeric", length=nrow(rawData))
WT_SE <- vector(mode="numeric", length=nrow(rawData))
WT_PAOx_FPKM <- vector(mode="numeric", length=nrow(rawData))
WT_PAOx_SE <- vector(mode="numeric", length=nrow(rawData))
Cyp79A2_FPKM <- vector(mode="numeric", length=nrow(rawData))
Cyp79A2_SE <- vector(mode="numeric", length=nrow(rawData))

for(i in 1:nrow(rawData)){
  WT_FPKM[i] <- mean(c(rawData$WT_rep1[i], rawData$WT_rep2[i], rawData$WT_rep3[i]))
  WT_SE[i] <- sd(c(rawData$WT_rep1[i], rawData$WT_rep2[i], rawData$WT_rep3[i]))/sqrt(3)
  
  WT_PAOx_FPKM[i] <- mean(c(rawData$WT.PAOx_rep1[i], rawData$WT.PAOx_rep2[i], rawData$WT.PAOx_rep3[i]))
  WT_PAOx_SE[i] <- sd(c(rawData$WT.PAOx_rep1[i], rawData$WT.PAOx_rep2[i], rawData$WT.PAOx_rep3[i]))/sqrt(3)
  
  Cyp79A2_FPKM[i] <- mean(c(rawData$Cyp79A2_rep1[i], rawData$Cyp79A2_rep2[i], rawData$Cyp79A2_rep3[i]))
  Cyp79A2_SE[i] <- sd(c(rawData$Cyp79A2_rep1[i], rawData$Cyp79A2_rep2[i], rawData$Cyp79A2_rep3[i]))/sqrt(3)
}

set_dec_places <- function(num, dec_places){
  if(is.na(num)){
    return(NA)
  }else{
    return(as.numeric(trimws(format(round(num, dec_places), nsmall=dec_places))))
  }
}
  
allData <- data.frame(
  locus=rawData$locus,
  short_name=rawData$short_name,
  name=rawData$name,
  aliases=rawData$aliases,
  
  WT_FPKM=sapply(WT_FPKM, set_dec_places, dec_places=3),
  WT_SE=sapply(WT_SE, set_dec_places, dec_places=3),
  
  WT_PAOx_FPKM=sapply(WT_PAOx_FPKM, set_dec_places, dec_places=3),
  WT_PAOx_SE=sapply(WT_PAOx_SE, set_dec_places, dec_places=3),
  
  Cyp79A2_FPKM=sapply(Cyp79A2_FPKM, set_dec_places, dec_places=3),
  Cyp79A2_SE=sapply(Cyp79A2_SE, set_dec_places, dec_places=3)
)

save(allData, file="data/allData.rda")