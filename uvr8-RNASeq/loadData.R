# Script to load rpkm and gene description data from .csv files, merge them, calculate rpkm mean
# and standard error for each set of samples, and save the prepared data as a .RData file.

# Install and load packages, and set working directory
setwd("C:/Users/Bryce/Documents/interactive-data/uvr8-RNASeq/")

# Load rpkm and gene description data from .csv files, and merge into a single object
rpkm <- read.csv("data/rpkm.csv")
colnames(rpkm)[1] <- "locus"
araport11 <- read.csv("data/Araport11.csv", header=FALSE)
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")
rawData <- merge(araport11, rpkm, by="locus", all.y=TRUE)

# Mutate/transmute functions from dplyr don't work with mean/sd functions
# Mean and sd need to be calculated with for loop
Col0_WL_RPKM <- vector(mode="numeric", length=nrow(rawData))
Col0_WL_SE <- vector(mode="numeric", length=nrow(rawData))
Col0_WL_UVB6h_RPKM <- vector(mode="numeric", length=nrow(rawData))
Col0_WL_UVB6h_SE <- vector(mode="numeric", length=nrow(rawData))
uvr8_WL_RPKM <- vector(mode="numeric", length=nrow(rawData))
uvr8_WL_SE <- vector(mode="numeric", length=nrow(rawData))
uvr8_WL_UVB6h_RPKM <- vector(mode="numeric", length=nrow(rawData))
uvr8_WL_UVB6h_SE <- vector(mode="numeric", length=nrow(rawData))

for(i in 1:nrow(rawData)){
  Col0_WL_RPKM[i] <- mean(c(rawData$SRR9200658[i], rawData$SRR9200661[i], rawData$SRR9200664[i]))
  Col0_WL_SE[i] <- sd(c(rawData$SRR9200658[i], rawData$SRR9200661[i], rawData$SRR9200664[i]))/sqrt(3)
  
  Col0_WL_UVB6h_RPKM[i] <- mean(c(rawData$SRR9200667[i], rawData$SRR9200670[i], rawData$SRR9200674[i]))
  Col0_WL_UVB6h_SE[i] <- sd(c(rawData$SRR9200667[i], rawData$SRR9200670[i], rawData$SRR9200674[i]))/sqrt(3)
  
  uvr8_WL_RPKM[i] <- mean(c(rawData$SRR9200677[i], rawData$SRR9200680[i], rawData$SRR9200683[i]))
  uvr8_WL_SE[i] <- sd(c(rawData$SRR9200677[i], rawData$SRR9200680[i], rawData$SRR9200683[i]))/sqrt(3)
  
  uvr8_WL_UVB6h_RPKM[i] <- mean(c(rawData$SRR9200686[i], rawData$SRR9200690[i], rawData$SRR9200693[i]))
  uvr8_WL_UVB6h_SE[i] <- sd(c(rawData$SRR9200686[i], rawData$SRR9200690[i], rawData$SRR9200693[i]))/sqrt(3)
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
  
  Col0_WL_RPKM=sapply(Col0_WL_RPKM, set_dec_places, dec_places=3),
  Col0_WL_SE=sapply(Col0_WL_SE, set_dec_places, dec_places=3),
  
  Col0_WL_UVB6h_RPKM=sapply(Col0_WL_UVB6h_RPKM, set_dec_places, dec_places=3),
  Col0_WL_UVB6h_SE=sapply(Col0_WL_UVB6h_SE, set_dec_places, dec_places=3),
  
  uvr8_WL_RPKM=sapply(uvr8_WL_RPKM, set_dec_places, dec_places=3),
  uvr8_WL_SE=sapply(uvr8_WL_SE, set_dec_places, dec_places=3),
  
  uvr8_WL_UVB6h_RPKM=sapply(uvr8_WL_UVB6h_RPKM, set_dec_places, dec_places=3),
  uvr8_WL_UVB6h_SE=sapply(uvr8_WL_UVB6h_SE, set_dec_places, dec_places=3)
)

save(allData, file="data/allData.rda")