# Script to load rpkm and gene description data from .csv files, merge them, calculate rpkm mean
# and standard error for each set of samples, and save the prepared data as a .RData file.

# Install and load packages, and set working directory
setwd("C:/Users/Bryce/Documents/interactive-data/hy5-RNASeq/")

# Load rpkm and gene description data from .csv files, and merge into a single object
rpkm <- read.csv("data/rpkm.csv")
colnames(rpkm)[1] <- "locus"
araport11 <- read.csv("data/Araport11.csv", header=FALSE)
araport11 <- araport11[, c(1, 3, 13)]
colnames(araport11) <- c("locus", "name", "aliases")
rawData <- merge(araport11, rpkm, by="locus", all=TRUE)

# Mutate/transmute functions from dplyr don't work with mean/sd functions
# Mean and sd need to be calculated with for loop
Col0_d3d_RPKM <- vector(mode="numeric", length=nrow(rawData))
Col0_d3d_SE <- vector(mode="numeric", length=nrow(rawData))
hy5_d3d_RPKM <- vector(mode="numeric", length=nrow(rawData))
hy5_d3d_SE <- vector(mode="numeric", length=nrow(rawData))
Col0_d3d_l1.5h_RPKM <- vector(mode="numeric", length=nrow(rawData))
Col0_d3d_l1.5h_SE <- vector(mode="numeric", length=nrow(rawData))
hy5_d3d_l1.5h_RPKM <- vector(mode="numeric", length=nrow(rawData))
hy5_d3d_l1.5h_SE <- vector(mode="numeric", length=nrow(rawData))

for(i in 1:nrow(rawData)){
  Col0_d3d_RPKM[i] <- mean(c(rawData$SRR9313209[i], rawData$SRR9313210[i], rawData$SRR9313211[i]))
  Col0_d3d_SE[i] <- sd(c(rawData$SRR9313209[i], rawData$SRR9313210[i], rawData$SRR9313211[i]))/sqrt(3)
  
  hy5_d3d_RPKM[i] <- mean(c(rawData$SRR9313212[i], rawData$SRR9313213[i], rawData$SRR9313214[i]))
  hy5_d3d_SE[i] <- sd(c(rawData$SRR9313212[i], rawData$SRR9313213[i], rawData$SRR9313214[i]))/sqrt(3)
  
  Col0_d3d_l1.5h_RPKM[i] <- mean(c(rawData$SRR9313223[i], rawData$SRR9313224[i], rawData$SRR9313225[i]))
  Col0_d3d_l1.5h_SE[i] <- sd(c(rawData$SRR9313223[i], rawData$SRR9313224[i], rawData$SRR9313225[i]))/sqrt(3)
  
  hy5_d3d_l1.5h_RPKM[i] <- mean(c(rawData$SRR9313226[i], rawData$SRR9313227[i], rawData$SRR9313228[i]))
  hy5_d3d_l1.5h_SE[i] <- sd(c(rawData$SRR9313226[i], rawData$SRR9313227[i], rawData$SRR9313228[i]))/sqrt(3)
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
  name=rawData$name,
  aliases=rawData$aliases,
  
  Col0_d3d_RPKM=sapply(Col0_d3d_RPKM, set_dec_places, dec_places=3),
  Col0_d3d_SE=sapply(Col0_d3d_SE, set_dec_places, dec_places=3),
  
  hy5_d3d_RPKM=sapply(hy5_d3d_RPKM, set_dec_places, dec_places=3),
  hy5_d3d_SE=sapply(hy5_d3d_SE, set_dec_places, dec_places=3),
  
  Col0_d3d_l1.5h_RPKM=sapply(Col0_d3d_l1.5h_RPKM, set_dec_places, dec_places=3),
  Col0_d3d_l1.5h_SE=sapply(Col0_d3d_l1.5h_SE, set_dec_places, dec_places=3),
  
  hy5_d3d_l1.5h_RPKM=sapply(hy5_d3d_l1.5h_RPKM, set_dec_places, dec_places=3),
  hy5_d3d_l1.5h_SE=sapply(hy5_d3d_l1.5h_SE, set_dec_places, dec_places=3)
)

save(allData, file="data/allData.rda")