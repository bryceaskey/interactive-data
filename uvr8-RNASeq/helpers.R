# Helper functions for main app

# Function to create dataframe containing Gene_ID, Gene_name, Gene_family, mean FPKM, and se 
# FKPM values for only selected genotypes from full RNAseq dataframe
subsetGenotype <- function(data, selectedSamples){
  subsetData <- data[, 1:4]
  for(genotype in selectedSamples){
    meanColName <- paste(genotype, "_RPKM", sep="")
    meanCol <- data[, colnames(data)==meanColName]
    seColName <- paste(genotype, "_SE", sep="")
    seCol <- data[, colnames(data)==seColName]
    subsetData[, meanColName] <- meanCol
    subsetData[, seColName] <- seCol
  }
  return(subsetData)
}

# Function to restore genotype data to selectedData (in case a genotype is deselected, then
# reselected)
restoreGenotype <- function(selectedData, allData, selectedSamples){
  selectedRows <- rownames(selectedData)
  selectedData <- NULL
  rowCount <- 0
  for(rowName in selectedRows){
    rowCount <- rowCount + 1
    selectedData <- rbind(selectedData, allData[rownames(allData)==rowName, ])
    rownames(selectedData)[rowCount] <- rowName
  }
  return(selectedData)
}

# Function to convert data from wide format (with RPKM and SE values spread out across multiple
# columns) to tidy format for plotting with ggplot
convertToTidy <- function(data, selectedSamples){
  tidyData <- NULL
  for(sample in selectedSamples){
    sampleData <- data[, 1:4]
    meanColName <- paste(sample, "_RPKM", sep="")
    meanCol <- data[, colnames(data)==meanColName]
    seColName <- paste(sample, "_SE", sep="")
    seCol <- data[, colnames(data)==seColName]
    sampleData$RPKM <- meanCol
    sampleData$SE <- seCol
    sampleData$sample <- sample
    tidyData <- rbind(tidyData, sampleData)
  }
  return(tidyData)
}

# Without standard error - for display in app
getGenotypeCols <- function(selectedSamples){
  genotypeCols <- vector(mode="numeric", length=0)
  if(sum(selectedSamples=="Col0_WL")==1){genotypeCols <- append(genotypeCols, 5)}
  if(sum(selectedSamples=="uvr8_WL")==1){genotypeCols <- append(genotypeCols, 9)}
  if(sum(selectedSamples=="Col0_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, 7)}
  if(sum(selectedSamples=="uvr8_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, 11)}
  return(genotypeCols)
}

# With standard error - for download
getGenotypeCols2 <- function(selectedSamples){
  genotypeCols <- vector(mode="numeric", length=0)
  if(sum(selectedSamples=="Col0_WL")==1){genotypeCols <- append(genotypeCols, c(5, 6))}
  if(sum(selectedSamples=="uvr8_WL")==1){genotypeCols <- append(genotypeCols, c(9, 10))}
  if(sum(selectedSamples=="Col0_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, c(7, 8))}
  if(sum(selectedSamples=="uvr8_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, c(11, 12))}
  return(genotypeCols)
}