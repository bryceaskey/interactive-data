# Helper functions for main app

# Function to create dataframe containing Gene_ID, Gene_name, Gene_family, mean FPKM, and se 
# FKPM values for only selected genotypes from full RNAseq dataframe
subsetGenotype <- function(data, selectedGenotypes){
  subsetData <- data[, 1:3]
  for(genotype in selectedGenotypes){
    meanColName <- paste(genotype, "_FPKM", sep="")
    meanCol <- data[, colnames(data)==meanColName]
    seColName <- paste(genotype, "_FPKM_SE", sep="")
    seCol <- data[, colnames(data)==seColName]
    subsetData[, meanColName] <- meanCol
    subsetData[, seColName] <- seCol
  }
  return(subsetData)
}

# Function to restore genotype data to selectedData (in case a genotype is deselected, then
# reselected)
restoreGenotype <- function(selectedData, allData, selectedGenotypes){
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

# Function to convert data from wide format (with FPKM and SE values spread out across multiple
# columns) to tidy format for plotting with ggplot
convertToTidy <- function(data, selectedGenotypes){
  tidyData <- NULL
  for(genotype in selectedGenotypes){
    genoData <- data[, 1:3]
    meanColName <- paste(genotype, "_FPKM", sep="")
    meanCol <- meanCol <- data[, colnames(data)==meanColName]
    seColName <- paste(genotype, "_FPKM_SE", sep="")
    seCol <- seCol <- data[, colnames(data)==seColName]
    genoData$FPKM <- meanCol
    genoData$SE <- seCol
    genoData$Genotype <- genotype
    tidyData <- rbind(tidyData, genoData)
  }
  return(tidyData)
}

getGenotypeCols <- function(selectedGenotypes){
  genotypeCols <- vector(mode="numeric", length=0)
  if(sum(selectedGenotypes=="WT")==1){genotypeCols <- append(genotypeCols, 4)}
  if(sum(selectedGenotypes=="Ref5")==1){genotypeCols <- append(genotypeCols, 6)}
  if(sum(selectedGenotypes=="Ref2")==1){genotypeCols <- append(genotypeCols, 8)}
  if(sum(selectedGenotypes=="Med5")==1){genotypeCols <- append(genotypeCols, 10)}
  if(sum(selectedGenotypes=="Ref5Med5")==1){genotypeCols <- append(genotypeCols, 12)}
  if(sum(selectedGenotypes=="Ref2Med5")==1){genotypeCols <- append(genotypeCols, 14)}
  return(genotypeCols)
}

getGenotypeCols2 <- function(selectedGenotypes){
  genotypeCols <- vector(mode="numeric", length=0)
  if(sum(selectedGenotypes=="WT")==1){genotypeCols <- append(genotypeCols, c(4, 5))}
  if(sum(selectedGenotypes=="Ref5")==1){genotypeCols <- append(genotypeCols, c(6, 7))}
  if(sum(selectedGenotypes=="Ref2")==1){genotypeCols <- append(genotypeCols, c(8, 9))}
  if(sum(selectedGenotypes=="Med5")==1){genotypeCols <- append(genotypeCols, c(10, 11))}
  if(sum(selectedGenotypes=="Ref5Med5")==1){genotypeCols <- append(genotypeCols, c(12, 13))}
  if(sum(selectedGenotypes=="Ref2Med5")==1){genotypeCols <- append(genotypeCols, c(14, 15))}
  return(genotypeCols)
}