# Helper functions for main app

# Function to create dataframe containing Gene_ID, Gene_name, Gene_family, mean FPKM, and se 
# FKPM values for only selected genotypes from full RNAseq dataframe
subsetGenotype <- function(data, selectedSamples, PRJNA){
  if(PRJNA==PRJNA549285){
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
}

# Function to restore genotype data to selectedData (in case a genotype is deselected, then
# reselected)
restoreGenotype <- function(selectedData, allData, selectedSamples, PRJNA){
  if(PRJNA==PRJNA549285){
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
}

# Function to convert data from wide format (with RPKM and SE values spread out across multiple
# columns) to tidy format for plotting with ggplot
convertToTidy <- function(data, selectedSamples, PRJNA, unit=NULL){
  tidyData <- NULL
  if(PRJNA=="PRJNA549285" | PRJNA=="PRJNA546251"){
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
  }else if(PRJNA=="PRJNA230565" | PRJNA=="PRJNA388948" | PRJNA=="Cyp79A2"){
    for(sample in selectedSamples){
      sampleData <- data[, 1:4]
      meanColName <- paste(unit, "_", sample, sep="")
      meanCol <- data[, colnames(data)==meanColName]
      seColName <- paste(unit, "SE_", sample, sep="")
      seCol <- data[, colnames(data)==seColName]
      pColName <- paste("pValue_", sample, sep="")
      pCol <- data[, colnames(data)==pColName]
      if(unit=="tmm" & sample != "control" & sample != "WT"){
        signifCol <- NULL
        for(pVal in pCol){
          if(pVal < 0.001){
            signifCol <- c(signifCol, "***")
          }else if(pVal < 0.01){
            signifCol <- c(signifCol, "**")
          }else if(pVal < 0.05){
            signifCol <- c(signifCol, "*")
          }else{
            signifCol <- c(signifCol, NA)
          }
        }
      }else{
        signifCol <- rep(NA, length(meanCol))
        
      }
      sampleData$sample <- sample
      sampleData$exp <- meanCol
      sampleData$maxExp <- max(sampleData$exp)
      sampleData$se <- seCol
      sampleData$signif <- signifCol
      tidyData <- rbind(tidyData, sampleData)        
    }
    return(tidyData)
  }
}

# Without standard error - for display in app
getGenotypeCols <- function(selectedSamples, PRJNA, unit=NULL){
  if(PRJNA=="PRJNA549285"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(sum(selectedSamples=="Col0_d3d")==1){genotypeCols <- append(genotypeCols, 5)}
    if(sum(selectedSamples=="hy5_d3d")==1){genotypeCols <- append(genotypeCols, 7)}
    if(sum(selectedSamples=="Col0_d3d_l1.5h")==1){genotypeCols <- append(genotypeCols, 9)}
    if(sum(selectedSamples=="hy5_d3d_l1.5h")==1){genotypeCols <- append(genotypeCols, 11)}
    return(genotypeCols)
  }else if(PRJNA=="PRJNA546251"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(sum(selectedSamples=="Col0_WL")==1){genotypeCols <- append(genotypeCols, 5)}
    if(sum(selectedSamples=="uvr8_WL")==1){genotypeCols <- append(genotypeCols, 9)}
    if(sum(selectedSamples=="Col0_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, 7)}
    if(sum(selectedSamples=="uvr8_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, 11)}
    return(genotypeCols)
  }else if(PRJNA=="PRJNA230565"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(unit=="fpkm"){
      if(sum(selectedSamples=="control")==1){genotypeCols <- append(genotypeCols, 6)}
      if(sum(selectedSamples=="NAA")==1){genotypeCols <- append(genotypeCols, 10)}
    }else if(unit=="tmm"){
      if(sum(selectedSamples=="control")==1){genotypeCols <- append(genotypeCols, 8)}
      if(sum(selectedSamples=="NAA")==1){genotypeCols <- append(genotypeCols, c(12,14))}
    }
    return(genotypeCols)
  }else if(PRJNA=="PRJNA388948"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(unit=="fpkm"){
      if(sum(selectedSamples=="WT")==1){genotypeCols <- append(genotypeCols, 6)}
      if(sum(selectedSamples=="ref5")==1){genotypeCols <- append(genotypeCols, 10)}
      if(sum(selectedSamples=="ref2")==1){genotypeCols <- append(genotypeCols, 16)}
      if(sum(selectedSamples=="med5")==1){genotypeCols <- append(genotypeCols, 22)}
      if(sum(selectedSamples=="ref5med5")==1){genotypeCols <- append(genotypeCols, 28)}
      if(sum(selectedSamples=="ref2med5")==1){genotypeCols <- append(genotypeCols, 34)}
    }else if(unit=="tmm"){
      if(sum(selectedSamples=="WT")==1){genotypeCols <- append(genotypeCols, 8)}
      if(sum(selectedSamples=="ref5")==1){genotypeCols <- append(genotypeCols, 12)}
      if(sum(selectedSamples=="ref5")==1){genotypeCols <- append(genotypeCols, 18)}
      if(sum(selectedSamples=="med5")==1){genotypeCols <- append(genotypeCols, 24)}
      if(sum(selectedSamples=="ref5med5")==1){genotypeCols <- append(genotypeCols, 30)}
      if(sum(selectedSamples=="ref2med5")==1){genotypeCols <- append(genotypeCols, 36)}
    }
    return(genotypeCols)
  }else if(PRJNA=="Cyp79A2"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(unit=="fpkm"){
      if(sum(selectedSamples=="WT")==1){genotypeCols <- append(genotypeCols, 6)}
      if(sum(selectedSamples=="PAOx")==1){genotypeCols <- append(genotypeCols, 10)}
      if(sum(selectedSamples=="Cyp79A2")==1){genotypeCols <- append(genotypeCols, 16)}
    }else if(unit=="tmm"){
      if(sum(selectedSamples=="WT")==1){genotypeCols <- append(genotypeCols, 8)}
      if(sum(selectedSamples=="PAOx")==1){genotypeCols <- append(genotypeCols, c(12,14))}
      if(sum(selectedSamples=="Cyp79A2")==1){genotypeCols <- append(genotypeCols, c(18,20))}
    }
    return(genotypeCols)
  }
}

# With standard error - for download
getGenotypeCols2 <- function(selectedSamples, PRJNA, unit=NULL){
  if(PRJNA=="PRJNA549285"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(sum(selectedSamples=="Col0_d3d")==1){genotypeCols <- append(genotypeCols, c(5, 6))}
    if(sum(selectedSamples=="hy5_d3d")==1){genotypeCols <- append(genotypeCols, c(7, 8))}
    if(sum(selectedSamples=="Col0_d3d_l1.5h")==1){genotypeCols <- append(genotypeCols, c(9, 10))}
    if(sum(selectedSamples=="hy5_d3d_l1.5h")==1){genotypeCols <- append(genotypeCols, c(11, 12))}
    return(genotypeCols)
  }else if(PRJNA=="PRJNA546251"){
    genotypeCols <- vector(mode="numeric", length=0)
    if(sum(selectedSamples=="Col0_WL")==1){genotypeCols <- append(genotypeCols, c(5, 6))}
    if(sum(selectedSamples=="uvr8_WL")==1){genotypeCols <- append(genotypeCols, c(9, 10))}
    if(sum(selectedSamples=="Col0_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, c(7, 8))}
    if(sum(selectedSamples=="uvr8_WL_UVB6h")==1){genotypeCols <- append(genotypeCols, c(11, 12))}
    return(genotypeCols)
  }else if(PRJNA=="PRJNA230565"){
    genotypeCols <- c(6:15)
    return(genotypeCols)
  }else if(PRJNA=="PRJNA388948"){
    genotypeCols <- c(6:39)
    return(genotypeCols)
  }else if(PRJNA=="Cyp79A2"){
    genotypeCols <- c(6:21)
    return(genotypeCols)
  }
}