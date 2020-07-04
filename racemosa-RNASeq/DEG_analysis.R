#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(tibble)

setwd("C:/Users/Bryce/Documents/interactive-data/racemosa-RNASeq/")

# Read data from .csv file
rawData <- read.csv("read_counts_annotated.csv", row.names=1)
colnames(rawData)[1:6] <- c("leaf1", "leaf2", "leaf3", "root1", "root2", "root3")

readData <- rawData[, 1:6]

# Create metadata for samples
metadata <- data.frame(
  row.names=colnames(readData),
  condition=factor(c(rep("leaf", 3), rep("root", 3))),
  type=factor(c(rep("paired-end", 6)))
)

# Prefilter data by dropping contigs having at least 1 sample with read count = 0
readDataFiltered <- readData %>%
  rownames_to_column("contig") %>%
  filter(leaf1!=0 & leaf2!=0 & leaf3!=0 & root1!=0 & root2!=0 &root3!=0) %>%
  column_to_rownames("contig")

# Construct DESeqDataSet object from filtered read count data and metadata
DESeqData <- DESeqDataSetFromMatrix(
  countData=readDataFiltered,
  colData=metadata,
  design=~condition
)

DESeqData <- DESeq(DESeqData)
DESeqResults <- results(DESeqData, alpha=0.05)

# Construct final dataframe of outputs with contig names, BLAST results, FPKM values, 
# log2 fold change, log2 standard error, and adjusted p-values