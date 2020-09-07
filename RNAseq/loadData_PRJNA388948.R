# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data

# Install and load packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

library(edgeR)
library(tibble)

# Specify paths to data files
justCounts <- "C:/Users/Bryce/Documents/Cyp79A2-RNAseq/data/PRJNA388948/expression_counts/"
geneLengths <- "C:/Users/Bryce/Documents/Cyp79A2-RNAseq/data/PRJNA388948/gene_lengths.txt"
araport11 <- read.csv("C:/Users/Bryce/Documents/Cyp79A2-RNAseq/data/Araport11.csv", header=FALSE)

# Load gene expression and gene length data
countFiles <- paste(justCounts, dir(justCounts), sep="")
countLabels <- c("WT_rep1", "WT_rep2", "WT_rep3",
                 "ref5_rep1", "ref5_rep2", "ref5_rep3",
                 "ref2_rep1", "ref2_rep2", "ref2_rep3",
                 "med5_rep1", "med5_rep2", "med5_rep3",
                 "ref5med5_rep1", "ref5med5_rep2", "ref5med5_rep3",
                 "ref2med5_rep1", "ref2med5_rep2", "ref2med5_rep3")
expressionData <- readDGE(countFiles,
                          labels=countLabels,
                          group=c("A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E", "G", "G", "G"))
expressionData$genes <- read.delim(geneLengths, row.names=1)

# Only keep genes with significant level of expression
keep <- filterByExpr(expressionData)
expressionData <- expressionData[keep, , keep.lib.sizes=FALSE]

# Normalize for library size
expressionData <- calcNormFactors(expressionData, method="TMM")

# Create MDS plot to verify integrity of RNAseq data
plotMDS(expressionData, method="bcv")

# Perform DEG analysis to compare gene expression of ref5 and ref2 to WT
expDesign <- model.matrix(~0+group, data=expressionData$samples)
expressionData <- estimateDisp(expressionData, design=expDesign)
colnames(expDesign) <- levels(expressionData$samples$group)
modelFit <- glmQLFit(expressionData, expDesign)

contrast_ref5_vs_WT <- makeContrasts(B-A, levels=expDesign)
ref5_vs_WT <- glmQLFTest(modelFit, contrast=contrast_ref5_vs_WT)
ref5_DEGs <- topTags(ref5_vs_WT, n=nrow(expressionData), p.value=1)$table
ref5_DEGs <- rownames_to_column(ref5_DEGs, var="locus")

contrast_ref2_vs_WT <- makeContrasts(C-A, levels=expDesign)
ref2_vs_WT <- glmQLFTest(modelFit, contrast=contrast_ref2_vs_WT)
ref2_DEGs <- topTags(ref2_vs_WT, n=nrow(expressionData), p.value=1)$table
ref2_DEGs <- rownames_to_column(ref2_DEGs, var="locus")

contrast_med5_vs_WT <- makeContrasts(D-A, levels=expDesign)
med5_vs_WT <- glmQLFTest(modelFit, contrast=contrast_med5_vs_WT)
med5_DEGs <- topTags(med5_vs_WT, n=nrow(expressionData), p.value=1)$table
med5_DEGs <- rownames_to_column(med5_DEGs, var="locus")

contrast_ref5med5_vs_WT <- makeContrasts(E-A, levels=expDesign)
ref5med5_vs_WT <- glmQLFTest(modelFit, contrast=contrast_ref5med5_vs_WT)
ref5med5_DEGs <- topTags(ref5med5_vs_WT, n=nrow(expressionData), p.value=1)$table
ref5med5_DEGs <- rownames_to_column(ref5med5_DEGs, var="locus")

contrast_ref2med5_vs_WT <- makeContrasts(G-A, levels=expDesign)
ref2med5_vs_WT <- glmQLFTest(modelFit, contrast=contrast_ref2med5_vs_WT)
ref2med5_DEGs <- topTags(ref2med5_vs_WT, n=nrow(expressionData), p.value=1)$table
ref2med5_DEGs <- rownames_to_column(ref2med5_DEGs, var="locus")

# Load gene annotations
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")

# Add gene annotations to DEG dataframes
ref5_DEGs <- merge(araport11, ref5_DEGs, by="locus", all.y=TRUE)
ref2_DEGs <- merge(araport11, ref2_DEGs, by="locus", all.y=TRUE)
med5_DEGs <- merge(araport11, med5_DEGs, by="locus", all.y=TRUE)
ref5med5_DEGs <- merge(araport11, ref5med5_DEGs, by="locus", all.y=TRUE)
ref2med5_DEGs <- merge(araport11, ref2med5_DEGs, by="locus", all.y=TRUE)

# Calculate fpkm for all genes
fpkmAllGenes <- as.data.frame(rpkm(expressionData$counts, gene.length=expressionData$genes$Length, normalized.lib.sizes=TRUE, log=FALSE))
fpkmAllGenes <- rownames_to_column(fpkmAllGenes, var="locus")

# Calculate tmm for all genes
tmmAllGenes <- as.data.frame(cpm(expressionData, log=FALSE))
tmmAllGenes <- rownames_to_column(tmmAllGenes, var="locus")

# Add tmm and fpkm to ref5 DEG dataframe ----
fpkm_WT <- vector(mode="numeric", length=nrow(ref5_DEGs))
fpkmSE_WT <- vector(mode="numeric", length=nrow(ref5_DEGs))
tmm_WT <- vector(mode="numeric", length=nrow(ref5_DEGs))
tmmSE_WT <- vector(mode="numeric", length=nrow(ref5_DEGs))
fpkm_ref5 <- vector(mode="numeric", length=nrow(ref5_DEGs))
fpkmSE_ref5 <- vector(mode="numeric", length=nrow(ref5_DEGs))
tmm_ref5 <- vector(mode="numeric", length=nrow(ref5_DEGs))
tmmSE_ref5 <- vector(mode="numeric", length=nrow(ref5_DEGs))

for(i in 1:nrow(ref5_DEGs)){
  geneLocus <- ref5_DEGs$locus[i]
  
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
  fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  fpkm_ref5[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 5:7]))
  fpkmSE_ref5[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 5:7])))/sqrt(3)
  
  tmm_ref5[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 5:7]))
  tmmSE_ref5[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 5:7])))/sqrt(3)
}

ref5_DEGs$fpkm_WT <- fpkm_WT
ref5_DEGs$fpkmSE_WT <- fpkmSE_WT

ref5_DEGs$tmm_WT <- tmm_WT
ref5_DEGs$tmmSE_WT <- tmmSE_WT

ref5_DEGs$fpkm_ref5 <- fpkm_ref5
ref5_DEGs$fpkmSE_ref5 <- fpkmSE_ref5

ref5_DEGs$tmm_ref5 <- tmm_ref5
ref5_DEGs$tmmSE_ref5 <- tmmSE_ref5

ref5_DEGs <- ref5_DEGs[c(1:5,11:18,6,9)]
colnames(ref5_DEGs) <- c("locus", "short_name", "name", "aliases", "length",
                        "fpkm_WT", "fpkmSE_WT", "tmm_WT", "tmmSE_WT",
                        "fpkm_ref5", "fpkmSE_ref5", "tmm_ref5", "tmmSE_ref5", "log2FC_ref5", "pValue_ref5")

# Add tmm and fpkm to ref2 DEG dataframe ----
fpkm_WT <- vector(mode="numeric", length=nrow(ref2_DEGs))
fpkmSE_WT <- vector(mode="numeric", length=nrow(ref2_DEGs))
tmm_WT <- vector(mode="numeric", length=nrow(ref2_DEGs))
tmmSE_WT <- vector(mode="numeric", length=nrow(ref2_DEGs))
fpkm_ref2 <- vector(mode="numeric", length=nrow(ref2_DEGs))
fpkmSE_ref2 <- vector(mode="numeric", length=nrow(ref2_DEGs))
tmm_ref2 <- vector(mode="numeric", length=nrow(ref2_DEGs))
tmmSE_ref2 <- vector(mode="numeric", length=nrow(ref2_DEGs))

for(i in 1:nrow(ref2_DEGs)){
  geneLocus <- ref2_DEGs$locus[i]
  
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
  fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  fpkm_ref2[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10]))
  fpkmSE_ref2[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
  
  tmm_ref2[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10]))
  tmmSE_ref2[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10])))/sqrt(3)
}

ref2_DEGs$fpkm_WT <- fpkm_WT
ref2_DEGs$fpkmSE_WT <- fpkmSE_WT

ref2_DEGs$tmm_WT <- tmm_WT
ref2_DEGs$tmmSE_WT <- tmmSE_WT

ref2_DEGs$fpkm_ref2 <- fpkm_ref2
ref2_DEGs$fpkmSE_ref2 <- fpkmSE_ref2

ref2_DEGs$tmm_ref2 <- tmm_ref2
ref2_DEGs$tmmSE_ref2 <- tmmSE_ref2

ref2_DEGs <- ref2_DEGs[c(1,15:18,6,9)]
colnames(ref2_DEGs) <- c("locus",
                         "fpkm_ref2", "fpkmSE_ref2", "tmm_ref2", "tmmSE_ref2", "log2FC_ref2", "pValue_ref2")


# Add tmm and fpkm to med5 DEG dataframe ----
fpkm_WT <- vector(mode="numeric", length=nrow(med5_DEGs))
fpkmSE_WT <- vector(mode="numeric", length=nrow(med5_DEGs))
tmm_WT <- vector(mode="numeric", length=nrow(med5_DEGs))
tmmSE_WT <- vector(mode="numeric", length=nrow(med5_DEGs))
fpkm_med5 <- vector(mode="numeric", length=nrow(med5_DEGs))
fpkmSE_med5 <- vector(mode="numeric", length=nrow(med5_DEGs))
tmm_med5 <- vector(mode="numeric", length=nrow(med5_DEGs))
tmmSE_med5 <- vector(mode="numeric", length=nrow(med5_DEGs))

for(i in 1:nrow(med5_DEGs)){
  geneLocus <- med5_DEGs$locus[i]
  
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
  fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  fpkm_med5[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 11:13]))
  fpkmSE_med5[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 11:13])))/sqrt(3)
  
  tmm_med5[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 11:13]))
  tmmSE_med5[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 11:13])))/sqrt(3)
}

med5_DEGs$fpkm_WT <- fpkm_WT
med5_DEGs$fpkmSE_WT <- fpkmSE_WT

med5_DEGs$tmm_WT <- tmm_WT
med5_DEGs$tmmSE_WT <- tmmSE_WT

med5_DEGs$fpkm_med5 <- fpkm_med5
med5_DEGs$fpkmSE_med5 <- fpkmSE_med5

med5_DEGs$tmm_med5 <- tmm_med5
med5_DEGs$tmmSE_med5 <- tmmSE_med5

med5_DEGs <- med5_DEGs[c(1,15:18,6,9)]
colnames(med5_DEGs) <- c("locus",
                         "fpkm_med5", "fpkmSE_med5", "tmm_med5", "tmmSE_med5", "log2FC_med5", "pValue_med5")

# Add tmm and fpkm to ref5med5 DEG dataframe ----
fpkm_WT <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
fpkmSE_WT <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
tmm_WT <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
tmmSE_WT <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
fpkm_ref5med5 <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
fpkmSE_ref5med5 <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
tmm_ref5med5 <- vector(mode="numeric", length=nrow(ref5med5_DEGs))
tmmSE_ref5med5 <- vector(mode="numeric", length=nrow(ref5med5_DEGs))

for(i in 1:nrow(ref5med5_DEGs)){
  geneLocus <- ref5med5_DEGs$locus[i]
  
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
  fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  fpkm_ref5med5[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 14:16]))
  fpkmSE_ref5med5[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 14:16])))/sqrt(3)
  
  tmm_ref5med5[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 14:16]))
  tmmSE_ref5med5[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 14:16])))/sqrt(3)
}

ref5med5_DEGs$fpkm_WT <- fpkm_WT
ref5med5_DEGs$fpkmSE_WT <- fpkmSE_WT

ref5med5_DEGs$tmm_WT <- tmm_WT
ref5med5_DEGs$tmmSE_WT <- tmmSE_WT

ref5med5_DEGs$fpkm_ref5med5 <- fpkm_ref5med5
ref5med5_DEGs$fpkmSE_ref5med5 <- fpkmSE_ref5med5

ref5med5_DEGs$tmm_ref5med5 <- tmm_ref5med5
ref5med5_DEGs$tmmSE_ref5med5 <- tmmSE_ref5med5

ref5med5_DEGs <- ref5med5_DEGs[c(1,15:18,6,9)]
colnames(ref5med5_DEGs) <- c("locus",
                             "fpkm_ref5med5", "fpkmSE_ref5med5", "tmm_ref5med5", "tmmSE_ref5med5", "log2FC_ref5med5", "pValue_ref5med5")

# Add tmm and fpkm to ref2med5 DEG dataframe ----
fpkm_WT <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
fpkmSE_WT <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
tmm_WT <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
tmmSE_WT <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
fpkm_ref2med5 <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
fpkmSE_ref2med5 <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
tmm_ref2med5 <- vector(mode="numeric", length=nrow(ref2med5_DEGs))
tmmSE_ref2med5 <- vector(mode="numeric", length=nrow(ref2med5_DEGs))

for(i in 1:nrow(ref2med5_DEGs)){
  geneLocus <- ref2med5_DEGs$locus[i]
  
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
  fpkmSE_WT[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  tmmSE_WT[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4])))/sqrt(3)
  
  fpkm_ref2med5[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 17:19]))
  fpkmSE_ref2med5[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 17:19])))/sqrt(3)
  
  tmm_ref2med5[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 17:19]))
  tmmSE_ref2med5[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 17:19])))/sqrt(3)
}

ref2med5_DEGs$fpkm_WT <- fpkm_WT
ref2med5_DEGs$fpkmSE_WT <- fpkmSE_WT

ref2med5_DEGs$tmm_WT <- tmm_WT
ref2med5_DEGs$tmmSE_WT <- tmmSE_WT

ref2med5_DEGs$fpkm_ref2med5 <- fpkm_ref2med5
ref2med5_DEGs$fpkmSE_ref2med5 <- fpkmSE_ref2med5

ref2med5_DEGs$tmm_ref2med5 <- tmm_ref2med5
ref2med5_DEGs$tmmSE_ref2med5 <- tmmSE_ref2med5

ref2med5_DEGs <- ref2med5_DEGs[c(1,15:18,6,9)]
colnames(ref2med5_DEGs) <- c("locus",
                             "fpkm_ref2med5", "fpkmSE_ref2med5", "tmm_ref2med5", "tmmSE_ref2med5", "log2FC_ref2med5", "pValue_ref2med5")

# Merge all DEG dataframes, and limit decimal places displayed ----

#write.csv(ref5_DEGs, file="C:/Users/Bryce/Documents/Cyp79A2-RNAseq/data/ref5_DEGs.csv", row.names=FALSE)
#write.csv(ref2_DEGs, file="C:/Users/Bryce/Documents/Cyp79A2-RNAseq/data/ref2_DEGs.csv", row.names=FALSE)

DEGs_PRJNA388948_med5 <- merge(merge(merge(merge(ref5_DEGs, ref2_DEGs, by="locus"), med5_DEGs, by="locus"), ref5med5_DEGs, by="locus"), ref2med5_DEGs, by="locus")
for(col in 6:ncol(DEGs_PRJNA388948_med5)){
  DEGs_PRJNA388948_med5[,col] <- as.numeric(formatC(DEGs_PRJNA388948_med5[,col], width=5, format="G"))
}
saveRDS(DEGs_PRJNA388948_med5, file="C:/Users/Bryce/Documents/interactive-data/RNAseq/data/PRJNA388948_med5.rds")