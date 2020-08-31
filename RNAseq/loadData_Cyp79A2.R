# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data

# Install and load packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

library(edgeR)
library(tibble)

# Specify paths to data files
justCounts <- "C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/Cyp79A2/expression_counts/"
geneLengths <- "C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/Cyp79A2/gene_lengths.txt"
araport11 <- read.csv("C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/Araport11.csv")

# Load gene expression and gene length data
countFiles <- paste(justCounts, dir(justCounts), sep="")
countLabels <- c("Cyp79A2_rep1", "Cyp79A2_rep2", "Cyp79A2_rep3", "WT-PAOx_rep1", "WT-PAOx_rep2", "WT-PAOx_rep3", "WT_rep1", "WT_rep2", "WT_rep3")
expressionData <- readDGE(countFiles,
                          labels=countLabels,
                          group=c("A", "A", "A", "B", "B", "B", "C", "C", "C"))
expressionData$genes <- read.delim(geneLengths, row.names=1)

# Only keep genes with significant level of expression
keep <- filterByExpr(expressionData)
expressionData <- expressionData[keep, , keep.lib.sizes=FALSE]

# Normalize for library size
expressionData <- calcNormFactors(expressionData, method="TMM")

# Create MDS plot to verify integrity of RNAseq data
plotMDS(expressionData, method="bcv")

# Perform DEG analysis to compare gene expression of NAA treated plants to control
expDesign <- model.matrix(~0+group, data=expressionData$samples)
expressionData <- estimateDisp(expressionData, design=expDesign)
colnames(expDesign) <- levels(expressionData$samples$group)
modelFit <- glmQLFit(expressionData, expDesign)

contrast_Cyp79A2_vs_WT <- makeContrasts(A-C, levels=expDesign)
Cyp79A2_vs_WT <- glmQLFTest(modelFit, contrast=contrast_Cyp79A2_vs_WT)
Cyp79A2_DEGs <- topTags(Cyp79A2_vs_WT, n=nrow(expressionData), p.value=1)$table
Cyp79A2_DEGs <- rownames_to_column(Cyp79A2_DEGs, var="locus")

contrast_PAOx_vs_WT <- makeContrasts(B-C, levels=expDesign)
PAOx_vs_WT <- glmQLFTest(modelFit, contrast=contrast_PAOx_vs_WT)
PAOx_DEGs <- topTags(PAOx_vs_WT, n=nrow(expressionData), p.value=1)$table
PAOx_DEGs <- rownames_to_column(PAOx_DEGs, var="locus")

# Load gene annotations
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")

# Add gene annotations to DEG dataframes
Cyp79A2_DEGs <- merge(araport11, Cyp79A2_DEGs, by="locus", all.y=TRUE)
PAOx_DEGs <- merge(araport11, PAOx_DEGs, by="locus", all.y=TRUE)

# Calculate fpkm for all genes
fpkmAllGenes <- as.data.frame(rpkm(expressionData$counts, gene.length=expressionData$genes$Length, normalized.lib.sizes=TRUE, log=FALSE))
fpkmAllGenes <- rownames_to_column(fpkmAllGenes, var="locus")

# Calculate tmm for all genes
tmmAllGenes <- as.data.frame(cpm(expressionData, log=FALSE))
tmmAllGenes <- rownames_to_column(tmmAllGenes, var="locus")

# Add tmm and fpkm to Cyp79A2 DEG dataframe ----
tmm_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
tmm_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
fpkm_WT <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
fpkm_Cyp79A2 <- vector(mode="numeric", length=nrow(Cyp79A2_DEGs))
for(i in 1:nrow(Cyp79A2_DEGs)){
  geneLocus <- Cyp79A2_DEGs$locus[i]
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10]))
  tmm_Cyp79A2[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:4]))
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10]))
  fpkm_Cyp79A2[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:4]))
}
Cyp79A2_DEGs$tmm_WT <- tmm_WT
Cyp79A2_DEGs$tmm_Cyp79A2 <- tmm_Cyp79A2
Cyp79A2_DEGs$fpkm_WT <- fpkm_WT
Cyp79A2_DEGs$fpkm_Cyp79A2 <- fpkm_Cyp79A2

Cyp79A2_DEGs <- Cyp79A2_DEGs[c(1:5,13,11,14,12,6,9)]
colnames(Cyp79A2_DEGs) <- c("locus", "short_name", "name", "aliases", "length",
                            "fpkm_WT", "tmm_WT", "fpkm_Cyp79A2", "tmm_Cyp79A2",
                            "log2FC_Cyp79A2", "pValue_Cyp79A2")

# Add tmm and fpkm to PAOx DEG dataframe ----
tmm_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
tmm_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
fpkm_WT <- vector(mode="numeric", length=nrow(PAOx_DEGs))
fpkm_PAOx <- vector(mode="numeric", length=nrow(PAOx_DEGs))
for(i in 1:nrow(PAOx_DEGs)){
  geneLocus <- PAOx_DEGs$locus[i]
  tmm_WT[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 8:10]))
  tmm_PAOx[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 5:7]))
  fpkm_WT[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 8:10]))
  fpkm_PAOx[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 5:7]))
}
PAOx_DEGs$tmm_WT <- tmm_WT
PAOx_DEGs$tmm_PAOx <- tmm_PAOx
PAOx_DEGs$fpkm_WT <- fpkm_WT
PAOx_DEGs$fpkm_PAOx <- fpkm_PAOx

PAOx_DEGs <- PAOx_DEGs[c(1,14,12,6,9)]
colnames(PAOx_DEGs) <- c("locus", "fpkm_PAOx", "tmm_PAOx", "log2FC_PAOx", "pValue_PAOx")

# Merge DEG dataframes and save as an RDS ----
DEGs_Cyp79A2 <- merge(Cyp79A2_DEGs, PAOx_DEGs, by="locus")
saveRDS(DEGs_Cyp79A2, file="C:/Users/bca08_000/Documents/interactive-data/RNAseq/data/DEGs_Cyp79A2.rds")