# Differential gene expression analysis using edgeR package
# Goal - perform pairwise comparisons of expression data

# Install and load packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

library(edgeR)
library(tibble)

# Specify paths to data files
justCounts <- "C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/PRJNA230565/expression_counts/"
geneLengths <- "C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/PRJNA230565/gene_lengths.txt"
araport11 <- read.csv("C:/Users/bca08_000/Documents/Cyp79A2-RNAseq/data/Araport11.csv", header=FALSE)

# Load gene expression and gene length data
countFiles <- paste(justCounts, dir(justCounts), sep="")
countLabels <- c("control_rep1", "control_rep2", "NAA_rep1", "NAA_rep2")
expressionData <- readDGE(countFiles,
                          labels=countLabels,
                          group=c("A", "A", "B", "B"))
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

contrast_NAA_vs_control <- makeContrasts(B-A, levels=expDesign)
NAA_vs_control <- glmQLFTest(modelFit, contrast=contrast_NAA_vs_control)
NAA_DEGs <- topTags(NAA_vs_control, n=nrow(expressionData), p.value=1)$table
NAA_DEGs <- rownames_to_column(NAA_DEGs, var="locus")

# Load gene annotations
araport11 <- araport11[, c(1, 3, 4, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "aliases")

# Add gene annotations to DEG dataframes
NAA_DEGs <- merge(araport11, NAA_DEGs, by="locus", all.y=TRUE)

# Calculate fpkm for all genes
fpkmAllGenes <- as.data.frame(rpkm(expressionData$counts, gene.length=expressionData$genes$Length, normalized.lib.sizes=TRUE, log=FALSE))
fpkmAllGenes <- rownames_to_column(fpkmAllGenes, var="locus")

# Calculate tmm for all genes
tmmAllGenes <- as.data.frame(cpm(expressionData, log=FALSE))
tmmAllGenes <- rownames_to_column(tmmAllGenes, var="locus")

# Add tmm and fpkm to NAA treated DEG dataframe
fpkm_control <- vector(mode="numeric", length=nrow(NAA_DEGs))
fpkmSE_control <- vector(mode="numeric", length=nrow(NAA_DEGs))
tmm_control <- vector(mode="numeric", length=nrow(NAA_DEGs))
tmmSE_control <- vector(mode="numeric", length=nrow(NAA_DEGs))
fpkm_NAA <- vector(mode="numeric", length=nrow(NAA_DEGs))
fpkmSE_NAA <- vector(mode="numeric", length=nrow(NAA_DEGs))
tmm_NAA <- vector(mode="numeric", length=nrow(NAA_DEGs))
tmmSE_NAA <- vector(mode="numeric", length=nrow(NAA_DEGs))

for(i in 1:nrow(NAA_DEGs)){
  geneLocus <- NAA_DEGs$locus[i]
  
  fpkm_control[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:3]))
  fpkmSE_control[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 2:3])))/sqrt(2)
  
  tmm_control[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:3]))
  tmmSE_control[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 2:3])))/sqrt(2)
  
  fpkm_NAA[i] <- unname(rowMeans(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 4:5]))
  fpkmSE_NAA[i] <- sd(c(unlist(fpkmAllGenes[fpkmAllGenes$locus==geneLocus, 4:5])))/sqrt(2)
  
  tmm_NAA[i] <- unname(rowMeans(tmmAllGenes[tmmAllGenes$locus==geneLocus, 4:5]))
  tmmSE_NAA[i] <- sd(c(unlist(tmmAllGenes[tmmAllGenes$locus==geneLocus, 4:5])))/sqrt(2)
}

NAA_DEGs$fpkm_control <- fpkm_control
NAA_DEGs$fpkmSE_control <- fpkmSE_control

NAA_DEGs$tmm_control <- tmm_control
NAA_DEGs$tmmSE_control <- tmmSE_control

NAA_DEGs$fpkm_NAA <- fpkm_NAA
NAA_DEGs$fpkmSE_NAA <- fpkmSE_NAA

NAA_DEGs$tmm_NAA <- tmm_NAA
NAA_DEGs$tmmSE_NAA <- tmmSE_NAA



NAA_DEGs <- NAA_DEGs[c(1:5,11:18,6,9)]
colnames(NAA_DEGs) <- c("locus", "short_name", "name", "aliases", "length",
                        "fpkm_control", "fpkmSE_control", "tmm_control", "tmmSE_control",
                        "fpkm_NAA", "fpkmSE_NAA", "tmm_NAA", "tmmSE_NAA", "log2FC_NAA", "pValue_NAA")

# Limit decimal places written
for(col in 6:ncol(NAA_DEGs)){
  NAA_DEGs[,col] <- as.numeric(formatC(NAA_DEGs[,col], width=5, format="G"))
}
saveRDS(NAA_DEGs, file="C:/Users/bca08_000/Documents/interactive-data/RNAseq/data/DEGs_PRJNA230565.rds")