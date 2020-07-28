# Script to retrieve FPKM data for list of IAOx specific genes
setwd("C:/Users/Bryce/Documents/interactive-data/arabidopsis-RNASeq/")

load("data/RNAseq.rda")
drop_cols <- c("gene_name", "gene_family")
RNAseq <- RNAseq[ , !(names(RNAseq) %in% drop_cols)]

IAOx_specific <- read.csv("data/Copy of IAOx specific.csv", header=FALSE)
colnames(IAOx_specific) <- c("gene_id", "gene_short_name", "gene_name", "X4", "chromosome",
                             "start", "end", "X8", "X9", "X10", "X11")
IAOx_specific$gene_id[1] <- strsplit(IAOx_specific$gene_id[1], split="¿")[[1]][2]

merged <- merge(IAOx_specific, RNAseq, by="gene_id", all.x=TRUE)
write.csv(merged, file="data/IAOx_specific.csv")