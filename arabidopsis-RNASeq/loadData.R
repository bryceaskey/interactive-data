library(stringr)

# Script to load data from .csv file and save as a .rda file
setwd("C:/Users/Bryce/Documents/interactive-data/arabidopsis-RNAseq")

RNAseq_raw <- read.csv("data/RNAseq_REF5_BB.csv")[ , 2:25]
geneFamilies <- read.csv("data/gene_families_sep_20_09_update.csv")
geneFamilies$Genomic_Locus_Tag <- str_to_upper(geneFamilies$Genomic_Locus_Tag)
colnames(geneFamilies)[1] <- "Gene_Family"

# Function to retrieve gene family based on gene id
getGeneFamily <- function(geneID, geneFamilyData=geneFamilies){
  rowNum <- grep(geneID, geneFamilies$Genomic_Locus_Tag)
  if(length(rowNum)==0){
    return(NA)
  }else{
    geneFamily <- geneFamilyData$Gene_Family[rowNum]
    return(geneFamily)
  }
}

# Mutate/transmute functions from dplyr don't work with mean/sd functions
# Mean and sd need to be calculated with for loop
gene_family <- vector(mode="character", length=nrow(RNAseq_raw))
WT_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
WT_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref2_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref2_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Med5_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Med5_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5Med5_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5Med5_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref2Med5_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref2Med5_se_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))

for(i in 1:nrow(RNAseq_raw)){
  gene_family[i] <- getGeneFamily(RNAseq_raw$gene_id[i])
  
  WT_mean_FPKM[i] <- mean(c(RNAseq_raw$X011351_WT.1_FPKM[i], RNAseq_raw$X011352_WT.2_FPKM[i], RNAseq_raw$X011353_WT.3_FPKM[i]))
  WT_se_FPKM[i] <- sd(c(RNAseq_raw$X011351_WT.1_FPKM[i], RNAseq_raw$X011352_WT.2_FPKM[i], RNAseq_raw$X011353_WT.3_FPKM[i]))/sqrt(3)
  
  Ref5_mean_FPKM[i] <- mean(c(RNAseq_raw$X011354_Ref5.1_FPKM[i], RNAseq_raw$X011355_Ref5.2_FPKM[i], RNAseq_raw$X011356_Ref5.3_FPKM[i]))
  Ref5_se_FPKM[i] <- sd(c(RNAseq_raw$X011354_Ref5.1_FPKM[i], RNAseq_raw$X011355_Ref5.2_FPKM[i], RNAseq_raw$X011356_Ref5.3_FPKM[i]))/sqrt(3)
  
  Ref2_mean_FPKM[i] <- mean(c(RNAseq_raw$X011357_Ref2.1_FPKM[i], RNAseq_raw$X011358_Ref2.2_FPKM[i], RNAseq_raw$X011359_Ref2.3_FPKM[i]))
  Ref2_se_FPKM[i] <- sd(c(RNAseq_raw$X011357_Ref2.1_FPKM[i], RNAseq_raw$X011358_Ref2.2_FPKM[i], RNAseq_raw$X011359_Ref2.3_FPKM[i]))/sqrt(3)
  
  Med5_mean_FPKM[i] <- mean(c(RNAseq_raw$X011360_Med5.1_FPKM[i], RNAseq_raw$X011361_Med5.2_FPKM[i], RNAseq_raw$X011362_Med5.3_FPKM[i]))
  Med5_se_FPKM[i] <- sd(c(RNAseq_raw$X011360_Med5.1_FPKM[i], RNAseq_raw$X011361_Med5.2_FPKM[i], RNAseq_raw$X011362_Med5.3_FPKM[i]))/sqrt(3)
  
  Ref5Med5_mean_FPKM[i] <- mean(c(RNAseq_raw$X011363_Ref5Med5.1_FPKM[i], RNAseq_raw$X011364_Ref5Med5.2_FPKM[i], RNAseq_raw$X011365_Ref5Med5.3_FPKM[i]))
  Ref5Med5_se_FPKM[i] <- sd(c(RNAseq_raw$X011363_Ref5Med5.1_FPKM[i], RNAseq_raw$X011364_Ref5Med5.2_FPKM[i], RNAseq_raw$X011365_Ref5Med5.3_FPKM[i]))/sqrt(3)
  
  Ref2Med5_mean_FPKM[i] <- mean(c(RNAseq_raw$X011366_Ref2Med5.1_FPKM[i], RNAseq_raw$X011367_Ref2Med5.2_FPKM[i], RNAseq_raw$X011368_Ref2Med5.3_FPKM[i]))
  Ref2Med5_se_FPKM[i] <- sd(c(RNAseq_raw$X011366_Ref2Med5.1_FPKM[i], RNAseq_raw$X011367_Ref2Med5.2_FPKM[i], RNAseq_raw$X011368_Ref2Med5.3_FPKM[i]))/sqrt(3)
}

set_dec_places <- function(num, dec_places){
  if(is.na(num)){
    return(NA)
  }else{
    return(as.numeric(trimws(format(round(num, dec_places), nsmall=dec_places))))
  }
}
  
RNAseq <- data.frame(
  gene_id=RNAseq_raw$gene_id,
  gene_name=RNAseq_raw$gene_short_name,
  gene_family=factor(gene_family),
  
  WT_mean_FPKM=sapply(WT_mean_FPKM, set_dec_places, dec_places=3),
  WT_se_FPKM=sapply(WT_se_FPKM, set_dec_places, dec_places=3),
  
  Ref5_mean_FPKM=sapply(Ref5_mean_FPKM, set_dec_places, dec_places=3),
  Ref5_se_FPKM=sapply(Ref5_se_FPKM, set_dec_places, dec_places=3),
  
  Ref2_mean_FPKM=sapply(Ref2_mean_FPKM, set_dec_places, dec_places=3),
  Ref2_se_FPKM=sapply(Ref2_se_FPKM, set_dec_places, dec_places=3),
  
  Med5_mean_FPKM=sapply(Med5_mean_FPKM, set_dec_places, dec_places=3),
  Med5_se_FPKM=sapply(Med5_se_FPKM, set_dec_places, dec_places=3),
  
  Ref5Med5_mean_FPKM=sapply(Ref5Med5_mean_FPKM, set_dec_places, dec_places=3),
  Ref5Med5_se_FPKM=sapply(Ref5Med5_se_FPKM, set_dec_places, dec_places=3),
  
  Ref2Med5_mean_FPKM=sapply(Ref2Med5_mean_FPKM, set_dec_places, dec_places=3),
  Ref2Med5_se_FPKM=sapply(Ref2Med5_se_FPKM, set_dec_places, dec_places=3)
)

save(RNAseq, file="data/RNAseq.rda")