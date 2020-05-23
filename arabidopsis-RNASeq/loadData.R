library(stringr)

# Script to load data from .csv file and save as a .rda file
setwd("C:/Users/Bryce/Documents/interactive-data/arabidopsis-RNAseq")

RNAseq_raw <- read.csv("data/RNAseq_REF5_BB.csv")[ , 2:11]
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

#gene_family <- vector(mode="character", length=nrow(RNAseq_raw))
#for(i in 1:nrow(RNAseq_raw)){
#  gene_family[i] <- getGeneFamily(RNAseq_raw$gene_id[i])
#}
#RNAseq_raw$gene_family <- gene_family
#RNAseq_raw$gene_family <- as.factor(RNAseq_raw$gene_family)


# Mutate/transmute functions from dplyr don't work with mean/sd functions
# Mean and sd need to be calculated with for loop
WT_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
WT_sd_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5_mean_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
Ref5_sd_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))
ratio_FPKM <- vector(mode="numeric", length=nrow(RNAseq_raw))

for(i in 1:nrow(RNAseq_raw)){
  WT_mean_FPKM[i] <- mean(c(RNAseq_raw$X011351_WT.1_FPKM[i], RNAseq_raw$X011352_WT.2_FPKM[i], RNAseq_raw$X011353_WT.3_FPKM[i]))
  WT_sd_FPKM[i] <- sd(c(RNAseq_raw$X011351_WT.1_FPKM[i], RNAseq_raw$X011352_WT.2_FPKM[i], RNAseq_raw$X011353_WT.3_FPKM[i]))
  Ref5_mean_FPKM[i] <- mean(c(RNAseq_raw$X011354_Ref5.1_FPKM[i], RNAseq_raw$X011355_Ref5.2_FPKM[i], RNAseq_raw$X011356_Ref5.3_FPKM[i]))
  Ref5_sd_FPKM[i] <- sd(c(RNAseq_raw$X011354_Ref5.1_FPKM[i], RNAseq_raw$X011355_Ref5.2_FPKM[i], RNAseq_raw$X011356_Ref5.3_FPKM[i]))
  if(Ref5_mean_FPKM[i]==0){
    ratio_FPKM[i] <- NA
  }else{
    ratio_FPKM[i] <- log10(WT_mean_FPKM[i]/Ref5_mean_FPKM[i])
  }
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
  WT_mean_FPKM=sapply(WT_mean_FPKM, set_dec_places, dec_places=5),
  WT_sd_FPKM=sapply(WT_sd_FPKM, set_dec_places, dec_places=5),
  Ref5_mean_FPKM=sapply(Ref5_mean_FPKM, set_dec_places, dec_places=5),
  Ref5_sd_FPKM=sapply(Ref5_sd_FPKM, set_dec_places, dec_places=5),
  ratio_FPKM=sapply(ratio_FPKM, set_dec_places, dec_places=5)
)

save(RNAseq, file="data/RNAseq.rda")