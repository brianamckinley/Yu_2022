## Calculate pairwise Pearson’s correlation coefficient for differentially expressed genes.
#########################################################################
setwd("...")

cor.test.mk.R <- function(expvalues, cormethod = "pearson", corcutoff = 0.1, evalcutoff = 0.05) {


  correlation <- matrix(0, nrow = 1, ncol = 4)

  iminus1 <- nrow(expvalues) - 1
  corvalues <- numeric()

  for (i in 1:iminus1) {
    print(i)
    inext <- i + 1
    if (inext == nrow(expvalues)) {
      tmp.corvalues <- mycortest(expvalues[inext,], expvalues[i,], cormethod = cormethod)
      tmp.corvalues <- cbind(rownames(expvalues)[inext],
                             c(rownames(expvalues)[i]), t(tmp.corvalues))
    } else {
      tmp.corvalues <- apply(expvalues[inext:nrow(expvalues),], 1, mycortest,
                             expvalues[i, ], cormethod)
      tmp.corvalues <- cbind(rownames(expvalues)[inext:nrow(expvalues)],
                             c(rownames(expvalues)[i]), t(tmp.corvalues))
    }

    significantcorvalues <- which(as.numeric(tmp.corvalues[, 4]) <= evalcutoff)
    if (length(significantcorvalues) >= 1) {
      tmp.corvalues <- tmp.corvalues[significantcorvalues,]
      if (class(tmp.corvalues) != "matrix") {
        tmp.corvalues <- matrix(tmp.corvalues, ncol = 4, byrow = T)
      }
      poscorvalues <- which(as.numeric(tmp.corvalues[, 3]) >= corcutoff)
      if (length(poscorvalues) >= 1) {
        tmp.corvaluesup <- tmp.corvalues[poscorvalues,]
        if (class(tmp.corvaluesup) != "matrix") {
          tmp.corvaluesup <- matrix(tmp.corvaluesup, ncol = 4, byrow = T)
        }
        correlation <- rbind(correlation, tmp.corvaluesup)
        #write.table(tmp.corvaluesup, file = fileout, quote = F, row.names = F, append = T, col.names = F, sep = ",")
      }

      negcorvalues <- which(as.numeric(tmp.corvalues[, 3]) <= (-1 * corcutoff))
      if (length(negcorvalues) >= 1) {
        tmp.corvaluesdown <- tmp.corvalues[negcorvalues,]
        if (class(tmp.corvaluesdown) != "matrix") {
          tmp.corvaluesdown <- matrix(tmp.corvaluesdown, ncol = 4, byrow = T)
        }
        correlation <- rbind(correlation, tmp.corvaluesdown)
        #write.table(tmp.corvaluesdown, file = fileout, quote = F, row.names = F, append = T, col.names = F, sep = ",")
      }
    }
  }
  correlation <- correlation[-1, ]
  rownames(correlation) <- c()

  return(correlation)
}

#########################################################################

mycortest<-function(vector1, vector2, cormethod="pearson") {

  vector1=as.numeric(vector1)
  vector2=as.numeric(vector2)
  tmp.cor=corr.test(vector1, vector2, method=cormethod, adjust="fdr")
  tmp.result=c(tmp.cor$r, tmp.cor$p)
  as.numeric(tmp.result)

}

#########################################################################

library(psych)
cond <- c("_np_vs_wb")

for (i in 1:2){
  print(i)
  output <- paste("5.1_Corr_PCC", cond, ".cor", sep="")
  input <- paste("4.1_5TPM_5DE", cond, ".txt", sep="")
  # cpm values of all genes for allcondirions
  Data <- read.table(input, header = T, sep = "\t")
  row.names(Data) <- Data[, 1]
  Cor_Network <- cor.test.mk.R(Data)
  write.table(
    Cor_Network,
    file = output,
    quote = F,
    sep = ",",
    row.names = F
  )
}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

## Calculate mutual rank from pairwise PCC generated above
#########################################################################

mutual_Rank <- function (cor_data) {

  library(data.table)
  if (class(cor_data) == "character") {
    cor_data <- fread(cor_data, header = FALSE)
  }
  cor_data <- setDT(cor_data)

  setkey(cor_data)
  cor_data <- unique(cor_data)
  cor_data <- cor_data[!duplicated(cor_data[, c("V1", "V2")])]

  Neg_cor <- cor_data[cor_data$V3 < 0, ]
  setkey(Neg_cor)
  Neg_cor <- unique(Neg_cor)

  Neg_cor_R1 <- data.frame(Neg_cor$V1)
  Neg_cor_R1 <- data.frame(Neg_cor_R1[!duplicated(Neg_cor_R1), ])
  colnames(Neg_cor_R1)[1] <- "Genes"

  Neg_cor_R2 <- data.frame(Neg_cor$V2)
  Neg_cor_R2 <- data.frame(Neg_cor_R2[!duplicated(Neg_cor_R2), ])
  colnames(Neg_cor_R2)[1] <- "Genes"

  Neg_cor_rank <- data.frame()

  for (i in Neg_cor_R1$Genes)
  {
    Samp1 <- Neg_cor[Neg_cor$V1 == i]
    Samp1 <- Samp1[order(V3), ]
    Samp1$PCC_rank_G2toG1 <- seq.int(nrow(Samp1))
    Neg_cor_rank <- rbind(Neg_cor_rank, Samp1)
  }

  Neg_cor_rank_new <- data.frame()
  for (i in Neg_cor_R2$Genes)
  {
    Samp1 <- Neg_cor_rank[Neg_cor_rank$V2 == i]
    Samp1 <- Samp1[order(V3), ]
    Samp1$PCC_rank_G1toG2 <- seq.int(nrow(Samp1))
    Neg_cor_rank_new <- rbind(Neg_cor_rank_new, Samp1)
  }

  Pos_cor <- cor_data[cor_data$V3 >= 0, ]
  setkey(Pos_cor)
  Pos_cor <- unique(Pos_cor)

  Pos_cor_R1 <- data.frame(Pos_cor$V1)
  Pos_cor_R1 <- data.frame(Pos_cor_R1[!duplicated(Pos_cor_R1), ])
  colnames(Pos_cor_R1)[1] <- "Genes"

  Pos_cor_R2 <- data.frame(Pos_cor$V2)
  Pos_cor_R2 <- data.frame(Pos_cor_R2[!duplicated(Pos_cor_R2), ])
  colnames(Pos_cor_R2)[1] <- "Genes"

  Pos_cor_rank <- data.frame()

  for (i in Pos_cor_R1$Genes)
  {
    Samp1 <- Pos_cor[Pos_cor$V1 == i]
    Samp1 <- Samp1[order(-V3), ]
    Samp1$PCC_rank_G2toG1 <- seq.int(nrow(Samp1))
    Pos_cor_rank <- rbind(Pos_cor_rank, Samp1)
  }

  Pos_cor_rank_new <- data.frame()
  for (i in Pos_cor_R2$Genes)
  {
    Samp1 <- Pos_cor_rank[Pos_cor_rank$V2 == i]
    Samp1 <- Samp1[order(-V3), ]
    Samp1$PCC_rank_G1toG2 <- seq.int(nrow(Samp1))
    Pos_cor_rank_new <- rbind(Pos_cor_rank_new, Samp1)
  }

  CorNet_ranks <- rbind(Neg_cor_rank_new, Pos_cor_rank_new)
  CorNet_ranks$Mutual_rank <- (CorNet_ranks$PCC_rank_G2toG1 *
                                 CorNet_ranks$PCC_rank_G1toG2) ^ (1 / 2)

  CorNet_ranks <- data.frame(CorNet_ranks)
  CorNet_ranks <- CorNet_ranks[, c("V1", "V2", "Mutual_rank", "V3", "V4",
                                   "PCC_rank_G2toG1", "PCC_rank_G1toG2")]
  CorNet_ranks$New <- 0
  for (i in 1:nrow(CorNet_ranks)){
    CorNet_ranks[i,8] <- (1-(CorNet_ranks[i,3]/(max(CorNet_ranks[,3])-min(CorNet_ranks[,3]))))
  }

  return(CorNet_ranks)
}

#########################################################################

require(data.table)
cond <- c("_np_vs_wb")

for (x in 1:2){
  output <- paste("5.2_Corr_MRRank_", cond, ".cor", sep="")
  input <-paste("5.1_Corr_PCC_", cond, ".cor", sep="")
  MR_Network<- mutual_Rank(input)
  write.table(MR_Network,file=output,quote=F,sep=",", row.names = FALSE)
}

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

## Calculate TF binding site enrichment for each promoter of the gene set
#########################################################################

enrich_TFBS_Family <- function(SelValsPos, Gene, AveragesFile, FoldChange){

  Enrich <- data.frame(matrix(ncol = 3))
  colnames(Enrich) <- c("Gene_ID", "TFBS_ID", "FC")

  SelValsPos_Gene <- subset(SelValsPos, SelValsPos$Sequence.ID == Gene)
  TF <- count(SelValsPos_Gene, SelValsPos_Gene$TFBS.ID)
  colnames(TF) <- c("TFBS_ID", "Count")

  for (y in 1:nrow(TF)){
    FC <- (TF[y,2])/(AveragesFile[AveragesFile$TFBS_ID == TF[y,1],][,2])
    Temp_n <- cbind(as.character(Gene),as.character(TF[y,1]),as.numeric(FC))
    colnames(Temp_n) <- c("Gene_ID", "TFBS_ID", "FC")
    Enrich <- rbind(Enrich, Temp_n)
  }

  t <-subset(Enrich, Enrich$FC > FoldChange)
  return(t)

}

#########################################################################

BiocManager::install("qvalue")
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(qvalue)

# Provide Parameters
SelValsPos <- read.table("6.1_PlantPAN1KBGenesTFIDs_bothStrand_PlantTFDBFamilyUpdated.txt",header= TRUE, sep="\t")
AveragesFile <- read.table("6.2_AvgFrequencyTFOccurencePlantPan1KbUpstream_bothStrand.txt", header= TRUE, sep="\t")
TFBS_Family <- read.csv("6.3_TFBS_Family.csv", header = T, as.is = T)
TFBS_Family <- TFBS_Family[!duplicated(TFBS_Family$Matrix_ID),]
AllDEG <- read.csv("4.1_atlas_GRN_15XDE_20TPM_genes.csv", header = T, as.is = T)
GeneList <- unique(AllDEG$Gene_ID)

# Enriched CREs for all genes in the genome based on FC > 1.1
enrich_CRE <- data.frame(matrix(data = NA, ncol = 3))
colnames(enrich_CRE) <- c("Gene_ID", "TFBS_ID", "FC")
for (i in 1:length(GeneList)){
  print(i)
  t <- enrich_TFBS_Family(SelValsPos, GeneList[i], AveragesFile, 1.1)
  enrich_CRE <- rbind(enrich_CRE, t)
}

# CRE TF Family
enrich_CRE_Family <- merge(enrich_CRE, TFBS_Family, by.x = "TFBS_ID", by.y = "Matrix_ID")
# Write output
write.csv(enrich_CRE_Family, "7.1_CRE_Enrichment_FC_1.1_All_DEG.csv", row.names = F)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Connect TFs to promoters of the dataset
#########################################################################

# TFBS of all genes
input1 <- "6.1_PlantPAN1KBGenesTFIDs_bothStrand_PlantTFDBFamilyUpdated.txt"
input2 <- "7.1_CRE_Enrichment_FC_1.1_All_DEG.csv"
input3 <- "7.2_Sbi_TF_all.txt"
input4 <- "4.1_5TPM_5DE_np_vs_wb.csv"
output <- "8.1_All_DEG_reg_intr_CRE_Enrichment_FC1.1.txt"

Motif <- read.table(input1, header=T, sep="\t")
Enriched_CRE <- read.csv(input2, header = T, as.is = T)
TF_Family <- read.table(input3, header = T, sep= "\t", as.is = T)
DEG <- read.csv(input4, header = T, as.is = T)

uniq <- unique(DEG$Gene_ID)
network <- data.frame()

for (i in 1:length(uniq)){
  print(i)
  Enriched_CRE_sub <- subset(Enriched_CRE, Enriched_CRE$Gene_ID == uniq[i])
  Motif_sub <- subset(Motif, Motif$Sequence.ID == uniq[i] & Motif$Family %in% Enriched_CRE_sub$Family)
  TF <- unique(Motif_sub[Motif_sub$Sequence.ID == uniq[i],]$Family)
  tt <- TF_Family[TF_Family$Family %in% TF,]
  df <- cbind(rep(uniq[i], times= nrow(tt)), tt)
  network <- rbind(network,df)
}

final_reg <- as.data.frame(cbind(as.character(network[,2]), rep("interacts", times= nrow(network)), as.character(network[,1])))
write.table(final_reg, output, sep="\t", quote=F, row.names = F)

#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################
#########################################################################

# Integrate PCC, MR, and TF-promoter binding to create the final GRN
#########################################################################

# Correlation files
input1 <- data.frame()
input1 <- list.files(path = paste(getwd(), sep=""), pattern = "5.2_Corr_MRRank_*", full.names = T)

# regulatory interaction files
input2 <- "8.1_All_DEG_reg_intr_CRE_Enrichment_FC1.1.txt"

# Output files
output1 <- data.frame()
output1 <- gsub("Coexpression", "GRN", input1, fixed=TRUE)
output1 <- gsub(".cor", ".txt", output1, fixed=TRUE)
output1 <- gsub("5.2_Corr_MRRank", "9_SYM_GRN_PCC_0.0_pValue_0.05_MR_0.0_BAM", output1, fixed=TRUE)

final_reg <- read.csv(file = input2, header = T, sep= "\t", as.is = T)

#########################################################################

for (i in 1:length(input1)){

  Cor1 <- read.table(input1[i], header = T, sep= ",", as.is = T)
  colnames(Cor1) <- c("Gene1", "Gene2", "MR", "Corr", "P-Value", "PCC1", "PCC2", "MR_Scale")
  Cor_TS_All <- subset(Cor1, Cor1$Corr >= 0.0 & Cor1$`P-Value` <= 0.05 & Cor1$MR_Scale >= 0.0)

  colnames(final_reg) <- c("Gene1", "rg", "Gene2")
  temp1 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene1" = "Gene1", "Gene2" = "Gene2"))
  temp2 <- dplyr::inner_join(final_reg, Cor_TS_All, by=c("Gene2" = "Gene1", "Gene1" = "Gene2"))
  Cor_TS_Final <- rbind(temp1,temp2)
  Cor_TS_Final <- Cor_TS_Final[!(duplicated(Cor_TS_Final)| duplicated(Cor_TS_Final, fromLast=TRUE)),]
  Cor_TS_Final <- Cor_TS_Final[which(Cor_TS_Final$Gene1 != Cor_TS_Final$Gene2),]
  Cor_TS_Final <- Cor_TS_Final[complete.cases(Cor_TS_Final), ]

  Cor_TS_Final$rg <- "rg"
  Cor_TS_Final$rg <- paste(Cor_TS_Final$rg, round(as.numeric(Cor_TS_Final$Corr),digits = 1), sep = "")
  Cor_TS_Final <- Cor_TS_Final[,c(1,2,3,5,6,9)]

  write.table(Cor_TS_Final, output1[i], sep="\t", quote=F, row.names = F)
}
