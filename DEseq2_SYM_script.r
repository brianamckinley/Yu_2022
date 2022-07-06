if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
library(DESeq2)

# Read in the raw read counts ** Don't lable row names... row names are just numbers for this code
rawCounts1 <- read.csv("~/R_SYM Read Counts_Deanna.csv")
rawCounts1 <- na.omit(rawCounts1)

# Read in the sample mappings
sampleData1 <- read.delim("~/SYM_data.txt", row.names= NULL)

# Also save a copy for later
sampleData_v2 <- sampleData1

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts1$transcriptIDV3
sampleIndex <- grepl("6.", colnames(rawCounts1))
rawCounts <- as.matrix(rawCounts1[,sampleIndex]) #creates matrix by removing row names
rownames(rawCounts) <- geneID #adds Gene IDs as row.names to the matrix

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData1) <- sampleData1$Run #Assigns "Run" column from sampleData to rownames
keep <- c("Phytomer","Tissue") # combines elements "" to form a vector "keep" **** CHANGE THIS ****
sampleData <- sampleData1[,keep] # matches elements from keep vector to column names in sampleData, and creates new sampleDate with just those columns
colnames(sampleData) <- c("phytomer","tissue") # changes column names
sampleData$tissue <- factor(sampleData$tissue) # encodes "tissue" vector as a factor
include_list <- c(colnames(rawCounts)) #creates list of column names, to subset only the tissues of interest
sampleData <- subset(sampleData, rownames(sampleData) %in% include_list ) #subsets sample data to only include information about tissues of interest

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))] #reorders column names of count data to be same as sampleData
all(colnames(rawCounts) == rownames(sampleData)) # check if column  names are the same
sampleData$tissue <- factor(sampleData$tissue, levels=c("internode_lower","internode_upper", "NP", "pulvinus", "WB"))

# Create the DEseq2DataSet object ***** formula follows: design ~ specifies formula and + separates factors
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ tissue)

################################ Pre-filtering of Data
# reduces size of DEseq2DataSet object, speeds up run time

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 3, ])

# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 3, ]

################################ Set up multi-cores (optional)

# Install and load the library
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel", version = "3.8")

# Register the number of cores to use
library(BiocParallel)
register(MulticoreParam(2))

################################ Differential Expression Analysis choose one of the three

# 1. Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)

SYM_P6_data <- deseq2Data

################################ Extracting Results

# Extract differential expression results
# For "tissueType" perform primary vs normal comparison
SYM_6.3v6.4_results <- results(SYM_P6_data, contrast=c("tissue","NP","internode_upper")) #***CHANGE***

# View summary of results
summary(SYM_6.3v6.4_results)

# Writes output file for results as a csv file
write.csv(SYM_6.3v6.4_results ,file="SYM_6.3v6.4_results.csv")

################################ Pipeline to pull out relevant TPM data
library(dplyr)
deseq2results <- read.csv("~/SYM P6 DEseq/SYM_6.3v6.4_results.csv", stringsAsFactors=FALSE)
## To subset DE output file with FC cut-off
FCcull1 <- mutate(deseq2results, FC= (2^(abs(log2FoldChange))))
FCcull <- subset(FCcull1, FC > 5)
## To subset DE output file with FDR and p-value cut offs
## Downregulated at high planting density
## log2FoldChange >= 1 means upregulated, <= 1 means downregulated in 6.3
FDRcull <- subset(FCcull, padj < 0.05 & pvalue < 0.05 & log2FoldChange > 0)
rownames(FDRcull) <- FDRcull$X
write.csv(FDRcull, file="FDRcutoff_6.3v6.4_up.csv") #*Use this if you want to make a list of the transcripts of interest

## To subset TPM file using outputs from DE
total_dataD <- read.csv("~/R_SYM_Full_annotated_dataset_TPM (Copy) (1).csv", row.names=1)
## Compares the TranscriptIDs from the FDR cutoff file to the TPM file,
#then creates a new file of TPM values of transcripts that are DE and up or down regulated depending on the log2FoldChange
include_list <- c(row.names(FDRcull))
total_data_cull <- subset(total_dataD, rownames(total_dataD) %in% include_list )
colnames(total_data_cull)

## To subset TPM >5
total_data_cullup <- subset(total_data_cull,
                            SYM_Top_1cm_of_Internode_P6.3_TPMmean > 3
)
write.csv(total_data_cullup, file="totaldatacull_SYM_6.3v6.4_results_enrichedin6.3.csv")

################################ Pipeline to pull out relevant TPM data
library(dplyr)
deseq2results <- read.csv("~/SYM P6 DEseq/SYM_6.3v6.4_results.csv", stringsAsFactors=FALSE)
## To subset DE output file with FC cut-off
FCcull1 <- mutate(deseq2results, FC= (2^(abs(log2FoldChange))))
FCcull <- subset(FCcull1, FC > 5)
## To subset DE output file with FDR and p-value cut offs
## Downregulated at high planting density
## log2FoldChange >= 1 means upregulated, <= 1 means downregulated in 6.3
FDRcull <- subset(FCcull, padj < 0.05 & pvalue < 0.05 & log2FoldChange < 0)
rownames(FDRcull) <- FDRcull$X
write.csv(FDRcull, file="FDRcutoff_6.3v6.4_down.csv") #*Use this if you want to make a list of the transcripts of interest

## To subset TPM file using outputs from DE
# SAS_total_dataD <- read.csv("~/Desktop/SAS_data analysis/17_TX08001_SAS_means_reps_counts_BAM.csv", row.names=1)
total_dataD <- read.csv("~/R_SYM_Full_annotated_dataset_TPM (Copy) (1).csv", row.names=1)
## Compares the TranscriptIDs from the FDR cutoff file to the TPM file,
#then creates a new file of TPM values of transcripts that are DE and up or down regulated depending on the log2FoldChange
include_list <- c(row.names(FDRcull))
total_data_cull <- subset(total_dataD, rownames(total_dataD) %in% include_list )
colnames(total_data_cull)

## To subset TPM >5
total_data_culldown <- subset(total_data_cull,
                              SYM_Upper_White_Band_of_P6.4_TPMmean > 3
)

write.csv(total_data_culldown, file="totaldatacull_SYM_6.3v6.4_results_enrichedin6.4.csv")
