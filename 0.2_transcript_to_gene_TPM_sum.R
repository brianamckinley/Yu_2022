library(dplyr)
setwd("/Users/brianmckinley/desktop/sorghum/72_SYM/")
transcript_tpm <- read.csv("3.1_SYM_read_counts_total_dataset_input for summing.csv", header = T, as.is = T)
gene_tpm <- group_by(transcript_tpm, Gene_ID) %>% summarise(across(everything(), sum))
write.csv(gene_tpm, "3.2_SYM_read_counts_total_dataset_DE_genes.csv", row.names = FALSE)








