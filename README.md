# Yu_et_al._2022_Bionergy-sorghum-stem-intercalary-meristem-dynamics
This is a repository of the scripts that were used to conduct gene regulatory network analysis for the manuscript entitled: Bioenergy sorghum stem intercalary meristem localization, imaging, and developmental transcriptome dynamics.

The script works with gene level data only, not transcript level data. If your data is in transcripts, I have also included another R script for summing TPM values or count values in a dataset. 

Next, to slim down the set of genes to a more manageable size, differential expression analysis with edgeR was performed. The edgeR script used is included in the repository.

'0.0_transcript_to_gene_TPM_sum.R' - R script to sum transcript-level TPMs and counts to gene-level TPMs and counts.

'0.1_SYM_edgeR.R' - R script to calculate differential expression using the edgeR algorithm.

'0.3_SYM_GRN_workflow.R' - R script that contains the code to create the sorghum gene regulatory network.

'3.2_SYM_read_counts_total_dataset_DE_genes.csv' - dataset of reads counts for all sorghum genes across the samples of the SYM dataset.

'6.2_AvgFrequencyTFOccurencePlantPan1KbUpstream_bothStrand.txt' - text file that contains the average frequency of occurrence of each annotated transcription factor binding site. 

'6.3_TFBS_Family.csv' - .csv file that contains a matrix ID associated with TF family.

'7.2_Sbi_TF_all.txt' - list of sorghum transcription factors.

'SYM_np_vs_wb_DE_GRN.cys' - Cytoscape file that where the network can be visualized and queried.

