library("edgeR")
library("gtools")
setwd("...")
# import data
raw.data <- read.csv("3.2_SYM_read_counts_total_dataset_DE_genes.csv", header = TRUE, sep = ",")
ncol(raw.data)
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,14,14,14,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22))
y <- DGEList(counts=raw.data[,2:58], genes=raw.data[, 1], group=group)
levels(y$samples$group)
keep <- rowSums(cpm(y) > .5) > 2
summary(keep)
y <- y[keep, ,keep.lib.sizes=FALSE]
plotMDS(y, top=25000, main = "MDS of Count Data", cex=.75, labels = group, pch = 19)
y <- calcNormFactors(y, norm.method = "tmm", test.method = "edger", na.rm = TRUE)
barplot(y$samples$lib.size*1e-6, xpd = TRUE, axes = TRUE, col = 491, border = 300, xlab = "RNA-seq Libraires", ylab="Reads aligned to exonic sequence (millions)")
design <- model.matrix(~0+group)
y <- estimateDisp(y, design, robust=TRUE, tagwise=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, winsor.tail.p=c(0.05, 0.1), robust=TRUE)
plotQLDisp(fit, cex = 0.2)# Perform qausi-likelihood (QL) F-test to test for differentail expression


qlf1 <- glmQLFTest(fit, contrast = c(0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0))
qlf1$table$FC <- logratio2foldchange(qlf1$table$logFC, base=2)
A <- topTags(qlf1, n=47206, sort.by="PValue", p.value=0.05)
write.table(A, file = ".csv", sep = ",")
