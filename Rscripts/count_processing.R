
library(DESeq2)
library(tidyverse)
library(airway)
library(apeglm)
library(pheatmap)
library(RColorBrewer)
library(stringr) 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gage)
library(pathview)
library(gageData)
library(enrichplot)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]


# ____________ This block of code is to load and save variables ____________
# setwd("/Users/ngoc/ADSC/results") # set working directory to read and store files
# load("/Users/ngoc/ADSC/run.RData")
# save.image("/Users/ngoc/ADSC/run.RData")

# read the count file produced by featureCounts
counts_data <- read.csv('./featureCounts/featureCounts_merged.txt', sep="\t", row.names = 1)

# read sample info
colData <- read.csv('./sample_info.csv', row.names = 1)

# making sure the rownames in colData match column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
# are they in the same order?
all(colnames(counts_data) == rownames(colData))

# construct DESeq object
dds0 <- DESeqDataSetFromMatrix(countData = counts_data,
                               colData = colData,
                               design = ~ condition)

dds0 # 62710

# keeping rows that have at least 10 reads per sample (60 reads in total)
keep <- rowSums(counts(dds0)) >= 60
dds <- dds0[keep,]
dds # 19643

# contrasts:
resultsNames(dds)

# set the control sample
dds$condition <- relevel(dds$condition, ref="control")
dds$condition

# get statistics
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res)
plotMA(res)
dev.off()

# getting differentially expressed genes (DEGs)
deseq_result <- as.data.frame(res)
deseq_result <- deseq_result[order(deseq_result$pvalue),]
filtered <- deseq_result %>% 
              filter((padj < 0.05) & (abs(log2FoldChange) > 1))

# getting gene symbol, gene name and entrez ID
filtered$symbol <- mapIds(org.Hs.eg.db, keys=rownames(filtered),keytype = "ENSEMBL", column = 'SYMBOL')
filtered$name <- mapIds(org.Hs.eg.db, keys=rownames(filtered), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
filtered$entrez <- mapIds(org.Hs.eg.db, keys=rownames(filtered), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")

deseq_result$symbol <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result),keytype = "ENSEMBL", column = 'SYMBOL')
deseq_result$name <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
deseq_result$entrez <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")

resLFC$symbol <- mapIds(org.Hs.eg.db, keys=rownames(resLFC),keytype = "ENSEMBL", column = 'SYMBOL')
resLFC$name <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
resLFC$entrez <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")

# getting normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# dispersion plot
plotDispEsts(dds)
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup="condition")

# remove noise
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC <- as.data.frame(resLFC)
plotMA(resLFC, ylim=c(-2,2))


