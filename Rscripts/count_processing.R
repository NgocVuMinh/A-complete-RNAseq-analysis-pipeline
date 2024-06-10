
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
library(scales) 

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

dds0 # 60675

# keeping rows that have at least 10 reads per sample (60 reads in total)
# keeping genes that meet the following minimum count (at least 3 samples have count > 50)
dds <- dds0[rowSums(counts(dds0) >= 20) >= 3,]
dds 

# contrasts:
resultsNames(dds)

# set the control sample
dds$condition <- relevel(dds$condition, ref="control")
dds$condition

# getting normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# get statistics
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05) # setting alpha threshold = 0.05
summary(res)

# dispersion plot
dev.off()
plotMA(res)
plotDispEsts(dds)
vsd <- vst(dds) # , blind = FALSE
plotPCA(vsd, intgroup="condition") +
  geom_point(size = 5) +
  scale_color_manual(values = c("control" = "black", "co-cultured" = "red")) +
  theme_minimal() +
  theme(
    axis.text.y   = element_text(size=12, colour = "black"),
    axis.text.x   = element_text(size=12, colour = "black"),
    axis.title.y  = element_text(size=12, colour = "black"),
    axis.title.x  = element_text(size=12, colour = "black"),
    legend.text = element_text(size=11),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5)
  )

##### (NO SHRINKAGE) getting differentially expression results 
deseq_result <- as.data.frame(res)
deseq_result <- deseq_result[order(deseq_result$pvalue),]
# getting gene symbol, gene name and entrez ID
deseq_result$symbol <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result),keytype = "ENSEMBL", column = 'SYMBOL')
deseq_result$name <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
deseq_result$entrez <- mapIds(org.Hs.eg.db, keys=rownames(deseq_result), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")
# removing all uncharacterized genes (no identified GENE SYMBOL)
dim(deseq_result)
deseq_result_clean <- deseq_result[is.na(deseq_result$symbol) == F, ]
dim(deseq_result_clean) # 14669     9


##### (SHRINKAGE) getting differentially expression results 
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2])
resLFC <- as.data.frame(resLFC)
# getting gene symbol, gene name and entrez ID
resLFC$symbol <- mapIds(org.Hs.eg.db, keys=rownames(resLFC),keytype = "ENSEMBL", column = 'SYMBOL')
resLFC$name <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
resLFC$entrez <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")
# removing all uncharacterized genes (no identified GENE SYMBOL)
dim(resFLC)
resLFC_clean <- resLFC[is.na(resLFC$symbol) == F, ]
dim(resLFC_clean)

# get uncharacterized genes for further inspection
unchar_resLFC <- resLFC[is.na(resLFC$symbol) == T, ]
unchar_resLFC_DE <- dplyr::filter(unchar_resLFC, (unchar_resLFC$padj < 0.05) & (abs(unchar_resLFC$log2FoldChange) > 1))
write.csv(unchar_resLFC_DE, file="uncharacterized_genes.csv")


##### Get differentially expressed genes (DEGs)
filtered <- resLFC_clean %>% 
  dplyr::filter((resLFC_clean$padj < 0.05) & (abs(resLFC_clean$log2FoldChange) > 1))
filtered.noshrink <- deseq_result_clean %>% 
  dplyr::filter((deseq_result_clean$padj < 0.05) & (abs(deseq_result_clean$log2FoldChange) > 1))


dim(filtered)
dim(filtered[filtered$log2FoldChange>0,]) # upregulated DEGs
dim(filtered[filtered$log2FoldChange<0,]) # downregulated DEGs

