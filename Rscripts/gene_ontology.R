
library(stringr) 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(DOSE)
library(circlize)
library(gage)
library(pathview)
library(gageData)
library(enrichplot)


# GO enrichment
DE_genes <- rownames(filtered)
GO_results_BP <- enrichGO(gene = DE_genes, 
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENSEMBL",
                          ont = "BP", readable = T)
GO_df_BP <- as.data.frame(GO_results_BP)


plot(barplot(ego, showCategory = 15, width = 0.5)) +
  theme(
    axis.text.y=element_text(size=rel(0.9), vjust = 0),
    plot.margin = margin(5, 5, 5, 90, "pt")) +
  scale_x_continuous(expand = expansion(mult = c(0, .1)))

plot(dotplot(ego, showCategory = 20, font.size = 10))


ego <- GO_results_BP
# select several pathways to visualize using cnetplot
ego@result <- ego@result[c(2,4,7,11,12,13,14,15,16,17,18,19,20),]
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=filtered$log2FoldChange, 
         vertex.label.font=6,
         node_label = 'gene',
         colorEdge = T,
         color_category = "indianred2")