
# Ref tutorial Heatmap: https://www.youtube.com/watch?v=ht1r34-ifVI&list=PLi1VnGoeDGjvHvl83QySD2oAQYFHPRYso&index=7
# Ref tutoral Volcano: https://youtu.be/vRr78s37CI4?si=ePEzagdpSYtroa6s 
# Ref tutorial GSEA plot: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html 


#BiocManager::install("clusterProfiler")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("EnhancedVolcano")
# BiocManager::install("ComplexHeatmap")
#install.packages("stringr") 

library(stringr) 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)



# ___________ Volcano plots ___________ 

dev.off()
options(repr.plot.width = 2, repr.plot.height =3) 
# help("EnhancedVolcano")
EnhancedVolcano(resLFC, 
                x = "log2FoldChange", 
                y = "padj",
                lab = resLFC$symbol,
                title = '',
                subtitle = '',
                caption = '', #paste0(nrow(resLFC), " genes in total"),
                titleLabSize = 16, subtitleLabSize = 16, captionLabSize = 16,
                pointSize = 6,
                labSize = 4, axisLabSize = 16,
                xlim = c(-4, 10), ylim = c(-3, 170),
                legendLabSize = 14,
                legendIconSize = 7,
                pCutoff = 1e-2,
                cutoffLineWidth = 0.4,
                shape = 20,
                gridlines.major = F,
                gridlines.minor = T, 
                #drawConnectors = T,
                #legendLabSize = 12,
                #legendIconSize = 4,
                legendPosition = "right",
                legendLabels = c("Not significant", "|Log2FC| > 1", "P-value < 0.05", "P-value < 0.05 and |Log2FC| > 1")
                )


# ___________ GSEA plot for GO biological processes (bp) subontology ___________ 

gsea_stat <- res[order(-res$stat),]
gene_list <- gsea_stat$stat
names(gene_list) <- rownames(gsea_stat)

gse_bp <- gseGO(gene_list,
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)
View(as.data.frame(gse_bp)) # examine this dataframe to select GO pathway for gseaplot
gseaplot2(gse_bp, 
         geneSetID = c("GO:0001666"), # select a GO pathway
         color = "red",
         pvalue_table = TRUE)



# ___________ GSEA plot for KEGG ___________ 


kegg_gene_list <- resLFC$log2FoldChange
names(kegg_gene_list) <- resLFC$entrez
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

gse_kegg <- gseKEGG(kegg_gene_list,
                    organism = "hsa",
                    keyType = "ncbi-geneid",
                    nPermSimple = 10000,
                    eps = 0,
                    minGSSize    = 3,
                    maxGSSize    = 800,
                    pvalueCutoff = 0.9,
                    pAdjustMethod = "BH")

gse_kegg_df <- gse_kegg@result # examine this dataframe to select KEGG pathway for gseaplot
View(gse_kegg_df)
gseaplot2(gse_kegg, 
          geneSetID = c("hsa04066", # select KEGG pathways
                        "hsa00010",
                        "hsa00190"), 
          color = c("royalblue2", "red3", "darkgoldenrod2"),
          pvalue_table = F,
          base_size = 14)







