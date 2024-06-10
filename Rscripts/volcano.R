
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


