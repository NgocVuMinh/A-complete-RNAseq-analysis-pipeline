
# Reference:
# https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/ 


# BiocManager::install("gage")
# BiocManager::install("gageData")
# BiocManager::install("pathview")

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
library(org.Hs.eg.db)
data(kegg.sets.hs) # is a named list of 229 elements, each is a character vector of member gene Entrez IDs for a single KEGG pathway
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]



# __________ GAGE: Getting top-regulated KEGG pathways ___________

foldchanges <- resLFC$log2FoldChange
names(foldchanges) <- resLFC$entrez
kegg <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Get pathways related to specific genes only
# emt_ID <- which(grepl("^(VIM|COL6A|FYN|ANXA6|OLFM2|ACVR1|PMP22|PCOLCE|TIMP2|COL5A1|COL1A|COL3A|VCAN|POSTN|SPARC|FN1)", resLFC$symbol, fixed = F))
# foldchanges = resLFC[emt_ID,]$log2FoldChange

# Get the top 10 upregulated pathways
keggup <- data.frame(id=rownames(kegg$greater), kegg$greater) %>% 
  as_tibble() %>%
  filter(row_number()<=50) %>% 
  .$id %>% 
  as.character()

# Get the top 10 downregulated pathways
keggdown <- data.frame(id=rownames(kegg$less), kegg$less) %>% 
  as_tibble() %>%
  filter(row_number()<=50) %>% 
  .$id %>% 
  as.character()

# Get the IDs for plotting upregulated
keggupids <- substr(keggup, start=1, stop=8)
keggdownids <- substr(keggdown, start=1, stop=8)

# define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, 
                                      species="hsa", new.signature=FALSE)

# plot and save multiple pathways
tmp_up <- sapply(keggupids, function(pid) pathview(gene.data=foldchanges, 
                                                  pathway.id=pid, species="hsa"))
tmp_down <- sapply(keggdownids, function(pid) pathview(gene.data=foldchanges, 
                                                      pathway.id=pid, species="hsa"))





# __________ KEGG pathways ___________

DE_genes_kegg <- unlist(mget(DE_genes, envir=org.Hs.egENSEMBL2EG, ifnotfound = NA))
kegg <- enrichKEGG(gene = DE_genes_kegg,
                   organism = 'hsa',
                   #universe = geneUniverse,
                   pvalueCutoff = 0.05,
                   qvalueCutoff)

kegg_df <- as.data.frame(kegg@result)
kegg_df$gene_symbol <- fetch_gene_symbols(kegg_df$geneID, "ENTREZID")
kegg_pairism <- pairwise_termsim(kegg)
emapplot(kegg_pairism)  # this produce a dot plot but no indication of up/down regulated
cnetplot(kegg)

# had to use keggdownids and keggresids as index to filter them out and
# visualize them ndividually
kegg_df <- as.data.frame(kegg@result)
View(kegg_df[keggupids, ])

kegg_down <- kegg
kegg_up <- kegg
kegg_down@result <- kegg_df[keggdownids, ]
kegg_up@result <- kegg_df[keggupids, ]
dotplot(kegg, showCategory=30)
dotplot(kegg_up, showCategory=20)
# dotplot(kegg_down, showCategory=20)

# Function to fetch gene symbols for Gene IDs
fetch_gene_symbols <- function(my_ids, id_type) { # id_type is either ENTREZID or ENSEMBL
  id_list <- unlist(strsplit(my_ids, "/"))   # split multiple NCBI Gene IDs in a cell separated by '/'
  gene_symbols <- mapIds(org.Hs.eg.db, keys = id_list, column = "SYMBOL", keytype = id_type)
  symbol_combined <- sapply(strsplit(my_ids, "/"), function(ids) paste(gene_symbols[ids], collapse = "/"))
  return(symbol_combined)
}


# _____________ Visualizing upregulated KEGG using geom_dotplot 
# _____________ instead of dotplot for more flexibility

dummy <- kegg_df[keggupids, ]
dummy$GeneRatiof <- sapply(strsplit(dummy$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
dummy <- dummy[order(dummy$Count, decreasing=TRUE), ]
dummy$gene_symbol <- fetch_gene_symbols(dummy$geneID, "ENTREZID")

# Create the dot plot using geom_dotplot
ggplot(filter(dummy, (dummy$p.adjust < 0.16) & (dummy$Count >4)), aes(x = GeneRatiof, y = reorder(Description, Count), size = Count, fill = p.adjust)) +
  geom_point(alpha = 0.9, shape = 21) +
  scale_size(range = c(3, 12), breaks = c(5, 10, 15, 20)) +  # Adjust bubble size range
  scale_fill_gradient(low = "red2", high = "white", breaks = c(0,0.01,0.025, 0.05, 0.1, 0.15)) +
  labs(x = "GeneRatiof", y = "Description", size = "Count", fill = "p.adjust") +
  theme_minimal() +
  labs(x = "Gene ratio") +
  theme(
    #panel.grid.major = element_line(color = "gray", linewidth = 0.5),  # Change grid lines color and size
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Change axis line color
    axis.text = element_text(size = 12),  # Change axis text size
    axis.title = element_text(size = 12),#face = "bold"),  # Change axis title size and style
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10),  # Change legend text size
    legend.title = element_text(size = 12, face = "bold")  # Change legend title size and style
  )

