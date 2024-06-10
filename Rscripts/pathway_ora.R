# tutorial GO enrichment: https://www.youtube.com/watch?v=JPwdqdo_tRg&list=PLi1VnGoeDGjvHvl83QySD2oAQYFHPRYso&index=11

library(stringr) 
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(topGO)


######################## Overrepresentation GO analysis using only DEGs ########################

DE_genes <- rownames(filtered)

# biological process
ora_go_bp <- enrichGO(gene = DE_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
ora_go_bp_df <- as.data.frame(ora_go_bp)
ora_go_bp_df$genes <- fetch_gene_symbols(ora_go_bp_df$geneID, "ENSEMBL")
ora_go_bp_df_sel <- ora_go_bp_df[c("GO:0001666", "GO:0036293", "GO:0030198"), ]
ora_go_bp_df_sel <- ora_go_bp_df_sel[order(ora_go_bp_df_sel$p.adjust, decreasing=F), ]

# molecular function
ora_go_mf <- enrichGO(gene = DE_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
ora_go_mf_df <- as.data.frame(ora_go_mf)
ora_go_mf_df$genes <- fetch_gene_symbols(ora_go_mf_df$geneID, "ENSEMBL")
ora_go_mf_df_sel <- ora_go_mf_df[c("GO:0005201", "GO:0048029", "GO:0005539"),]
ora_go_mf_df_sel <- ora_go_mf_df_sel[order(ora_go_mf_df_sel$p.adjust, decreasing=F), ]

# cellular function
ora_go_cc <- enrichGO(gene = DE_genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")
ora_go_cc_df <- as.data.frame(ora_go_cc)
ora_go_cc_df$genes <- fetch_gene_symbols(ora_go_cc_df$geneID, "ENSEMBL")
ora_go_cc_df_sel <- ora_go_cc_df[c("GO:0062023", "GO:0005604", "GO:0031983"), ]
ora_go_cc_df_sel <- ora_go_cc_df_sel[order(ora_go_cc_df_sel$p.adjust, decreasing=F), ]


ggplot(ora_go_bp_df_sel, # save 1050x380
       aes(x = reorder(Description, p.adjust, decreasing=T), y = factor(Count), fill=p.adjust)) +
  geom_col() + #color="black"
  coord_flip() + 
  #scale_fill_gradient(low = "brown3", high = "#FFC3B0") +
  #scale_fill_gradient(low = "deepskyblue4", high = "skyblue2") +
  scale_fill_gradient(low = "seagreen", high = "#aae3b9") + # also change x_discrete to mult = c(.15, .15)
  scale_y_discrete(expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(mult = c(.1, .1))) +
  labs(x = "Pathway", y = "Gene hit", fill = "Adjusted p-value") +
  theme_minimal() +
  theme(
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"), 
    #panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 12, color = "black", face = "bold"), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10, color = "black"),  
    legend.title = element_text(size = 12, face = "bold", color = "black")  
  )


plotGOgraph(ora_go_cc)




######################## Overrepresentation KEGG analysis using only DEGs ########################


DE_genes_entrez <- filtered$entrez
#DE_genes_entrez_up <- dplyr::filter(filtered, filtered$log2FoldChange>0)$entrez
#DE_genes_entrez_down <- dplyr::filter(filtered, filtered$log2FoldChange<0)$entrez

ora_kegg <- enrichKEGG(gene = DE_genes_entrez, 
                       organism = "hsa", 
                       keyType = "ncbi-geneid")
ora_kegg_df <- as.data.frame(ora_kegg)
ora_kegg_df$genes <- fetch_gene_symbols(ora_kegg_df$geneID, "ENTREZID")
ora_kegg_df_sel <- ora_kegg_df[c("hsa04066", "hsa04512", "hsa04820"),]
#write_csv(ora_kegg_df_sel, file="ora_kegge_df_selected.csv")


ggplot(ora_kegg_df_sel, # save 1050x380
       aes(x = reorder(Description, p.adjust, decreasing=T), y = factor(Count), fill=p.adjust)) +
  geom_col() +
  coord_flip() + 
  scale_fill_gradient(low = "#004b63", high = "deepskyblue3") +
  scale_y_discrete(expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(mult = c(.1, .1))) +
  labs(x = "Pathway", y = "Gene count", fill = "Adjusted p-value") +
  theme_minimal() +
  theme(
    #axis.text.y = element_text(color = colors_common),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"),  
    #panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 12, color = "black"), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10, color = "black"), 
    legend.title = element_text(size = 12, face = "bold", color = "black") 
  )




