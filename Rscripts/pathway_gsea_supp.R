
#### A few more ways to visualize GSEA results aside from lollipop plots


################## GO ##################

#### Lollipop SWITCH TITLE LABELS 
ggplot(gse_bp_df_sel, aes(x = NES, y = reorder(Description, NES, decreasing = F))) +
  geom_segment(aes(xend = 0, yend = Description), color = "black") +
  geom_point(aes(size = GeneRatio, fill = p.adjust), shape = 21, color = "black") +
  scale_fill_gradient(low = "brown3", high = "white") +
  scale_size(range = c(6, 14), breaks = c(0.3, 0.5, 0.8)) + # 
  #scale_x_discrete(breaks = NULL) +
  theme_minimal() +
  #opts(legend.position = "none") + 
  labs(x = "", y = "") + 
  geom_text(aes(x = 0, label = Description, hjust = ifelse(NES < 0, 0, 1))) +
  #labs(x = "Normalized enrichment score (NES)", y = NULL, size = "Gene ratio", fill = "Adjusted p-value") +
  theme(
    #axis.text.y = element_text(color = colors_common_gse),
    axis.line.x = element_blank(), #element_line(color="black"),
    axis.line.y = element_blank(),
    panel.grid.major = element_line(color="grey97"),  # Change grid lines color and size
    panel.grid.minor = element_blank(),
    #axis.line = element_blank(),
    #axis.line = element_line(color = "black"),  # Change axis line color
    #axis.text = element_text(size = 12, color = "black"),  # Change axis text size
    axis.text.y = element_blank(),
    axis.title = element_text(size = 12, color = "black", face = "bold"), #face = "bold"),  # Change axis title size and style
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10, color = "black"),  # Change legend text size
    legend.title = element_text(size = 12, face = "bold", color = "black")  # Change legend title size and style
  )



### bubble plot
ggplot(gse_bp_df_sel, aes(x = NES, y = reorder(Description, NES, decreasing = F))) +
  geom_point(aes(size = GeneRatio, fill = p.adjust), alpha = 0.9, shape = 21) +
  scale_size(range = c(6, 14), breaks = c(0.3, 0.5, 0.8)) + # 
  scale_fill_gradient(low = "brown3", high = "white") +
  scale_y_discrete(labels = label_wrap(50)) +
  labs(x = "Gene Count", y = "Description", size = "Count", fill = "p.adjust") +
  theme_minimal() +
  theme_minimal() +
  theme(
    #axis.text.y = element_text(color = colors_common_gse),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"),  # Change grid lines color and size
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),  # Change axis line color
    axis.text = element_text(size = 14, color = "black"),  # Change axis text size
    axis.title = element_text(size = 14, color = "black", face = "bold"), #face = "bold"),  # Change axis title size and style
    axis.title.y = element_blank(),
    legend.text = element_text(size = 14, color = "black"),  # Change legend text size
    legend.title = element_text(size = 14, face = "bold", color = "black")  # Change legend title size and style
  ) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, size = "Gene ratio", fill = "Adjusted p-value")

# Normal bar chart, drop gene count
ggplot(gse_bp_df_sel, # save 1050x380
       aes(x = reorder(Description, NES, decreasing=F), y = NES, fill=p.adjust)) +
  geom_col() +
  coord_flip()

ggplot(gse_bp_df_sel, # save 1050x380
       aes(x = reorder(Description, p.adjust, decreasing=T), y = NES, fill=p.adjust)) +
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
    panel.grid.major = element_line(color="grey97"),  # Change grid lines color and size
    #panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Change axis line color
    axis.text = element_text(size = 12, color = "black"),  # Change axis text size
    axis.title = element_text(size = 12, color = "black"), #face = "bold"),  # Change axis title size and style
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10, color = "black"),  # Change legend text size
    legend.title = element_text(size = 12, face = "bold", color = "black")  # Change legend title size and style
  )



################## KEGG ##################

############ BUBBLE PLOT
ggplot(gse_kegg_df_sel, aes(x = setSize, y = reorder(Description, NES, decreasing = F), size = NES, fill = p.adjust)) +
  geom_point(alpha = 0.9, shape = 21) +
  scale_size(range = c(3, 12), breaks = c(5, 10, 15, 20)) +  # Adjust bubble size range
  scale_fill_gradient(low = "red2", high = "white", breaks = c(0,0.01,0.025, 0.05, 0.1)) +
  labs(x = "Gene Count", y = "Description", size = "Count", fill = "p.adjust") +
  theme_minimal() +
  #labs(x = "Gene ratio", fill = "Adjusted p-value") +
  theme(
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.text = element_text(size = 12, color = "black"),  
    axis.title = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10, color = "black"),  
    legend.title = element_text(size = 12, face = "bold", color = "black")  
  )


############ BAR PLOT
ggplot(gse_kegg_df_sel, # save 1000x700
       aes(x = reorder(Description, NES, decreasing = F), y = NES, fill = setSize)) +
  geom_col() +
  coord_flip() + 
  scale_fill_continuous(low="deepskyblue2", high="grey10") +
  #scale_y_discrete(expand = expansion(mult = c(0, .1))) +
  #scale_x_discrete(expand = expansion(mult = c(0, .1))) +
  labs(x = "KEGG pathway", y = "Normalized enrichment score (NES)", fill = "Gene count") +
  theme_minimal() +
  theme(
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    #panel.grid.major = element_line(color="grey"), 
    #panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 10),  
    legend.title = element_text(size = 12, face = "bold")  
  )


################## HALLMARKS ##################

############ BAR PLOT
ggplot(fgseaRes_sel, 
       aes(x = reorder(pathway, NES, decreasing = F), y = NES, fill = padj)) +
  geom_col() +
  coord_flip() + 
  scale_fill_continuous(low="darkslateblue", high="#bebad9") + #breaks=c(4e-25,1e-3)
  #scale_y_discrete(expand = expansion(mult = c(0, .1))) +
  #scale_x_discrete(expand = expansion(mult = c(0, .1))) +
  labs(x = "Hallmark gene sets", y = "Normalized enrichment score (NES)", fill = "Adjusted p-value") +
  theme_minimal() +
  theme(
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color="black", size = 12), 
    axis.title = element_text(color="black", size = 12, face = "bold"), 
    axis.title.y = element_blank(),
    legend.text = element_text(color="black", size = 10), 
    legend.title = element_text(color="black", size = 12, face = "bold")  
  )