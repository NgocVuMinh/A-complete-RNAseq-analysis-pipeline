# Reference: https://youtu.be/S2_FTg9kaZU?si=VkRntFVoMnL8FMBX 


# ________ Process data for HEATMAP
# ________ The following heatmap is for visualizing top DEGs

# only taking DEGs using the filtered dataframe
rlog_out <- rlog(dds, blind=FALSE) # get normalized counts
mat <- assay(rlog_out)[rownames(filtered), rownames(colData)]
colnames(mat) <- rownames(colData)
mat.scaled <- t(apply(mat, 1, scale)) # center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- colnames(mat)

# taking all genes using the deseq_result:
mat_all <- assay(rlog_out)[rownames(deseq_result), rownames(colData)]
colnames(mat_all) <- rownames(colData)
mat_all.scaled <- t(apply(mat_all, 1, scale)) # center and scale each column (Z-score) then transpose
colnames(mat_all.scaled) <- colnames(mat_all)

# Get 50 top DE:
num_keep <- 50 # keep 25 highest and lowest log2FoldChange -> total 50
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled) - num_keep), 
                                    nrow(mat.scaled) ))

# get log2FC for each gene
l2_val <- as.matrix(filtered[rows_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"
# get baseMean for each gene
mean <- as.matrix(filtered[rows_keep,]$baseMean)
colnames(mean) <- "AveExpr"

col_logFC <- colorRamp2(c(min(l2_val), -0.8, max(l2_val)), c("mediumseagreen", "white", "red3"))
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red2"))
col_Z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# the heatmap:
ha <- columnAnnotation(Condition = colData[,]$condition,
                       col = list(Condition = c("co-cultured" = "grey35", "control" = "grey")),
                       show_annotation_name = F,
                       annotation_legend_param = list(
                         labels_gp = gpar(fontsize = 10),
                         title_gp = gpar(fontsize = 10, fontface = "bold"),
                         #at = c("co-cultured", "control"),
                         #labels = c("low", "zero", "high"),
                         #title = "Some values",
                         legend_height = unit(4, "cm"),
                         grid_width = unit(0.8, "cm")
                         #title_position = "leftbot-rot"
                       )) #+ Legend(at = colData[,]$condition, nrow = 1)

h1 <- Heatmap(mat.scaled[rows_keep,], 
              col = col_Z,
              cluster_rows = T,
              cluster_columns = T, 
              column_labels = colnames(mat.scaled), 
              name="Z-score",
              show_row_names = F,
              show_column_names = F,
              #left_annotation = ha_rows,
              bottom_annotation = ha,
              row_title_gp = gpar(fontsize = 3.5),
              column_names_gp = gpar(fontsize = 12),
              column_title_gp = gpar(fontsize = 12),
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 10),
                title_gp = gpar(fontsize = 10, fontface = "bold"),
               # at = c(-2, -1, 0, 1, 2),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop"
              )
) 


h2 <- Heatmap(l2_val,
              row_labels = filtered$symbol[rows_keep],
              row_names_gp = gpar(fontsize = 9),
              cluster_rows = T,
              name="Log2FC", 
              #top_annotation = ha,
              col = col_logFC,
              #cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid
              #  grid.text(round(l2_val[i, j], 2), x, y,
              #            gp=gpar(col="grey90"))
              #},
              width = unit(5, "mm"),
              column_names_rot = 90,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 10),
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                at = seq(round(min(l2_val)), round(max(l2_val)), by = 2),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                # legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop-rot"
              )
)


h3 <- Heatmap(mean, row_labels = filtered$symbol[rows_keep],
              cluster_rows = F,
              name = "AveExpr",
              col = col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid
                grid.text(round(mean[i, j], 2), x, y)
              },
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 12),
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                #at = seq(0, round(max(mean)), by = 5000),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop-rot"
              )
)

dev.off()
h1+h2+h3 # save at width 900 x height 1600
h1 + h2
