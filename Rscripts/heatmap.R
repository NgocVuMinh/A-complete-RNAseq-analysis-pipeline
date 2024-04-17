
# Reference: https://youtu.be/S2_FTg9kaZU?si=VkRntFVoMnL8FMBX 


# ________ Process data for HEATMAP


# only taking DEGs using the filtered dataframe
mat <- assay(rlog_out)[rownames(filtered), rownames(colData)]
colnames(mat) <- rownames(colData)
mat.scaled <- t(apply(mat, 1, scale)) # center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- colnames(mat)

# taking all genes using the deseq_result:
mat_all <- assay(rlog_out)[rownames(deseq_result), rownames(colData)]
colnames(mat_all) <- rownames(colData)
mat_all.scaled <- t(apply(mat_all, 1, scale)) # center and scale each column (Z-score) then transpose
colnames(mat_all.scaled) <- colnames(mat_all)



# ________ Getting row index that contain selected GENES:

# e.g.: chemokine
# gene_fam <- which(grepl("^(ACKR|CCR|CXC|CCL)", filtered_p$symbol, fixed = F))

# breast cancer prognosis reporters (taken from Supp table 2 DOI:10.1038/415530a )
prog_reporters <- c('ACADS', 'AD024', 'ADM', 'AKAP2', 'ALDH4', 'AP2B1', 'ASNS', 'BBC3', 'BIRC5', 'BM037', 'BM039', 'BNIP3', 'BTG2', 'BUB1', 'CA9', 'CCNB2', 'CCNE2', 'CDC25B', 'CEGP1', 'CENPA', 'CENPF', 'CFFM4', 'CIRBP', 'CKS2', 'COL4A2', 'CP', 'CTPS', 'CTSL2', 'DC13', 'DCK', 'DEGS', 'DKFZP434B168', 'DKFZP434I114', 'DKFZP564D0462', 'DKFZP586A011', 'DKFZP586F1018', 'DKFZP761H171', 'DKFZP761L0424', 'DKFZp762E1312', 'ECT2', 'ERP70', 'ESM1', 'EXT1', 'EZH2', 'FBP1', 'FBXO5', 'FGF18', 'FLJ10134', 'FLJ10156', 'FLJ10461', 'FLJ10474', 'FLJ10511', 'FLJ10549', 'FLJ10901', 'FLJ11190', 'FLJ11354', 'FLJ12150', 'FLJ12443', 'FLJ20354', 'FLJ21924', 'FLJ22341', 'FLJ22477', 'FLJ23468', 'FLT1', 'FUT8', 'GBE1', 'GCN1L1', 'GGH', 'GMPS', 'GNAZ', 'GSTM3', 'HEC', 'HMG4', 'HRB', 'HS1119D91', 'HSA250839', 'HSU54999', 'IGFBP5', 'IGFBP5', 'KIAA0008', 'KIAA0042', 'KIAA0175', 'KIAA0882', 'KIAA1104', 'KIAA1181', 'KIAA1324', 'KIAA1442', 'KIAA1683', 'KIF3B', 'KRT18', 'L13', 'L2DTL', 'LAP18', 'LCHN', 'LOC51203', 'LOC56901', 'LOC56994', 'LOC57110', 'MAD2L1', 'MAPRE2', 'MATN3', 'MCCC1', 'MCM6', 'MGAT4A', 'MGC2827', 'MIR', 'MMP9', 'MMSDH', 'MMSDH', 'MP1', 'MTMR2', 'NDRG1', 'NMB', 'NMU', 'NS1-BP', 'NSAP1', 'ORC6L', 'OXCT', 'PACE4', 'PCTK1', 'PECI', 'PEX12', 'PFKP', 'PGK1', 'PIB5PA', 'PIR', 'PK428', 'PRAME', 'PRC1', 'PRO2000', 'PSMD2', 'PSMD7', 'PTDSS1', 'QDPR', 'RAB27B', 'RAB6B', 'RAD21', 'RAI2', 'RBP3', 'RFC4', 'RNB6', 'RPS4X', 'RRM2', 'SACS', 'SEC14L2', 'SERF1A', 'SLC2A3', 'SM-20', 'SMC4L1', 'ST7', 'STK15', 'STK3', 'STK6', 'STX1A', 'TBX3-iso', 'TFRC', 'TGFB3', 'TK1', 'TMEFF1', 'TRIP13', 'UCH37', 'VEGF', 'WISP1')


# Matching exact gene symbols in the 'symbol' column 
gene_fam <- which(filtered$symbol %in% prog_reporters) # replace filtered with deseq_result if needed
gene_fam <- gene_fam[order(filtered$log2FoldChange[gene_fam])]

# get log2FC for each gene
l2_val <- as.matrix(filtered[gene_fam,]$log2FoldChange)
colnames(l2_val) <- "log2FC"

# get baseMean for each gene
mean <- as.matrix(filtered[gene_fam,]$baseMean)
colnames(mean) <- "AveExpr"

# maps values between b/w/r for min and max l2_val values
col_logFC <- colorRamp2(c(min(l2_val), max(l2_val)), c("blue3", "red2"))

# maps between 0% quantile and 75%
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red2"))
col_Z = colorRamp2(c(-2, 0, 2), c("blue3", "white", "red2"))



# __________ A heatmap for selected genes

ha <- columnAnnotation(Condition = colData[,]$condition,
                       col = list(Condition = c("co-cultured" = "grey35", "control" = "grey")),
                       show_annotation_name = F,
                       annotation_legend_param = list(
                         labels_gp = gpar(fontsize = 12),
                         title_gp = gpar(fontsize = 12, fontface = "bold"),
                         #title = "Some values",
                         legend_height = unit(4, "cm"),
                         grid_width = unit(0.8, "cm")
                         #title_position = "lefttop-rot"
                       )
                       ) #+ Legend(at = colData[,]$condition, nrow = 1)

h1 <- Heatmap(mat.scaled[gene_fam,], 
              col = col_Z,
              cluster_rows = F,
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
                labels_gp = gpar(fontsize = 12),
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                at = c(-2, -1, 0, 1, 2),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm")
                #direction = "horizontal"
                #title_position = "lefttop-rot"
              )
              ) 

h2 <- Heatmap(l2_val,
              row_labels = filtered$symbol[gene_fam],
              cluster_rows = F,
              name="Log2FC", 
              #top_annotation = ha,
              col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid
                grid.text(round(l2_val[i, j], 2), x, y,
                          gp=gpar(col="grey90"))
              },
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 12),
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                at = seq(round(min(l2_val)), round(max(l2_val)), by = 1.5),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm")
                #title_position = "lefttop-rot"
              )
              )

h3 <- Heatmap(mean, row_labels = filtered$symbol[gene_fam],
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
                grid_width = unit(0.8, "cm")
                #title_position = "lefttop-rot"
              )
              )

dev.off()
h <- h1+h2+h3
h
draw(h, heatmap_legend_side = "right", annotation_legend_side = "bottom")


