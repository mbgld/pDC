# set-up
library(pheatmap); library(Seurat); library(dplyr)
outdir = "./Results/5_PrettyHeatmap/"
YS_Itg.sfj <- readRDS("./Results/4_Assign/Sim_2_9/YS_Itg.sfj")

# Take out marker genes
YS_Itg_All.mk <- FindAllMarkers(YS_Itg.sfj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
YS_Itg_Top3.mk <- YS_Itg_All.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# get average matrix
YS_Itg1.avg <- AverageExpression(YS_Itg.sfj)

# generate scaled RNA expression matrix and get primitive graph 
YS_Itg2.avg <- as.data.frame(YS_Itg1.avg[["RNA"]][YS_Itg_Top3.mk$gene, ]); mode(YS_Itg2.avg)

# To exaggerate 
YS_Itg3.avg  <- t(scale(t(YS_Itg2.avg)))

########################################################
# End of number matrix
########################################################

# 1. Legends for cluster
column.lgd <- as.data.frame(colnames(YS_Itg3.avg), row.names = colnames(YS_Itg3.avg))
colnames(column.lgd) <- 'Cell type'

# 2. Row legends
df_row.lgd <- as.data.frame(row.names(YS_Itg3.avg), row.names = row.names(YS_Itg3.avg))
colnames(df_row.lgd) <- 'Marker gene'

# 2.1 Marker genes
MY.gen <- c("IFI30", "C1QA", "LYZ")
NKT.gen <- c("CCL5", "NKG7", "IL32")
MA.gen <- c("TPSB2", "TPSAB1", "CPA3")
BC.gen <- c("IGKC", "CD79A", "IGLC2")
EC.gen <- c("SPARCL1", "RAMP2", "CALCRL")
FB.gen <- c("DCN", "MGP", "SPARCL1")
EP.gen <- c("SLPI", "WFDC2", "CAPS")
CA.gen <- c("NAPSA", "SFTPB", "KRT8")

# 2.2 MARKER GENE GROUPING FOR ANNOTATION
df_row.lgd$'Marker genes' <- 
  ifelse(df_row.lgd$`Marker gene` %in% MY.gen,"MY",
    ifelse(df_row.lgd$`Marker gene` %in% MA.gen,"MA",
      ifelse(df_row.lgd$`Marker gene` %in% BC.gen,"BC",
        ifelse(df_row.lgd$`Marker gene` %in% NKT.gen,"NK/T",
          ifelse(df_row.lgd$`Marker gene` %in% EC.gen,"EC",
            ifelse(df_row.lgd$`Marker gene` %in% FB.gen,"FB",
              ifelse(df_row.lgd$`Marker gene` %in% EP.gen,"EP",
                 "CA")))))))

# Extraction for annotation because options require data.frame format! Stupid Work!
row.lgd <- df_row.lgd %>% dplyr::select('Marker genes')

# 3. Set color #'DC'='#00008B',
my_color <- list(
  'Marker genes' = c('MY'='#0EBE43', 'NK/T'='#06BBE3', 'EP'='#90A7FF',
                     'FB'='#F87C74', 'MA'='#92A900', 'EC'='#D79D19',
                     'BC'='#00C19F', 'CA'='#FF7ACC'),
  'Cell type' =    c('MY'='#0EBE43', 'NK/T'='#06BBE3', 'EP'='#90A7FF',
                     'FB'='#F87C74', 'MA'='#92A900', 'EC'='#D79D19',
                     'BC'='#00C19F', 'CA'='#FF7ACC'))

# End products
# Figure 1
jpeg(filename = paste0(outdir, "Heatmap_Major_cells.jpeg"), width = 3500, height =4000, res = 600)
pheatmap(YS_Itg3.avg, 
         annotation_col = column.lgd,
         annotation_row = row.lgd,
         annotation_colors = my_color,
         cutree_cols = F, cluster_rows = F, cluster_cols = F,
         gaps_col = c(1,2,3,4,5,6,7),
         gaps_row = c(3,6,9,12,15,18,21),
         # labels_row = as.expression(row.lgd),
         annotation_legend = T
)
dev.off()
