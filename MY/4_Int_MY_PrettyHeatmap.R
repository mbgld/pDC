### 1.Set-up
rm(list = ls()); gc(); outdir = './Results/'
library(pheatmap); library(Seurat); library(dplyr)
MY.sfj <- readRDS('./Results/MY.sfj')

### 2. Generate marker genes
# 2.1 Take out marker genes
MY_All.mk <- read.csv("./Results/MD_markers.tsv", sep = "\t", row.names = NULL)
# MY_All.mk <- MY_All.mk[order(MY_All.mk$avg_log2FC), ]
MY_All.mk <- MY_All.mk[-which(duplicated(MY_All.mk$gene)), ]
MY_All.mk$subcell.id <- gsub("-", "_", MY_All.mk$subcell.id)
MY_All.mk$subcell.id <- gsub("/", "_", MY_All.mk$subcell.id)
MY_Top3.mk <- MY_All.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

# 2.2 Assign marker genes
subcell.id <- unique(MY_Top3.mk$subcell.id)
for (idx in subcell.id) {
  a <- MY_Top3.mk[MY_Top3.mk$subcell.id==idx, ]$gene
  assign(paste0(idx, ".gen"), a)
  rm(a)
}

### 3. Generate number matrix
# 3.1 Get average matrix
MY1.avg <- AverageExpression(MY.sfj)

# 3.2 Scaled RNA expression matrix
MY2.avg <- as.data.frame(MY1.avg[["RNA"]][MY_Top3.mk$gene, ]); mode(MY2.avg)

# 3.3 Exaggeration
MY3.avg  <- t(scale(t(MY2.avg)))

# 3.4 Reorder dataframe
row_order<- c("Mo-lineage", "FCN1-mono", "Mono-Mc", "CD163/LGMN", 
              "Alv-Mc3", "Alv-Mc2", "Alv-Mc1", "Alv-Mc5", "Alv-Mc4",
              "Prol-Mc", "Alv-Mc6", "Alv-Mc8", "Alv-Mc7",
              "mo-DC", "cDC1", "cDC2", "pDC")
row_order <-gsub("-", "_", row_order); row_order <-gsub("/", "_", row_order)

MY4.avg = data.frame()

for (idx in row_order){
  a = MY3.avg[row.names(MY3.avg) %in% get(paste0(idx, ".gen")), ]
  MY4.avg = rbind(MY4.avg, a)
}

### 4. Decorate
# 4.1. Legends for cluster
column.lgd <- as.data.frame(colnames(MY4.avg), row.names = colnames(MY4.avg))
colnames(column.lgd) <- 'Cell type'

# 4.2. Row legends
df_row.lgd <- as.data.frame(row.names(MY4.avg), row.names = row.names(MY4.avg))
colnames(df_row.lgd) <- 'Marker gene'

### 5. MARKER GENE GROUPING FOR ANNOTATION
df_row.lgd$'Marker genes' <- 
  ifelse(df_row.lgd$`Marker gene` %in% Mo_lineage.gen,"Mo_lineage",
   ifelse(df_row.lgd$`Marker gene` %in% FCN1_mono.gen,"FCN1_mono",
    ifelse(df_row.lgd$`Marker gene` %in% Mono_Mc.gen,"Mono_Mc",
     ifelse(df_row.lgd$`Marker gene` %in% CD163_LGMN.gen,"CD163/LGMN",
      ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc3.gen,"Alv_Mc3",       
       ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc2.gen,"Alv_Mc2",
        ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc1.gen,"Alv_Mc1",
         ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc5.gen,"Alv_Mc5",
          ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc4.gen,"Alv_Mc4",
           ifelse(df_row.lgd$`Marker gene` %in% Prol_Mc.gen ,"Prol_Mc",
            ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc6.gen,"Alv_Mc6",
             ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc8.gen,"Alv_Mc8",
              ifelse(df_row.lgd$`Marker gene` %in% Alv_Mc7.gen,"Alv_Mc7",
               ifelse(df_row.lgd$`Marker gene` %in% mo_DC.gen,"mo_DC",
                ifelse(df_row.lgd$`Marker gene` %in% cDC1.gen,"cDC1",
                 ifelse(df_row.lgd$`Marker gene` %in% cDC2.gen,"cDC2",
                 "pDC"))))))))))))))))

# Extraction for annotation because options require data.frame format! Stupid Work!
row.lgd <- df_row.lgd %>% dplyr::select('Marker genes')


# End products
# Figure 1
jpeg(filename = paste0(outdir, "Heatmap_MY_cells.jpeg"), width = 4500, height =8000, res = 600)
pheatmap(MY4.avg, 
         annotation_col = column.lgd,
         annotation_row = row.lgd,
         # annotation_colors = my_color,
         cutree_cols = F, cluster_rows = F, cluster_cols = F,
         gaps_col = 1:16,
         gaps_row = seq(3, 48, 3),
         # labels_row = as.expression(row.lgd),
         annotation_legend = T
)
dev.off()
