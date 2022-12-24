# SET-UP Codes
rm(list = ls()); gc(); library(Seurat)
outdir = "./Results/3_Markers/Sim_2_9/"
YS_Itg.suj <- readRDS('./Results/2_9_Itg/YS_Itg.suj')

# Define function
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

### 1. Check previous clustering: Hiracheal clustering and ClusterTree
DimPlot(YS_Itg.suj, reduction = "umap", label = T)
PlotClusterTree(BuildClusterTree(YS_Itg.suj))

### 2. Set default assay as RNA for DEG
DefaultAssay(YS_Itg.suj) <- 'RNA'

### 3. FindAllmarkers
YS_Itg_All.mk <- FindAllMarkers(YS_Itg.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(YS_Itg_All.mk, file = paste0(outdir, "/YS_Itg_All_mk.tsv"), quote = F, row.names = FALSE, sep = "\t")

YS_Itg_All.pw <- FindAllMarkers(YS_Itg.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
YS_Itg_All <- merge(YS_Itg_All.mk, YS_Itg_All.pw, by="row.names", all.x = F, all.y = F)
YS_Itg_All <- subset(YS_Itg_All, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
YS_Itg_All <- YS_Itg_All[rev(order(YS_Itg_All$cluster.x, YS_Itg_All$avg_log2FC.x, YS_Itg_All$power)), ]
write.table(YS_Itg_All, file = paste0(outdir, "/YS_Itg_All_mp.tsv"), quote = F, sep = "\t")

### 4. Find each cluster markers. 
for (idx in levels(YS_Itg.suj)) {
  a <- FindMarkers(YS_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(YS_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("YS_Itg_",idx,".mp"), c)
  write.table(c, file = paste0(outdir, "YS_Itg_", idx, "mp.tsv"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}

### 5. Save