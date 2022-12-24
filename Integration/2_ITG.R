# SET-UP Codes
rm(list = ls()); gc(); library(Seurat)
outdir = "./Results/2_9_Itg/"

### 1~5 Use 2_ITG.R integrated ~ jackstraw data
YS_Itg.sjj <- readRDS('./Results/2_Itg/YS_Itg.sjj')

### 6. Determine the number of cluster
YS_Itg.snj <- FindNeighbors(YS_Itg.sjj, dims = 1:30) 
YS_Itg.snj <- FindClusters(YS_Itg.snj, resolution = 2.0)

### 7. Determine distance between clusters using UMAP
YS_Itg.suj <- RunUMAP(YS_Itg.snj, reduction = 'pca',  dims = 1:30); rm(YS_Itg.snj)

### 8. Visualize
# plot 1
jpeg(filename = paste0(outdir, "UMAP_Itg_origin.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Itg.suj, group.by= 'orig.ident', label = T)
dev.off()

# plot 2
jpeg(filename = paste0(outdir, "UMAP_Itg_cluster.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Itg.suj, label = T)
dev.off()

# plot 3
jpeg(filename = paste0(outdir, "PCA_Plot_Itg.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Itg.suj, reduction = 'pca')
dev.off()

### 9. save
saveRDS(YS_Itg.suj, paste0(outdir, 'YS_Itg.suj'))