### SET-UP Codes
rm(list = ls()); gc(); library(Seurat);
outdir= "./Results/2_Markers/"
Case_01_Tu.sdj <- readRDS('./Results/1_QC/Case_01_Tu.sdj')
Case_01_Tu.sdj$tissue.id <- "Tu"
# Define function
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}
# colors
EP.color <- c('#CCCCCC', '#66FFFF', '#6699FF')
# Lung EP markers
Lamb_AL.mk <- c("CLDN18", "FOLR1", "AQP4", "PEBP4")
Lamb_BR.mk <- c("CAPS", "TMEM190", "PIFO", "SNTN")
Kims_EP.mk <- c("EPCAM", "KRT19", "CDH1", "KRT18")

### 1. LogNormalization
Case_01_Tu.snj = NormalizeData(Case_01_Tu.sdj, normalization.method = "LogNormalize", scale.factor = 10000); rm(Case_01_Tu.sdj) 

### 2. Scaling
# 1. Find variable features 
Case_01_Tu.snj <- FindVariableFeatures(Case_01_Tu.snj, selection.method = "vst", nfeatures = 2000)
v.genes <- VariableFeatures(Case_01_Tu.snj)
LabelPoints(VariableFeaturePlot(Case_01_Tu.snj), points = v.genes[1:15],repel = T)
#	2. Scaling
Case_01_Tu.ssj <- ScaleData(Case_01_Tu.snj, features = v.genes); rm(Case_01_Tu.snj)

### 3. Linear Dimensional Reduction by PCA. 
Case_01_Tu.spj <- RunPCA(Case_01_Tu.ssj, features = v.genes); rm(Case_01_Tu.ssj)
# 1. Print PCA components
print(Case_01_Tu.spj[["pca"]], dims = 1:5, nfeatures=5) 
VizDimLoadings(Case_01_Tu.spj, dims = 1:2, reduction ="pca")
DimPlot(Case_01_Tu.spj, reduction = 'pca')
DimHeatmap(Case_01_Tu.spj, dims = 1, cells = 500, balanced =TRUE)
DimHeatmap(Case_01_Tu.spj, dims = 1:15, cells = 500, balanced =TRUE)

### 4. Non-linear dimension reduction 
ElbowPlot(Case_01_Tu.spj)
Case_01_Tu.spj <- JackStraw(Case_01_Tu.spj, num.replicate = 100)
Case_01_Tu.spj <- ScoreJackStraw(Case_01_Tu.spj, dims = 1:20)
JackStrawPlot(Case_01_Tu.spj, dims = 1:20)

### 5. Clustering analysis
# 1 Determine the number
Case_01_Tu.spj <- FindNeighbors(Case_01_Tu.spj, dims = 1:10)
Case_01_Tu.spj <- FindClusters(Case_01_Tu.spj, resolution = 0.2)
# 2 Determine distance
Case_01_Tu.suj <- RunUMAP(Case_01_Tu.spj, dims = 1:10); rm(Case_01_Tu.spj)

### 6. Visualize
jpeg(filename = paste0(outdir, "UMAP_Case_01_Tu.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_Tu.suj, reduction = "umap", label = T)
dev.off()

### 7. Check EP cluster identity
jpeg(filename = paste0(outdir, "AL_Case_01_Tu.jpeg"), width = 4000, height=3000, res = 600)
FeaturePlot(Case_01_Tu.suj, cols = EP.color, features = Lamb_AL.mk, label = T)
dev.off()

jpeg(filename = paste0(outdir, "BR_Case_01_Tu.jpeg"), width = 4000, height=3000, res = 600)
FeaturePlot(Case_01_Tu.suj, cols = EP.color, features = Lamb_BR.mk, label = T)
dev.off()

jpeg(filename = paste0(outdir, "EP_Case_01_Tu.jpeg"), width = 4000, height=3000, res = 600)
FeaturePlot(Case_01_Tu.suj, cols = EP.color, features = Kims_EP.mk, label = T)
dev.off()

### 8. FindAllmarkers
# 8.1 Find AllMarker
Case_01_Tu.mk <- FindAllMarkers(Case_01_Tu.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(Case_01_Tu.mk, file = paste0(outdir, "Case_01_Tu.mk"), quote = F, row.names = FALSE, sep = "\t")
# 8.2 Get power and merge
Case_01_Tu.pw <- FindAllMarkers(Case_01_Tu.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
# 8.3 Merge
Case_01_Tu <- merge(Case_01_Tu.mk, Case_01_Tu.pw, by="row.names", all.x = F, all.y = F)
Case_01_Tu <- subset(Case_01_Tu, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
Case_01_Tu <- Case_01_Tu[rev(order(Case_01_Tu$cluster.x, Case_01_Tu$avg_log2FC.x, Case_01_Tu$power)), ]
write.table(Case_01_Tu, paste0(outdir, "Case_01_Tu.mp"), quote = F, sep = "\t")

### 9. Findmarkers
for (idx in levels(Case_01_Tu.suj)) {
  a <- FindMarkers(Case_01_Tu.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(Case_01_Tu.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("Case_01_Tu_",idx,".mp"), c)
  write.table(c, paste0(outdir, "Case_01_Tu_", idx, ".mp"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}

### 10. Take-out Case_01_Tu epithelial cell
Case_01_Tu_EP.sbj <- subset(Case_01_Tu.suj, idents = c(3, 8, 10))

### 11. Save End product
saveRDS(Case_01_Tu.suj, paste0(outdir, 'Case_01_Tu.suj'))
write.table(colnames(Case_01_Tu_EP.sbj), paste0(outdir, "Case_01_Tu_EP.bc"))
