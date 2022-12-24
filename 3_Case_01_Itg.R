### SET-UP Codes
rm(list = ls()); gc(); library(Seurat);
outdir= "./Results/3_Itg/"
Case_01_Tu.sdj <- readRDS('/media/yschang/V/Case_001/Tu/Results/1_QC/Case_01_Tu.sdj')
Case_01_Tu.sdj$tissue.id <- "Tu"
Case_01_NL.sdj <- readRDS('/media/yschang/V/Case_001/NL/Results/1_QC/Case_01_NL.sdj')
Case_01_NL.sdj$tissue.id <- "NL"

### 1. Merge
Case_01_Combi.sdj <- merge(Case_01_NL.sdj, y = Case_01_Tu.sdj); rm(Case_01_NL.sdj, Case_01_Tu.sdj)

### 1. Dataset Preprocessing 
Case_01_dataset.list <- SplitObject(Case_01_Combi.sdj, split.by = 'orig.ident')
Case_01_dataset.list <- lapply(Case_01_dataset.list, FUN = function(x){
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})
Case_01.anchors <- FindIntegrationAnchors(object.list = Case_01_dataset.list, dims=1:30)

### 2. Integration
Case_01_Itg.sbj <- IntegrateData(anchorset=Case_01.anchors, dims=1:30)
DefaultAssay(Case_01_Itg.sbj) <- 'integrated'; rm(Case_01_dataset.list, Case_01.anchors)

### 3. Scaling: .sbj starting point
Case_01_Itg.ssj <- ScaleData(Case_01_Itg.sbj, features = rownames(Case_01_Itg.sbj))

### 4. Linear dimensional reduction
Case_01_Itg.spj <- RunPCA(Case_01_Itg.ssj, npcs = 30, verbose = FALSE); rm(Case_01_Itg.ssj); print(Case_01_Itg.spj[["pca"]], dims = 1:5, nfeatures=5)
# Visualize data
VizDimLoadings(Case_01_Itg.spj, dims = 1:2, reduction ="pca"); DimPlot(Case_01_Itg.spj, reduction = 'pca'); DimHeatmap(Case_01_Itg.spj, dims = 1, cells = 500, balanced =TRUE); DimHeatmap(Case_01_Itg.spj, dims = 1:15, cells = 500, balanced =TRUE)

### 5.  Non-linear dimension reduction
ElbowPlot(Case_01_Itg.spj)
Case_01_Itg.spj <- JackStraw(Case_01_Itg.spj, num.replicate = 100)
Case_01_Itg.spj <- ScoreJackStraw(Case_01_Itg.spj, dims = 1:20)
JackStrawPlot(Case_01_Itg.spj, dims = 1:20)

### 6. Determine number
Case_01_Itg.spj <- FindNeighbors(Case_01_Itg.spj, dims = 1:10) 
Case_01_Itg.spj <- FindClusters(Case_01_Itg.spj, resolution = 0.2)

### 7. Determine distance
Case_01_Itg.suj <- RunUMAP(Case_01_Itg.spj, reduction = 'pca',  dims = 1:10)

### 8. Visualize
# plot 1
jpeg(filename = paste0(outdir, "UMAP_Case_01_Itg_origin.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_Itg.suj, group.by= 'orig.ident', label = T)
dev.off()

# plot 2
jpeg(filename = paste0(outdir, "PCA_Plot_Case_01_Itg.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_Itg.suj, reduction = 'pca')
dev.off()

### 9. Save
saveRDS(Case_01_Itg.suj, paste0(outdir, 'Case_01_Itg.suj'))
