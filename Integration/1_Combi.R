# SET-UP Codes
rm(list = ls()); gc();library(Seurat); library(dplyr)
outdir = './Results/1_Combi/'

### 1. Call datasets
# Metadata
Metadata<-read.csv('./metadata.csv', header = T, row.names = 1)
# dataset
Case_01_NL.sdj <- readRDS('../Case_001/NL/Results/1_QC/Case_01_NL.sdj')
Case_01_Tu.sdj <- readRDS('../Case_001/Tu/Results/1_QC/Case_01_Tu.sdj')

Case_02_NL.sdj <- readRDS('../Case_002/NL/Results/1_QC/Case_02_NL.sdj')
Case_02_Tu.sdj <- readRDS('../Case_002/Tu/Results/1_QC/Case_02_Tu.sdj')

Case_03_NL.sdj <- readRDS('../Case_003/NL/Results/1_QC/Case_03_NL.sdj')
Case_03_Tu.sdj <- readRDS('../Case_003/Tu/Results/1_QC/Case_03_Tu.sdj')

Case_04_NL.sdj <- readRDS('../Case_004/NL/Results/1_QC/Case_04_NL.sdj')
Case_04_Tu.sdj <- readRDS('../Case_004/Tu/Results/1_QC/Case_04_Tu.sdj')

Case_05_NL.sdj <- readRDS('../Case_005/NL/Results/1_QC/Case_05_NL.sdj')
Case_05_Tu.sdj <- readRDS('../Case_005/Tu/Results/1_QC/Case_05_Tu.sdj')

Case_07_NL.sdj <- readRDS('../Case_007/NL/Results/1_QC/Case_07_NL.sdj')
Case_07_Tu.sdj <- readRDS('../Case_007/Tu/Results/1_QC/Case_07_Tu.sdj')

Case_08_NL.sdj <- readRDS('../Case_008/NL/Results/1_QC/Case_08_NL.sdj')
Case_08_Tu.sdj <- readRDS('../Case_008/Tu/Results/1_QC/Case_08_Tu.sdj')

Case_09_NL.sdj <- readRDS('../Case_009/NL/Results/1_QC/Case_09_NL.sdj')
Case_09_Tu.sdj <- readRDS('../Case_009/Tu/Results/1_QC/Case_09_Tu.sdj')

Case_10_NL.sdj <- readRDS('../Case_010/NL/Results/1_QC/Case_10_NL.sdj')
Case_10_Tu.sdj <- readRDS('../Case_010/Tu/Results/1_QC/Case_10_Tu.sdj')

Case_11_NL.sdj <- readRDS('../Case_011/NL/Results/1_QC/Case_11_NL.sdj')
Case_11_Tu.sdj <- readRDS('../Case_011/Tu/Results/1_QC/Case_11_Tu.sdj')

Case_12_NL.sdj <- readRDS('../Case_012/NL/Results/1_QC/Case_12_NL.sdj')
Case_12_Tu.sdj <- readRDS('../Case_012/Tu/Results/1_QC/Case_12_Tu.sdj')

### 2. Merge
YS_Combi.sdj <- merge(Case_01_NL.sdj, y=c(Case_01_Tu.sdj,
                                          Case_02_NL.sdj, Case_02_Tu.sdj,
                                          Case_03_NL.sdj, Case_03_Tu.sdj, 
                                          Case_04_NL.sdj, Case_04_Tu.sdj,
                                          Case_05_NL.sdj, Case_05_Tu.sdj,
                                          Case_07_NL.sdj, Case_07_Tu.sdj, 
                                          Case_08_NL.sdj, Case_08_Tu.sdj, 
                                          Case_09_NL.sdj, Case_09_Tu.sdj,
                                          Case_10_NL.sdj, Case_10_Tu.sdj,
                                          Case_11_NL.sdj, Case_11_Tu.sdj, 
                                          Case_12_NL.sdj, Case_12_Tu.sdj
                                         ))
rm(Case_01_NL.sdj, Case_01_Tu.sdj, Case_02_NL.sdj, Case_02_Tu.sdj,Case_03_NL.sdj, Case_03_Tu.sdj, Case_04_NL.sdj, Case_04_Tu.sdj, Case_05_NL.sdj, Case_05_Tu.sdj, Case_07_NL.sdj, Case_07_Tu.sdj, Case_08_NL.sdj, Case_08_Tu.sdj, Case_09_NL.sdj, Case_09_Tu.sdj, Case_10_NL.sdj, Case_10_Tu.sdj,Case_11_NL.sdj, Case_11_Tu.sdj, Case_12_NL.sdj, Case_12_Tu.sdj)

### 3. Add metadata
# gender.id
gender <- unique(Metadata$gender.id)
for (idx in gender){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(gender.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
Idents(YS_Combi.sdj); rm(a, idx, gender)
YS_Combi.sdj$gender.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# patient.id
patient <- unique(Metadata$patient.id)
for (idx in patient){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(patient.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, patient)
YS_Combi.sdj$patient.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# tissue.id
tissue <- unique(Metadata$tissue.id)
for (idx in tissue){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(tissue.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, tissue)
YS_Combi.sdj$tissue.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# pathology.id
pathology <- unique(Metadata$pathology.id)
for (idx in pathology){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(pathology.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, pathology)
YS_Combi.sdj$pathology.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# smoking.id
smoking <- unique(Metadata$smoking.id)
for (idx in smoking){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(smoking.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, smoking)
YS_Combi.sdj$smoking.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# EGFR.id
EGFR <- unique(Metadata$EGFR.id)
for (idx in EGFR){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(EGFR.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, EGFR)
YS_Combi.sdj$EGFR.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident

# stage.id
stage <- unique(Metadata$stage.id)
for (idx in stage){
  a <- WhichCells(YS_Combi.sdj, idents = row.names(Metadata %>% filter(stage.id == idx))); Idents(YS_Combi.sdj, cells = a) <- idx
}
unique(Idents(YS_Combi.sdj)); rm(a, idx, stage)
YS_Combi.sdj$stage.id <- Idents(YS_Combi.sdj)
Idents(YS_Combi.sdj) <- YS_Combi.sdj$orig.ident


### 4. QC check
VlnPlot(YS_Combi.sdj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),group.by = "orig.ident", ncol = 3) 
plot1 <- FeatureScatter(YS_Combi.sdj, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(YS_Combi.sdj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots=list(plot1, plot2)); saveRDS(YS_Combi.sdj, paste0(outdir, 'YS_Combi.sdj'))

# subset and save
YS_Combi.sdj <- subset(YS_Combi.sdj, subset=
                           orig.ident=="Case_01_NL"|orig.ident=="Case_01_Tu"|
                           orig.ident=="Case_02_NL"|orig.ident=="Case_02_Tu"|
                           orig.ident=="Case_03_NL"|orig.ident=="Case_03_Tu"|
                           orig.ident=="Case_04_NL"|orig.ident=="Case_04_Tu"|
                           orig.ident=="Case_05_NL"|orig.ident=="Case_05_Tu"|
                           orig.ident=="Case_07_NL"|orig.ident=="Case_07_Tu"|
                           orig.ident=="Case_08_NL"|orig.ident=="Case_08_Tu"|
                           orig.ident=="Case_09_NL"|orig.ident=="Case_09_Tu"|
                           orig.ident=="Case_10_NL"|orig.ident=="Case_10_Tu"|
                           orig.ident=="Case_12_NL"|orig.ident=="Case_12_Tu")

saveRDS(YS_Combi.sdj, paste0(outdir, 'YS_Combi.sdj'))

### 5. Normalization and scaling
## 5.1 Normalization
YS_Combi.snj = NormalizeData(object = YS_Combi.sdj, normalization.method = "LogNormalize", scale.factor = 10000); rm(YS_Combi.sdj) 

## 5.2 Variable feature selection and visualize
YS_Combi.snj <- FindVariableFeatures(YS_Combi.snj, selection.method = "vst", nfeatures = 2000)
v.genes <- VariableFeatures(YS_Combi.snj)
# Visualize
LabelPoints(plot = VariableFeaturePlot(YS_Combi.snj), points = v.genes[1:15],repel = T)

## 5.3 Scaling
YS_Combi.ssj <- ScaleData(object = YS_Combi.snj, features = rownames(YS_Combi.snj)); rm(YS_Combi.snj)

### 6. Linear Dimensional Reduction by PCA. 
YS_Combi.spj <- RunPCA(YS_Combi.ssj, features = v.genes); rm(YS_Combi.ssj)
# Print PCA components
print(YS_Combi.spj[["pca"]], dims = 1:5, nfeatures=5) 
# Visualize data
VizDimLoadings(YS_Combi.spj, dims = 1:2, reduction ="pca")
DimPlot(YS_Combi.spj, reduction = 'pca')
DimHeatmap(YS_Combi.spj, dims = 1, cells = 500, balanced =TRUE)
DimHeatmap(YS_Combi.spj, dims = 1:15, cells = 500, balanced =TRUE)

### 7. NON-LINEAR DIMENSION REDUCTION 
ElbowPlot(YS_Combi.spj)
YS_Combi.spj <- JackStraw(YS_Combi.spj, num.replicate = 100)
YS_Combi.spj <- ScoreJackStraw(YS_Combi.spj, dims = 1:20)
JackStrawPlot(YS_Combi.spj, dims = 1:20)

## 8. Clustering analysis
### 8.1 Determine the number of cluster
YS_Combi.spj <- FindNeighbors(YS_Combi.spj, dims = 1:10)
YS_Combi.spj <- FindClusters(YS_Combi.spj, resolution = 0.2)

### 8.2 Determine distance between cluster
YS_Combi.suj <- RunUMAP(YS_Combi.spj, dims = 1:10)

## 9. Visualize and save
# plot1
jpeg(filename = paste0(outdir, "UMAP_Combi.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Combi.suj, reduction = "umap", label = T)
dev.off()

# plot2
jpeg(filename = paste0(outdir, "UMAP_Combi_by_origin.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Combi.suj, group.by= 'orig.ident', label = T)
dev.off()

# plot3
jpeg(filename = paste0(outdir, "PCA_Plot_Combi.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(YS_Combi.suj, reduction = 'pca')
dev.off()

# save
saveRDS(YS_Combi.suj, paste0(outdir, 'YS_Combi.suj'))
