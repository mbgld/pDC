### 1.Set-up
rm(list = ls()); gc()
library(Seurat); library(dplyr); outdir = './Results/'

# 1.1 get previous work
YS_Itg.sfj <- readRDS("/media/yschang/V/Combi/Results/4_Assign/Sim_2_9/YS_Itg.sfj")
MY.sbj <- readRDS('/media/yschang/V/Combi/Results/6_Split/Sim_2_9/MY.sbj')
# 1.2 Color and markers
MY.color <- c('#CCCCCC', '#99FF00', '#009900')
# 1.3 Markers
My_MY.mk <- c("LYZ", "MARCO", "CD68", "FCGR3A", "PTPRC")

# 1.4 Check MY identity
DimPlot(YS_Itg.sfj, cells.highlight = WhichCells(MY.sbj), cols.highlight = '#009900', reduction = "umap", label = T)
FeaturePlot(YS_Itg.sfj, cols = MY.color, features = My_MY.mk)

### 2. Re-clustering
DefaultAssay(MY.sbj) <- 'RNA'
MY.sbj <- FindVariableFeatures(MY.sbj, selection.method = "vst", nfeatures = 2000, assay = 'RNA'); v.genes <- VariableFeatures(MY.sbj)
DefaultAssay(MY.sbj) <- 'integrated'
MY.ssj <- ScaleData(MY.sbj, features = rownames(MY.sbj)); rm(MY.sbj)
MY.spj <- RunPCA(MY.ssj, features = v.genes)

# 2.1 Dimension reduction
ElbowPlot(MY.spj)
MY.spj <- JackStraw(MY.spj, num.replicate = 100); MY.spj <- ScoreJackStraw(MY.spj, dims = 1:20); JackStrawPlot(MY.spj, dims = 1:20)
MY.suj <- RunUMAP(MY.spj, dims = 1:30, seed.use = 0)

# 2.2 Clustering
MY.suj <- FindNeighbors(MY.suj, dims = 1:30)
MY.suj <- FindClusters(MY.suj, resolution = 2)
DimPlot(MY.suj, reduction = "umap", label = T, group.by = 'orig.ident', repel = T)
DimPlot(MY.suj, reduction = "umap", label = T, pt.size = 0.1)
DimPlot(MY.suj, reduction = "umap", label = F, group.by = 'tissue.id', pt.size = 0.1)
MY.suj <- SetIdent(MY.suj, value = 'subcell.id')

# 2.3 DimPlot
jpeg(filename = paste0(outdir, "MD_suj_DimPlot.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(MY.suj, reduction = "umap",  label = TRUE, label.size = 2, label.color = "black", repel = TRUE)
dev.off()

### 4. Save
saveRDS(MY.suj, paste0(outdir,"MY.suj"))