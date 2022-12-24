### 1.Set-up
rm(list = ls()); gc()
library(Seurat);library(slingshot); library(tidyverse); library(BUSpaRse); library(tidymodels); library(scales); library(viridis); library(Matrix); library(SingleCellExperiment)

outdir = './Results/'

# 1.1 get previous work
MY.suj <- readRDS('./Results/MY.suj')

# 1.2 set color palettes
# custom palettes (color customizing)
custom_colors <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
                   '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
                   '#EE5A24','#009432','#0652DD','#9980FA','#833471',
                   '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
                   '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
                   '#2c2c54','#474787','#aaa69d','#227093','#218c74',
                   '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
                   '#b33939','#cd6133','#84817a','#cc8e35','#ccae62')

custom_colors <- c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
                   '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
                   '#EE5A24','#009432','#0652DD','#9980FA','#833471',
                   '#EA2027','#006266')

default_colors <- c('#F8766D', '#E68613', '#7CAE00', '#0CB702',
                   '#00BFC4', '#00B8E7', '#C77CFF', '#ED68ED')

# 1.4 confirm previous setting and data
MY.suj <- SetIdent(MY.suj, value = "subcell.id")

# Entire
jpeg(filename = paste0(outdir, "MY_DimPlot.jpeg"), width = 6000, height=5000, res = 600)
DimPlot(MY.suj, reduction = "umap", label = TRUE, label.size = 5, 
        pt.size = 1, repel = TRUE, seed = 0)+
  guides(color = guide_legend(override.aes = list(size=10), ncol=1) )
dev.off()

# By smoking
jpeg(filename = paste0(outdir, "MY_DimPlot_by_smoking.jpeg"), width = 2500, height=1500, res = 600)
DimPlot(MY.suj, reduction = "umap", split.by = 'smoking.id', label = FALSE,
        pt.size = 0.01, repel = TRUE) +  NoLegend()
dev.off()

# By tissue
jpeg(filename = paste0(outdir, "MY_DimPlot_by_tissue.jpeg"), width = 2500, height=1500, res = 600)
DimPlot(MY.suj, reduction = "umap", split.by = 'tissue.id', label = FALSE,
        pt.size = 0.01, repel = TRUE) +  NoLegend()
dev.off()

### 2. Modify Seurat to sce format
MY.sce <- as.SingleCellExperiment(MY.suj)

# 2.1 transport UMAP/PCA information (coordinate)
reducedDim(MY.sce, "UMAP", withDimnames=TRUE) <- MY.suj[['umap']]@cell.embeddings
reducedDim(MY.sce, "PCA", withDimnames=TRUE) <- MY.suj[['pca']]@cell.embeddings

### 3. create a slingshot data
# Labeling by metadata in seurat object (MY.suj@meta.data$subcluster.id)
sds <- slingshot(MY.sce, clusterLabels = MY.suj$subcluster.id, reducedDim = "UMAP", start.clus ="Mo-lineage", stretch = 2)

### 4. Plot
# 4.1 umap
plot(reducedDim(sds, "UMAP"), col = custom_colors[sds$subcluster.id], pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd = 1, col = 'black', show.constraints = TRUE)
lines(SlingshotDataSet(sds), type = 'l',lwd = 2, col = 'black')

# 4.2 line
lin1 <- getLineages(SlingshotDataSet(sds), start.clus = "Mo-lineage")
plot(reducedDim(sds, "UMAP"), col = custom_colors[sds$subcluster.id],  pch = 16, cex = 0.5)
lines(SlingshotDataSet(lin1), lwd = 2, col = 'black', show.constraints = TRUE)

# 4.2.1 Save line plot
jpeg(filename = paste0(outdir, "MY_Trj_line.jpeg"), width = 6000, height =4500, res = 600)
lin1 <- getLineages(SlingshotDataSet(sds), start.clus = "Mo-lineage")
plot(reducedDim(sds, "UMAP"), col = custom_colors[sds$subcluster.id],  pch = 16, cex = 0.5)
DimPlot(MY.suj, reduction = "umap", group.by = "subcell.id",
        pt.size = 0.5, label = FALSE, repel = TRUE) + NoLegend()
lines(SlingshotDataSet(lin1), lwd = 2, col = 'black', show.constraints = TRUE)
dev.off()

# 4.3 curve
lin1 <- getLineages(SlingshotDataSet(sds), start.clus = "Mo-lineage")
crv1 <- getCurves(lin1)
plot(reducedDim(sds, "UMAP"), col = custom_colors[sds$subcluster.id], pch = 16, cex =0.5)
lines(SlingshotDataSet(crv1), lwd = 2, col = 'black', show.constraints = TRUE)

# 4.3.1 Save curve plot
jpeg(filename = paste0(outdir, "MY_Trj_curv.jpeg"), width =9000, height =4500, res = 600)
lin1 <- getLineages(SlingshotDataSet(sds), start.clus = "Mo-lineage")
crv1 <- getCurves(lin1)
plot(reducedDim(sds, "UMAP"), col = custom_colors[sds$subcluster.id], pch = 16, cex =0.5)
DimPlot(MY.suj, reduction = "umap", group.by = "subcell.id",
        pt.size = 0.5, label = TRUE, repel = TRUE)  + NoLegend()
lines(SlingshotDataSet(crv1), lwd = 2, col = 'black', show.constraints = TRUE)
dev.off()