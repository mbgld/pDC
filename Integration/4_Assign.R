### 0. Set-up
rm(list = ls()); gc(); library(Seurat); library(ggplot2)
outdir = "./Results/4_Assign/Sim_2_9/"
YS_Itg.suj <- readRDS('./Results/2_Itg/Sim_2_9/YS_Itg.suj')
Cell_id <- read.csv('./Results/4_Assign/Sim_2_9/Cell_id.csv')

### 1. Assign celltype id
YS_Itg.sfj <- YS_Itg.suj
# 1.1. Keep current active.id. as cluster.id
YS_Itg.sfj$cluster.id <- Idents(YS_Itg.sfj)

### 2. CA
# 2.1 Take out barcode
case <- c('01', '02', '03', '04', '05', '07', '08', '09', '10', '12')
for (idx in case){
  assign('CA_barcode', paste0('/media/yschang/T/Case_0', idx,'/Combi/Results/Case_',idx,'_CA.bc'))
  a <- readLines(CA_barcode)
  a <- a[-1]
  a <-gsub('\\.', '-', a)
  assign(paste0("Case_0", idx, "_CA.bc"), a)
  }
# 2.2 Cancer cell barcodes
YS_CA.bc<- c(Case_001_CA.bc, Case_002_CA.bc, Case_003_CA.bc, Case_004_CA.bc, Case_005_CA.bc, Case_007_CA.bc, Case_008_CA.bc, Case_009_CA.bc, Case_010_CA.bc, Case_012_CA.bc); rm(Case_001_CA.bc, Case_002_CA.bc, Case_003_CA.bc, Case_004_CA.bc, Case_005_CA.bc, Case_007_CA.bc, Case_008_CA.bc, Case_009_CA.bc, Case_010_CA.bc, Case_012_CA.bc, a, CA_barcode, case, idx)

### 3. Naming celltype.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, cells = YS_CA.bc, value = "CA")
celltype.id <- Cell_id$celltype.id
names(celltype.id) <- levels(YS_Itg.sfj)[-1]
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, celltype.id)
DimPlot(YS_Itg.sfj, reduction = "umap", pt.size = 0.1, label = T)
# 3.1 Keep current active.id. as subcluster.id and reset active.id as cluster.id
YS_Itg.sfj$celltype.id <- Idents(YS_Itg.sfj); Idents(YS_Itg.sfj) <- YS_Itg.sfj$cluster.id

### 4. Naming subcluster.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, cells = YS_CA.bc, value = "CA")
subcluster.id <- Cell_id$subcluster.id
names(subcluster.id) <- levels(YS_Itg.sfj)[-1]
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, subcluster.id)
DimPlot(YS_Itg.sfj, reduction = "umap", pt.size = 0.1, label = T, repel = T)
# 4.1 Keep current active.id. as subcluster.id and reset active.id as subcluster.id
YS_Itg.sfj$subcluster.id <- Idents(YS_Itg.sfj); Idents(YS_Itg.sfj) <- YS_Itg.sfj$cluster.id

### 5. Naming subcell.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, cells = YS_CA.bc, value = "CA")
subcell.id <- Cell_id$subcell.id
names(subcell.id) <- levels(YS_Itg.sfj)[-1]
YS_Itg.sfj <- RenameIdents(YS_Itg.sfj, subcell.id)
DimPlot(YS_Itg.sfj, reduction = "umap", pt.size = 0.1, label = T, repel = T)
# 5.1 Keep current active.id. as subcell.id and reset active.id as subcell.id
YS_Itg.sfj$subcell.id <- Idents(YS_Itg.sfj); Idents(YS_Itg.sfj) <- YS_Itg.sfj$cluster.id

### 6. Visualize by celltype.id
Idents(YS_Itg.sfj) <- YS_Itg.sfj$celltype.id
# 6.1 CA
jpeg(filename = paste0(outdir, "YS_Itg_CA.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c('CA')), cols.highlight = "#FF66CC", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_CA.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c('CA')), cols.highlight = "#FF66CC", pt.size = 0.1, label = F)
dev.off()

# 6.2 MY
jpeg(filename = paste0(outdir, "YS_Itg_MY.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MY")), cols.highlight = "#009900", pt.size = 0.1,label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_MY.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MY")), cols.highlight = "#009900", pt.size = 0.1,label = F)
dev.off()

# 6.3 MA
jpeg(filename = paste0(outdir, "YS_Itg_MA.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MA")), cols.highlight = "#999900", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "MA.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("MA")), cols.highlight = "#999900", pt.size = 0.1, label = F)
dev.off()

# 6.4 BC
jpeg(filename = paste0(outdir, "YS_Itg_BC.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("BC")), cols.highlight = "#66CC99", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_BC.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("BC")), cols.highlight = "#66CC99", pt.size = 0.1, label = F)
dev.off()

# 6.5 NK/T
jpeg(filename = paste0(outdir, "YS_Itg_NK_T.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("NK/T")), cols.highlight = "#00CCFF", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_NK_T.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("NK/T")), cols.highlight = "#00CCFF", pt.size = 0.1, label = F)
dev.off()

# 6.6 EC
jpeg(filename = paste0(outdir, "YS_Itg_EC.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EC")), cols.highlight = "#CC9933", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_EC.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EC")), cols.highlight = "#CC9933", pt.size = 0.1, label = F)
dev.off()

# 6.7 FB
jpeg(filename = paste0(outdir, "YS_Itg_FB.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("FB")), cols.highlight = "#FF6600", pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_FB.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("FB")), cols.highlight = "#FF6600", pt.size = 0.1, label = F)
dev.off()

# 6.8 EP
jpeg(filename = paste0(outdir, "YS_Itg_EP.jpeg"), width =4000, height=2500, res =600)
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EP")), cols.highlight = "#6699FF",  pt.size = 0.1, label = F)
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_EP.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", cells.highlight = WhichCells(YS_Itg.sfj, idents = c("EP")), cols.highlight = "#6699FF",  pt.size = 0.1, label = F)
dev.off()

### 7. Final DimPlot
jpeg(filename = paste0(outdir, "YS_Itg_DimPlot.jpeg"), width=5000, height=3500, res = 1200)
DimPlot(YS_Itg.sfj, reduction = "umap", label = T, label.size = 2.5) +
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_DimPlot.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", label = T, label.size = 0.5)
dev.off()

### 8. DimPlot Split by tissue.id
jpeg(filename = paste0(outdir, "YS_Itg_DimPlot_by_tissue.jpeg"), width=8000, height=10000, res = 1200)
DimPlot(YS_Itg.sfj, reduction = "umap", label = F, label.size = 0.5, split.by = 'tissue.id', ncol = 1) +
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
dev.off()
# Save as EPS
postscript(paste0(outdir, "YS_Itg_DimPlot_by_tissue.eps"))
DimPlot(YS_Itg.sfj, reduction = "umap", label = F, label.size = 0.5, split.by = 'tissue.id', ncol = 2)
dev.off()

### 9. Save using celltype.id as active id
saveRDS(YS_Itg.sfj, paste0(outdir,"YS_Itg.sfj"))
