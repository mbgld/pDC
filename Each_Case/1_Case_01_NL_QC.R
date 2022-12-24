#############################################
# This code is for Supp table 2.
# 1. to get original cell count no
# 2. to estimate cut-off % mt
# 3. to estimate cut-off % ribo
# 4. to get expected count of doublets (nExp)
# 5. to get pKa value
# 6. And get QC Seurat object; sdj
############################################

# SET-UP Codes
rm(list = ls()); gc(); library(Seurat); library(DoubletFinder)
outdir = './Results/1_QC/'

### Part A. Quality Control
## 1. Get dataset and apply default detection-based filtering
Case_01_NL.rdt <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
Case_01_NL.sbj <- CreateSeuratObject(Case_01_NL.rdt, project = "Case_01_NL", min.cells = 3, min.features = 200); rm(Case_01_NL.rdt)
Case_01_NL.sbj <- RenameCells(Case_01_NL.sbj, add.cell.id = 'Case_01_NL')

## 2. Preprocessing of QC
# 2.1 Get mito percentage
Case_01_NL.sbj[["percent_mt"]] <- PercentageFeatureSet(Case_01_NL.sbj, "^MT-"); summary(Case_01_NL.sbj[["percent_mt"]])
# 2.2 Get ribosomal protein percentage
Case_01_NL.sbj[["percent_ribo"]] <- PercentageFeatureSet(Case_01_NL.sbj, "^RP[SL]"); summary(Case_01_NL.sbj[["percent_ribo"]])
# 2.3 Get hemoglobin genes percentage
Case_01_NL.sbj[["percent_hb"]] <- PercentageFeatureSet(Case_01_NL.sbj, "^HB[^(P)]"); summary(Case_01_NL.sbj[["percent_hb"]])
# 2.4 Get platelet genes percentage
Case_01_NL.sbj[["percent_plat"]] <- PercentageFeatureSet(Case_01_NL.sbj, "PECAM1|PF4"); summary(Case_01_NL.sbj[["percent_plat"]])
# 2.5 Plot before_QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(Case_01_NL.sbj, features = feats, pt.size = 0.1, ncol = 3)
FeatureScatter(Case_01_NL.sbj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
# 2.6 Save important parameter
# plot 1
jpeg(filename = paste0(outdir, "Bef_QC_nFeature_percent_mt.jpeg"), width = 4000, height=3000, res = 600)
VlnPlot(Case_01_NL.sbj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0.1, ncol = 3)
dev.off()
# plot 2
jpeg(filename = paste0(outdir, "Bef_QC_FeatureScatter.jpeg"), width = 4000, height=3000, res = 600)
FeatureScatter(Case_01_NL.sbj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()

## 3. Apply basic QC parameter (percent_mt <= 10 & percent_ribo > 5)
Case_01_NL.sqj <- subset(Case_01_NL.sbj, subset = percent_mt <= 10 & percent_ribo >= 5)
# 3.1 Plot after_QC
VlnPlot(Case_01_NL.sqj, features = feats, pt.size = 0.1, ncol = 3)
FeatureScatter(Case_01_NL.sqj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
# 3.2 Remove unnecessary ^MT and ^RP genes
Case_01_NL.sqj <- Case_01_NL.sqj[!grepl("MALAT1", rownames(Case_01_NL.sqj)), ]
Case_01_NL.sqj <- Case_01_NL.sqj[!grepl("^MT-", rownames(Case_01_NL.sqj)), ]

## 4. Dimension reduction
# 4.1 Normalization
Case_01_NL.snj  = NormalizeData(Case_01_NL.sqj , normalization.method = "LogNormalize", scale.factor = 10000)
# 4.2 Feature selection & Plot
Case_01_NL.snj <- FindVariableFeatures(Case_01_NL.snj , selection.method = "vst", nfeatures = 2000); v.genes <- VariableFeatures(Case_01_NL.snj)
LabelPoints(plot = VariableFeaturePlot(Case_01_NL.snj), points = v.genes[1:15],repel = T)
# 4.3 Scaling
Case_01_NL.ssj <- ScaleData(Case_01_NL.snj, features = v.genes, vars.to.regress = c("nFeature_RNA", "percent_mt"))
# 4.4 Linear Dimensional Reduction by PCA
Case_01_NL.spj  <- RunPCA(object = Case_01_NL.ssj , features = v.genes, npcs = 20)
# 4.5 Determine distance
Case_01_NL.suj  <- RunUMAP(Case_01_NL.spj , dims = 1:10)
# 4.6 Visualize
jpeg(filename = paste0(outdir, "Bef_QC_DimPlot.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_NL.suj, reduction = "umap", label = T)
dev.off()
# 4.7 Count recovered cells
Recoverd_cell_count <- ncol(Case_01_NL.suj)

## 5. DoubletFinder
# 5.1 pK Parameter adjustment
sweep.res <- paramSweep_v3(Case_01_NL.suj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)

# 5.2 visualization by linear curve
jpeg(filename = paste0(outdir, "Optimal_pK.jpeg"), width = 4000, height=3000, res = 600)
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b", col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
dev.off()

# 5.4 define the expected number of doublet cells: nExp = No_recovered_cells X 0.0007589 + 0.0527214
Multiplet_rate <- (Recoverd_cell_count*0.0007589 + 0.0527214)/100
nExp <- round(Recoverd_cell_count * Multiplet_rate)
Case_01_NL.suj <- doubletFinder_v3(Case_01_NL.suj, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10)

# 5.5 Extract the correct column name and visualize
DF.name = colnames(Case_01_NL.suj@meta.data)[grepl("DF.classification", colnames(Case_01_NL.suj@meta.data))]
# Plot 1
jpeg(filename = paste0(outdir, "Localization_doublet.jpeg"), width = 4000, height=3000, res = 600)
cowplot::plot_grid(ncol = 2, DimPlot(Case_01_NL.suj, group.by = "orig.ident") + NoAxes(), DimPlot(Case_01_NL.suj, group.by = DF.name) + NoAxes())
dev.off()
# Plot 2
VlnPlot(Case_01_NL.suj, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

## 6. get final product and save
Case_01_NL.sdj = Case_01_NL.suj[, Case_01_NL.suj@meta.data[, DF.name] == "Singlet"]
cat("Remained cell count after QC is", ncol(Case_01_NL.sdj))
# Plot 1
jpeg(filename = paste0(outdir, "Aft_QC_FeatureScatter.jpeg"), width = 4000, height=3000, res = 600)
FeatureScatter(Case_01_NL.sdj, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
dev.off()
# Plot 2
jpeg(filename = paste0(outdir, "Aft_QC_DimPlot.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_NL.sdj, group.by = "orig.ident")
dev.off()
# Save final product
saveRDS(Case_01_NL.sdj, paste0(outdir, 'Case_01_NL.sdj'))

### End and print ###

cat("Mitochondrial gene percentage is", summary(Case_01_NL.sbj[["percent_mt"]]), "\n")
cat("Ribosomal protein percentage is", summary(Case_01_NL.sbj[["percent_ribo"]]), "\n")

cat("Recovered cell count is", Recoverd_cell_count, "\n")
cat("Remained cell count after QC is", ncol(Case_01_NL.sdj), "\n")

cat("Estimated multiplet rate is", Multiplet_rate, "\n")
cat("Selected pKa is", pK_choose, "\n")
