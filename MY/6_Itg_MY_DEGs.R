### 0. Call assigned MY.sfj dataset
library(Seurat); library(dplyr); outdir = './Results/'
YS_Itg.sfj <- readRDS("/media/yschang/V1/Combi/Results/4_Assign/Sim_2_9/YS_Itg.sfj")
MY.sfj <- readRDS('./Results/MY.sfj')
DimPlot(MY.sfj, reduction = "umap",  label = TRUE, label.size = 4, label.color = "black", repel = T)
levels(MY.sfj)
levels(YS_Itg.sfj)

# 1.0  Get cell barcodes
# 1.1 Cells identified by tissue.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "tissue.id")
YS_Tu.bc <- WhichCells(YS_Itg.sfj, idents = 'Tu')
YS_NL.bc <- WhichCells(YS_Itg.sfj, idents = 'NL')

# 1.2 Cells identified from smoking.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "smoking.id")
YS_never_smoker.bc <- WhichCells(YS_Itg.sfj, idents = 'never')
YS_cur_smoker.bc <- WhichCells(YS_Itg.sfj, idents = 'cur')

# 1.3 Take out barcode of MY
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "celltype.id")
MY.bc <- WhichCells(YS_Itg.sfj, idents = 'MY')
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "subcluster.id")

# 1.4 Take out barcode of pDC
pDC.bc <- WhichCells(MY.sfj, idents="pDC")
DimPlot(YS_Itg.sfj, cells.highlight = pDC.bc, split.by = 'tissue.id', sizes.highlight = 0.001)
DimPlot(YS_Itg.sfj, cells.highlight = pDC.bc, split.by = 'smoking.id', sizes.highlight = 0.001)

### Part A
# 2.1 Current smoker's vs Never smoker's pDC entire 
pDC_cur.bc <- intersect(pDC.bc, YS_cur_smoker.bc)
pDC_never.bc <- intersect(pDC.bc, YS_never_smoker.bc)
a <- FindMarkers(MY.sfj, ident.1 = pDC_cur.bc, ident.2 = pDC_never.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "pDC_Cur_vs_Never_Entire.tsv"), quote = F, row.names = TRUE, sep = "\t")

# 2.2 Current smoker's pDC vs. Never smoker's pDC in Tu
pDC_cur_Tu.bc <- intersect(pDC.bc, intersect(YS_cur_smoker.bc, YS_Tu.bc))
pDC_never_Tu.bc <- intersect(pDC.bc, intersect(YS_never_smoker.bc, YS_Tu.bc))
a <- FindMarkers(MY.sfj, ident.1 = pDC_cur_Tu.bc, ident.2 = pDC_never_Tu.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "pDC_Cur_vs_Never_Tu.tsv"), quote = F, row.names = TRUE, sep = "\t")

# 2.3 Current smoker's pDC in NL vs Tu
pDC_cur_Tu.bc <- intersect(pDC.bc, intersect(YS_cur_smoker.bc, YS_Tu.bc))
pDC_cur_NL.bc <- intersect(pDC.bc, intersect(YS_cur_smoker.bc, YS_NL.bc))
a <- FindMarkers(MY.sfj, ident.1 = pDC_cur_Tu.bc, ident.2 = pDC_cur_NL.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "pDC_Tu_vs_NL_in_cur.tsv"), quote = F, row.names = TRUE, sep = "\t")

# get ModuleScore
pDC_cur_Tu.sfj <- subset(MY.sfj, cells = pDC_cur_Tu.bc)
a <- read.csv('./pDC_GO.txt', header = TRUE, sep='\t')
AddModuleScore(pDC_cur_Tu.sfj, features = )

### Part B
# 3.1 Current smoker's vs Never smoker's MY entire 
MY_cur.bc <- intersect(MY.bc, YS_cur_smoker.bc)
MY_never.bc <- intersect(MY.bc, YS_never_smoker.bc)
a <- FindMarkers(MY.sfj, ident.1 = MY_cur.bc, ident.2 = MY_never.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "MY_Cur_vs_Never_Entire.tsv"), quote = F, row.names = TRUE, sep = "\t")

# 3.2 Current smoker's MY vs. Never smoker's MY in Tu
MY_cur_Tu.bc <- intersect(MY.bc, intersect(YS_cur_smoker.bc, YS_Tu.bc))
MY_never_Tu.bc <- intersect(MY.bc, intersect(YS_never_smoker.bc, YS_Tu.bc))
a <- FindMarkers(MY.sfj, ident.1 = MY_cur_Tu.bc, ident.2 = MY_never_Tu.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "MY_Cur_vs_Never_Tu.tsv"), quote = F, row.names = TRUE, sep = "\t")

# 3.3 Current smoker's MY in NL vs Tu
MY_cur_Tu.bc <- intersect(MY.bc, intersect(YS_cur_smoker.bc, YS_Tu.bc))
MY_cur_NL.bc <- intersect(MY.bc, intersect(YS_cur_smoker.bc, YS_NL.bc))
a <- FindMarkers(MY.sfj, ident.1 = MY_cur_Tu.bc, ident.2 = MY_cur_NL.bc)
b <- a[a$p_val_adj < 0.05, ]
c <- b[order((b$avg_log2FC), decreasing = T), ]
write.table(c, file = paste0(outdir, "MY_Tu_vs_NL_in_cur.tsv"), quote = F, row.names = TRUE, sep = "\t")
