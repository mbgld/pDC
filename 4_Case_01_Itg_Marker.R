### Set-up
rm(list = ls()); gc(); library(Seurat)
outdir = "./Results/4_Markers/"
Case_01_Itg.suj <- readRDS('./Results/3_Itg/Case_01_Itg.suj')
# Define function
MK_merge <- function(x, y){
  z <- merge(x, y, by="row.names", all.x = F, all.y = F)
  z <- z[rev(order(z$avg_log2FC.x, z$power)), ]
  z <- subset(z, select = c('Row.names', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj'))
}

### 1. FindAllmarkers
# 1.1 Find AllMarker
Case_01_Itg.mk <- FindAllMarkers(Case_01_Itg.suj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
write.table(Case_01_Itg.mk, file = paste0(outdir, "Case_01_Itg.mk"), quote = F, row.names = FALSE, sep = "\t")
# 1.2 Get power and merge
Case_01_Itg.pw <- FindAllMarkers(Case_01_Itg.suj, min.pct = 0.25, only.pos = TRUE , logfc.threshold = 0.25, assay = 'RNA', test.use = "roc")
# 1.3 Merge
Case_01_Itg <- merge(Case_01_Itg.mk, Case_01_Itg.pw, by="row.names", all.x = F, all.y = F)
Case_01_Itg <- subset(Case_01_Itg, select = c('Row.names', 'cluster.x', 'avg_log2FC.x', 'power', 'pct.1.x', 'pct.2.x', 'p_val', 'p_val_adj', 'gene.y'))
Case_01_Itg <- Case_01_Itg[rev(order(Case_01_Itg$cluster.x, Case_01_Itg$avg_log2FC.x, Case_01_Itg$power)), ]
write.table(Case_01_Itg, paste0(outdir, "Case_01_Itg.mp"), quote = F, sep = "\t")

### 2. Findmarkers
for (idx in levels(Case_01_Itg.suj)) {
  a <- FindMarkers(Case_01_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE)
  b <- FindMarkers(Case_01_Itg.suj, ident.1 = idx, logfc.threshold = 0.25, assay = 'RNA', min.pct = 0.25, only.pos = TRUE, test.use = "roc")
  c <- MK_merge(a, b)
  assign(paste0("Case_01_Itg_",idx,".mp"), c)
  write.table(c, paste0(outdir, "Case_01_Itg_", idx, ".mp"), quote = F, row.names = FALSE, sep = "\t")
  rm(a, b, c)
}

### 3. Save End product
saveRDS(Case_01_Itg.suj, paste0(outdir, 'Case_01_Itg.suj'))