### 1.Set-up
rm(list = ls()); gc()
library(Seurat);library(ggplot2)
outdir = './Results/'

### 1. Prep
# 1.1 get previous work
MY.suj <- readRDS('./Results/MY.suj')
MY.sfj <- MY.suj

# 1.2 Reorder clusters
MY.sfj@active.ident <- factor(MY.sfj@active.ident, 
               levels = c("Mo-lineage", "FCN1-mono", "Mono-Mc", "CD163/LGMN",
                          "Alv-Mc3", "Alv-Mc2", "Alv-Mc1", "Alv-Mc5",
                          "Alv-Mc4", "Prol-Mc", "Alv-Mc6", "Alv-Mc8", "Alv-Mc7",
                          "mo-DC", "cDC1", "cDC2", "pDC"))

### 2. Select genes
DefaultAssay(MY.sfj) <- 'RNA'
# 2.1 Take out marker genes
# MY_All.mk <- FindAllMarkers(MY.sfj, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25, assay = 'RNA')
# MY_Top3.mk <- MY_All.mk %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
# MY_gene_set_1 <- unique(MY_Top3.mk$gene)

# 2.2 Custom gene sets_2
Mon <- c("FCN1", "CXCL8", "IL1B", "S100A8", "S100A9", "S100A12")
ECM <- c("TGFB1", "TGFBI", "INHBA", "LGMN","CD163")
M2 <- c("MRC1")
M1 <- c("CD14", "CD68", "CD86", "FCGR2A","FCGR3A",  "TLR4")
Alv <- c("MARCO", "FABP4", "MCEMP1")
# TAM <-c("TREM1") # IL10", "CCL2", "CCL5") -> all very bad
Lip <- c("APOC1", "APOE", "LIPA")
Prol <- c("STMN1", "MKI67", "TOP2A")
Cdc <- c("FCER1A", "CLEC10A", "PKIB", "CD1E", "BTLA", "CLEC9A", "CADM1", "XCR1", "LY75",  "CCR7", "LAMP3")
Mdc <- c("MRC1", "CD209", "SIRPA", "ITGAM", "CD1A")
Pdc <- c("CLEC4C", "IL3RA", "PTCRA")

MY_gene_set_2 <- c(Mon, M1, ECM, M2, Alv, Lip, Prol, Cdc, Pdc)

# 1.2 Final plot
jpeg(filename = paste0(outdir, "MY_DotPlot.jpeg"), width = 8000, height =3200, res = 600)
DotPlot(MY.sfj, features = MY_gene_set_2, cols=c('black', 'red'), col.min=-2.5, col.max=2.5, dot.scale = 7.5) &
  theme(axis.text.x = element_text(size=10)) +
  RotatedAxis()
dev.off()

saveRDS(MY.sfj, paste0(outdir,"MY.sfj"))
