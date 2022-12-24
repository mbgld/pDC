# set up
rm(list = ls()); gc();
library(Seurat); library(dplyr)
outdir = "./Results/6_Split/Sim_2_9/"

# load
YS_Itg.sfj <- readRDS("./Results/4_Assign/Sim_2_9/YS_Itg.sfj")

# FB: #FF6600
FB.sbj <- subset(YS_Itg.sfj, idents = 'FB')
FB.sbj <- SetIdent(FB.sbj, value = "subcluster.id" )
DimPlot(FB.sbj, label=T, repel = TRUE)
saveRDS(FB.sbj, paste0(outdir,"FB.sbj"))

# EC: #CC9933
EC.sbj <- subset(YS_Itg.sfj, idents = 'EC')
EC.sbj <- SetIdent(EC.sbj, value = "subcluster.id" )
DimPlot(EC.sbj, label=T, repel = TRUE)
saveRDS(EC.sbj, paste0(outdir,"EC.sbj"))

# MY: #009900
MY.sbj <- subset(YS_Itg.sfj, idents = 'MY')
MY.sbj <- SetIdent(MY.sbj, value = "subcluster.id" )
DimPlot(MY.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(MY.sbj, paste0(outdir,"MY.sbj"))

# MA: #999900
MA.sbj <- subset(YS_Itg.sfj, idents = 'MA')
MA.sbj <- SetIdent(MA.sbj, value = "subcluster.id" )
DimPlot(MA.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(MA.sbj, paste0(outdir,"MA.sbj"))

# NKT: #00CCFF
NKT.sbj <- subset(YS_Itg.sfj, idents = 'NK/T')
NKT.sbj <- SetIdent(NKT.sbj, value = "subcluster.id" )
DimPlot(NKT.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(NKT.sbj, paste0(outdir,"NK_T.sbj"))

# BC: ##66CC99
BC.sbj <- subset(YS_Itg.sfj, idents = 'BC')
BC.sbj <- SetIdent(BC.sbj, value = "subcluster.id" )
DimPlot(BC.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(BC.sbj, paste0(outdir,"BC.sbj"))

# EP: #6699FF
EP.sbj <- subset(YS_Itg.sfj, idents = 'EP')
EP.sbj <- SetIdent(EP.sbj, value = "subcluster.id" )
DimPlot(EP.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(EP.sbj, paste0(outdir,"EP.sbj"))

# CA: #FF66CC
CA.sbj <- subset(YS_Itg.sfj, idents = 'CA')
CA.sbj <- SetIdent(CA.sbj, value = "subcluster.id" )
DimPlot(CA.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(CA.sbj, paste0(outdir,"CA.sbj"))

# EPIC: EP + CA
EPIC.sbj <- subset(YS_Itg.sfj, idents = c('CA', 'EP'))
EPIC.sbj <- SetIdent(EPIC.sbj, value = "subcluster.id" )
DimPlot(EPIC.sbj, label=T, label.size = 2, repel = TRUE)
saveRDS(EPIC.sbj, paste0(outdir,"EPIC.sbj"))
