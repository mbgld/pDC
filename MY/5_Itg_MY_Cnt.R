### 1.Set-up
rm(list = ls()); gc()
library(Seurat); library(dplyr); outdir = './Results/'

# 1.1 get previous work
YS_Itg.sfj <- readRDS("/media/yschang/V1/Combi/Results/4_Assign/Sim_2_9/YS_Itg.sfj")
MY.sfj <- readRDS('./Results/MY.sfj')

# 1.2 Cells identified by tissue.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "tissue.id")
YS_Tu.bc <- WhichCells(YS_Itg.sfj, idents = 'Tu')
YS_NL.bc <- WhichCells(YS_Itg.sfj, idents = 'NL')

# 1.3 Cells identified from smoking.id
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "smoking.id")
YS_never_smoker.bc <- WhichCells(YS_Itg.sfj, idents = 'never')
YS_cur_smoker.bc <- WhichCells(YS_Itg.sfj, idents = 'cur')

# 1.4 MY cell counts
## MY from never_smoker's NL 
length(intersect(colnames(MY.sfj), intersect(YS_NL.bc, YS_never_smoker.bc)))
## MY from cur_smoker's NL
length(intersect(colnames(MY.sfj), intersect(YS_NL.bc, YS_cur_smoker.bc)))
## MY from never_smoker's Tu
length(intersect(colnames(MY.sfj), intersect(YS_Tu.bc, YS_never_smoker.bc)))
## MY from cur_smoker's Tu
length(intersect(colnames(MY.sfj), intersect(YS_Tu.bc, YS_cur_smoker.bc)))

### 2 Assign active.id as subcell.id
# 2.1 Set id
MY.sfj <- SetIdent(MY.sfj, value = 'subcell.id')

# 2.2 Plots
DimPlot(MY.sfj, reduction = "umap", label = T, pt.size = 1.0)

### 3. Count MY's subcells
sink(paste0(outdir, 'Subcell_count_of_MY.txt'), split = TRUE)
for (i in levels(MY.sfj)){
  cat(i, "count:", length(WhichCells(MY.sfj, idents = i)), "\n")
  rm(i)
}
sink()

### 4. count by tissue.id and smoking.id
sink(paste0(outdir,'MY_subcell_by_tissue_and_smoking.txt'), split = TRUE)
for(a in unique(MY.sfj$smoking.id)){
  for (i in levels(MY.sfj)){
    cat(i, "in", a, "_smoker's NL tissue:", length(intersect(get(paste0("YS_", a, "_smoker.bc")), intersect(YS_NL.bc, WhichCells(MY.sfj, idents = i)))), "\n")
    cat(i, "in", a, "_smoker's Tu Tissue:", length(intersect(get(paste0("YS_", a, "_smoker.bc")), intersect(YS_Tu.bc, WhichCells(MY.sfj, idents = i)))), "\n")
    rm(i)
  }
  rm(a)
}
sink()

### 5. count by orig.ident
YS_Itg.sfj <- SetIdent(YS_Itg.sfj, value = "orig.ident")

sink(paste0(outdir, 'MY_count_by_original_ident.txt'), split = TRUE)
for(a in levels(YS_Itg.sfj)){
  for (b in levels(MY.sfj)){
    cat("MY:", a, ":", b, ":", length(intersect(WhichCells(YS_Itg.sfj, idents = a), WhichCells(MY.sfj, idents = b))), "\n")
    rm(b)}
  rm(a)
}
sink()

