# Set-up codes
rm(list = ls()); gc(); library(Seurat); library(infercnv)
outdir = "./Results/5_infercnv/"
Case_01_Itg.suj <- readRDS('./Results/4_Markers/Case_01_Itg.suj')


###	1. Set input files
### 1.1 Create read count matrix file
Case_01_Itg_infercnv.matrix.file <- paste0(outdir, "inferCNV.matrix")
write.table(Case_01_Itg.suj@assays$RNA@counts, file = Case_01_Itg_infercnv.matrix.file, quote = F, sep = "\t")

# 1.2 Tu EP barcodes and visualize
# Get Tu EP barcodes
Case_01_Tu_EP.bc <- read.table('/media/yschang/V/Case_001/Tu/Results/2_Markers/Case_01_Tu_EP.bc')
Case_01_Tu_EP.bc <- Case_01_Tu_EP.bc$x
Idents(Case_01_Itg.suj, cells = Case_01_Tu_EP.bc) <- 'Tu_EP'
# Plot Tu_EP
jpeg(filename = paste0(outdir, "Tu_EP_Case_01_Itg.jpeg"), width = 4000, height=3000, res = 600)
DimPlot(Case_01_Itg.suj, cells.highlight= Case_01_Tu_EP.bc, reduction = "umap", label = T, label.size = 4)
dev.off()
# Check level
levels(Case_01_Itg.suj)

# 1.3 Create ann.file from Seurat.obj
Case_01_Itg_infercnv.ann.file <- paste0(outdir, "inferCNV.annotation")
infercnv.ann.df <- data.frame(barcode = rownames(Case_01_Itg.suj@meta.data), clusters = Case_01_Itg.suj@active.ident)
infercnv.ann.df$barcode <- gsub("-", ".", infercnv.ann.df$barcode)
write.table(infercnv.ann.df, file = Case_01_Itg_infercnv.ann.file, quote = F, sep = "\t", row.names = F, col.names = F)

### 2. Run inferCNV
# 2.1.Indicate gtf file location
gene.order.file <- "/media/yschang/V/hg19.inferCNV.gtf"

# 2.2 Check DimPlot and indicate ref.group
DimPlot(Case_01_Itg.suj, reduction = "umap", label = T)
levels(Case_01_Itg.suj@active.ident)
ref.group <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")

# 2.3. Create infer cnv objects
infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste0(outdir,'inferCNV.matrix'), gene_order_file = gene.order.file, annotations_file = paste0(outdir,'inferCNV.annotation'), ref_group_names = ref.group, delim = "\t", chr_exclude = c("X", "Y"))

# 2.4. Run infercnv
infercnv.obj <- infercnv::run(infercnv_obj = infercnv.obj, num_threads = 8, cutoff = 0.1, cluster_by_groups = T, denoise = T, HMM = F, out_dir = outdir)

# 2.5 Plot infercnv
plot_cnv(infercnv_obj = infercnv.obj, out_dir = outdir, title = "Case_01_Itg_infercnv_Plot", obs_title = "Lung cancer", ref_title = "Normal cells", output_format = "png", png_res = 1200)

### 3. Reanalyze with adjusted 'k_obs_groups = 7'
# 3.1 Create New.outdir
New.outdir = './Results/5_Infercnv/New/'

# 3.2 Generate new infercnv.obj with new_reference
new_ref.group <-  c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")

infercnv.obj <- CreateInfercnvObject(raw_counts_matrix = paste0(outdir,'inferCNV.matrix'), gene_order_file = gene.order.file, annotations_file = paste0(outdir,'inferCNV.annotation'), ref_group_names = new_ref.group, delim = "\t", chr_exclude = c("X", "Y"))

# 3.3 run infercnv with fine tuned parameters; analysis_mode
infercnv_run = infercnv::run(infercnv_obj = infercnv.obj, num_threads = 8, cutoff=0.1, out_dir= New.outdir, analysis_mode="subclusters", cluster_by_groups=F, plot_steps=T, denoise=T, no_prelim_plot=F, k_obs_groups = 4, HMM=T)

# 3.4 Plot infercnv
plot_cnv(infercnv_obj = infercnv.obj, out_dir = New.outdir, title = "Case_01_Itg_New_InferCNV_Plot", obs_title = "Lung cancer", ref_title = "Normal cells", output_format = "png", png_res = 1200)

### Save