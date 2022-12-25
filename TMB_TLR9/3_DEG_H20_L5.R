### 1. Set-ups
library(DESeq2); library(dplyr); library(ggplot2)
output='./Results/Unmatched_H20_L5/'
gene_list <- read.csv('../Results/RNA/gene_list.tsv', sep='\t')
Qvalue = 0.05; log2FC = 1

### 2. Call dataset
# 2.1 Generate cnt data
df_RNASeq_tmb = read.csv("./Results/RNA_seq_with_Low_High_TMB.csv", row.names = 1)
df_cnt= as.matrix(df_RNASeq_tmb[1:60660, 2:317])
# rm(df_RNASeq_tmb)

# 2.2 Generate col data
df_col= read.csv("./Results/Colnames_RNA_seq_with_Low_High_TMB.csv", row.names = NULL, sep='\t')

### 3. Generate DEG matrix
DDS_1 <- DESeqDataSetFromMatrix(countData = df_cnt, colData = df_col, design = ~ condition); DDS_1

### 4. Normalization and call normalized table
DDS_1.norm <- DESeq(DDS_1)
DDS_1.norm.table <- counts(DDS_1.norm, normalized=T)

### 5. Get results
res.DDS_1.norm <- results(DDS_1.norm, pAdjustMethod="BH")
mcols(res.DDS_1.norm, use.names=TRUE)

### 6. Sample QC and check the expression table
DDS_1.norm.table.filter = DDS_1.norm.table[!is.na(res.DDS_1.norm[,2]),]
dim(DDS_1.norm.table)
dim(DDS_1.norm.table.filter)

### 7. Gene count distribution
par(mfcol=c(1,1),mar=c(2,2,2,1), oma=c(2,2,.5,.5), lwd=.5)
plot(density(log(DDS_1.norm.table.filter[,1:305])), col="red", ylim=c(0, 0.25), main="Gene count distribution")
lines(density(log(DDS_1.norm.table.filter[,306:316])),col="blue")	
mtext("Expression (log10)",side=1,line=2, cex=1.5)
mtext("Density",side=2,line=2.5, cex=1.5)
legend("topright",legend=c("TMB_Low","TMB_High"),col=c("red","blue"),pch=15,cex=1.8)

### 8. boxplot
par(mfcol=c(1,1), mar=c(2,2,1,1), oma=c(3,3,1,1), font=1)
boxplot(log10(DDS_1.norm.table.filter+0.001), cex.axis=1.2, main="", ylab="", xlab="", range=2, outline=F, las=2, outcol="#004C9990", outpch=16, outcex=1, lwd=1.3, axes=T, border="#004C99")
mtext("Expression (log10)",side=2.5,line=2, cex=1.5)

### 9. PCA analysis
pca <- prcomp(DDS_1.norm.table.filter, center=TRUE, scale=TRUE) 
par(mfcol=c(1,1), mar=c(4,4.5,.5,.5), oma=c(.5,.5,0,0))
plot(pca$rotation[,1], pca$rotation[,2], main="", xlab="PC1",ylab="PC2", pch=20, cex=5, cex.axis=1.5, cex.lab=1.5, lwd=.5, col="yellow", ylim=c(-0.5, 0.5), xlim=c(-0.2, 0.2))
text(pca$rotation[,1], pca$rotation[,2], labels=row.names(pca$rotation), cex=.7, col="blue", font=2)

### 10. DEG
res.DDS_1.norm.exp <- cbind(res.DDS_1.norm, DDS_1.norm.table)
res.DDS_1.norm.exp <- res.DDS_1.norm.exp[order(res.DDS_1.norm.exp$pvalue),]
res.DDS_1.norm.exp.filt <- res.DDS_1.norm.exp[!is.na(res.DDS_1.norm.exp[,6]),]
res.DDS_1.norm.exp.filt$ensembl <- sapply(strsplit(rownames(res.DDS_1.norm.exp.filt), split="\\+" ), "[", 1)
res.DDS_1.norm.exp.filt <- merge(res.DDS_1.norm.exp.filt, gene_list, by.x="ensembl", by.y="gene_id")
res.DDS_1.norm.exp.filt <- res.DDS_1.norm.exp.filt[, c(1:7, 324:325)]
res.DDS_1.norm.exp.filt <- res.DDS_1.norm.exp.filt[order(res.DDS_1.norm.exp.filt$log2FoldChange, decreasing = TRUE), ]

## 10.1 Save entire DEG list
write.table(res.DDS_1.norm.exp.filt, file=paste0(output, "Unmatched_DEG_table_entire_gene.tsv"), row.names=T, col.names=T, sep="\t", quote=F)

## 10.2 Generate gene table
up.regulated.genes <- res.DDS_1.norm.exp.filt[res.DDS_1.norm.exp.filt$padj <= Qvalue & res.DDS_1.norm.exp.filt$log2FoldChange > log2FC ,]
down.regulated.genes <- res.DDS_1.norm.exp.filt[res.DDS_1.norm.exp.filt$padj <= Qvalue & res.DDS_1.norm.exp.filt$log2FoldChange < -log2FC ,]
dim(up.regulated.genes)
dim(down.regulated.genes)

# Save lists of up-regulated and down-regulated
write.table(up.regulated.genes, file=paste0(output, "Unmatched_Up_regulated_genes.tsv"), row.names=T, col.names=T,sep="\t",quote=F)
write.table(down.regulated.genes, file=paste0(output, "Unmatched_Down_regulated_genes.tsv"), row.names=T, col.names=T,sep="\t",quote=F)

## 10.3 Select genes with significant meaning
Sig.DEGs <- res.DDS_1.norm.exp.filt[res.DDS_1.norm.exp.filt$padj <= Qvalue & (res.DDS_1.norm.exp.filt$log2FoldChange > log2FC|res.DDS_1.norm.exp.filt$log2FoldChange < -log2FC),]
write.table(Sig.DEGs, file=paste0(output,"Unmatched_Sig_DEGs.tsv"), row.names=T, col.names=T, sep="\t", quote=F)


## 11 Volcano plot (1) with annotation
# In the preset Q value <= 0.05 and fold change > 2x or fold change < 0.5
tiff(paste0(output, "Volcano_DEG_1.tiff"), width = 8, height = 5, units = 'in', res = 300)
par(mfcol=c(1,1), mar=c(3,3,2,2), oma=c(2,2,0, 0), lwd=1)
plot(res.DDS_1.norm.exp.filt$log2FoldChange, -log10(res.DDS_1.norm.exp.filt$padj), col="#004C9950", pch=16, cex=1, xlab="", ylab="")
abline(h=1.3, lty=3, lwd=1.5, col="red")
abline(v=1, lty=3, lwd=1.5, col="red")
abline(v=-1, lty=3, lwd=1.5, col="red") 
text(res.DDS_1.norm.exp.filt$log2FoldChange[1:20], -log10(res.DDS_1.norm.exp.filt$padj[1:20]), res.DDS_1.norm.exp.filt$gene_name[1:20], col="black", cex=0.5, adj=c(0,0))
mtext("Fold change (log2)", side=1, line=2.5, cex=1.2)
mtext("Significance (-log10)", side=2, line=2.5, cex=1.2)
dev.off()

## 12 Volcano plot (2) with coloring
tiff(paste0(output, "Volcano_DEG_2.tiff"), width = 8, height = 5, units = 'in', res = 300)
cut_lfc <- 1; cut_pvalue <- 0.01 # Set-up cut off value
par(mar=c(5,5,5,5), cex=0.75, cex.main=2.0, cex.axis=1.0, cex.lab=1.4)
topT <- as.data.frame(res.DDS_1.norm.exp.filt)
# Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)
text(res.DDS_1.norm.exp.filt$log2FoldChange[1:20], -log10(res.DDS_1.norm.exp.filt$padj[1:20]), res.DDS_1.norm.exp.filt$gene_name[1:20], col="black", cex=0.5)
dev.off()

## 13. Additional plot
# 13.1 Call geneset from reactome and extract from DEG table
TLR9_cascade <- read.table('./REACTOME_TLR9.csv', sep = ',', header = TRUE)
TLR9_cascade <- TLR9_cascade[19,]
df_TLR9_set <- res.DDS_1.norm.exp.filt[res.DDS_1.norm.exp.filt$gene_name %in% TLR9_cascade, ]

# 13.2 Draw_A
tiff(paste0(output, "Volcano_TLR9_A.tiff"), width = 8, height = 5, units = 'in', res = 300)
par(mfcol=c(1,1), mar=c(3,3,2,2), oma=c(2,2,0, 0), lwd=1)
plot(res.DDS_1.norm.exp.filt$log2FoldChange, -log10(res.DDS_1.norm.exp.filt$padj), col="#004C9950", pch=16, cex=1, xlab="", ylab="")
abline(h=1.3, lty=3, lwd=1.5, col="red")
abline(v=1, lty=3, lwd=1.5, col="red")
abline(v=-1, lty=3, lwd=1.5, col="red") 
text(df_TLR9_set$log2FoldChange[1:10], -log10(df_TLR9_set$padj[1:10]), df_TLR9_set$gene_name[1:10], col="black", cex=0.5, adj=c(0,0))
mtext("Fold change (log2)", side=1, line=2.5, cex=1.2)
mtext("Significance (-log10)", side=2, line=2.5, cex=1.2)
dev.off()

# 13.2 Draw_B
tiff(paste0(output, "Volcano_TLR9_B.tiff"), width = 8, height = 8, units = 'in', res = 300)
cut_lfc <- 1; cut_pvalue <- 0.01 # Set-up cut off value
par(mar=c(5,5,5,5), cex=0.75, cex.main=2.0, cex.axis=1.0, cex.lab=1.4)
topT <- as.data.frame(res.DDS_1.norm.exp.filt)
# Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", col='grey', cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), points(log2FoldChange, -log10(padj), pch=20, col='red', cex=1.5))
with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col='blue', cex=1.5))
## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=3, lwd=1.0)
abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)
text(df_TLR9_set$log2FoldChange[1:5], -log10(df_TLR9_set$padj[1:5]), df_TLR9_set$gene_name[1:5], col="black", cex=0.75, adj=c(0,0))
dev.off()
