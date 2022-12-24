### GSEA; fold change가 큰 순으로, 비슷한 발현양상을 보이는 유전자들의 SNR 방법으로 순위를 매겨 정렬하고, 관심있는 molecular profile상관관계를 이용해 ES 및 P-value를 계산한 뒤, 관심있는 GO 정보를 제공하는 방법

### 0. Get library
library(clusterProfiler); library(enrichplot); library(dplyr)
library(ggplot2) # we use ggplot2 to add x axis labels (ex: ridgeplot)
organism = 'org.Hs.eg.db'
library(organism, character.only = TRUE)
outdir = "./Results/"

### 1. Prep input dataset: Get data from deseq2: 자료의 형태는 반드시 names()를 사용하여 named number의 형태로 지정되어야 한다
df_1 = read.csv('./Results/pDC_Cur_vs_Never_Tu.tsv', header = TRUE, sep='\t', row.names = 1)
gene_list = df_1$avg_log2FC # log2FoldChange의 값이 그림을 그리는데 핵심자료이므로 이것만 불러와 vector로 지정한다.
names(gene_list) = rownames(df_1)# names 함수를 이용하여 위에서 불러온 vector에다가 유전자명의 열을 이용해 이름을 가져다 붙인다..
gene_list = na.omit(gene_list) # 결측치를 제거하고서
gene_list = sort(gene_list, decreasing = TRUE) # 내림차순으로 정렬한다. 안하면 clusterProfiler에서 에러발생한다. 

### 2.INPUT 에서 사용되어 있는 유전자의 keytype을 확인
keytypes(org.Hs.eg.db) # org.Hs.eg.db에서 제공하는 keytype의 list를 확인한다.
columns(org.Hs.eg.db) # org.Hs.eg.db에서 제공하는 column의 이름을 확인한다. 
id_exam = c('EEF1G')
select(org.Hs.eg.db, keys = id_exam, columns= c('ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENENAME'), keytype = 'SYMBOL')

### 3. Perform gse. Watchout parameters
gse <- gseGO(geneList = gene_list, ont = 'ALL', keyType = "SYMBOL", minGSSize=3, maxGSSize=800, pvalueCutoff = 0.05, verbose=TRUE, OrgDb=org.Hs.eg.db, pAdjustMethod='BH')

### 4. Plots
# 4.1 DOTPLOT을 그려본다. 
require(DOSE)
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_dotplot.tiff"), width = 600, height = 800, type="cairo") 
dotplot(gse, showCategory=10, split='.sign') + facet_grid(.~.sign)
dev.off()

# 4.2 Enrichment MAP을 그린다.
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_enrichment_map.tiff"), width = 1000, height = 800, units = "px", type="cairo") 
emapplot(gse, showCategory = 10)
dev.off()

# 4.3 Category Netplot을 그리고 싶다
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_cat_netplot.tiff"), width = 900, height = 600, units = "px", type="cairo")
cnetplot(gse, categorySize='pvalue', foldChange = gene_list, showCategory = 3)
dev.off()

# 4.4 RIDGEPLOT으로 표현하고 싶다.
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_ridgeplot.tiff"), width = 600, height = 1000, units = "px", type="cairo") 
ridgeplot(gse, showCategory = 10) + labs(x='enrichment distribution')
dev.off()

# 4.5  GSEA plot
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_GSEAplot.tiff"), width = 1400, height = 1800, units = "px", type="cairo") 
gseaplot(gse, by = 'all', title=gse$Description[1], geneSetID = 1)
dev.off()

# 4.6 PUBMED trend of enriched terms
terms <- gse$Description[1:3]
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Tu_pubmed_terms.tiff"), width = 1400, height = 800, units = "px", type="cairo")
pmcplot(terms, 2018:2022, proportion = FALSE)
dev.off()

# Part B (entire pDC in Cur vs Never)
### 1. Prep input dataset
df_2 = read.csv('./Results/pDC_Cur_vs_Never_Entire.tsv', header = TRUE, sep='\t', row.names = 1)
gene_list = df_2$avg_log2FC 
names(gene_list) = rownames(df_2)
gene_list = na.omit(gene_list) 
gene_list = sort(gene_list, decreasing = TRUE)

### 2. Perform gse. Watch out parameters
gse <- gseGO(geneList = gene_list, ont = 'ALL', keyType = "SYMBOL", minGSSize=3, maxGSSize=800, pvalueCutoff = 0.05, verbose=TRUE, OrgDb=org.Hs.eg.db, pAdjustMethod='BH')

### 4. Plots
# 4.1 DOTPLOT
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_dotplot.tiff"), width = 600, height = 800, type="cairo") 
dotplot(gse, showCategory=7, split='.sign') + facet_grid(.~.sign)
dev.off()

# 4.2 Enrichment MAP
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_enrichment_map.tiff"), width = 1000, height = 800, units = "px", type="cairo") 
emapplot(gse, showCategory = 10)
dev.off()

# 4.3 Category Netplot
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_cat_netplot.tiff"), width = 600, height = 300, units = "px", type="cairo")
cnetplot(gse, categorySize='pvalue', foldChange = gene_list, showCategory = 3)
dev.off()

# 4.4 RIDGEPLOT
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_ridgeplot.tiff"), width = 600, height = 1000, units = "px", type="cairo") 
ridgeplot(gse, showCategory = 10) + labs(x='enrichment distribution')
dev.off()

# 4.5  GSEA plot
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_GSEAplot.tiff"), width = 1400, height = 1800, units = "px", type="cairo") 
gseaplot(gse, by = 'all', title=gse$Description[1], geneSetID = 1)
dev.off()

# 4.6 PUBMED trend of enriched terms
terms <- gse$Description[1:3]
tiff(filename = paste0(outdir, "pDC_Smoker_vs_Never_in_Entire_pubmed_terms.tiff"), width = 1400, height = 800, units = "px", type="cairo")
pmcplot(terms, 2010:2018, proportion = FALSE)
dev.off()


# Part C (entire pDC in Cur vs Never)
### 1. Prep input dataset
df_3 = read.csv('./Results/pDC_Tu_vs_NL_in_cur.tsv', header = TRUE, sep='\t', row.names = 1)
gene_list = df_3$avg_log2FC 
names(gene_list) = rownames(df_3)
gene_list = na.omit(gene_list) 
gene_list = sort(gene_list, decreasing = TRUE)

### 2. Perform gse. Watch out parameters
gse <- gseGO(geneList = gene_list, ont = 'ALL', keyType = "SYMBOL", minGSSize=3, maxGSSize=800, pvalueCutoff = 0.05, verbose=TRUE, OrgDb=org.Hs.eg.db, pAdjustMethod='BH')

### 4. Plots
# 4.1 DOTPLOT
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_dotplot.tiff"), width = 600, height = 800, type="cairo") 
dotplot(gse, showCategory=7, split='.sign') + facet_grid(.~.sign)
dev.off()

# 4.2 Enrichment MAP
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_in_Entire_enrichment_map.tiff"), width = 1000, height = 800, units = "px", type="cairo") 
emapplot(gse, showCategory = 10)
dev.off()

# 4.3 Category Netplot
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_Catnetplot.tiff"), width = 600, height = 300, units = "px", type="cairo")
cnetplot(gse, categorySize='pvalue', foldChange = gene_list, showCategory = 3)
dev.off()

# 4.4 RIDGEPLOT
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_ridgeplot.tiff"), width = 600, height = 1000, units = "px", type="cairo") 
ridgeplot(gse, showCategory = 10) + labs(x='enrichment distribution')
dev.off()

# 4.5  GSEA plot
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_GSEAplot.tiff"), width = 1400, height = 1800, units = "px", type="cairo") 
gseaplot(gse, by = 'all', title=gse$Description[1], geneSetID = 1)
dev.off()

# 4.6 PUBMED trend of enriched terms
terms <- gse$Description[1:3]
tiff(filename = paste0(outdir, "pDC_NL_vs_Tu_in_Smoker_pubmed_terms.tiff"), width = 1400, height = 800, units = "px", type="cairo")
pmcplot(terms, 2010:2018, proportion = FALSE)
dev.off()

