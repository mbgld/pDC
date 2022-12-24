# Symbol to ensembl
### 1. Set-up
library("Seurat"); library("biomaRt"); outdir = './Results/CPDB/'
MY.sfj <- readRDS('./Results/MY.sfj')

# 1. Ann.file from Seurat.obj
annotation.file <- paste0(outdir, "MY_annotation.txt")
annotation <- data.frame(Cell=rownames(MY.sfj@meta.data), cell_type=MY.sfj@active.ident)
annotation$Cell <- gsub("-1", ".1", annotation$Cell)
write.table(annotation, file = annotation.file, quote = F, sep = "\t", row.names = F, col.names = T)

### 2. get count dataset
row.count <-data.frame(MY.sfj@assays$RNA@counts)
count_matrix <- apply(row.count, 2, function(x) (x/sum(x))*10000)

### 3.Change Symbol to Ensembl
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
Symbols <- row.names(count_matrix)
df_gene <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol"), values = Symbols, mart= mart)

### 4. Merging count_matrix with df_gene (row.name vs hgnc_symbol)
# 4.1 Change count_matrix row.name to 1st column
count_matrix <-cbind(hgnc_symbol=rownames(row.count), count_matrix)
rownames(count_matrix) <- 1:nrow(count_matrix)
count_matrix <-as.data.frame(count_matrix)

# 4.2 Merge
df_count_matrix <-merge(df_gene, count_matrix, by= 'hgnc_symbol')

# 4.3 Check and remove duplication
df_count_matrix$hgnc_symbol[duplicated(df_count_matrix$hgnc_symbol)]  
df_count_matrix = df_count_matrix[-which(duplicated(df_count_matrix$hgnc_symbol)), ]

# 4.4 Remove Symbol and set row.names
df_count_matrix <-df_count_matrix[,-1]
df_count_matrix <- data.frame(df_count_matrix, row.names = T)

# 4.5 Export and save
count_matrix.file <- paste0(outdir, "MY_count_matrix.txt")
write.table(df_count_matrix, file = count_matrix.file, quote = F, sep = "\t", col.names = T, row.names = T)
