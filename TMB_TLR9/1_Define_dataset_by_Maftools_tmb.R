### 0. Set-up
library(maftools)
Output='./Results/'

### 1. Call dataset
# 1.0 Generate MAF file and generate new barcodes
TCGA_LUAD_1=read.maf('./Maf_files.tsv')
df_MAF_01<-tmb(maf = TCGA_LUAD_1)
df_MAF_01$submitter.id<-substr(df_MAF_01$Tumor_Sample_Barcode, 1, 12)
df_MAF_02<-df_MAF_01[!duplicated(df_MAF_01[ , c('submitter.id')]),]

# 1.1 Take-out available folder list
RNASeq_folder_list <- read.csv('../RNASeq/Folder_list.csv')
RNASeq_folder_list <- sort(RNASeq_folder_list$X0)
RNASeq_folder_list <- RNASeq_folder_list[-1:-85]

# 1.2 Select case only with RNASeq
df_MAF_02 <-df_MAF_02[df_MAF_02$submitter.id %in% RNASeq_folder_list, ]

### 2. Take out Q1 TMB and Q3 TMB
# 2.1 Set Q1 and Q3
TMB_q1<-quantile(df_MAF_02$total_perMB, 0.25)
TMB_q3<-quantile(df_MAF_02$total_perMB, 0.75)
# 2.2 Take out dataframe
df_TMB_q1<-df_MAF_02[df_MAF_02$total_perMB<=TMB_q1]
df_TMB_q3<-df_MAF_02[df_MAF_02$total_perMB>=TMB_q3]
# 2.3 Get list of Q1 and Q3
TCGA_LUAD_q1_TMB<-df_TMB_q1$submitter.id
TCGA_LUAD_q3_TMB<-df_TMB_q3$submitter.id
# 2.4 Save
write.table(df_MAF_02, file=paste0(Output, 'TCGA_LUAD_MAF_tmb.tsv'), sep = '\t')
write.table(TCGA_LUAD_q1_TMB, file=paste0(Output, 'TCGA_LUAD_q1_tmb.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE)
write.table(TCGA_LUAD_q3_TMB, file=paste0(Output, 'TCGA_LUAD_q3_tmb.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE)

### 3. Take out Low, intermediate and high TMB dafaframe
# 3.1 Take out dataframe
df_Low_TMB <- df_MAF_02[df_MAF_02$total_perMB<=5]
df_Int_TMB <- df_MAF_02[(df_MAF_02$total_perMB>5)&(df_MAF_02$total_perMB<20)  ]
df_High_TMB <- df_MAF_02[df_MAF_02$total_perMB>=20]
# 3.2 Get list of Low, intermediate and high TMB
TCGA_LUAD_Low_TMB<-df_Low_TMB$submitter.id
TCGA_LUAD_Int_TMB<-df_Int_TMB$submitter.id
TCGA_LUAD_High_TMB<- df_High_TMB$submitter.id
# 3.3 Save
write.table(TCGA_LUAD_Low_TMB, file=paste0(Output, 'TCGA_LUAD_Low_TMB.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE)
write.table(TCGA_LUAD_Int_TMB, file=paste0(Output, 'TCGA_LUAD_Int_TMB.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE)
write.table(TCGA_LUAD_High_TMB, file=paste0(Output, 'TCGA_LUAD_High_TMB.tsv'), sep = '\t', col.names = FALSE, row.names = FALSE)

