# 1. Set-up
library(ggplot2); library(dplyr); library(ggpubr); library(stringr)
PATH = "/media/yschang/V/Specific_cluster/"
outdir = './Results/'
Cases = c("Case_01", "Case_02", "Case_03", "Case_04", "Case_05", 
          "Case_07", "Case_08", "Case_09", "Case_10", "Case_12")
Tissues = c("NL", "Tu")
Smoking_Hx = c("Never", "Cur")
Never_smoker = c("Case_01", "Case_02", "Case_03", "Case_04", "Case_05")
Cur_smoker = c("Case_07", "Case_08", "Case_09", "Case_10", "Case_12")

# 1.1 Call datasets and merge
df_MY = read.csv(paste0(PATH, '/1_MY/Results/MY_count_by_original_ident.txt'), header = FALSE)
df_MA = read.csv(paste0(PATH, '/2_MA/Results/MA_count_by_original_ident.txt'), header = FALSE)
df_BC = read.csv(paste0(PATH, '/3_BC/Results/BC_count_by_original_ident.txt'), header = FALSE)
df_NKT = read.csv(paste0(PATH, '/4_NKT/Results/NK_T_count_by_original_ident.txt'), header = FALSE)
df_FB = read.csv(paste0(PATH, '/5_FB/Results/FB_count_by_original_ident.txt'), header = FALSE)
df_EC = read.csv(paste0(PATH, '/6_EC/Results/EC_count_by_original_ident.txt'), header = FALSE)
df_EP = read.csv(paste0(PATH, '/7_EP/Results/EP_count_by_original_ident.txt'), header = FALSE)
df_CA = read.csv(paste0(PATH, '/8_CA/Results/CA_count_by_original_ident.txt'), header = FALSE)
df_01 = rbind(df_MY, df_MA, df_BC, df_NKT, df_FB, df_EC, df_EP, df_CA)
rm(df_MY, df_MA, df_BC, df_NKT, df_FB, df_EC, df_EP, df_CA)

# 1.2 Organize dataset, Add colnames, change attribute of data
df_02 <- data.frame(do.call(rbind, strsplit(df_01$V1, split = ':')))
df_02$Case <- substr(df_02$X2, 2, 8); df_02$Origin <- substr(df_02$X2, 10, 11); df_02 <- df_02[ ,-c(2)]
names(df_02) <- c("Cell_name", "Annotation", "Count",  "Case", "Origin")
df_02 = df_02 %>% mutate(Smoking_Hx = (ifelse(Case == "Case_01"|Case == "Case_02"|Case == "Case_03"|Case == "Case_04"|Case == "Case_05", "Never", "Cur")))
# str_trim(df_02$Cell_name); str_trim(df_02$Annotation);str_trim(df_02$Count);str_trim(df_02$Case);str_trim(df_02$Origin)
df_02$Count <- as.numeric(df_02$Count)

# 1.3 Major cell counts
Total_cell.cnt = sum(df_02$Count)
Cell_name <- unique(df_02$Cell_name)

# Major cell cluster count
sink(paste0(outdir, 'Each_cluster_count.txt'), split = TRUE)
for (idx in Cell_name) {
  a = subset(df_02, Cell_name==idx)
  assign(paste0(idx, ".cnt"), sum(a$Count))
  cat(paste0(idx, ".cnt"), ':', sum(a$Count), '\n')
  rm(a)
}
sink()

# Major cell cluster count by tissue origin
sink(paste0(outdir, 'Each_cluster_count_by_origin.txt'), split = TRUE)
for (idx in Cell_name) {
  for (idy in Tissues) {
    a = subset(df_02, Cell_name == idx)
    b = subset(a, Origin == idy)
    assign(paste0(idy, "_", idx, ".cnt"), sum(b$Count))
    cat(paste0(idy, "_", idx, ".cnt"), ':', sum(b$Count), '\n')
    rm(a, b)
  }
}
sink()

# Major cell cluster count by smoking history
sink(paste0(outdir, 'Each_cluster_count_by_smoking.txt'), split = TRUE)
for (idx in Cell_name) {
  for (idy in Smoking_Hx) {
    a = subset(df_02, Cell_name == idx)
    b = subset(a, Smoking_Hx == idy)
    assign(paste0(idy, "_", idx, ".cnt"), sum(b$Count))
    cat(paste0(idy, "_", idx, ".cnt"), ':', sum(b$Count), '\n')
    rm(a, b)
  }
}
sink()

# Major cell cluster count by smoking history and tissue origin
sink(paste0(outdir, 'Each_cluster_count_by_smoking_tissue_origin.txt'), split = TRUE)
for (idx in Cell_name) {
  for (idy in Smoking_Hx) {
    for (idz in Tissues) {
    a = subset(df_02, Cell_name == idx)
    b = subset(a, Smoking_Hx == idy)
    c = subset(b, Origin == idz)
    assign(paste0(idy, "_", idz, "_", idx, ".cnt"), sum(c$Count))
    cat(paste0(idy, "_", idz, "_", idx, ".cnt"), ':', sum(c$Count), '\n')
    rm(a, b)
    }
  }
}
sink()

# 1.4 Remove EP and CA cells from further counting
df_03 <- df_02[!(df_02$Cell_name=="CA"|df_02$Cell_name=="EP"), ] # CA and EP cells were removed

# 1.5 Add Smoking_Tissue Group
df_03$Group[df_03$Origin == "NL" & df_03$Smoking_Hx == "Never"] <- "Never_NL"
df_03$Group[df_03$Origin == "Tu" & df_03$Smoking_Hx == "Never"] <- "Never_Tu"
df_03$Group[df_03$Origin == "NL" & df_03$Smoking_Hx == "Cur"] <- "Cur_NL"
df_03$Group[df_03$Origin == "Tu" & df_03$Smoking_Hx == "Cur"] <- "Cur_Tu"

# 2. Basic counting job
# 2.1 Total cell count(NL vs. Tu/ Never vs. Cur)
df_NL = subset(df_03, Origin=="NL"); NL.cnt = sum(df_NL$Count)
df_Tu = subset(df_03, Origin=='Tu'); Tu.cnt = sum(df_Tu$Count)
df_Never = subset(df_03, Smoking_Hx=="Never"); Never.cnt = sum(df_Never$Count)
df_Cur = subset(df_03, Smoking_Hx=='Cur'); Cur.cnt = sum(df_Cur$Count)

# 2.2 Cell count of each case
sink(paste0(outdir, 'Total_cells_each_case.txt'), split = TRUE)
for (idx in Cases) {
  a = subset(df_03, Case==idx)
  assign(paste0(idx, ".cnt"), sum(a$Count))
  cat(paste0(idx, ".cnt"), ':', sum(a$Count), '\n')
  rm(a)
}
sink()

# 2.3 NL cell count of each case
sink(paste0(outdir, 'NL_cells_each_case.txt'), split = TRUE)
for (idx in Cases){
  a = subset(df_03, Case==idx & Origin =="NL")
  assign(paste0(idx, "_NL", ".cnt"), sum(a$Count))
  cat(paste0(idx, "_NL.cnt"), ':', sum(a$Count), '\n')
}
sink()

# 2.4 Tu cell count of each case
sink(paste0(outdir, 'Tu_cells_each_case.txt'), split = TRUE)
for (idx in Cases){
  a = subset(df_03, Case==idx & Origin =="Tu")
  assign(paste0(idx, "_Tu", ".cnt"), sum(a$Count))
  cat(paste0(idx, "_Tu.cnt"), ':', sum(a$Count), '\n')
}
sink()

# 3.1 Calculate percentile of each cluster according to tissue origin
df_04A=data.frame()

for (idx in Tissues){
  for (idy in unique(df_03$Cell_name)){
    a=subset(df_03, Cell_name==idy & Origin ==idx)
    b=sum(a$Count)
    c=b/(get(paste0(idx, ".cnt")))*100
    d=data.frame(Origin=idx, Cell_name=idy, Percentile=c)
    df_04A=rbind(df_04A, d)
  }
}

# 3.2 Prop.test
prop.test(c(NL_MY.cnt, Tu_MY.cnt), c(NL.cnt, Tu.cnt))
prop.test(c(NL_MA.cnt, Tu_MA.cnt), c(NL.cnt, Tu.cnt))
prop.test(c(NL_BC.cnt, Tu_BC.cnt), c(NL.cnt, Tu.cnt))
prop.test(c(NL_NK_T.cnt, Tu_NK_T.cnt), c(NL.cnt, Tu.cnt))
prop.test(c(NL_FB.cnt, Tu_FB.cnt), c(NL.cnt, Tu.cnt))
prop.test(c(NL_EC.cnt, Tu_EC.cnt), c(NL.cnt, Tu.cnt))

# 3.3 Make stacked bar graph of Major cluster by tissue orgin
library(ggplot2)
ColorPalate = c("#C77CFF",  "#7CAE00",  "#0CB702", "#00B8E7",
                "#F8766D", "#CD9600") #"#00A9FF",

ggplot(df_04A, aes(Origin, Percentile, fill = Cell_name))+
  geom_bar(position = 'stack', stat='identity')+
  scale_fill_manual(values = ColorPalate)+
  theme(axis.title=element_text(size=1), legend.title = element_text(size=1), #change legend title font size
  legend.text = element_text(size=1))+
  theme_light()
ggsave(paste0(outdir, "Fraction_Major_Cell_by_tissue.tiff"), width = 6, height = 6, units = "cm", dpi=600)

# 4.1 Calculate percentile of each cluster according to smoking_Hx
df_04B=data.frame()

for (idx in Smoking_Hx){
  for (idy in unique(df_03$Cell_name)){
    a=subset(df_03, Smoking_Hx == idx)    # Cur.cnt form df_03
    b=sum(a$Count)
    c=subset(a, Cell_name==idy)
    d=sum(c$Count)
    e=d/b*100
    f=data.frame(Smoking_Hx= idx, Cell_name=idy, Percentile=e)
    df_04B=rbind(df_04B, f)
  }
}

# 4.2 Prop.test
prop.test(c(Never_MY.cnt, Cur_MY.cnt), c(Never.cnt, Cur.cnt))
prop.test(c(Never_MA.cnt, Cur_MA.cnt), c(Never.cnt, Cur.cnt))
prop.test(c(Never_BC.cnt, Cur_BC.cnt), c(Never.cnt, Cur.cnt))
prop.test(c(Never_NK_T.cnt, Cur_NK_T.cnt), c(Never.cnt, Cur.cnt))
prop.test(c(Never_FB.cnt, Cur_FB.cnt), c(Never.cnt, Cur.cnt))
prop.test(c(Never_EC.cnt, Cur_EC.cnt), c(Never.cnt, Cur.cnt))

# 4.3 Make stacked bar graph of Major cluster
ggplot(df_04B, aes(Smoking_Hx, Percentile, fill = Cell_name))+
  geom_bar(position = 'stack', stat='identity')+
  scale_fill_manual(values = ColorPalate)+
  scale_x_discrete(limit = c("Never", "Cur"))+
  theme(axis.title=element_text(size=1), legend.title = element_text(size=1), #change legend title font size
        legend.text = element_text(size=1))+
  theme_light()
ggsave(paste0(outdir, "Fraction_Major_Cell_by_Smoking_Hx.tiff"), width = 6, height = 6, units = "cm", dpi=600)

# 5. Calculate percentile of each cluster according to smoking_Hx and Tissue_origin
df_04C = data.frame()

for (idx in Cases) {
  for (idy in Tissues) {
    for (idz in unique(df_03$Cell_name)){
      a= subset(df_03, Case==idx & Origin ==idy)
      b= sum(a$Count)
      c= subset(a, Cell_name== idz)
      d= sum(c$Count)
      e= d/b*100
      f=data.frame(Case= idx, Origin =idy , Cell_name=idz, Percentile=e)
      df_04C=rbind(df_04C, f)
    }
  }
}

df_04C$Group[df_04C$Origin == "NL" & df_04C$Case %in% Never_smoker] <- "Never_NL"
df_04C$Group[df_04C$Origin == "Tu" & df_04C$Case %in% Never_smoker] <- "Never_Tu"
df_04C$Group[df_04C$Origin == "NL" & df_04C$Case %in% Cur_smoker] <- "Cur_NL"
df_04C$Group[df_04C$Origin == "Tu" & df_04C$Case %in% Cur_smoker] <- "Cur_Tu"

for(i in unique(df_03$Cell_name)){
  print(paste0('--------', i, '-----------'))
  a = subset(df_04C, Cell_name==i)
  ggboxplot(a, x="Group", y="Percentile", color = "black", fill="Group",
            xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
            palette = c("cyan", "darkblue", "gray70", "gray20"), add="jitter", legend="none") +
    stat_compare_means(method = "kruskal.test")
  ggsave(paste0(outdir, "Part_A_JPEG]", i, ".jpeg"), width = 5, height = 4)
}

# 3.3.1 Add p-value on above MY plot
a = subset(df_04C, Cell_name=="MY")
ggboxplot(a, x="Group", y="Percentile", color = "black", fill="Group",
          xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
          palette = c("#0000FF", "#FF0000", "#FFFF00", "#808080"), add="jitter", legend="none") +
  stat_compare_means(method = "kruskal.test", label.y = 90)
ggsave(paste0(outdir, "Part_A_JPEG]", "MY", ".jpeg"), width = 5, height = 4)

library(FSA)
dunnTest(Percentile ~ Group, data=a, method="bonferroni")

# 4. Calculate percentile of major cluster according to each case
df_05=data.frame()

for(idx in Cases){
  for (idy in df_03$Cell_name) {
    a=subset(df_03, Case==idx & Cell_name==idy)
    b=sum(a$Count)
    c=b/get(paste0(idx, ".cnt"))*100
    d=data.frame(Case=idx, Cell_name=idy, Percentile=c)
    df_05=rbind(df_05, d)
  }
}
df_05= unique(df_05)
df_05 = df_05%>% mutate(Smoking_Hx = ifelse(Case %in% Never_smoker, "never", "cur"))

# Plot
library(ggplot2)
ggplot(df_05, aes(Case, Percentile, fill = Cell_name))+
  geom_bar(position = 'stack', stat='identity')+
  scale_fill_manual(values = ColorPalate)+
  # facet_grid(.~Smoking_Hx)+
  # facet_wrap(.~Smoking_Hx)+
  theme_light()+
  coord_flip()
ggsave(paste0(outdir, "Pctl_by_case.tiff"), width = 8, height = 12, units = "cm", dpi=900)

# save.image(paste0(outdir, "Subcluster_comp_13Sep.RData"))
