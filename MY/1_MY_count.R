# 3. MY (NL vs Tu)
MY.outdir = './Results/1_MY/'
df_MY = subset(df_03, Cell_name=='MY')
df_MY$Annotation <- gsub("CD163/LGMN", "CD163-LGMN", df_MY$Annotation)
MY.ann <- unique(df_MY$Annotation)

df_NL_MY = subset(df_MY, Origin =="NL"); NL_MY.cnt = sum(df_NL_MY$Count)
df_Tu_MY = subset(df_MY, Origin =="Tu"); Tu_MY.cnt = sum(df_Tu_MY$Count)

# 3.1 Basic count

### Never smokers MY in NL and Tu vs Cur smokers' MY in NL and Tu.
# Ratio of NL_MY in Never smoker
Never_smoker_NL_MY.cnt = sum(subset(df_NL_MY, Smoking_Hx=="Never")$Count)
Never_smoker_NL.cnt = sum(subset(df_NL, Smoking_Hx=="Never")$Count)
# Ratio of Tu_MY in Never smoker
Never_smoker_Tu_MY.cnt = sum(subset(df_Tu_MY, Smoking_Hx=="Never")$Count)
Never_smoker_Tu.cnt = sum(subset(df_Tu, Smoking_Hx=="Never")$Count)
# Final value
cat("The ratio of MY in NL tissue vs in Tu tissue in Never smoker is :",
(Never_smoker_NL_MY.cnt / Never_smoker_NL.cnt *100)/(Never_smoker_Tu_MY.cnt / Never_smoker_Tu.cnt *100))

# Ratio of NL_MY in Cur smoker
Cur_smoker_NL_MY.cnt = sum(subset(df_NL_MY, Smoking_Hx=="Cur")$Count)
Cur_smoker_NL.cnt = sum(subset(df_NL, Smoking_Hx=="Cur")$Count)
# Ratio of Tu_MY in Cur smoker
Cur_smoker_Tu_MY.cnt = sum(subset(df_Tu_MY, Smoking_Hx=="Cur")$Count)
Cur_smoker_Tu.cnt = sum(subset(df_Tu, Smoking_Hx=="Cur")$Count)
# Final value
cat("The ratio of MY in NL tissue vs in Tu tissue in Cur smoker is : ", (Cur_smoker_NL_MY.cnt / Cur_smoker_NL.cnt *100)/(Cur_smoker_Tu_MY.cnt / Cur_smoker_Tu.cnt *100))

##
cat("The ratio of NL_MY in Never vs Current smoker is : ",(Never_smoker_NL_MY.cnt / Never_smoker_NL.cnt *100)/(Cur_smoker_NL_MY.cnt / Cur_smoker_NL.cnt *100))
cat("The ratio of Tu_MY in Never vs Current smoker is : ",(Never_smoker_Tu_MY.cnt / Never_smoker_Tu.cnt *100)/(Cur_smoker_Tu_MY.cnt / Cur_smoker_Tu.cnt *100))


# 3.1 Count MY cells of each case in NL or Tu
for (idx in Cases){
  assign(paste0(idx, "_NL_MY.cnt"), sum((df_NL_MY %>% filter(Case==idx))$Count))
  assign(paste0(idx, "_Tu_MY.cnt"), sum((df_Tu_MY %>% filter(Case==idx))$Count))}


# 3.2 Add Pct_by_case_origin by each case in NL or Tu
df_06 = data.frame()

for (idx in Cases){
  for (idy in Tissues) {
    a = subset(df_MY, Case == idx & Origin == idy)
    b = sum(a$Count)
    c = a %>% mutate(Pct_by_case_origin = Count/b*100)
    df_06 = rbind(df_06, c)
    }
}
# mondify dataset
df_06$Group[df_06$Origin == "NL" & df_06$Smoking_Hx == "Never"] <- "Never_NL"
df_06$Group[df_06$Origin == "Tu" & df_06$Smoking_Hx == "Never"] <- "Never_Tu"
df_06$Group[df_06$Origin == "NL" & df_06$Smoking_Hx == "Cur"] <- "Cur_NL"
df_06$Group[df_06$Origin == "Tu" & df_06$Smoking_Hx == "Cur"] <- "Cur_Tu"

# 3.3 Comparison of Pct_by_case_origin according to Smoking status and Tissue
for(i in MY.ann){
  print(paste0('--------', i, '-----------'))
  a = subset(df_06, Annotation==i)
  ggboxplot(a, x="Group", y="Pct_by_case_origin", color = "black", fill="Group",
    xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
    palette = c("#0000FF", "#FF0000", "#FFFF00", "#808080"), add="jitter", legend="none") +
    stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_A/Experimental/JPEG]", i, ".jpeg"), width = 5, height = 3)
}

# 4.1 Proportion test
# Both tissue "EP and CA" removed: comparison of MY in NL vs Tu
prop.test(c(NL_MY.cnt, Tu_MY.cnt), c(NL.cnt, Tu.cnt))

# Part A
# 1. In both tissue: comparison between NL vs Tu
sink('./Results/1_MY/Part_A/Both/[Both] MYs (NL vs Tu).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('--------', i, '-----------'))
  a = subset(df_04, Annotation==i)
  # Distribution test
  print(shapiro.test(a$Pct_by_case_origin))
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin ~ Origin, data = a, paired = TRUE, exact = NULL))
  
  # EPS
  postscript(paste0(MY.outdir, "Part_A/Both/EPS [Both]", i, " (NL vs Tu).eps" ))
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_A/Both/JPEG [Both]", i, " (NL vs Tu).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# 2. Comparison between NL vs Tu in cur smokers.
sink('./Results/1_MY/Part_A/Cur/[Cur] MYs (NL vs Tu).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('-------', i, '----------'))
  a = subset(df_04, Annotation==i & Smoking_Hx=="Cur")
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin ~ Origin, data = a, paired = TRUE, exact = NULL))
  
  # EPS 
  postscript(paste0(MY.outdir, "Part_A/Cur/EPS [Cur]", i, " (NL vs Tu).eps" ))
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_A/Cur/JPEG [Cur]", i, " (NL vs Tu).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# 3. Comparison between NL vs Tu in never smokers.
sink('./Results/1_MY/Part_A/Never/[Never] MYs (NL vs Tu).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('-------', i, '----------'))
  a = subset(df_04, Annotation==i & Smoking_Hx=="Never")
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin ~ Origin, data = a, paired = TRUE, exact = NULL))
  
  # EPS 
  postscript(paste0(MY.outdir, "Part_A/Never/EPS [Never]", i, " (NL vs Tu).eps" ))
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Origin", y="Pct_by_case_origin", color = "black", fill="Origin", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_A/Never/JPEG [Never]", i, " (NL vs Tu).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# Part B
# 1. Both tissue: comparison between cur vs Never smoker
sink('./Results/1_MY/Part_B/Both/[Both] MYs (Cur vs Never).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('-------', i, '----------'))
  a = subset(df_04, Annotation==i)
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin ~ Smoking_Hx, data = a, paired = FALSE, exact = NULL))
  
  # EPS 
  postscript(paste0(MY.outdir, "Part_B/Both/EPS [Both]", i, " (Cur vs Never).eps" ))
  ggboxplot(a, x="Smoking_Hx", y="Pct_by_case_origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Smoking_Hx", y="Pct_by_case_origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_B/Both/JPEG [Both]", i, " (Cur vs Never).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# 2. NL tissue: Comparison of proportion between Cur vs Never
sink('./Results/1_MY/Part_B/NL/[NL] MYs (Cur vs Never).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('--------', i, '-----------'))
  a = subset(df_04, Annotation==i & Origin=="NL")
  # Distribution test
  print(shapiro.test(a$Pct_by_case_origin))
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin~Smoking_Hx, data = a, paired = FALSE, exact = NULL))
  
  # EPS
  postscript(paste0(MY.outdir, "Part_B/NL/EPS [NL]", i, " (Cur vs Never).eps" ))
  ggboxplot(a, x="Smoking_Hx", y="Pct_by_case_origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Smoking_Hx", y="Pct_by_case_origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_B/NL/JPEG [NL]", i, " (Cur vs Never).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# 3.Tu tissue: comparison of proportion between Cur vs Never
sink('./Results/1_MY/Part_B/Tu/[Tu] MYs (Cur vs Never).txt', split=TRUE)
for(i in MY.ann){
  print(paste0('--------', i, '-----------'))
  a = subset(df_04, Annotation==i & Origin=="Tu")
  # Distribution test
  #print(shapiro.test(a$Pct_by_case_origin))
  # Test of mean by Wilcox test
  print(wilcox.test(Pct_by_case_origin~Smoking_Hx, data = a, paired = FALSE, exact = NULL))
  
  # EPS 
  postscript(paste0(MY.outdir, "Part_B/Tu/EPS [Tu]", i, " (Cur vs Never).eps" ))
  ggboxplot(a, x="Smoking_Hx", y="Origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  dev.off()
  
  # JPEG
  ggboxplot(a, x="Smoking_Hx", y="Pct_by_case_origin", color = "black", fill="Smoking_Hx", palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") + stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
  ggsave(paste0(MY.outdir, "Part_B/Tu/JPEG [Tu]", i, "(Cur vs Never).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()

# Part C. Calculate TN.ratio of specific cluster
sink('./Results/1_MY/Part_C/[Ratio] MYs (Tu vs NL).txt', split=TRUE)
df_05 = data.frame()
for (idx in Cases){
  for (idy in MY.ann){
    a = subset(df_MY, Case == idx & Origin == "Tu" & Annotation == idy)
    b = subset(df_MY, Case == idx & Origin == "NL" & Annotation == idy)
    c = a$Count/b$Count*100
    d = a %>% mutate(TN.ratio = c)
    df_05 = rbind(df_05, d)
  }
}

for(i in MY.ann){
  print(paste0('--------', i, '-----------'))
  a = subset(df_05, Annotation==i)
  # Distribution test
  #print(shapiro.test(a$Pct_by_case_origin))
  # Test of mean by Wilcox test
  print(wilcox.test(TN.ratio~Smoking_Hx, data = df_05, paired = FALSE, exact = NULL))
  
  # JPEG
  ggboxplot(a, x="Smoking_Hx", y="TN.ratio", color = "black", fill="Smoking_Hx",
    palette = c("#0000FF", "#FF0000"), add="jitter", legend="none") +
    stat_compare_means(aes(label = sprintf("p = %5.4f", as.numeric(..p.format..))))
    ggsave(paste0(MY.outdir, "Part_C/[Ratio]", i, "(Cur vs Never).jpeg"), width = 2.5, height = 3)
  rm(a, i)
}
sink()
