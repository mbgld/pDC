library(FSA)
my_comparison_A <- list(c("Never_NL", "Never_Tu"))
my_comparison_B <- list(c("Cur_NL", "Cur_Tu"))

# 1. BC
a = subset(df_04C, Cell_name=='BC')
ggboxplot(a, x="Group", y="Percentile", color = "black", fill="Group",
          xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
          palette = c("#0000FF", "#FF0000", "#FFFF00", "#808080"), add="jitter", legend="none") +
  stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 10) +
  stat_compare_means(comparisons = my_comparison_B, method = "wilcox.test", label.y = 15) +
  stat_compare_means(method = "kruskal.test", label.y = 15)
  ggsave(paste0(outdir, "BC", ".jpeg"), width = 5, height = 4)

# 1.1 Post-Hoc
dunnTest(Percentile ~ Group, data=a, method="bonferroni")

# 2. DC
a = subset(df_04C, Cell_name=='DC')
ggboxplot(a, x="Group", y="Percentile", color = "black", fill="Group",
          xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
          palette = c("#0000FF", "#FF0000", "#FFFF00", "#808080"), add="jitter", legend="none") +
  stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 10) +
  stat_compare_means(comparisons = my_comparison_B, method = "wilcox.test", paired = TRUE, label.y = 10) +
  stat_compare_means(method = "kruskal.test", label.y = 12)
ggsave(paste0(outdir, "DC", ".jpeg"), width = 5, height = 4)

# 2.1 Post-Hoc
dunnTest(Percentile ~ Group, data=a, method="bonferroni")

# 3. NK_T
a = subset(df_04C, Cell_name=='NK_T')
ggboxplot(a, x="Group", y="Percentile", color = "black", fill="Group",
          xlab = "Group", ylab = 'Percentile', font.label = list(size = 1, face = "plain"), 
          palette = c("#0000FF", "#FF0000", "#FFFF00", "#808080"), add="jitter", legend="none") +
  stat_compare_means(comparisons = my_comparison_A, method = "wilcox.test", paired = TRUE, label.y = 80) +
  stat_compare_means(comparisons = my_comparison_B, method = "wilcox.test", paired = FALSE, label.y = 75) +
  stat_compare_means(method = "kruskal.test", label.y=90)
ggsave(paste0(outdir, "NK_T", ".jpeg"), width = 5, height = 4)

# 3.1 Post-Hoc
dunnTest(Percentile ~ Group, data=a, method="bonferroni")