setwd("/boot3/znh/8月病毒TB细胞比对/BCR/分析")
library(scRepertoire)
library(immunarch)
library(Seurat)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggplot2)

# Load B cell Seurat object and combined BCR data
Seucount <- readRDS("/boot3/znh/8月病毒TB细胞比对/BCR/分析/B_Plasma_dob.rds") 
combined <- readRDS("/boot3/znh/8月病毒TB细胞比对/BCR/分析/combined-BCR.rds")

# Combine BCR clonotype information into Seurat object
seurat <- combineExpression(combined, Seucount,
                           cloneCall = "gene",
                           # group.by = "Sample",
                           proportion = TRUE,
                           cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1))

# Color palettes for clone sizes
col <- c("#7DA26D", "#743D8B", "#C079A7")
col_deep <-  c("#5A9E3D", "#7A28A0", "#C05088")
col <- c("#7DA26D","#743D8B","#C079A7")

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

# t-SNE plot colored by clone size
p222 <- DimPlot(seurat, group.by = "cloneSize", label = F, reduction = "tsne") + #NoLegend() +
  scale_color_manual(values = col_deep) + 
  theme(plot.title = element_blank())
ggsave("sc_clontype.pdf", plot = p222, width = 8, height = 5)

# Reclassify clinical groups
seurat$clin <- seurat$illness
seurat$clin[which(seurat$clin == "Nomal")] <- "Normal"
seurat$clin[which(seurat$clin == "Severe")] <- "Patient"
seurat$clin[which(seurat$clin == "Mild")] <- "Patient"

# Clonal occupancy by clinical status (proportion)
p222 <- clonalOccupy(seurat, x.axis = "clin", proportion = T) +
  scale_fill_manual(values = colorblind_vector(4)[-1])  #
ggsave("sc克隆分布-Clin-percent.pdf", plot = p222, width = 4, height = 5)

# Per clinical group, clonal occupancy by cell type (recluster)
for (i in unique(seurat$clin)){
  seurat_T <- subset(seurat, clin == i)
  p222 <- clonalOccupy(seurat_T, x.axis = "recluster", proportion = T) +
    scale_fill_manual(values = col) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(paste0(i, "sc克隆分布.pdf"), plot = p222, width = 8, height = 6)
} 

# ---------------------------
# Add module scores for gene signatures
genes <- list()
genes$activation <- c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", "CD47", "EOMES",
                       "IFNG", "CD8A", "CD8B", "FASLG", "LAMP1",
                       "HLA-DRA", "TNFRSF4", "ICOS", "TNFRSF9", "TNFRSF18","CD38")

# Add "cyto" branch to genes list
genes$cyto <- c("GZMA", "GZMB", "GZMM", "NKG7", "GNLY", "PRF1")

# Add "EX" branch to genes list
genes$EX <- c("HAVCR2", "LAG3", "TIGIT", "CTLA4", "PDCD1", "LAYN")

seurat <- AddModuleScore(object = seurat, 
                         features = genes, 
                         name = names(genes))

dim(seurat@meta.data)
colnames(seurat@meta.data)[25:27] <- names(genes)

library(ggsignif)
library(ggpubr)

# Activation signature vs clone type
df <- data.frame(
  cloneType = seurat$cloneType,
  signature_score = seurat$activation) 
colors <- colorblind_vector(5)
df <- na.omit(df)

p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(x = "Clone Type", y = "Signature Score", fill = "Clone Type", title = "Activation signature Score by Clone Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(
  c("Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)"),
  c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)"),
  c("Small (1e-04 < X <= 0.001)", "Large (0.01 < X <= 0.1)")
))
ggsave("Activation.pdf", plot = p, width = 5, height = 5.5)

# EX signature vs clone type
df <- data.frame(
  cloneType = seurat$cloneType,
  signature_score = seurat$EX) 
colors <- colorblind_vector(5)
dim(df)
df <- na.omit(df)
p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(x = "Clone Type", y = "Signature Score", fill = "Clone Type", title = "EX signature Score by Clone Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(
  c("Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)"),
  c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)"),
  c("Small (1e-04 < X <= 0.001)", "Large (0.01 < X <= 0.1)")
))
ggsave("EX.pdf", plot = p, width = 5, height = 5.5)

# Cyto signature vs clone type
df <- data.frame(
  cloneType = seurat$cloneType,
  signature_score = seurat$cyto) 
colors <- colorblind_vector(5)
df <- na.omit(df)
p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(x = "Clone Type", y = "Signature Score", fill = "Clone Type", title = "Cyto signature Score by Clone Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = list(
  c("Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)"),
  c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)"),
  c("Small (1e-04 < X <= 0.001)", "Large (0.01 < X <= 0.1)")
))
ggsave("Cyto.pdf", plot = p, width = 5, height = 5.5)

# ---------------------------
# BCR diversity analysis using immunarch style
contig_annotations_list <- readRDS("/boot3/znh/8月病毒TB细胞比对/BCR/分析/bcr_contig_annotations_list.rds")
combined <- combineBCR(contig_annotations_list,
                       ID = c("CXJ3", "JH4", "LTY11", "HPJ1", "GLM6", "JN12", "ZJL2", "PYL10"),
                       samples = c("Normal", "Normal", "Normal", "Mild", "Mild","Mild","Severe","Severe"))
data_tcr <- addVariable(combined, variable.name = "lnc",   # corrected parameter name
                       variables = c("Normal", "Normal", "Normal", "illness", "illness","illness", "illness", "illness")) 

# Diversity analysis (bootstrap)
p10 <- clonalDiversity(data_tcr, 
                       cloneCall = "strict", 
                       n.boots = 1000, 
                       x.axis = 'lnc', # grouping variable (here 'lnc' corresponds to clinical status)
                       group = 'ID') + # ID is each sample
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
setwd("/boot3/znh/8月病毒TB细胞比对/BCR/新分析/新结果/多样性")
ggsave("多样性分析-N-P.pdf", plot = p10, width = 8, height = 6)