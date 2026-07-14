####
# TCR analysis and visualization
library(scRepertoire)
library(immunarch)
library(Seurat)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggplot2)
setwd("/home/xudahua/TCR")

### T cell data
Seucount <- readRDS("T.rds") 
combined <- readRDS("combined-tcr.rds")

# Combine TCR clonotype information into Seurat object
seurat <- combineExpression(combined, Seucount,
                           cloneCall = "gene",
                           proportion = TRUE,
                           cloneSize = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1))

# Color palette for clone sizes
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
col <- c("#7DA26D", "#743D8B", "#C079A7")
col_deep <-  c("#5A9E3D", "#7A28A0", "#C05088")
col <- c("#7DA26D","#743D8B","#C079A7")

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
ggsave("sc克隆分布-Clin-percent.pdf", plot = p222, width = 4, height = 5) # file name remains

# Per clinical group, clonal occupancy by cell type (recluster)
for (i in unique(seurat$clin)){
  seurat_T <- subset(seurat, clin == i)
  p222 <- clonalOccupy(seurat_T, x.axis = "recluster", proportion = T) +
    scale_fill_manual(values = col) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ggsave(paste0(i, "sc克隆分布.pdf"), plot = p222, width = 8, height = 6) # file name with Chinese
} 

# ---------------------------
# Add module scores for gene signatures
genes <- list()
genes$activation <- c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", "CD47", "EOMES",
                       "IFNG", "CD8A", "CD8B", "FASLG", "LAMP1",
                       "HLA-DRA", "TNFRSF4", "ICOS", "TNFRSF9", "TNFRSF18","CD38")

genes$cyto <- c("GZMA", "GZMB", "GZMM", "NKG7", "GNLY", "PRF1")

genes$EX <- c("HAVCR2", "LAG3", "TIGIT", "CTLA4", "PDCD1", "LAYN")

seurat <- AddModuleScore(object = seurat, 
                         features = genes, 
                         name = names(genes))

dim(seurat@meta.data)
colnames(seurat@meta.data)[25:27] <- names(genes)

# ---------------------------
# Violin/boxplot for Activation signature vs clone size
library(ggsignif)
library(ggpubr)

df <- data.frame(
  cloneType = seurat$cloneSize,
  signature_score = seurat$activation) 
colors <- colorblind_vector(5)[-1]
colors <- colors[-3]
df <- na.omit(df)

p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  labs(x = "Clone Type", y = "Signature Score", fill = "Clone Type", title = "Activation signature Score by Clone Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Add pairwise Wilcoxon significance annotations with adjusted y positions
y_max <- max(df$signature_score, na.rm = TRUE)
y_min <- min(df$signature_score, na.rm = TRUE)
step <- (y_max - y_min) * 0.08
y_pos1 <- y_max + step
y_pos2 <- y_max + step * 2.5
y_pos3 <- y_max + step * 4

p <- p + geom_signif(
    comparisons = list(c("Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos1,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
) + geom_signif(
    comparisons = list(c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos2,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
) + geom_signif(
    comparisons = list(c("Small (1e-04 < X <= 0.001)", "Large (0.01 < X <= 0.1)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos3,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
)

print(p)
ggsave("Activation.pdf", plot = p, width = 5, height = 5.5)

# EX signature vs clone size
df <- data.frame(
  cloneType = seurat$cloneSize,
  signature_score = seurat$EX) 
dim(df)
df <- na.omit(df)
p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = col) +
  theme_minimal() +
  labs(x = "Clone Type", y = "Signature Score", fill = "Clone Type", title = "EX signature Score by Clone Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

y_max <- max(df$signature_score, na.rm = TRUE)
y_min <- min(df$signature_score, na.rm = TRUE)
step <- (y_max - y_min) * 0.08
y_pos1 <- y_max + step
y_pos2 <- y_max + step * 2.5
y_pos3 <- y_max + step * 4

p <- p + geom_signif(
    comparisons = list(c("Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos1,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
) + geom_signif(
    comparisons = list(c("Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos2,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
) + geom_signif(
    comparisons = list(c("Small (1e-04 < X <= 0.001)", "Large (0.01 < X <= 0.1)")),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1),
    y_position = y_pos3,
    tip_length = 0.01,
    textsize = 3.5,
    vjust = 0.2
)
ggsave("EX.pdf", plot = p, width = 5, height = 5.5)

# Cyto signature vs clone size
df <- data.frame(
  cloneType = seurat$cloneSize,
  signature_score = seurat$cyto) 
df <- na.omit(df)
p <- ggplot(df, aes(x = cloneType, y = signature_score, fill = cloneType)) +
  geom_violin() +
  geom_boxplot(width = 0.1, position = position_dodge(0.75), fill = "white") +
  scale_fill_manual(values = col) +
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
# TCR diversity analysis using immunarch style
list <- readRDS("contig_annotations_list.rds")
contig_annotations_list <- list
combined <- combineTCR(contig_annotations_list,
                       ID = c("CXJ3", "JH4", "LTY11", "HPJ1", "GLM6", "JN12", "ZJL2", "PYL10"),
                       samples = c("Normal", "Normal", "Normal", "Mild", "Mild","Mild","Severe","Severe"))
data_tcr <- addVariable(combined, variable.name = "lnc", 
                       variables = c("Normal", "Normal", "Normal", "illness", "illness","illness", "illness", "illness")) 

# Diversity analysis (bootstrap)
p10 <- clonalDiversity(data_tcr, 
                       cloneCall = "strict", 
                       n.boots = 1000, 
                       x.axis = 'lnc', 
                       group = 'ID') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("多样性分析-N-P.pdf", plot = p10, width = 8, height = 6) # file name with Chinese

# ---------------------------
# Feature plots for T cell marker genes
genes <- c("CD3D", "CD3E", "CD3G", "CD247")
pic <- list()
for(i in genes){
  pic[[i]] <- FeaturePlot(Seucount, features = i, reduction = "tsne") + 
             scale_color_gradientn(colors = c('#cacaca30','#cacaca30','#f6a11a','red','red'))
}
ggsave(paste0(getwd(), '/marker.tiff'), cowplot::plot_grid(plotlist = pic, ncol = 4), width = 12, height = 3)
ggsave(paste0(getwd(), '/marker.pdf'), cowplot::plot_grid(plotlist = pic, ncol = 4), width = 12, height = 3)

# ---------------------------
# Clustering at different resolutions
seurat <- FindClusters(Seucount, resolution = 0.8)
seurat <- FindClusters(Seucount, resolution = 1.5)

DimPlot(seurat, reduction = "tsne", label = TRUE, group.by = "seurat_clusters", pt.size = 0.5) +
    theme_minimal() +
    labs(title = "tsne - Clusters (res=1)")
ggsave("tsne_clusters_res15.pdf", width = 8, height = 6)

DimPlot(Seucount, reduction = "umap", label = TRUE, group.by = "seurat_clusters", pt.size = 0.5) +
    theme_minimal() +
    labs(title = "UMAP - Clusters (res=1)")
ggsave("umap_clusters_res1.pdf", width = 8, height = 6)