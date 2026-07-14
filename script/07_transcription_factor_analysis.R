library(SCENIC)
library(SCopeLoomR)
library(data.table)
library(Seurat)
library(circlize)

setwd("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/B/Pyscienic")

scenicLoomPath <- 'sample_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)
close_loom(loom)

rownames(regulonAUC)
names(regulons)

seucount <- readRDS("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/B/B_Plasma_dob.rds")

# Distinguish normal and disease groups:
seucount$sample_znh <- seucount$illness
seucount$sample_znh[which(seucount$illness == "Nomal")] <- "Normal"
seucount$sample_znh[which(seucount$illness == "Mild" | seucount$illness == "Severe")] <- "Patient"

cellTypes <- seucount$recluster
cellTypes <- as.data.frame(cellTypes)
colnames(cellTypes) <- "celltype"
rownames(cellTypes) <- names(seucount$recluster)

##############################
selectedResolution <- "celltype"  # select resolution
cellsPerGroup <- split(rownames(cellTypes),
                       cellTypes[, selectedResolution])

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(regulonAUC)[, cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale = T))
# Scale each regulon across different clusters
dim(regulonActivity_byGroup_Scaled)

regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

library(colorRamp2)
library(ComplexHeatmap)
library(pheatmap)
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from = -2, to = 2, length = 11), rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows       = "ward.D2",
  clustering_method_columns    = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

#############
library(dplyr)
rss <- regulonActivity_byGroup_Scaled
head(rss)

df <- do.call(rbind,
              lapply(1:ncol(rss), function(i) {
                dat <- data.frame(
                  path    = rownames(rss),
                  cluster = colnames(rss)[i],
                  sd.1    = rss[, i],
                  sd.2    = apply(rss[, -i], 1, median)
                )
              }))
df$fc <- df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn <- data.frame(path = top5$cluster)
n <- rss[unique(top5$path), ]

# If duplicates exist, they are removed (commented out)
# n <- rss[top5$path, ]
# duplicates <- duplicated(rownames(n))
# duplicate_indices <- which(duplicates)
# n <- n[-duplicate_indices, ]

rownames(rowcn) <- rownames(n)
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

pdf("celltype.pdf")
pheatmap(n, show_colnames = T, cluster_cols = F, annotation_row = rowcn,
         scale = "row", border = "white",
         color = c(colorRampPalette(colors = c("#4459CB", "white"))(length(bk) / 2),
                   colorRampPalette(colors = c("white", "#B1001F"))(length(bk) / 2)),
         legend_breaks = seq(-2, 2, 0.5),
         breaks = bk)
dev.off()

pdf("celltype-ncluster.pdf", height = 8, width = 7)
pheatmap(n, show_colnames = T, cluster_cols = F, annotation_row = rowcn,
         # annotation_colors = colors,
         scale = "row", border = "white",
         color = c(colorRampPalette(colors = c("#4459CB", "white"))(length(bk) / 2),
                   colorRampPalette(colors = c("white", "#B1001F"))(length(bk) / 2)),
         legend_breaks = seq(-2, 2, 0.5),
         breaks = bk,
         cluster_rows = FALSE)
dev.off()

######################
# Extract cell types combining sample status and recluster
cellTypes <- paste(seucount$sample_znh, seucount$recluster, sep = "_")
cellTypes <- as.data.frame(cellTypes)
colnames(cellTypes) <- "celltype"
rownames(cellTypes) <- names(seucount$recluster)

selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes),
                       cellTypes[, selectedResolution])

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(regulonAUC)[, cells]))

# Scale expression.
# Scale function normalizes columns; transpose so cells are rows and genes are columns
# Reference: https://www.jianshu.com/p/115d07af3029
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale = T))
dim(regulonActivity_byGroup_Scaled)  # [1] 209   9
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

############## pheatmap
rss <- regulonActivity_byGroup_Scaled
n <- rss[unique(top5$path), ]
n <- rss[top5$path, ]
# Remove duplicate rows if needed (commented)
# duplicates <- duplicated(rownames(n))
# duplicate_indices <- which(duplicates)
# n <- n[-duplicate_indices, ]

# Reorder columns as desired
id <- c("Normal_B_Memory", "Patient_B_Memory", "Normal_B_Naive", "Patient_B_Naive",
        "Normal_B_Plasma", "Patient_B_Plasma", "Normal_B_Plasmablast", "Patient_B_Plasmablast")
index <- match(id, colnames(n))
n <- n[, index]

rowcn <- data.frame(path = top5$cluster)
rownames(rowcn) <- rownames(n)
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

pdf("cell-tf.pdf")
pheatmap(n, show_colnames = T, cluster_cols = F, annotation_row = rowcn,
         scale = "row", border = "white",
         color = c(colorRampPalette(colors = c("#4459CB", "white"))(length(bk) / 2),
                   colorRampPalette(colors = c("white", "#B1001F"))(length(bk) / 2)),
         legend_breaks = seq(-2, 2, 0.5),
         breaks = bk)
dev.off()

pdf("cell-tf-cluster_rows.pdf", height = 8, width = 7)
pheatmap(n, show_colnames = T, cluster_cols = F, annotation_row = rowcn,
         # annotation_colors = colors,
         scale = "row", border = "white",
         color = c(colorRampPalette(colors = c("#4459CB", "white"))(length(bk) / 2),
                   colorRampPalette(colors = c("white", "#B1001F"))(length(bk) / 2)),
         legend_breaks = seq(-2, 2, 0.5),
         breaks = bk,
         cluster_rows = F)
dev.off()