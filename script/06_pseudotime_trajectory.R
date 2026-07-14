#
seucount <- readRDS("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/B/B_Plasma_dob.rds")
setwd("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/B/monocle")

# Tumor-associated macrophages (all samples)
library(stringr); library(reshape2); library(plyr); library(Seurat); library(dplyr);
library(Matrix); library(ggplot2); library(edgeR); library(data.table); library(pheatmap);
library("ggsci"); library(monocle); library(parallel); library(magrittr); library(ggcorrplot)
library(ggpubr); library(beeswarm); library(presto);
library(GSEABase)

ls <- seucount

# Based on gene expression heatmap for CD8: clusters 1,5,12,11 (cluster 3 not obvious, 8 is NK)
table(seucount$recluster)
# Idents(object = seucount) <- seucount$recluster
Idents(object = seucount) <- seucount$recluster
seucount[['cell_type']] <- seucount@active.ident

expr_matrix <- as(as.matrix(seucount@assays$RNA@counts), 'sparseMatrix')
p_data <- seucount@meta.data
p_data$celltype <- seucount@active.ident  
f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))

pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)
expr_matrix <- expr_matrix[, rownames(pd)]
cds <- newCellDataSet(expr_matrix, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

# 2. Estimate size factors and dispersion
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 3. Filter low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)  
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))  # filter genes expressed in fewer than 10 cells

# 4. Trajectory gene selection, visualization, and trajectory construction
# Step 1: Select differentially expressed genes
diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~cell_type", cores = 1)
head(diff)

deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval, decreasing = F),]
head(deg)

####################################
write.table(deg, file = "train.monocle.DEG.xls", col.names = T, row.names = F, sep = "\t", quote = F)

# Visualize genes used for trajectory construction
# ordergene <- rownames(deg)
# If too many genes, one can select top genes
ordergene <- row.names(deg)[order(deg$qval)][1:500]

cds <- setOrderingFilter(cds, ordergene) 

###### Step 2: Dimensionality reduction
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

###### Step 3: Build pseudotime trajectory and order cells along pseudotime
cds <- orderCells(cds)
# saveRDS(cds, 'cds.rds')

plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5, show_backbone = F, show_branch_points = F)
ggsave("trajectory_State.pdf", width = 4, height = 4.5)
plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5, show_backbone = F, show_branch_points = F)
ggsave("tPseudotime.pdf", width = 4, height = 4.5)

# colors <- c("orangered1","deepskyblue")
# scale_fill_manual(values = c("Normal" = "#D0E7ED", "Mild" = "#9392BE","Severe"="#D5E4A8"))
colors <- c("#9392BE","#D0E7ED","#D5E4A8")
plot_cell_trajectory(cds, color_by = "illness", size = 1, show_backbone = TRUE) + scale_color_manual(values = colors)
ggsave("illness.pdf", width = 4, height = 4.5)

colors <- c("#E1C855","#E07B54","#51B1B7","#F1C89A","#A797DA","#97C8AF") #,"#E79397","#97C8AF")
plot_cell_trajectory(cds, color_by = "recluster", size = 1, show_backbone = TRUE) + scale_color_manual(values = colors)
ggsave("illness.pdf", width = 4, height = 4.5)

cds$znh <- paste(cds$illness, cds$recluster)
plot_cell_trajectory(cds, color_by = "znh", size = 1, show_backbone = TRUE)

## Adjust root state as needed
a <- cds
cds$State1 <- cds$State
cds$State <- cds$recluster
HSMM_1 <- orderCells(cds, root_state = "B_Naive")
cds <- HSMM_1
saveRDS(cds, 'cds1.rds')

# cds <- a

# Visualization
# 1. Color by pseudotime values
# plot_cell_trajectory(cds,color_by="State",cell_size=0.5,show_backbone=F,show_branch_points = F)
# ggsave("trajectory_State.pdf",width = 4,height = 4)
plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.5, show_backbone = F, show_branch_points = F)
ggsave("tPseudotime.pdf", width = 4, height = 4)
# colors <- c("orangered1","deepskyblue")
# scale_fill_manual(values = c("Normal" = "#D0E7ED", "Mild" = "#9392BE","Severe"="#D5E4A8"))
colors <- c("#9392BE","#D0E7ED","#D5E4A8")
plot_cell_trajectory(cds, color_by = "illness", size = 1, show_backbone = TRUE) + scale_color_manual(values = colors)
ggsave("illness.pdf", width = 4, height = 4)

colors <- c("#C6B3D3","#ED9F9B","#80BABA","#E1C855")
plot_cell_trajectory(cds, color_by = "recluster", size = 1, show_backbone = TRUE) + scale_color_manual(values = colors)
ggsave("cell.pdf", width = 4, height = 4)

ordergene <- row.names(deg)[order(deg$qval)]
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)] 
write.csv(Time_diff, "Time_diff_all.csv", row.names = F)

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# Time_diff <- read.csv("/boot2/xdh/NPC/znh/result1/T/monocle2/Time_diff_all.csv")
pdf("timdf_pseudotime_heatmap100.pdf")
plot_pseudotime_heatmap(cds[Time_diff$gene_short_name[order(Time_diff$qval)][1:100],],
                        num_clusters = 3, show_rownames = T, return_heatmap = T)
dev.off()
pdf("timdf_pseudotime_heatmap100gai.pdf", height = 7, width = 4)
plot_pseudotime_heatmap(cds[Time_diff$gene_short_name[order(Time_diff$qval)][1:100],],
                        num_clusters = 3, show_rownames = T, return_heatmap = T)
dev.off()
# ggsave("timdf_pseudotime_heatmap100.pdf")

# saveRDS(cds,"cds.rds")
getwd()

library("monocle")
library(dplyr)
library(Seurat)
library(patchwork) 
library(ggridges)
library(RColorBrewer)
library(scales)

unique(cds$seurat_clusters)
cds$znh <- paste(cds$seurat_clusters, cds$recluster)
unique(cds$znh)
test <- cds
plotdf <- pData(test)
pdf("tmp2.pdf", height = 3)
ggplot(plotdf, aes(x = Pseudotime, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 1) +
  geom_vline(xintercept = c(5,10), linetype = 2) +
  scale_y_discrete("") +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  ) + scale_color_manual(values = colors)
dev.off()

# a <- readRDS("/boot2/xdh/NPC/znh/result2/myeloid/MM-pt-lnm-mo/cds1.rds")
# a$recluster_znh <- paste(a$seurat_clusters, a$recluster)
# plot_cell_trajectory(a,color_by="recluster_znh",size=1,show_backbone=F)
# ggsave("c_rec_tra.pdf")

a <- cds
plotdf <- pData(a)
a$znh <- paste(a$seurat_clusters, a$recluster)
pdf("tmp3.pdf", height = 3)
ggplot(plotdf, aes(x = Pseudotime, y = celltype, fill = znh)) +
  geom_density_ridges(scale = 1) +
  geom_vline(xintercept = c(5,10), linetype = 2) +
  scale_y_discrete("") +
  theme_minimal() +
  theme(
    panel.grid = element_blank()
  ) + scale_color_manual(values = colors)
dev.off()

### Pathway enrichment analysis:
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)

# First gene list
gene <- read.table("ex_deg.txt") # import data
EntrezID <- bitr(gene$V1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
id_list <- EntrezID$ENTREZID
id_list <- EntrezID$SYMBOL

go <- enrichGO(gene = EntrezID$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
# ggplot(go)
go.res <- data.frame(go)   # convert GO results to data frame
goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))[1:10,]
go.df <- goBP
# Ensure GO term order matches input
go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))
# Plot
go_bar <- ggplot(data = go.df,
                 aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "GeneNumber", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
ggsave(go_bar, filename = "GO_Barplot.pdf", width = 9, height = 7)

##############################
# Second gene list
gene <- read.table("ex_deg2.txt")
EntrezID <- bitr(gene$V1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = EntrezID$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
go.res <- data.frame(go)
goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))[1:10,]
go.df <- goBP
go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))
go_bar <- ggplot(data = go.df,
                 aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "GeneNumber", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
go_bar
ggsave(go_bar, filename = "GO_Barplot2.pdf", width = 9, height = 7)

##############
# Third gene list
gene <- read.table("ex_deg3.txt")
EntrezID <- bitr(gene$V1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db")
go <- enrichGO(gene = EntrezID$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
go.res <- data.frame(go)
goBP <- subset(go.res, subset = (ONTOLOGY == "BP"))[1:10,]
go.df <- goBP
go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))
go_bar <- ggplot(data = go.df,
                 aes(x = Description, y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.9) +
  coord_flip() + theme_bw() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO terms", y = "GeneNumber", title = "Barplot of Enriched GO Terms") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
ggsave(go_bar, filename = "GO_Barplot3.pdf", width = 9, height = 7)