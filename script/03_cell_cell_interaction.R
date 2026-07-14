# =============================================
# Primary and lymphatic metastasis signaling pathways
# =============================================
cellchat_nomal <- readRDS("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/CellChat/cellchat_nomal.rds")
cellchat_Patient <- readRDS("//boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/CellChat/cellchat_Patient.rds")
cellchat.normal <- cellchat_nomal
cellchat.Patient <- cellchat_Patient

library(CellChat)
cellchat.normal <- updateCellChat(cellchat_nomal)
cellchat.Patient <- updateCellChat(cellchat_Patient)

levels(cellchat.normal@idents)
levels(cellchat.Patient@idents)
object.list <- list(normal = cellchat.normal, Patient = cellchat.Patient)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# netVisual_diffInteraction(cellchat)
setwd("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/CellChat2/N-P")
cellchat <- updateCellChat(cellchat)

pdf("netVisual_diffInteraction.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T, label.edge = T, measure = "count", arrow.size = 0.5)
dev.off()

pdf("netVisual_diffInteraction_weight.pdf")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", label.edge = T, arrow.size = 0.5)
dev.off()

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf("heat_c+wei.pdf", height = 4, width = 7)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 2,
                        title.name = "Increased signaling in Patient", angle.x = 45, remove.isolate = T)
pdf("Increased signaling in Patient.pdf", height = 6, width = 10)
gg1
dev.off()

# Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, comparison = c(1, 2), max.dataset = 1,
                        title.name = "Decreased signaling in Patient", angle.x = 45, remove.isolate = T)
pdf("Decreased signaling in Patient.pdf", height = 3.5, width = 8)
gg2
dev.off()

# =============================================
# Identify signaling groups based on functional similarity
# =============================================
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# > Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
# > Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional", nCores = 5, k = 2)
# > Classification learning of the signaling networks for datasets 1 2

# Visualization in 2D-space
pdf("function.pdf", height = 10, width = 9)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf("function1.pdf", height = 4, width = 4.5)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5, groupNum = 2)

pdf("function_distance.pdf", height = 3, width = 2.7)
rankSimilarity(cellchat, type = "functional")
dev.off()

# =============================================
# Compare overall information flow of each signaling pathway
# =============================================
pdf("function_信息流.pdf", height = 4, width = 3)   # filename kept as original (includes Chinese characters)
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
dev.off()

