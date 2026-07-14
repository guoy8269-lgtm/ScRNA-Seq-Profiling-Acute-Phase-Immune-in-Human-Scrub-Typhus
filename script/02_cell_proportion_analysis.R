library(ggh4x)
# ======================
# Stacked bar plot of cell proportions
# ======================
seucount <- readRDS("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/rec.rds") # rec.rds contains single-cell data for cell proportion analysis
colors <- c("#E1C855", "#E07B54", "#51B1B7", "#F1C89A", "#A797DA", "#97C8AF") # Adjust based on cell types to display
# Custom theme (replace with your own mytheme if already defined)
mytheme <- theme_bw() + 
  theme(plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(1.2)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = rel(1.5), color = 'black'),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = rel(1.5), color = 'black'),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(color = "black", size = 1),
        axis.line = element_line(color = "black", size = 0.5))

# Split metadata by illness status
N <- seucount@meta.data[which(seucount$illness == "Nomal"), ]
P <- seucount@meta.data[-which(seucount$illness == "Nomal"), ]

# Normal group counts
coun <- as.data.frame(as.matrix(table(N$recluster), , 2))
coun
for (i in 1:nrow(coun)) {
  coun$percent[i] <- coun[i, 1] / sum(coun$V1)
}
coun
coun$sample <- "normal"
normal <- coun

# Patient group counts
coun <- as.data.frame(as.matrix(table(P$recluster), , 2))
coun
for (i in 1:nrow(coun)) {
  coun$percent[i] <- coun[i, 1] / sum(coun$V1)
}
coun
coun$sample <- "patient"
patient <- coun

normal$cell <- rownames(normal)
patient$cell <- rownames(patient)

data <- cbind(normal, patient)
data
data <- rbind(normal, patient)
# data$cell <- factor(data$cell, level=c("B","Plasma","Myeloid","NK","T"))

p <- ggplot(data, aes(x = sample, weight = percent, fill = cell)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = colors) +
  mytheme +
  labs(y = "Percent")
p
ggsave("sample_percent_n_p.pdf", width = 3)


# ======================
# Boxplots of mean proportion changes per cell type
# ======================
setwd("/groups/g900008/home/zhangnihui/实验备份/羌虫病/重新分析(美化图片)/细胞比例6月")
seucount <- readRDS("/groups/g900008/home/zhangnihui/实验备份/羌虫病/SingleR部分/rec.rds")

N <- seucount@meta.data[which(seucount$illness == "Nomal"), ]
P <- seucount@meta.data[-which(seucount$illness == "Nomal"), ]

result <- c()
a <- unique(unique(N$orig.ident))

for (i in 1:length(a)) {
  pri1 <- N[which(N$orig.ident == a[i]), ]
  number <- nrow(pri1)
  result1 <- table(pri1$recluster) / number
  result1 <- cbind(as.numeric(result1), names(result1), "N")
  result <- rbind(result, result1)
}

a <- unique(unique(P$orig.ident))
for (i in 1:length(a)) {
  pri1 <- P[which(P$orig.ident == a[i]), ]
  number <- nrow(pri1)
  result1 <- table(pri1$recluster) / number
  result1 <- cbind(as.numeric(result1), names(result1), "P")
  result <- rbind(result, result1)
}


library(ggplot2)
library(ggpubr)   # for stat_compare_means
library(patchwork)

colnames(result) <- c("p", "celltype", "sample")

df <- as.data.frame(result)
df$p <- as.numeric(df$p)   # convert p column from character to numeric


sample_colors <- c("N" = "#62B197", "P" = "#E18E6D")
celltype_colors <- colors

p <- ggplot(df, aes(x = sample, y = p)) +
  geom_boxplot(aes(color = sample), fill = NA,
               outlier.shape = 21, outlier.size = 1.5) +
  geom_jitter(aes(color = sample), width = 0.2, size = 1, alpha = 0.6) +
  scale_color_manual(values = sample_colors) +
  # Use ggh4x::facet_wrap2 with themed strips
  facet_wrap2(~ celltype, nrow = 1, scales = "free_y",
              strip = strip_themed(background_x = elem_list_rect(fill = celltype_colors))) +
  ylab("Proportion of cells") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.format",
                     comparisons = list(c("N", "P")),
                     tip.length = 0) +
  theme(legend.title = element_blank()) +
  mytheme

ggsave("All_cell_boxplots_per_celltype.pdf", plot = p, width = 10, height = 5)

# Calculate mean proportions per cell type and sample
library(dplyr)
result_df <- as.data.frame(df)

result_summary <- result_df %>%
  group_by(celltype, sample) %>%
  summarize(mean_p = mean(as.numeric(p)))

print(result_summary)

p1 <- ggplot(result_summary, aes(x = sample, y = mean_p, group = sample, color = celltype)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_line(aes(group = celltype), position = position_dodge(width = 0.7), size = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  labs(title = "Mean p-value by Sample and Cell Type", x = "Sample", y = "Mean p-value") +
  theme(legend.title = element_blank())

ggsave("all-box.pdf", plot = p1, height = 3, width = 3.5)