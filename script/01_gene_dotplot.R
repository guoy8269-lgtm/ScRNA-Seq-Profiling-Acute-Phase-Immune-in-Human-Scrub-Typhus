###############################
library(Seurat)
# ALL cells
seucount<-readRDS("/boot2/xdh/实验2/数据整合/修改标准化/数据整合/SingleR/rec.rds") 
genes = c('MS4A1','BANK1','CD79A','FCRL1', 'AFF3',  # B
          "JCHAIN","MZB1",'TNFRSF17','DERL3',      # Plasma (also B markers, merged into B)
          "MS4A6A",'CD14', 'CD68','S100A8','S100A9', # Myeloid
          "NKG7",'GNLY', 'KLRB1', 'KLRD1','KLRF1',  # NK (TRAC is also T marker, merged later)
          'CD3E', 'CD3G','CD3D','IL32','LCK'        # T
)
limiy = c("B","Plasma",'Myeloid','NK','T')
p <- DotPlot(seucount, features = genes) +
  scale_color_gradientn(colours = c('white','red')) + #, limits = c(0, 2.5),
  #na.value = "white")+
  #scale_x_discrete(limits = limix)+
  scale_y_discrete(limits = limiy) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_size(range = c(1, 10))
ggsave("qp-ALL2.pdf", p, width = 11, height = 5)
#ggsave("qp-ALL.pdf", p, width = 11, height = 5)

#############################
# Myeloid
seucount <- readRDS("/groups/g900008/home/zhangnihui/实验备份/羌虫病/重新分析(美化图片)/Myeloid/myeloid-gai.rds")
limiy <- c("mDC", "pDC", "M-LRMDA", "M-GZMA", "M-CSF3R", "M-CD16", "M-CD14")
# DotPlot with selected genes
gene_of_interest <- c(
  # pDC
  "LILRA4", "CLEC4C", "PLD4", "IRF7", "TLR9",
  # mDC
  "CD1C", "HLA-DRA", "HLA-DRB1", "FCER1A", "CLEC10A",
  # M-LRMDA
  "LRMDA", "PLXDC2", "ZEB2", "NEAT1",
  # M-GZMA
  "GZMA", "CCL5", "GZMH", "PRF1", "GNLY",
  # M-CSF3R
  "CSF3R", "FCGR3B", "S100A8", "S100A9", "MMP9", "IL1R2",
  # M-CD16
  "FCGR3A", "CSF1R", "C1QA", "C1QB", "C1QC", "CX3CR1", "LILRB1",
  # M-CD14
  "VCAN", "CD14", "CST3", "MS4A6A", "FCN1"
)
p <- DotPlot(seucount, features = gene_of_interest) +
  scale_color_gradientn(colours = c('white','red')) +
  scale_y_discrete(limits = limiy) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_size(range = c(1, 10))

ggsave("qp-mm2.pdf", p, width = 20, height = 5)

# B cells
seucount <- readRDS("/groups/g900008/home/zhangnihui/实验备份/羌虫病/SingleR部分/B/B_Plasma_dob.rds")
# B_Naive, B_Memory ('CD27', "CD19", "MS4A1"), B_Plasma, B_Plasmablast
genes = c('TCL1A','IGHD', 'FCER2', 'IL4R','BACH2',         # B_Naive
          'CD82','TFEC','CD27',"CD19","MS4A1", 'AIM2','CD79A', # B_Memory (CD27 also in multiple B subsets)
          "CD38",'MZB1', 'XBP1', 'TNFRSF17','SDC1','IRF4',   # B_Plasma
          "MKI67",'CDK1', 'TYMS', 'TOP2A', 'CCNB1')          # B_Plasmablast ("CD38",'MZB1','TNFRSF17' also markers for Plasma)
limiy = c("B_Naive","B_Memory","B_Plasma","B_Plasmablast")
p <- DotPlot(seucount, features = genes) +
  scale_color_gradientn(colours = c('white','red')) +
  scale_y_discrete(limits = limiy) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_size(range = c(1, 10))
ggsave("qp-B.pdf", p, width = 11, height = 5)

#############################
############# T cell major types
########
seucount <- readRDS("/groups/g900008/home/zhangnihui/实验备份/羌虫病/SingleR部分/T2/CD4CD8NKγδ.rds")
genes <- c('CD8A','CD8B','LAG3','GZMA','PLEK',   # CD8+ T
           'CD4','CCR7',"IL7R","TCF7",'LTB','CD40LG','IL6R', # CD4+ T
           'TRGV9','TRDV2','KLRB1','TRDC',       # γδ T (TRGV9, TRDV2, KLRB1, KLRG1)
           "GNLY",'KLRF1','KLRD1','TYROBP','AOAH','FCGR3A') # NK (CD56=NCAM1, CD16=FCGR3A)
p <- DotPlot(seucount, features = genes) +
  scale_color_gradientn(colours = c('white','red')) + #, limits = c(0, 2.5),
  #na.value = "white")+
  #scale_x_discrete(limits = limix)+
  #scale_y_discrete(limits = limiy)+ 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_size(range = c(1, 10))
ggsave("qp-ALL-T4.pdf", p, width = 9.5, height = 4.5)

### CD8+ T
seucount <- readRDS("/groups/g900008/home/zhangnihui/实验备份/羌虫病/SingleR部分/T2/CD83/cd8.rds")
Idents(seucount) <- seucount$recluster
genes <- c('CCR7', 'TCF7', 'LEF1', 'SELL','IL7R', 'NELL2',             # CD8_Naive
           'GZMK','CD27','LTB','CD44',                                  # CD8_Tcm (also high 'IL7R','CCR7','NELL2')
           'SLC4A10', 'KLRB1', 'RORC',                                  # CD8_TcEtr (also 'GZMK','LTB','IL7R')
           'CHI3L2', 'CD38', 'LAG3',                                    # CD8_eMem (high 'GZMK','GZMH','CD27')
           'GZMH','GNLY','FCRL6', 'CCL4', 'KLRG1', 'NKG7',              # CD8_Eff_GNLY
           'TYMS', 'MKI67', 'RRM2', 'STMN1', 'PCNA')                   # CD8_Pro
limiy <- c("CD8_Naive","CD8_Tcm","CD8_TcEtr","CD8_eMem","CD8_Eff_GNLY","CD8_Pro")
p <- DotPlot(seucount, features = genes) +
  scale_color_gradientn(colours = c('white','red')) + #, limits = c(0, 2.5),
  #na.value = "white")+
  #scale_x_discrete(limits = limix)+
  scale_y_discrete(limits = limiy) + 
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_size(range = c(1, 10))
ggsave("qp-CD8-T.pdf", p, width = 12, height = 4.5)