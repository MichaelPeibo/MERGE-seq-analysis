
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")
library(tidyverse)
library(Seurat)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")
barcode_genes <- c("barcode0","barcode1","barcode2","barcode3","barcode4")

load('./data_file/pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all <- mda; rm(mda)
pfc_all_meta <- read.csv('./data_file/pfc_seurat_merge_meta.csv',header = T,row.names = 1)
pfc_all <- subset(pfc_all, cells = rownames(pfc_all_meta))

pfc_all_bar <- subset(pfc_all,features = barcode_genes)
pfc_all_barmat <- as.data.frame(pfc_all_bar@assays$RNA@counts)
pfc_all_barmat <- as.data.frame(t(pfc_all_barmat))
pfc_all_barmat[,'EGFP']=rowSums(pfc_all_barmat[,1:5])
## find which cells have EGFP >0 ##
test <- pfc_all_barmat[pfc_all_barmat$EGFP > 0,]
pfc_all <- SetIdent(pfc_all, cells = colnames(pfc_all), value = 'EGFP-')
pfc_all <- SetIdent(pfc_all, cells = rownames(test), value = 'EGFP+')
table(Idents(pfc_all))
pfc_all@meta.data$egfp <- Idents(pfc_all)

temp <- pfc_all@meta.data %>% dplyr::select(sample,egfp)
df <- as.data.frame(table(temp$sample,temp$egfp))
df <- df %>% as_tibble() %>% group_by(Var1) %>% 
  mutate(ratio=Freq/sum(Freq)) %>% mutate(ratio=round(ratio,2))
df$ratio <- round(df$ratio*100,2)

ggplot(df, aes(x = Var1, y = ratio, fill = Var2)) +
  geom_bar(stat = "identity", position = ("stack")) +
  scale_fill_manual(values = c("#00c206", "grey")) +
  geom_text(aes(label = ratio), position = position_dodge(width = 0.9), vjust = -0.25, color = "black", size = 7, show.legend = FALSE) +
  theme_Publication() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20, face = "bold"),
    axis.ticks.length = unit(0.5, "cm")
  ) + labs(y = "Ratio")
ggsave('./results/figs1_scrnaseq_EGFP_cell distribution in sort and unsort.pdf',height = 5,width = 5)
