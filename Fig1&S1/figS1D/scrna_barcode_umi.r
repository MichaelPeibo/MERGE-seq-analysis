

library(tidyverse)
library(Seurat)
library(ggsci)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/scrna_vs_cdna/python_scrna_readtable")
barcode_genes <- c("barcode0","barcode1","barcode2","barcode3","barcode4")

pfc1 <- read.csv("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/pmacs_results/bar.table_pfc_1_v2_hd2.csv",row.names = 1)
pfc2 <- read.csv("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/pmacs_results/bar.table_pfc_2_v2_hd2.csv",row.names = 1)
pfc3 <- read.csv("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/pmacs_results/bar.table_pfc_3_v2_hd2.csv",row.names = 1)
pfc4 <- read.csv("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/pmacs_results/bar.table_pfc_4_v2_hd2.csv",row.names = 1)

add_prefix_to_rownames <- function(df, prefix) {
  df %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::mutate(rowname = paste(prefix, rowname, sep = "_")) %>%
    tibble::column_to_rownames("rowname")
}
pfc1 <- add_prefix_to_rownames(pfc1, "1")
pfc2 <- add_prefix_to_rownames(pfc2, "2")
pfc3 <- add_prefix_to_rownames(pfc3, "3")
pfc4 <- add_prefix_to_rownames(pfc4, "4")
pfc_cdna <- rbind(pfc1,pfc2,pfc3,pfc4)
pfc_cdna <- pfc_cdna[,c(1:6)]
colnames(pfc_cdna) <- c(barcode_genes,'nUMI')
pfc_cdna$lib <- "cdna"

load('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all <- mda; rm(mda)
pfc_all_meta <- read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/pfc_seurat_merge_meta.csv',header = T,row.names = 1)
pfc_all <- subset(pfc_all, cells = rownames(pfc_all_meta))

bar.mat <- read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
## filter pfc_all by bar.mat cells ##
pfc_all <- subset(pfc_all, cells = rownames(bar.mat))
## filter pfc_all_meta by pfc_all cells ##
pfc_all_meta <- pfc_all_meta[colnames(pfc_all),]

pfc_all_bar <- subset(pfc_all,features = barcode_genes)
pfc_all_barmat <- as.data.frame(pfc_all_bar@assays$RNA@counts)
pfc_all_barmat <- as.data.frame(t(pfc_all_barmat))
pfc_all_barmat[,'nUMI']=rowSums(pfc_all_barmat[,1:5])
pfc_all_barmat$lib <- "scrna"
pfc_combine <- rbind(pfc_all_barmat,pfc_cdna)


df_long <- pfc_combine %>%
  pivot_longer(
    cols = c(barcode_genes,'nUMI'),
    names_to = "barcode",
    values_to = "counts"
  )

ggplot(df_long, aes(x = lib, y = counts+1, fill=barcode)) +
  geom_boxplot() +  
  scale_y_continuous(trans = 'log10', labels = scales::label_log())+
  labs(x = " ", y = "Counts") +
  theme_Publication(base_size = 12) +
  scale_fill_aaas() +
  theme(axis.title.x = element_text(size = 20),  
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),  
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 20),
    axis.ticks.length = unit(0.5, "cm"))
ggsave("boxplot scrna counts vs cdna counts.pdf", 
       width = 8, height = 8)
