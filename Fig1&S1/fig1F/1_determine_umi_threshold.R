
### 1 ###
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")
library(tidyverse)
library(Seurat)
library(wesanderson)
library(ggpubr)
library(ggbeeswarm)
library(ggsci)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")

### read in raw barcode matrix
barcode_genes <- c("barcode0","barcode1","barcode2","barcode3","barcode4")
pfc1 <- read.csv("./pmacs_results/bar.table_pfc_1_v2_hd2.csv",row.names = 1)
pfc2 <- read.csv("./pmacs_results/bar.table_pfc_2_v2_hd2.csv",row.names = 1)
pfc3 <- read.csv("./pmacs_results/bar.table_pfc_3_v2_hd2.csv",row.names = 1)
pfc4 <- read.csv("./pmacs_results/bar.table_pfc_4_v2_hd2.csv",row.names = 1)


threshold <- quantile(pfc1$nUMI, probs = 0.95); threshold ## 1257, 403
pfc1_filter <- pfc1[pfc1$nUMI < 403,]

threshold <- quantile(pfc2$nUMI, probs = 0.95); threshold ## 486, 1488
pfc2_filter <- pfc2[pfc2$nUMI < 486,]

threshold <- quantile(pfc3$nUMI, probs = 0.95); threshold ## 529, 1699
pfc3_filter <- pfc3[pfc3$nUMI < 529,]

threshold <- quantile(pfc4$nUMI, probs = 0.95); threshold ## 1012, 1929
pfc4_filter <- pfc4[pfc4$nUMI < 1012,]

matrix_list <- list(pfc1_filter,pfc2_filter,pfc3_filter,pfc4_filter)
matrix_list_trim <- vector("list", 4)
for (i in 1:4){
  temp=matrix_list[[i]]
  temp=temp[,1:5]  
  rownames(temp) <- paste(sprintf("pfc_%s",i), rownames(temp), sep = "_")
  colnames(temp)=barcode_genes
  matrix_list_trim[[i]]=temp
}
bar.mat <- rbind(matrix_list_trim[[1]],matrix_list_trim[[2]],matrix_list_trim[[3]],matrix_list_trim[[4]])
write.csv(bar.mat, file = './results/bar_mat_filter_high_copy_quantile_005.csv')

### 3 ###
### non-neuron cells ###
non.neuron.cells <- read.csv('./data_file/non.neuron.cells.csv',row.names = 1);dim(non.neuron.cells)
## filter bar.mat by non.neuron.cells ## 
non.neuron.bar.mat <- bar.mat[rownames(bar.mat) %in% rownames(non.neuron.cells), ]
non.neuron.bar.mat[,'EGFP'] <- rowSums(non.neuron.bar.mat[,1:5]);dim(non.neuron.bar.mat)

pfc_4 <- Read10X(data.dir = "pfc_4")
pfc_4 <- CreateSeuratObject(counts = pfc_4, project = "pfc_4", min.cells = 3, min.features = 200)
pfc_4
pfc_4_bar <- subset(pfc_4,features = barcode_genes)
pfc4_barmat <- as.data.frame(pfc_4_bar@assays$RNA@counts)
pfc4_barmat <- as.data.frame(t(pfc4_barmat))
pfc4_barmat[,'EGFP']=rowSums(pfc4_barmat[,1:5])
sum(pfc4_barmat$EGFP > 0, na.rm = TRUE)/nrow(pfc4_barmat) # 0.7427466, we will use 0.7427466 for later analysis;

rownames(pfc4_barmat) <- gsub("-1", "", rownames(pfc4_barmat))
rownames(pfc4_barmat) <- paste(sprintf("pfc_4"), rownames(pfc4_barmat), sep = "_")


scrna_gfp_neg <- rownames(pfc4_barmat[pfc4_barmat$EGFP == 0,]); length(scrna_gfp_neg)
scrna_gfp_pos <- rownames(pfc4_barmat[pfc4_barmat$EGFP > 0,]); length(scrna_gfp_pos)

### filter bar.mat by scrna_gfp_pos ###
sort_pos_barmat <- bar.mat[rownames(bar.mat) %in% scrna_gfp_pos, ]; dim(sort_pos_barmat)

### filter bar.mat by scrna_gfp_neg ###
sort_neg_barmat <- bar.mat[rownames(bar.mat) %in% scrna_gfp_neg, ]; dim(sort_neg_barmat)


### 4 ###
for (i in 1:5){
    df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("Sorted",nrow(pfc4)))),
  umi = c(non.neuron.bar.mat[,i], pfc4[,i])
)
#Recreate ecdf data of non-neuron
nonneuron <- df %>% filter(celltype == 'Non-neuron')
dat_ecdf <- 
  data.frame(x = unique(nonneuron$umi),
             y = ecdf(nonneuron$umi)(unique(nonneuron$umi))*length(nonneuron$umi))
#rescale y to 0,1 range
dat_ecdf$y <- scale(dat_ecdf$y,center = min(dat_ecdf$y),scale = diff(range(dat_ecdf$y)))
idx <- which(abs(dat_ecdf$y - 0.999) == min(abs(dat_ecdf$y - 0.999))) 
print(dat_ecdf[idx,])## print threshold results ##
}

### sort-neg ###
for (i in 1:5){
    df <- data.frame(
  celltype = factor(c(rep("Sort-neg", nrow(sort_neg_barmat)),
                      rep("Sorted",nrow(pfc4)))),
  umi = c(sort_neg_barmat[,i], pfc4[,i])
)
#Recreate ecdf data of non-neuron
sort_neg <- df %>% filter(celltype == 'Sort-neg')
dat_ecdf <- 
  data.frame(x = unique(sort_neg$umi),
             y = ecdf(sort_neg$umi)(unique(sort_neg$umi))*length(sort_neg$umi))
#rescale y to 0,1 range
dat_ecdf$y <- scale(dat_ecdf$y,center = min(dat_ecdf$y),scale = diff(range(dat_ecdf$y)))
idx <- which(abs(dat_ecdf$y - 0.999) == min(abs(dat_ecdf$y - 0.999))) 
print(dat_ecdf[idx,])
}

#####################
#### plot density ######
#####################

i=1
df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("FAC-sorted",nrow(sort_pos_barmat)),
                      rep("EGFP-negative", nrow(sort_neg_barmat)))
                      ),
  umi = c(non.neuron.bar.mat[,i], sort_pos_barmat[,i], sort_neg_barmat[,i])
)

ggplot(df, aes(umi,color = celltype)) +
  scale_color_aaas() +
  geom_density()+
  scale_x_log10(breaks = c(1,2,3,28))+
  geom_vline(xintercept = 28, linetype="dashed", color = "black")+
  labs(
    x ="UMI-AI", y = "Density")+
  theme_Publication()+
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/density_plot_umi_ai.pdf", width = 4, height = 4)

i=2
df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("FAC-sorted",nrow(sort_pos_barmat)),
                      rep("EGFP-negative", nrow(sort_neg_barmat)))
  ),
  umi = c(non.neuron.bar.mat[,i], sort_pos_barmat[,i], sort_neg_barmat[,i])
)

ggplot(df, aes(umi,color = celltype)) +
  scale_color_aaas() +
  geom_density()+
  scale_x_log10(breaks = c(1,2,3,101))+
  geom_vline(xintercept = 101, linetype="dashed", color = "black")+
  labs(
    x ="UMI-DMS", y = "Density")+
  theme_Publication()+
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/density_plot_umi_dms.pdf", width = 4, height = 4)

i=3
df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("FAC-sorted",nrow(sort_pos_barmat)),
                      rep("EGFP-negative", nrow(sort_neg_barmat)))
  ),
  umi = c(non.neuron.bar.mat[,i], sort_pos_barmat[,i], sort_neg_barmat[,i])
)
ggplot(df, aes(umi,color = celltype)) +
  scale_color_aaas() +
  geom_density()+
  scale_x_log10(breaks = c(1,2,3,114))+
  geom_vline(xintercept = 114, linetype="dashed", color = "black")+
  labs(
    x ="UMI-MD", y = "Density")+
  theme_Publication()+
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/density_plot_umi_md.pdf", width = 4, height = 4)

i=4
df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("FAC-sorted",nrow(sort_pos_barmat)),
                      rep("EGFP-negative", nrow(sort_neg_barmat)))
  ),
  umi = c(non.neuron.bar.mat[,i], sort_pos_barmat[,i], sort_neg_barmat[,i])
)
ggplot(df, aes(umi,color = celltype)) +
  scale_color_aaas() +
  geom_density()+
  scale_x_log10(breaks = c(1,2,3,35))+
  geom_vline(xintercept = 35, linetype="dashed", color = "black")+
  labs(
    x ="UMI-BLA", y = "Density")+
  theme_Publication()+
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/density_plot_umi_bla.pdf", width = 4, height = 4)

i=5
df <- data.frame(
  celltype = factor(c(rep("Non-neuron", nrow(non.neuron.bar.mat)),
                      rep("FAC-sorted",nrow(sort_pos_barmat)),
                      rep("EGFP-negative", nrow(sort_neg_barmat)))
  ),
  umi = c(non.neuron.bar.mat[,i], sort_pos_barmat[,i], sort_neg_barmat[,i])
)
ggplot(df, aes(umi,color = celltype)) +
  scale_color_aaas() +
  geom_density()+
  scale_x_log10(breaks = c(1,2,3,103))+
  geom_vline(xintercept = 103, linetype="dashed", color = "black")+
  labs(
    x ="UMI-LH", y = "Density")+
  theme_Publication()+
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/density_plot_umi_lh.pdf", width = 4, height = 4)
