
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")
library(ggsci)
library(Seurat)
library(tidyverse)
library(rstatix)
library(ggpubr)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")

load('./data_file/sub_trim_seurat_pipe.Robj')
exn.sub <- mda; rm(mda)
exn.sub <- UpdateSeuratObject(exn.sub)
bar.mat <- read.csv('./results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)

AI_valid <- rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid <- rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid <- rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid <- rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid <- rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
length(AI_valid)
length(DMS_valid)
length(MD_valid)
length(BLA_valid)
length(LH_valid)
unique_valid_cells <- unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))
length(unique_valid_cells)

valid_cells_list=vector('list',5)
valid_cells_list[[1]]=AI_valid
valid_cells_list[[2]]=DMS_valid
valid_cells_list[[3]]=MD_valid
valid_cells_list[[4]]=BLA_valid
valid_cells_list[[5]]=LH_valid
valid_cells_names <- c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
barcode_names <- c('AI','DMS','MD','BLA','LH')
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)

### filter exn.sub by bar.mat cells ###
exn.sub <- subset(exn.sub, cells = rownames(bar.mat))
### filter exn_meta by exn.sub cells ###
exn_meta <- exn_meta[colnames(exn.sub),]

for (i in 1:5) {
  exn.sub <- SetIdent(exn.sub, cells = colnames(exn.sub), value = 'Others')
  exn.sub <- SetIdent(exn.sub, cells = valid_cells_list[[i]], value =valid_cells_names[i])
  temp <- as.data.frame(Idents(exn.sub))
  colnames(temp)=valid_cells_names[i]
  exn_meta=cbind(exn_meta,temp)
}
write.csv(exn_meta,file='./results/exn_meta_valid.csv')

exn_meta <- read.csv('./results/exn_meta_valid.csv',row.names = 1)
table(exn_meta$cluster)
exn_meta %>% group_by(cluster) %>% summarise(., count = n()) %>% mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))

table(exn_meta$AI_valid)
table(exn_meta$DMS_valid)
table(exn_meta$MD_valid)
table(exn_meta$BLA_valid)
table(exn_meta$LH_valid)

sub_exn_palette <- c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')
## define a null dataframe ##
df <- data.frame(Var1 = factor(), Freq = numeric(), Percentage = numeric())

for (i in 1:5) {
  temp <- exn_meta %>% as_tibble() %>% 
    dplyr::select(valid_cells_names[i],cluster) %>% 
    filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]])
  temp <- as.data.frame(table(temp$cluster))
  temp$Var1 <- factor(temp$Var1,levels=c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6', 
                                  'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6'))
  temp$Percentage <- temp$Freq / sum(temp$Freq) * 100
  ## add a column, named as valid_cells_names[i] ##
  temp <- cbind(temp, valid_cells_names[i])
  ### rbind the temp to df ###
  df <- rbind(df, temp)
}
colnames(df) <- c('cluster', 'Freq', 'Percentage', 'barcode')

### make a bar plot of df using ggplot2 ###
### change factor order of 'barcode' ###
df$barcode <- factor(df$barcode, levels = c('AI_valid', 'DMS_valid', 'BLA_valid', 'MD_valid','LH_valid'))
ggplot(df, aes(x = barcode, y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  ylab("Ratio") +
  scale_fill_manual(values = sub_exn_palette) +
  theme_Publication()+ 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave('./results/fig2_barplot_cluster_composition_in_barcoded_cells.pdf')

###########################
exn_sub_names <- c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6','L5-S100b', 'L6-Npy', 'L6-Syt6')#'L5-Htr2c'
## define a null dataframe ##
df <- data.frame(Var1 = factor(), Freq = numeric(), Percentage = numeric())

for (i in 1:6) {
  temp <- exn_meta %>% as_tibble() %>% filter(cluster==exn_sub_names[[i]]) %>%
    select(cluster,AI_valid,DMS_valid,BLA_valid,MD_valid,LH_valid)
  df_ai=as.data.frame(table(temp$AI_valid))
  df_dms=as.data.frame(table(temp$DMS_valid))
  df_md=as.data.frame(table(temp$MD_valid))
  df_bla=as.data.frame(table(temp$BLA_valid))
  df_lh=as.data.frame(table(temp$LH_valid))
  temp=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
  temp$Percentage <- temp$Freq / sum(temp$Freq) * 100
  ## add a column, named as valid_cells_names[i] ##
  temp <- cbind(temp, exn_sub_names[i])
  ### rbind the temp to df ###
  df <- rbind(df, temp)
}
colnames(df) <- c('barcode', 'Freq', 'Percentage', 'cluster')
## cacluate the sum of Freq grouped by each cluster
df %>% group_by(cluster) %>% summarise(., count = sum(Freq))

ggplot(df, aes(x = cluster, y = Freq, fill = barcode)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  ylab("Ratio") +
  scale_fill_npg() +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave('./results/fig2_stacked_barplot_region_composition_in_exn_subtype.pdf')


#### plot by mouse or facs ####
sub_exn_palette <- c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17','#666666')
exn_meta$cluster <- fct_relevel(exn_meta$cluster, 
                        "L2/3-Calb1", "L2/3-Rorb", "L5-Bcl6", "L5-Htr2c", "L5-S100b", "L6-Npy", "L6-Syt6")
all_clusters <- c("L2/3-Calb1", "L2/3-Rorb", "L5-Bcl6", "L5-Htr2c", "L5-S100b", "L6-Npy", "L6-Syt6")
# Function to map cluster names to colors
get_cluster_color <- function(cluster_name) {
  if (cluster_name %in% all_clusters) {
    return(sub_exn_palette[which(all_clusters == cluster_name)])
  } else {
    return("#999999")  # Default color for missing clusters
  }
}

for (i in 1:5) {
  df <- exn_meta %>% select(valid_cells_names[i],facs,cluster) %>% 
    filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
    group_by(facs,cluster) %>% summarise(., count = n()) %>% 
    mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
  df$cluster <- factor(df$cluster, levels = all_clusters)
  ggplot(df, aes(x = facs, y = ratio, color = cluster, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    xlab("") +
    ylab("Ratio") +
    scale_color_manual(values = sapply(levels(df$cluster), get_cluster_color)) +  # Use custom color palette
    scale_fill_manual(values = sapply(levels(df$cluster), get_cluster_color)) + 
    theme_Publication() + 
    theme(axis.ticks.length = unit(0.5, "cm"))
  ggsave(paste0('./results/fig2_', valid_cells_names[i],'_projection_motif_by_cluster.pdf'),height = 5,width = 3)
}

#######################
df <- exn_meta %>% 
  group_by(cluster,barcoded) %>% dplyr::count(name = "CellCount") %>% ungroup()

p <- ggbarplot(df, x = "cluster", y = "CellCount",
            color = "barcoded", fill = "barcoded",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",repel=T,lab.size = 7,
            xlab ="Cell type", 
            ylab = 'Cell Count', ggtheme=theme_pubclean()) + 
  rotate_x_text(45)+
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p
ggsave('./results/fig2_barcoded cell distribution.pdf',height = 10,width = 10)