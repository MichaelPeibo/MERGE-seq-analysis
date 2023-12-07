
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

library(tidyverse)
library(plyr)
library(UpSetR)
library(ggpubr)
library(EnhancedVolcano)
library(wesanderson)
library(Seurat)
library(readxl)
library(rstatix)
library(ggsci)
library(ComplexHeatmap)
library(dendextend)
library(ggalluvial)
#library(srvyr)
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
unique_valid_cells <- unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid));length(unique_valid_cells)
upset_list <- list(AI=AI_valid,DMS=DMS_valid,MD=MD_valid,BLA=BLA_valid,LH=LH_valid)

### find cells in bar.mat that AI, DMS, MD, BLA, LH > 0 ###
#bar.mat %>% filter(AI>0 & DMS>0 & LH>0 & MD>0 & BLA>0) %>% rownames() 

overlapGroups <- function (listInput, sort = TRUE) {
#https://github.com/hms-dbmi/UpSetR/issues/85
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = "+")]] <- myelements
    myelements
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}


li <- overlapGroups(upset_list)
df <- ldply (li, data.frame)
colnames(df) <- c('pattern','cellid')
temp_df <- as.data.frame(table(df$pattern))
rownames(temp_df) <- temp_df$Var1
temp_df[,1] <- NULL
  
targets <- c('AI','DMS','MD','BLA','LH')
df_pattern <- data.frame(matrix(ncol = 1, nrow = 2^length(targets)))
n=1
for (x in 1:5) {
  temp=combn(targets,x)
  print(temp)
  for (i in 1:ncol(temp)) {
    rownames(df_pattern)[n]=str_c(temp[,i],collapse = '+')
    n=n+1
    #print(n)
  }
}
rownames(df_pattern)
colnames(df_pattern) <- 'Freq'
df_pattern <- df_pattern %>% rownames_to_column
temp_df <- temp_df %>% rownames_to_column
final_df <- left_join(df_pattern,temp_df,by = 'rowname')
final_df[is.na(final_df)] <- 0
final_df[,'Freq.x'] <- NULL
colnames(final_df)[2] <- 'observed'
write.csv(final_df,file='./results/obeserved_counts.csv')

final_df <- read.csv('./results/obeserved_counts.csv',row.names = 1)
patterns.counter=vector('character',4)
patterns.counter[1]=sum(final_df$observed[1:5])
patterns.counter[2]=sum(final_df$observed[6:15])
patterns.counter[3]=sum(final_df$observed[16:25])
patterns.counter[4]=sum(final_df$observed[26:32])

patterns.counter=as.numeric(patterns.counter)
names(patterns.counter)=c('1_target','2_targets','3_targets','>3_targets')
df_counter=patterns.counter %>% as_tibble()
df_counter=cbind(df_counter,rep(sum(df_counter$value),4))
colnames(df_counter)=c('counts','all')
rownames(df_counter)=c('1_target','2_targets','3_targets','>3_targets')
df=df_counter %>% mutate(ratio=counts/all)
df$patterns=c('1 target','2 targets','3 targets','>3 targets')
df$pielabels <- paste(df$patterns, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")

levels(df$pielabels)=levels(df$patterns)
col_palette=c(wes_palette("BottleRocket2")[3],
              wes_palette("Darjeeling1")[2],
              'orange',wes_palette("IsleofDogs1")[1])

p=ggpie(
  df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
  lab.font =  c(6, "bold", "black"), fill = "patterns", color = "black", palette = col_palette
)+ theme(text=element_text(face='bold',size=10))

p
ggsave('./results/fig3_pie_targets_ratio.pdf')

### read in single neuron projectome data observation ###

compare_with_projectome_nn <- read_excel("./results/count_observed_expected_vs_nn.xls")
compare_with_projectome_nn
### transform to long format ###
compare_with_projectome_nn <- compare_with_projectome_nn %>%
  gather(key = "Variable", value = "Value", -pattern)
colnames(compare_with_projectome_nn) <- c("pattern", "study", "count")
## calculate the percentage of each category ##
compare_with_projectome_nn <- compare_with_projectome_nn %>% group_by(study) %>%
  mutate(percentage = count / sum(count) * 100)
compare_with_projectome_nn
### remove rows whose "precentage" by "pattern" in "merge_seq" or "nn_projectome" is 0 ###
compare_with_projectome_nn_filter <- compare_with_projectome_nn %>%
  group_by(pattern) %>%
  filter(!all(percentage == 0))
### reorder factor levels of "study" ###
compare_with_projectome_nn_filter$study <- factor(compare_with_projectome_nn_filter$study,levels = c("This study", "Gao et al., 2020"))

ggplot(compare_with_projectome_nn_filter, aes(x = fct_rev(fct_reorder(pattern, percentage)), y = percentage, fill = study)) +
  geom_bar(stat = "identity",position = "dodge2") +
  xlab("Projection pattern") +
  ylab("Percentage") +
  ggtitle("Percentage of projection pattern") +
  scale_fill_nejm() + 
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave('./results/fig3_barplot_nn_projectome.pdf',width = 13,height = 7)

### calculate the ratio by combining single projection, double projection, triple, and 3+ projection ###
compare_with_projectome_nn <- read_excel("./results/count_observed_expected_vs_nn.xls")
compare_with_projectome_nn$meta_pattern <- NA  # Create a new column with NA values
compare_with_projectome_nn[1:5, "meta_pattern"] <- "single"
compare_with_projectome_nn[6:15, "meta_pattern"] <- "double"
compare_with_projectome_nn[16:25, "meta_pattern"] <- "triple"
compare_with_projectome_nn[26:31, "meta_pattern"] <- "triple_more"
colnames(compare_with_projectome_nn) <- c("pattern", "this_study", "gao2020","meta_pattern")
compare_with_projectome_nn <- compare_with_projectome_nn[,-1]
### transform to long format ###
compare_with_projectome_nn <- compare_with_projectome_nn %>%
  gather(key = "Variable", value = "Value", -meta_pattern)
colnames(compare_with_projectome_nn) <- c("meta_pattern", "study", "count")
## calculate the percentage of each category ##
compare_with_projectome_nn <- compare_with_projectome_nn %>% group_by(study) %>%
  mutate(percentage = count / sum(count) * 100)
compare_with_projectome_nn
#write.csv(compare_with_projectome_nn,file='./results/compare_with_projectome_nn_metapattern.csv')
stat.test <- compare_with_projectome_nn %>% 
  group_by(meta_pattern) %>%
  rstatix::wilcox_test(percentage ~ study, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test,file='./results/fig3_meta_pattern_comparision_source data.csv')
### reorder factor levels of "study" ###
compare_with_projectome_nn$meta_pattern <- factor(compare_with_projectome_nn$meta_pattern,levels = c("single", "double", "triple", "triple_more"))
compare_with_projectome_nn$study <- factor(compare_with_projectome_nn$study,levels = c("this_study", "gao2020"))
ggplot(compare_with_projectome_nn, aes(x = meta_pattern, y = percentage, color = study)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, size = 0.5) +
  labs(title = " ",x = "",y = "percentage",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size=0.1) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "meta_pattern", label = "p.adj = {p.adj} ", y.position = 18) + 
  ylim(NA, 18) +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig3_meta_pattern_comparision.pdf", width = 4, height = 5)



### binary cluster in seurat object
pfc.bar <- subset(exn.sub,cells = unique_valid_cells)
temp.bar <- bar.mat[colnames(pfc.bar),]
colnames(temp.bar) <- c("bAI", "bDMS", "bMD",'bBLA','bLH')
# Transform table
# Add barcode data as a new assay independent from RNA
pfc.bar[["bar"]] <- CreateAssayObject(counts =as.data.frame(t(temp.bar)))
DefaultAssay(object = pfc.bar) <- "bar"

pfc.bar <- NormalizeData(object = pfc.bar, normalization.method = "CLR")
pfc.bar <- FindVariableFeatures(object = pfc.bar, selection.method = "vst", nfeatures = 5)
pfc.bar <- ScaleData(object = pfc.bar, vars.to.regress = c("nCount_RNA"))
pfc.bar <- RunPCA(object = pfc.bar, nfeatures.print = 10,npcs = 4)


AI <- WhichCells(pfc.bar, expression = bAI >0 & bDMS == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
BLA <- WhichCells(pfc.bar, expression = bBLA > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bLH == 0,
               slot='counts')
DMS <- WhichCells(pfc.bar, expression = bDMS > 0 & bAI == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
               slot='counts')
MD <- WhichCells(pfc.bar, expression = bMD > 0 & bAI == 0 & bDMS == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
LH <- WhichCells(pfc.bar, expression = bLH > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bBLA == 0,
              slot='counts')

AIDMS <- WhichCells(pfc.bar, expression = bAI > 0 & bDMS > 0 & bMD == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')
DMSMD <- WhichCells(pfc.bar, expression = bDMS > 0 & bMD > 0 & bAI == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')
DMSLH <- WhichCells(pfc.bar, expression = bDMS > 0 & bLH > 0 & bAI == 0 & bMD == 0 & bBLA == 0,
                 slot='counts')

MDAIDMS <- WhichCells(pfc.bar, expression = bAI > 0 & bDMS > 0 & bMD > 0 & bBLA == 0 & bLH == 0,
                 slot='counts')

LHMDAIDMS <- WhichCells(pfc.bar, expression = bDMS > 0 & bMD > 0 & bAI > 0 & bBLA == 0 & bLH > 0,
                 slot='counts')

pfc.bar <- SetIdent(object = pfc.bar, cells = colnames(pfc.bar), value = 'Other')
binary_cluster <- c('AI','DMS','MD','BLA','LH','DMS+AI','DMS+MD','DMS+LH','DMS+AI+MD','DMS+AI+MD+LH')
binary_cluster_cell_list <- list(AI,DMS,MD,BLA,LH,AIDMS,DMSMD,DMSLH,MDAIDMS,LHMDAIDMS)
### set Idents of pfc.bar by binary_cluster_cell_list ###
for(i in 1:length(binary_cluster)) {
  pfc.bar <- SetIdent(object = pfc.bar, cells = binary_cluster_cell_list[[i]], value = binary_cluster[i])
}
table(Idents(pfc.bar))
save(pfc.bar,file='./results/binary_seurat_exn.Robj')

pfc.bar <- subset(pfc.bar,idents = 'Other',invert = T)

load('./results/binary_seurat_exn.Robj')


pca_embed <- Embeddings(object = pfc.bar[["pca"]])
write.csv(pca_embed,file='./results/fig3_pca_embed.csv')

col_pal <- c("#FAD510","#ED0000FF","#42B540FF","#0099B4FF","#9986A5",
          "#FDAF91FF","#CCBA72","#ADB6B6FF","#5050FFFF",'#ABDDDE') #'grey20'

pdf('./results/figs3_pca plot of barcode by binary cluster only top10.pdf',width = 8)
DimPlot(pfc.bar, reduction = "pca",
        pt.size = 2,group.by = 'ident',label = T,label.size = 8,repel = T
) + scale_colour_manual(values=col_pal)
dev.off()


plot_list <-vector("list", 6)
for(i in 1:5) {
  plot_list[[i]] <- FeaturePlot(pfc.bar, features = rownames(pfc.bar)[i],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) & 
    theme(plot.title = element_text(size = 20,vjust = 1),legend.position = c(0.1,0.2))
  #plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}
plot_list[[6]]=DimPlot(pfc.bar, reduction = "pca",
                       pt.size = 1,group.by = 'ident',label = T,label.size = 5,repel = T) + 
  scale_colour_manual(values=col_pal)+ labs(title = "Projection Cluster") & NoLegend()

#plot_list[[6]]<- AugmentPlot(plot_list[[6]],dpi = 400)
pl <- cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=3, align="h")
cowplot::save_plot("./results/barcoded pc black only top10.pdf", pl, base_width=15,base_height = 10)



### complexheatmap binary projection ###
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
mda <- pfc.bar
temp <- exn_meta[colnames(pfc.bar),]
temp$binary <- Idents(pfc.bar)
head(temp)
df <- as.matrix(pfc.bar@assays$bar@data)
rownames(df) <- c("AI",'DMS','MD','BLA','LH')
mda@meta.data$cluster <- as.factor(temp$cluster)
pheno <- data.frame(Idents(mda),mda@meta.data$cluster)
colnames(pheno) <- c('projection_cluster','transcription_cluster')
rownames(pheno) <- colnames(mda)
sub_exn_palette <- c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')
ann_colors <- list(
  'projection_cluster'=col_pal,
  'transcription_cluster'=sub_exn_palette
)
names(ann_colors$'projection_cluster') <- levels(pheno$'projection_cluster')
names(ann_colors$'transcription_cluster') <- levels(pheno$'transcription_cluster')

dend <- hclust(dist(df))
dend <- color_branches(dend, k = 3)
col <- circlize::colorRamp2(c(0,1,2,3,4),viridis::magma(5))

anno <- HeatmapAnnotation(df=pheno,col=ann_colors,
                         annotation_legend_param = list(
                           'projection_cluster' = list(title='Projection',
                                                       legend_direction='horizontal',
                                                       nrow = 3,
                                                       title_gp = gpar(fontsize = 14)),
                           'transcription_cluster' = list(title='Transcription',
                                                          legend_direction='horizontal',
                                                          nrow = 3,
                                                          title_gp = gpar(fontsize = 14))))

p <- Heatmap(df,name = "Normalized Counts", 
          show_column_names=F,
          cluster_rows =F,row_title = NULL,
          #cluster_columns =hclust(dist(t(df))),
          column_split = 3,
          column_title = NULL,
          heatmap_legend_param = list(legend_direction = "horizontal"),
          top_annotation = anno,
          #right_annotation=rowAnnotation(df=pheno_row,col=ann_row_colors),
          col = col
)
pdf('./results/fig3 top 10 complexheatmap binary projection magma.pdf',
    width=15,height=6)
draw(p, heatmap_legend_side = "bottom", merge_legend = TRUE,
     annotation_legend_side = "bottom")
dev.off()

### alluvial plot ###
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
exn_barcode_meta <- exn_meta[colnames(pfc.bar),]
exn_barcode_meta <- exn_barcode_meta %>% select(layer_v1,cluster)
colnames(exn_barcode_meta) <- c('layer','cluster')

meta <- exn_barcode_meta %>% add_column(Binary_Projection = Idents(pfc.bar)) %>% 
  rownames_to_column() %>% as_tibble() %>%
  select(rowname,layer,cluster,Binary_Projection) %>% mutate_if(is.integer, as.factor) %>%
  group_by(layer,cluster,Binary_Projection) %>% srvyr::summarise(., count = n())
meta
meta$cluster <- factor(meta$cluster,levels=c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6', 
                                          'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6'))

ggplot(as.data.frame(meta),
       aes(axis1 = cluster, axis2 = Binary_Projection,
           y = count)) +
  scale_x_discrete(limits = c("Layer", "Projection"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = cluster),alpha=0.8,width =  0.3) +
  scale_fill_manual(values=c('#7fc97f', '#beaed4', '#fdc086', 
                             '#386cb0', '#f0027f', '#bf5b17', '#666666'))+
  #geom_flow()+
  geom_stratum(width = 0.3) + geom_text(stat = "stratum", infer.label = TRUE) +
  theme_void()+
  theme(legend.title=element_text(face='bold',size=10),
        legend.text=element_text(face = 'bold',size=15))+
  theme(legend.position = "bottom") +
  #theme_Publication() + 
  ggtitle("Binary projection pattern distribution in layer/subtype")
write.csv(meta,file='./results/fig3_source data file.csv')
ggsave('./results/fig3 alluvial plot cluster.pdf',height = 8,width = 6)


### pie plot ###
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
exn_barcode_meta <- exn_meta[colnames(pfc.bar),]
exn_barcode_meta <- exn_barcode_meta %>% select(layer_v1,cluster)
colnames(exn_barcode_meta) <- c('layer','cluster')
temp <- exn_barcode_meta %>% add_column(binary = Idents(pfc.bar)) %>% 
  rownames_to_column() %>% as_tibble()
temp$palette <- 'black'

temp$palette[temp$cluster == "L2/3-Calb1"] <- '#7fc97f'
temp$palette[temp$cluster == "L2/3-Rorb"] <- '#beaed4'
temp$palette[temp$cluster == "L5-Bcl6"] <- '#fdc086'
temp$palette[temp$cluster == "L5-Htr2c"] <- '#386cb0'
temp$palette[temp$cluster == "L5-S100b"] <- '#f0027f'
temp$palette[temp$cluster == "L6-Npy"] <- '#bf5b17'
temp$palette[temp$cluster == "L6-Syt6"] <- '#666666'

temp <- temp %>% mutate(palette = factor(palette, 
                      levels = c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                                 '#666666')))
temp <- temp %>% mutate(cluster = factor(cluster, levels = c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6',
                                                 'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6')))
write.csv(temp,file='./results/fig3G_source data file.csv')
diprojection_names <- c('DMS+AI','DMS+MD','DMS+LH')
plot_list <- list()
for (i in 1:3) {
  df= temp %>% select(binary,cluster) %>% 
    filter(binary==diprojection_names[[i]]) %>% 
    group_by(cluster) %>% srvyr::summarise(., count = n()) %>% 
    mutate(ratio=count/sum(count))
  df$pielabels <- paste(df$cluster, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")
  levels(df$pielabels)=levels(droplevels(df$cluster))
  palette_temp=temp %>% select(binary,palette) %>% filter(binary==diprojection_names[[i]])
  palette_temp$palette=factor(palette_temp$palette)
  exn_palette=levels(droplevels(palette_temp$palette))
  #print(exn_palette)
  p=ggpie(
    df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
    lab.font =  c(6, "bold", "black"), fill = "cluster", color = "black", palette = exn_palette
  )+ theme(text=element_text(face='bold',size=20))
  p=annotate_figure(p,top = text_grob(sprintf('%s',diprojection_names[i]), face = "bold", size = 20))
  pdf(sprintf('diprojection_piechart_cluster_distribution_%s.pdf',diprojection_names[i]))
  print(p)
  dev.off()
  plot_list[[i]]=p
}
pg <- cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=3, align="h")
cowplot::save_plot("./results/fig3_diprojection_piechart_plotlist.pdf", pg, base_width=9, base_height = 10)

#### dedicated projections ####
dediprojection_names <- c('AI','DMS','BLA','MD','LH')
plot_list = list()
for (i in 1:5) {
  df= temp %>% select(binary,cluster) %>% 
    filter(binary==dediprojection_names[[i]]) %>% 
    group_by(cluster) %>% dplyr::summarise(., count = n()) %>% 
    dplyr::mutate(ratio=count/sum(count))
  df$pielabels <- paste(df$cluster, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")
  levels(df$pielabels)=levels(droplevels(df$cluster))
  palette_temp=temp %>% select(binary,palette) %>% filter(binary==dediprojection_names[[i]])
  palette_temp$palette=factor(palette_temp$palette)
  exn_palette=levels(droplevels(palette_temp$palette))
  #print(exn_palette)
  p=ggpie(
    df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
    lab.font =  c(6, "bold", "black"), fill = "cluster", color = "black", palette = exn_palette
  )+ theme(text=element_text(face='bold',size=20))
  p=annotate_figure(p,top = text_grob(sprintf('%s',dediprojection_names[i]), face = "bold", size = 20))
  pdf(sprintf('dediprojection_piechart_cluster_distribution_%s.pdf',dediprojection_names[i]))
  print(p)
  dev.off()
  plot_list[[i]]=p
}

pg=cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=5, align="h")
cowplot::save_plot("./results/fig3_dediprojection_piechart_plotlist.pdf", pg,base_width=15,base_height = 10)


#### >= 3 projections ####
more3_projection_names <- c('DMS+AI+MD','DMS+AI+MD+LH')
plot_list <- list()
for (i in 1:2) {
  df= temp %>% select(binary,cluster) %>% 
    filter(binary==more3_projection_names[[i]]) %>% 
    group_by(cluster) %>% dplyr::summarise(., count = n()) %>% 
    dplyr::mutate(ratio=count/sum(count))
  df$pielabels <- paste(df$cluster, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")
  levels(df$pielabels)=levels(droplevels(df$cluster))
  palette_temp=temp %>% select(binary,palette) %>% filter(binary==more3_projection_names[[i]])
  palette_temp$palette=factor(palette_temp$palette)
  exn_palette=levels(droplevels(palette_temp$palette))
  #print(exn_palette)
  p=ggpie(
    df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
    lab.font =  c(6, "bold", "black"), fill = "cluster", color = "black", palette = exn_palette
  )+ theme(text=element_text(face='bold',size=20))
  p=annotate_figure(p,top = text_grob(sprintf('%s',more3_projection_names[i]), face = "bold", size = 20))
  pdf(sprintf('dediprojection_piechart_cluster_distribution_%s.pdf',more3_projection_names[i]))
  print(p)
  dev.off()
  plot_list[[i]]=p
}

pg <- cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=2, align="h")
cowplot::save_plot("./results/fig3_more3_projection_piechart_plotlist.pdf", pg,base_width=6, base_height = 10)


### convert to h5ad ###
library(SeuratDisk)
load('./data_file/sub_trim_seurat_pipe.Robj')
exn.sub <- mda; rm(mda)
exn.sub <- UpdateSeuratObject(exn.sub)
load('./results/binary_seurat_exn.Robj')
pfc.projection <- subset(exn.sub,cells = colnames(pfc.bar))
### keep genes that are only expressed at least in 3 cells ###
seurat_filter_genes <- function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
pfc.projection <- seurat_filter_genes(pfc.projection)
Idents(pfc.projection) <- Idents(pfc.bar)
new_names <- c('DMS+AI+MD+LH','DMS+AI+MD', 'DMS+LH','DMS+MD', 'DMS+AI', 'LH', 'BLA','MD','DMS','AI')
names(new_names) <- levels(pfc.projection)
pfc.projection <- RenameIdents(object = pfc.projection, new_names)
pfc.projection@meta.data$binary <- Idents(pfc.projection)
barcode_genes <- c('barcode0','barcode1','barcode2','barcode3','barcode4')
genes.use <- setdiff(rownames(pfc.projection),barcode_genes)
pfc.projection <- subset(pfc.projection, features = genes.use)
pfc_binary_projection_markers <- FindAllMarkers(pfc.projection,logfc.threshold = 0.25,
               test.use = "MAST", min.pct = 0.1,
               verbose = TRUE, only.pos = T)
pfc_binary_projection_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10$gene,file='./results/fig5_top10_binary_projection_markers.csv')
top10_degs <- read.csv('./results/fig5_top10_binary_projection_markers.csv',header = T)
SaveH5Seurat(pfc.projection, filename = "./results/pfc.projection.h5Seurat")
Convert("./results/pfc.projection.h5Seurat", dest = "h5ad")










