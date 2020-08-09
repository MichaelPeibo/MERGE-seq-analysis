

setwd('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue')

library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(ggsci)
library(wesanderson)
library(Seurat)
library(UpSetR)
library(ggrepel)
library(pheatmap)
library(ggforce)
library(grid)
library(webr)
library(moonBook)
library(EnhancedVolcano)
library(ggalluvial)
library(ComplexHeatmap)


### 1 ###

### read in raw barcode matrix
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
pfc1=read.csv('bar.table_pfc_1.csv',row.names = 1)
pfc2=read.csv('bar.table_pfc_2.csv',row.names = 1)
pfc3=read.csv('bar.table_pfc_3.csv',row.names = 1)
pfc4=read.csv('bar.table_pfc_4.csv',row.names = 1)

matrix_list=list(pfc1,pfc2,pfc3,pfc4)
matrix_list_trim <-vector("list", 4) 
for (i in 1:4){
  temp=matrix_list[[i]]
  temp=temp[,1:5]  
  rownames(temp) <- paste(sprintf("pfc_%s",i), rownames(temp), sep = "_")
  colnames(temp)=barcode_genes
  matrix_list_trim[[i]]=temp
}
bar.mat=rbind(matrix_list_trim[[1]],matrix_list_trim[[2]],matrix_list_trim[[3]],matrix_list_trim[[4]])

### 2 ###
### using ecdf to define each barcodes count threshold

### read in meta info
non.neuron.cells=read.csv('non.neuron.cells.csv',row.names = 1)
non.neuron.bar.mat=bar.mat[rownames(non.neuron.cells),]
neuron.bar.mat=read.csv('bar.table_pfc_4.csv',row.names = 1)

plotlist=vector('list',5)
## y=0.95
## for AI, threshold=24 
## for DMS, threshold=87 
## for MD, threshold=74  
## for BLA, threshold=9 
## for LH, threshold=53 
i=5
df = data.frame(
  celltype = factor(c(rep("EGFP+",length(rownames(neuron.bar.mat))), 
                      rep("Non-neuron",length(rownames(non.neuron.bar.mat))))),
  umi = c(neuron.bar.mat[,i],non.neuron.bar.mat[,i]))
p=ggplot(df, aes(umi,colour = celltype)) +
  scale_color_manual(values = wes_palette("BottleRocket1"))+
  stat_ecdf(geom = "step",pad = FALSE,size = 1) +
  scale_y_continuous(breaks=c(0,0.25,0.50,0.75,0.95,1))+
  scale_x_continuous(breaks=c(0,25,75,100,53),limits = c(0, 100))+
  geom_vline(xintercept=53, linetype="dashed", color = "black")+
  geom_hline(yintercept=0.95, linetype="dashed", color = "black")+
  labs(title="UMI empirical cumulative distribution of LH",
       x ="UMI", y = "Cumulative density")+
  theme_pubr() + labs_pubr()
p
plotlist[[i]]=p
plotlist
saveRDS(plotlist,file ='ecdf0.95_plotlist.rds')
# plotlist=readRDS('ecdf0.95_plotlist_modified.rds')
rw=5
pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=3, align="h",rel_widths=c(rw,rw,rw,rw,rw))

cowplot::save_plot("ecdf0.95_plotlist.pdf", pg,base_width=15,base_height = 10)

### 3 ###
### drop negative barcodes counts to zero
bar.mat=rbind(matrix_list_trim[[1]],matrix_list_trim[[2]],matrix_list_trim[[3]],matrix_list_trim[[4]])
colnames(bar.mat)=c('AI','DMS','MD','BLA','LH')
temp=bar.mat
temp = temp %>% mutate(AInew = case_when(AI <= 24 ~ 0,TRUE   ~ as.numeric(AI)),
                       DMSnew=case_when(DMS<=87 ~ 0, TRUE ~ as.numeric(DMS)),
                       MDnew=case_when(MD<=74 ~ 0,TRUE ~ as.numeric(MD)),
                       BLAnew=case_when(BLA<=9 ~ 0, TRUE ~ as.numeric(BLA)),
                       LHnew=case_when(LH<=53 ~ 0, TRUE ~ as.numeric(LH)))
head(temp,n=10)
rownames(temp)=rownames(bar.mat)
temp = temp[,6:10]
colnames(temp)=c('AI','DMS','MD','BLA','LH')

bar.mat.drop=temp
write.csv(bar.mat.drop,file='bar.mat.drop.csv')

df=gather(temp,key = 'barcodes',value='UMI')
head(df)

#pdf('test beeswarm violin log10 drop negative cells.pdf')
p=ggviolin(df,x='barcodes',y='UMI',outlier.shape = NA,
           xlab ="Projection targets", 
           ylab = parse(text = paste0(' ~ log[10] ~','Counts'))) + 
  yscale("log10", .format = TRUE) + geom_quasirandom(size=0.1) + 
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p
ggsave("beeswarm violin log10 drop negative cells.png", 
       width = 8, height = 8, dpi=1000)

### 4 ###
# some cells are filtered in seurat
# rownames(pfc1) <- paste("pfc_1", rownames(pfc1), sep = "_")
# unique_valid_cells[!unique_valid_cells %in% colnames(pfc_all)]
# 'pfc_1_AGTTCCCCAGCATTGT' %in% rownames(pfc1)
# 'pfc_1_AGTTCCCCAGCATTGT' %in% colnames(mda)

### barcoded cell distribution in all celltypes
load('./1_seurat_processing_all_cells/pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all=mda
rm(mda)
pfc_all_meta=read.csv('./2_scanpy_processing_all_cells/pfc_seurat_merge_meta.csv',header=T,row.names=1)
pfc_all=subset(pfc_all,cells=rownames(pfc_all_meta))
unique_valid_cells=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))
pfc_all <- SetIdent(pfc_all, cells = colnames(pfc_all), value = 'Unbarcoded')
pfc_all <- SetIdent(pfc_all, cells = unique_valid_cells, value = 'Barcoded')
pfc_all_meta$barcoded=Idents(pfc_all)
write.csv(pfc_all_meta,file='pfc_all_meta_barcoded.csv')
df = pfc_all_meta %>% select(facs,leiden_coarse,barcoded)
df=df %>% 
  group_by(leiden_coarse,barcoded) %>% count(name = "CellCount") %>% ungroup() %>%
  mutate(leiden_coarse = factor(leiden_coarse, 
                                levels = c('Excitatory','Microglia', 'Endo',
                                           'OPC','Oligo','Inhibitory',
                                           'Astro','Act.Microglia')))
df
# Create stacked bar graphs with labels
p=ggbarplot(df, x = "leiden_coarse", y = "CellCount",
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
ggsave('barcoded cell distribution.pdf',height = 10,width = 10)
pdf('barcoded cell distribution.pdf',width = 10)
p
dev.off()
###--------------
### layer distribution ###
exn_meta=read.csv('./3_scanpy_processing_exn/sub_exn_trim_meta.csv',row.names = 1,header = T)
exn_barcode_meta=pfc_all_meta[rownames(exn_meta),]
exn_barcode_meta=cbind(exn_barcode_meta,exn_meta$layer_v1,exn_meta$cluster)
colnames(exn_barcode_meta)[15:16]=c('layer','cluster')
df = exn_barcode_meta %>% as_tibble() %>% select(layer,cluster,barcoded)
df=df %>% 
  group_by(layer,barcoded) %>% count(name = "CellCount") %>% ungroup()
df
# Create stacked bar graphs with labels
p=ggbarplot(df, x = "layer", y = "CellCount",
            color = "barcoded", fill = "barcoded",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",
            xlab ="Layer",ylab = 'Cell Count')+theme_pubclean() + labs_pubr()
p

pdf('barcoded layer distribution.pdf')
p
dev.off()
###--------------
### subtype distribution ###
exn_meta=read.csv('./3_scanpy_processing_exn/sub_exn_trim_meta.csv',row.names = 1,header = T)
exn_barcode_meta=pfc_all_meta[rownames(exn_meta),]
exn_barcode_meta=cbind(exn_barcode_meta,exn_meta$layer_v1,exn_meta$cluster)
colnames(exn_barcode_meta)[15:16]=c('layer','cluster')
df = exn_barcode_meta %>% as_tibble() %>% select(layer,cluster,barcoded)
df=df %>% 
  group_by(cluster,barcoded) %>% count(name = "CellCount") %>% ungroup()
df
# Create stacked bar graphs with labels
p=ggbarplot(df, x = "cluster", y = "CellCount",
            color = "barcoded", fill = "barcoded",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",lab.size = 7,
            xlab ="Layer",ylab = 'Cell Count',ggtheme=theme_pubclean()) + 
  rotate_x_text(45)+
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p

pdf('barcoded exn distribution.pdf')
p
dev.off()
###--------------
### barcoded cell distribution in facs condition

df = pfc_all_meta %>% select(facs,barcoded) %>% group_by(facs,barcoded) %>% 
  summarise(., count = n()) %>% mutate(ratio=round(count/sum(count),2))

p=ggbarplot(df, "facs", "ratio",
            fill = "barcoded", color = "barcoded", palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.col = "black", lab.pos = "in",lab.size=6,
            xlab ="",ylab = 'Ratio') +
  theme_pubr() + labs_pubr()
p
pdf('barplot barcoded cell distribution in facs condition.pdf')
p
dev.off()


### 5 ###
### piechart distrubition ###
load('./3_scanpy_processing_exn/sub_trim_seurat_pipe.Robj')
bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)

AI_valid=rownames(bar.mat[bar.mat$AI>0,])
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])
MD_valid=rownames(bar.mat[bar.mat$MD>0,])
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])
LH_valid=rownames(bar.mat[bar.mat$LH>0,])

valid_cells_list=vector('list',5)
valid_cells_list[[1]]=AI_valid
valid_cells_list[[2]]=DMS_valid
valid_cells_list[[3]]=MD_valid
valid_cells_list[[4]]=BLA_valid
valid_cells_list[[5]]=LH_valid
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
barcode_names=c('AI','DMS','MD','BLA','LH')
plot_list = list()

exn_meta=read.csv('./3_scanpy_processing_exn/sub_exn_trim_meta.csv',row.names = 1)
for (i in 1:5) {
  mda <- SetIdent(mda, cells = colnames(mda), value = 'Others')
  mda <- SetIdent(mda, cells = valid_cells_list[[i]], value =valid_cells_names[i])
  temp=as.data.frame(Idents(mda))
  colnames(temp)=valid_cells_names[i]
  exn_meta=cbind(exn_meta,temp)
  df= exn_meta %>% select(valid_cells_names[i],layer_v1) %>% 
    filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
    group_by(layer_v1) %>% summarise(., count = n()) %>% 
    mutate(ratio=count/sum(count))
  print(df)
  df$pielabels <- paste(df$layer_v1, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")
  layer_palette = c('#1f77b4', '#8c564b', '#17becf')
  levels(df$pielabels)=levels(df$layer_v1)
  p=ggpie(
    df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
    lab.font =  c(6, "bold", "black"), fill = "layer_v1", color = "black", palette = layer_palette
  )+ theme(text=element_text(face='bold',size=20))
  p=annotate_figure(p,top = text_grob(sprintf('%s',barcode_names[i]), face = "bold", size = 20))
  pdf(sprintf('piechart_layer_distribution_%s.pdf',barcode_names[i]))
  print(p)
  dev.off()
  plot_list[[i]]=p
}
write.csv(exn_meta,file='exn_meta_valid.csv')
# 
# cowplot::plot_grid(plotlist=plot_list, labels = "auto",
#                    label_size = 20,
#                    ncol=3, align="h",rel_widths=c(rw,rw,rw,rw,rw))
# ggsave('five piechart.pdf',height = 12,width = 15)
# 
# ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
#           plot_list[[5]],ncol = 5, nrow = 1,
#           common.legend = TRUE)
# ggsave('five piechart.pdf',height = 5,width = 20)

sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')

source('PieDonut_modified.R')
for (i in 1:5) {
  df=exn_meta %>% as_tibble() %>% 
    select(valid_cells_names[i],layer_v1,cluster) %>% 
    filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]])
  levels(df$layer_v1)=levels(exn_meta$layer_v1)
  levels(df$cluster)=levels(exn_meta$cluster)
  pdf(sprintf('piedonut_layer_distribution_%s.pdf',barcode_names[i]))
  PieDonut_modified(df,aes(pies=layer_v1,donuts=cluster),center_label=barcode_names[i],
                    selected = 1,labelposition = 0,showRatioThreshold = 0.1,explode = c(1,2,3))  
  dev.off()
}



### ###
### facs sort cell distribution in all cell types
pfc_all_meta=read.csv('./2_scanpy_processing_all_cells/pfc_seurat_merge_meta.csv',header=T,row.names=1)
df=pfc_all_meta %>% 
  group_by(facs,leiden_coarse) %>%  
  summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% ungroup()
df = df %>% mutate(leiden_coarse = factor(leiden_coarse, 
                                          levels = c('Excitatory','Microglia', 'Endo',
                                                     'OPC','Oligo','Inhibitory',
                                                     'Astro','Act.Microglia')))
df=df %>% filter(facs=='sort')
df$pielabels <- paste(df$leiden_coarse, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")
celltype_pelette = c('#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c',
                     '#98df8a','#d62728','#ff9896')
levels(df$pielabels)=levels(df$leiden_coarse)
PieDonut(df,aes(leiden_coarse,count=count),color='white',
         r0=0.7,start=3*pi/2,labelpositionThreshold=0.2,showRatioThreshold=0.0004) + 
  scale_fill_manual(values =celltype_pelette )
ggsave('piedounout.pdf')



### ###
### upset plot of binary cluster

load('./3_scanpy_processing_exn/sub_trim_seurat_pipe.Robj')
exn.sub=mda;rm(mda)
bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)

AI_valid=rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid=rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid=rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]

upset_list=list(AI=AI_valid,DMS=DMS_valid,MD=MD_valid,BLA=BLA_valid,LH=LH_valid)

pdf('upset valid cell.pdf')
UpSetR::upset(fromList(upset_list), 
              order.by = "freq", sets.bar.color ='deep sky blue',
              point.size = 3.5, line.size = 2, 
              nintersects=10,
              mainbar.y.label = "Cell Count",
              text.scale = c(3, 1.3, 1, 1, 2, 1.5),
              queries = list(list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('DMS'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('LH'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('BLA'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI','DMS'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','BLA'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','DMS'), color='orange',active = T))
              
)
dev.off()
pdf('upset valid cell all groups.pdf')
UpSetR::upset(fromList(upset_list), 
              order.by = "freq", sets.bar.color ='deep sky blue',
              point.size = 3.5, line.size = 2, 
              #nintersects=10,
              mainbar.y.label = "Cell Count",
              text.scale = c(3, 1.3, 1, 1, 2, 1.5),
              queries = list(list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('DMS'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('LH'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('BLA'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI','DMS'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','BLA'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','DMS'), color='orange',active = T))
              
)
dev.off()
##
### ###
### binary cluster in seurat object

pfc.bar=subset(exn.sub,cells=unique_valid_cells)

seurat_filter_genes = function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
pfc.bar=subset(pfc.bar,features = setdiff(rownames(pfc.bar),barcode_genes))

pfc.bar=seurat_filter_genes(pfc.bar,min.cells=3)

temp.bar=bar.mat.drop[colnames(pfc.bar),]
#colnames(temp.bar)=c("barcode0", "barcode1", "barcode2",'barcode3','barcode4')
colnames(temp.bar)=c("bAI", "bDMS", "bMD",'bBLA','bLH')

# Transform table
# Add barcode data as a new assay independent from RNA
pfc.bar[["bar"]] <- CreateAssayObject(counts =as.data.frame(t(temp.bar)))
DefaultAssay(object = pfc.bar) <- "bar"

pfc.bar

AI=WhichCells(pfc.bar, expression = bAI >0 & bDMS == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
BLA=WhichCells(pfc.bar, expression = bBLA > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bLH == 0,
               slot='counts')
DMS=WhichCells(pfc.bar, expression = bDMS > 0 & bAI == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
               slot='counts')
MD=WhichCells(pfc.bar, expression = bMD > 0 & bAI == 0 & bDMS == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
LH=WhichCells(pfc.bar, expression = bLH > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bBLA == 0,
              slot='counts')

AIDMS=WhichCells(pfc.bar, expression = bAI > 0 & bDMS > 0 & bMD == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')
MDLH=WhichCells(pfc.bar, expression = bMD > 0 & bLH > 0 & bAI == 0 & bDMS == 0 & bBLA == 0,
                slot='counts')
DMSLH=WhichCells(pfc.bar, expression = bDMS > 0 & bLH > 0 & bAI == 0 & bMD == 0 & bBLA == 0,
                 slot='counts')
DMSBLA=WhichCells(pfc.bar, expression = bDMS > 0 & bBLA > 0 & bAI == 0 & bMD == 0 & bLH == 0,
                  slot='counts')
DMSMD=WhichCells(pfc.bar, expression = bDMS > 0 & bMD > 0 & bAI == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')


pfc.bar = SetIdent(object = pfc.bar, cells = colnames(pfc.bar), value = 'Other')
binary_cluster=c('AI','DMS','MD','BLA','LH','DMS+AI','MD+LH','DMS+LH','DMS+BLA','DMS+MD')
binary_cluster_cell_list=list(AI,DMS,MD,BLA,LH,AIDMS,MDLH,DMSLH,DMSBLA,DMSMD)

for (i in 1:10) {
  pfc.bar <- SetIdent(object = pfc.bar, cells = binary_cluster_cell_list[[i]], value = binary_cluster[i])
}
table(Idents(pfc.bar))

pfc.bar.sub=subset(pfc.bar,idents = 'Other',invert=T)
DefaultAssay(object = pfc.bar.sub) <- "RNA"
table(Idents(pfc.bar.sub))

### ###
### alluvial plot of binary cluster
meta = exn_barcode_meta[colnames(pfc.bar.sub),]
meta = meta %>% add_column(Binary_Projection=Idents(pfc.bar.sub)) %>% 
  rownames_to_column() %>% as_tibble() %>%
  select(rowname,layer,cluster,Binary_Projection) %>% mutate_if(is.integer, as.factor) %>%
  group_by(layer,cluster,Binary_Projection) %>% summarise(., count = n())
meta

ggplot(as.data.frame(meta),
       aes(axis1 = layer, axis2 = Binary_Projection,
           y= count)) +
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
  ggtitle("Binary projection pattern distribution in layer/subtype")
ggsave('alluvial plot.pdf',height = 8,width = 6)
### ###
#### calculate degs 
DefaultAssay(object = pfc.bar.sub) <- "bar"
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
bar_genes=c("bAI", "bDMS", "bMD",'bBLA','bLH')
n=1
top_deg_list_names <- c('AI','DMS','MD','BLA','LH')
top_deg_list <- vector("list", length(top_deg_list_names))
names(top_deg_list) <- top_deg_list_names

for (i in bar_genes){
  DefaultAssay(object = pfc.bar.sub) <- "bar"
  pfc.bar.sub = SetIdent(object = pfc.bar.sub, cells = colnames(pfc.bar.sub), value = 'Other')
  temp_expr <- FetchData(object = pfc.bar.sub, vars = i,slot='counts')
  temp_cells=colnames(pfc.bar.sub[, which(x = temp_expr > 0)])
  pfc.bar.sub <- SetIdent(object = pfc.bar.sub, cells = temp_cells, value = valid_cells_names[n])
  table(Idents(pfc.bar.sub))
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  temp_markers=FindMarkers(pfc.bar.sub,ident.1=valid_cells_names[n],
                           logfc.threshold = 0.25,
                           test.use = "MAST", min.pct = 0.1,
                           verbose = TRUE, only.pos = T)
  write.csv(temp_markers,file=paste(valid_cells_names[n],'.csv',sep = ''))
  top100 = rownames(temp_markers[1:100,])
  top_deg_list[[n]]=top100
  n=n+1
}


top_deg_list_names <- c('AI','DMS','MD','BLA','LH')
diprojection_names <- c('DMS+AI','MD+LH','DMS+LH','DMS+BLA','DMS+MD')
#volcano_deg_list <- vector("list", 10)
for (i in 1:5){
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  for (n in 1:5) {
    
    temp_markers=FindMarkers(pfc.bar.sub,ident.1=top_deg_list_names[i],ident.2=diprojection_names[n],
                             test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1,
                             verbose = TRUE, only.pos = F)
    write.csv(temp_markers,file=paste(top_deg_list_names[i],'vs',diprojection_names[n],'.csv',sep = ''))
    png(sprintf('EnhancedVolcano %s vs %s .png',top_deg_list_names[i],diprojection_names[n]),
        units="px", width=1600, height=1600, res=300)
    plot(EnhancedVolcano(temp_markers,
                         lab = rownames(temp_markers),
                         x = 'avg_logFC',
                         y = 'p_val_adj',
                         #xlim = c(-3, 3),
                         #ylim = c(0, -log10(10e-250)),
                         drawConnectors = TRUE,
                         widthConnectors = 0.5, 
                         colConnectors = 'grey30',
                         title = sprintf('%s vs %s projection',top_deg_list_names[i],diprojection_names[n]),
                         xlab = bquote(~Log[2]~ .(paste('FC(',top_deg_list_names[[i]],'/',diprojection_names[[n]],')',sep = ''))) ,
                         subtitle = NULL,
                         pCutoff = 10e-20,
                         FCcutoff = 1,
                         pointSize = 2,
                         col=c('black', 'black', 'black', '#cc163a'),
                         colAlpha=1,
                         labSize = 4.0)+ theme(legend.position = 'none')
    )
    dev.off()
  }
}

top_deg_list_names <- c('AI','DMS','MD','BLA','LH')

for (i in 1:5){
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  for (n in 1:5) {
    if(i==n){
      next
    }
    temp_markers=FindMarkers(pfc.bar.sub,ident.1=top_deg_list_names[i],ident.2=top_deg_list_names[n],
                             test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1,
                             verbose = TRUE, only.pos = F)
    write.csv(temp_markers,file=paste(top_deg_list_names[i],'vs',top_deg_list_names[n],'.csv',sep = ''))
    png(sprintf('EnhancedVolcano %s vs %s .png',top_deg_list_names[i],top_deg_list_names[n]),
        units="px", width=1600, height=1600, res=300)
    plot(EnhancedVolcano(temp_markers,
                         lab = rownames(temp_markers),
                         x = 'avg_logFC',
                         y = 'p_val_adj',
                         #xlim = c(-3, 3),
                         #ylim = c(0, -log10(10e-250)),
                         drawConnectors = TRUE,
                         widthConnectors = 0.5, 
                         colConnectors = 'grey30',
                         title = sprintf('%s vs %s dedicated projection',top_deg_list_names[i],top_deg_list_names[n]),
                         xlab = bquote(~Log[2]~ .(paste('FC(',top_deg_list_names[[i]],'/',top_deg_list_names[[n]],')',sep = ''))) ,
                         subtitle = NULL,
                         pCutoff = 10e-20,
                         FCcutoff = 1,
                         pointSize = 2,
                         col=c('black', 'black', 'black', '#cc163a'),
                         colAlpha=1,
                         labSize = 4.0)+ theme(legend.position = 'none')
    )
    dev.off()
  }
}
UpSetR::upset(UpSetR::fromList(top_deg_list), order.by = "freq")
unique_degs=unique(c(top_deg_list[[1]],top_deg_list[[2]],top_deg_list[[3]],
                     top_deg_list[[4]],top_deg_list[[5]]))

### ###
### heatmap of degs
bar.sub.mat=as.matrix(pfc.bar.sub@assays$RNA@scale.data)[unique_degs,]
col<- circlize::colorRamp2(c(-1,0,1,2,3),viridis::inferno(5))
bar.sub.mat=t(bar.sub.mat)

Heatmap(bar.sub.mat,name = "Scaled Expression", 
        #show_column_names=F,
        show_row_names=F,
        cluster_rows =hclust(dist(bar.sub.mat)),
        #row_order=as.character(unique(deg_genes)[drop=T]), 
        cluster_columns = hclust(dist(t(bar.sub.mat))),
        #top_annotation = HeatmapAnnotation(df=pheno,col=ann_colors),
        col = col
        #viridis::inferno(20)
        #col = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
)

pfc.bar.sub=subset(pfc.bar,idents = 'Other',invert=T)
DefaultAssay(object = pfc.bar.sub) <- "RNA"
table(Idents(pfc.bar.sub))
pfc.bar.sub <- ScaleData(object = pfc.bar.sub, 
                         features=unique_degs,
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt",'facs'))

binary_cluster_markers=FindAllMarkers(pfc.bar.sub,logfc.threshold = 0.25,
                                      test.use = "MAST", min.pct = 0.1,
                                      verbose = TRUE, only.pos = T)
write.csv(binary_cluster_markers,file='binary_cluster_markers.csv')
# Re-level object@ident
my_levels <- c('AI','DMS','BLA','DMS+AI','DMS+BLA','DMS+MD','DMS+LH','MD','LH','MD+LH')

pfc.bar.sub@active.ident <- factor(x = pfc.bar.sub@active.ident, levels = my_levels)
top10 <- binary_cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10$gene,file='top10_degs.txt', sep = "",row.names=F,col.names=F,quote = FALSE)

DoHeatmap(pfc.bar.sub, disp.min = -1,disp.max=1,
          features = top10$gene) + NoLegend() + viridis::scale_fill_viridis()


bar.sub.mat=as.matrix(pfc.bar.sub@assays$RNA@scale.data)[unique(top10$gene),]
col<- circlize::colorRamp2(c(-1,0,1,2,3),viridis::inferno(5))
#bar.sub.mat=t(bar.sub.mat)

p=Heatmap(bar.sub.mat,name = "Scaled Expression", 
          show_column_names=F,
          show_row_names=T,
          cluster_rows =hclust(dist(bar.sub.mat)),
          #row_order=as.character(unique(deg_genes)[drop=T]), 
          cluster_columns = hclust(dist(t(bar.sub.mat))),
          #top_annotation = HeatmapAnnotation(df=pheno,col=ann_colors),
          col = col
          #viridis::inferno(20)
          #col = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
)
p = draw(p)
row_order(p)
write.table(unique(top10$gene)[row_order(p)],file='top10_degs.txt', sep = "",row.names=F,col.names=F,quote = FALSE)

save(pfc.bar.sub,file='pfc.bar.sub.Robj')
### on cluster convert to loom file
load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/barcode_projection/pfc.bar.sub.Robj')
pfc.bar.sub@meta.data$binary=Idents(pfc.bar.sub)
head(pfc.bar.sub@meta.data)
meta_info=pfc.bar.sub@meta.data
table(is.na(meta_info))
pfc.bar.sub@graphs <- list()
pfc.bar.sub.loom=as.loom(x = pfc.bar.sub,filename='pfc.bar.sub.loom')
pfc.bar.sub.loom$close_all()

DefaultAssay(object = pfc.bar.sub) <- "bar"
pfc.bar.sub <- NormalizeData(object = pfc.bar.sub, normalization.method = "CLR")
pfc.bar.sub <- FindVariableFeatures(object = pfc.bar.sub, selection.method = "vst", nfeatures = 5)

pfc.bar.sub@meta.data$binary=Idents(pfc.bar.sub)
head(pfc.bar.sub@meta.data)
meta_info=pfc.bar.sub@meta.data
table(is.na(meta_info))
pfc.bar.sub@graphs <- list()
pfc.bar.sub.loom=as.loom(x = pfc.bar.sub,filename='pfc.bar.sub.bargenes.loom')
pfc.bar.sub.loom$close_all()

##################



