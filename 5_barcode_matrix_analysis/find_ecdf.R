

### 1 ###
setwd('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue')
setwd('/Users/peiboxu/Desktop/merge-seq analysis/')
library(tidyverse)
library(Seurat)

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
## y=0.95 ##wrong
## for AI, threshold=24
## for DMS, threshold=87 
## for MD, threshold=74  
## for BLA, threshold=9 
## for LH, threshold=53 

## y=0.9
## ai, 16
## dms, 226
## md, 148
## bla, 7
## lh, 71

## y=0.95
## ai, 39
## dms, 374
## md, 257
## bla, 10
## lh, 142

#### according to non-neuron
## y=0.999
## ai, 23
## dms, 71
## md, 134
## bla, 9
## lh, 125

i=1
df = data.frame(
  celltype = factor(c(rep("EGFP+",length(rownames(neuron.bar.mat))), 
                      rep("Non-neuron",length(rownames(non.neuron.bar.mat))))),
  umi = c(neuron.bar.mat[,i],non.neuron.bar.mat[,i]))
egfp=df %>% filter(celltype=='EGFP+')
non_neuron=df %>% filter(celltype=='Non-neuron')
#Recreate ecdf data of non-neuron
dat_ecdf <- 
  data.frame(x=unique(non_neuron$umi),
             y=ecdf(non_neuron$umi)(unique(non_neuron$umi))*length(non_neuron$umi))
#rescale y to 0,1 range
dat_ecdf$y <- 
  scale(dat_ecdf$y,center=min(dat_ecdf$y),scale=diff(range(dat_ecdf$y)))

idx=which(abs(dat_ecdf$y - 0.999) == min(abs(dat_ecdf$y - 0.999)))
dat_ecdf[idx,]

p=ggplot(df, aes(umi,colour = celltype)) +
  scale_color_manual(values = wes_palette("BottleRocket1"))+
  geom_density()+
  #stat_ecdf(geom = "step",pad = FALSE,size = 1) +
  #scale_y_continuous(breaks=c(0,0.25,0.50,0.75,0.95,1))+
  scale_x_log10(breaks=c(0,1,2,3,4,25,50,100))+
  geom_vline(xintercept=125, linetype="dashed", color = "black")+
  #geom_hline(yintercept=0.95, linetype="dashed", color = "black")+
  labs(
       x ="log10 UMI-LH", y = "Density")+
  theme_pubr() + labs_pubr()
p
ggsave('barcoded cell distribution LH.pdf',height = 10,width = 10)











