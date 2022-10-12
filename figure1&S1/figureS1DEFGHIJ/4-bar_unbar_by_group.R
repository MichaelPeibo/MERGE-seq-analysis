

setwd('/Users/peiboxu/Desktop/merge-seq analysis/')

library(tidyverse)
library(Seurat)
library(ggpubr)

bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)
AI_valid=rownames(bar.mat[bar.mat$AI>0,])
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])
MD_valid=rownames(bar.mat[bar.mat$MD>0,])
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])
LH_valid=rownames(bar.mat[bar.mat$LH>0,])

load('pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all=mda
rm(mda)
pfc_all_meta=read.csv('pfc_seurat_merge_meta.csv',header=T,row.names=1)

pfc_all=subset(pfc_all,cells=rownames(pfc_all_meta))
unique_valid_cells=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))
pfc_all <- SetIdent(pfc_all, cells = colnames(pfc_all), value = 'Unbarcoded')
pfc_all <- SetIdent(pfc_all, cells = unique_valid_cells, value = 'Barcoded')
pfc_all_meta$barcoded=Idents(pfc_all)
#write.csv(pfc_all_meta,file='pfc_all_meta_barcoded.csv')

### ### ### ### ### ### 
### start from here ###
### ### ### ### ### ### 
pfc_all_meta=read.csv('pfc_all_meta_barcoded.csv',row.names = 1)
temp = pfc_all_meta %>% select(facs,barcoded)
df=as.data.frame(table(temp$facs,temp$barcoded))
df=df %>% as_tibble() %>% group_by(Var1) %>% 
  mutate(ratio=Freq/sum(Freq)) %>% mutate(ratio=round(ratio,2))
p=ggbarplot(df, x = "Var1", y = "ratio",
            color = "Var2", fill = "Var2",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",
            repel=T,lab.size = 7,
            xlab ="", 
            ylab = 'Ratio', ggtheme=theme_pubclean()) + 
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p
ggsave('barcoded cell distribution in sort and unsort.pdf',height = 10,width = 10)
pdf('barcoded cell distribution.pdf',width = 10)
p
dev.off()




