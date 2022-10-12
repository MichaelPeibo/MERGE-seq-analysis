
setwd('/Users/peiboxu/Desktop/merge-seq analysis/')

library(tidyverse)
library(ggpubr)


### ### ### ### ### ### 
### start from here ###
### ### ### ### ### ### 

valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
#exn_meta=read.csv('exn_meta_valid.csv',row.names = 1)
exn_meta=exn_meta %>% mutate(cluster = factor(cluster, 
                                               levels = c('L2/3-Calb1','L2/3-Rorb', 
                                                          'L5-Bcl6','L5-Htr2c',
                                                          'L5-S100b','L6-Npy','L6-Syt6')))

mouse = rep(c('unsort_1-3','sort_4-6'), c(7765,1603))
exn_meta$mouse=mouse

sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')
### plot one by one due to color palette unmatches ###
i=1 # AI
df= exn_meta %>% select(valid_cells_names[i],cluster,mouse) %>% 
  filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
  group_by(mouse,cluster) %>% summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
df
p=ggbarplot(df, x = "mouse", y = "ratio",
            color = "cluster", fill = "cluster",
            palette = sub_exn_palette,
            label = TRUE, lab.pos = "in", lab.col = "black",
            repel=T,lab.size = 7,
            xlab ="", 
            ylab = 'Ratio', ggtheme=theme_pubclean()) + 
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20) 
AI_comparison=compare_means(ratio ~ mouse, data = df, 
                            group.by = "cluster")
AI_comparison
write.csv(AI_comparison,file='AI_comparison.csv')
p
ggsave('AI_projection_motif_by_mouse.pdf',height = 5,width = 3)
####
i=2 # DMS
df= exn_meta %>% select(valid_cells_names[i],cluster,mouse) %>% 
  filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
  group_by(mouse,cluster) %>% summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
df
p=ggbarplot(df, x = "mouse", y = "ratio",
            color = "cluster", fill = "cluster",
            palette = sub_exn_palette,
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
DMS_comparison=compare_means(ratio ~ mouse, data = df, 
                            group.by = "cluster")
DMS_comparison
write.csv(DMS_comparison,file='DMS_comparison.csv')

ggsave('DMS_projection_motif_by_mouse.pdf',height = 5,width = 3)

####
i=3 # MD
df= exn_meta %>% select(valid_cells_names[i],cluster,mouse) %>% 
  filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
  group_by(mouse,cluster) %>% summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
df
sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086',  '#f0027f', '#bf5b17',
                  '#666666')
p=ggbarplot(df, x = "mouse", y = "ratio",
            color = "cluster", fill = "cluster",
            palette = sub_exn_palette,
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
MD_comparison=compare_means(ratio ~ mouse, data = df, 
                             group.by = "cluster")
MD_comparison
write.csv(MD_comparison,file='MD_comparison.csv')
ggsave('MD_projection_motif_by_mouse.pdf',height = 5,width = 3)

####
i=4 #BLA
df= exn_meta %>% select(valid_cells_names[i],cluster,mouse) %>% 
  filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
  group_by(mouse,cluster) %>% summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
df
sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086',  '#f0027f', '#bf5b17',
                  '#666666')
p=ggbarplot(df, x = "mouse", y = "ratio",
            color = "cluster", fill = "cluster",
            palette = sub_exn_palette,
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
BLA_comparison=compare_means(ratio ~ mouse, data = df, 
                            group.by = "cluster")
BLA_comparison
write.csv(BLA_comparison,file='BLA_comparison.csv')
ggsave('BLA_projection_motif_by_mouse.pdf',height = 5,width = 3)

####
i=5 #LH
df= exn_meta %>% select(valid_cells_names[i],cluster,mouse) %>% 
  filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]]) %>% 
  group_by(mouse,cluster) %>% summarise(., count = n()) %>% 
  mutate(ratio=count/sum(count)) %>% mutate(ratio=round(ratio,3))
df
sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#f0027f', '#bf5b17',
                  '#666666')
p=ggbarplot(df, x = "mouse", y = "ratio",
            color = "cluster", fill = "cluster",
            palette = sub_exn_palette,
            label = TRUE, lab.pos = "in", lab.col = "black",
            repel=T,lab.size = 7,
            xlab ="", 
            ylab = 'Ratio', ggtheme=theme_pubclean()) +
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)

LH_comparison=compare_means(ratio ~ mouse, data = df, 
                             group.by = "cluster")
LH_comparison

write.csv(LH_comparison,file='LH_comparison.csv')
p
ggsave('LH_projection_motif_by_mouse.pdf',height = 5,width = 3)










