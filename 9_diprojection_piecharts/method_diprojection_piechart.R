

library(ggpubr)
library(tidyverse)
library(Seurat)
load('./5_barcode_matrix_analysis/pfc.bar.sub.Robj')


pfc.bar.sub
exn_meta=read.csv('./3_scanpy_processing_exn/sub_exn_trim_meta.csv',row.names = 1,header = T)
temp=exn_meta[colnames(pfc.bar.sub),]
temp$binary=Idents(pfc.bar.sub)
head(temp)
temp$palette <- 'black'

temp$palette[temp$cluster == "L2/3-Calb1"] <- '#7fc97f'
temp$palette[temp$cluster == "L2/3-Rorb"] <- '#beaed4'
temp$palette[temp$cluster == "L5-Bcl6"] <- '#fdc086'
temp$palette[temp$cluster == "L5-Htr2c"] <- '#386cb0'
temp$palette[temp$cluster == "L5-S100b"] <- '#f0027f'
temp$palette[temp$cluster == "L6-Npy"] <- '#bf5b17'
temp$palette[temp$cluster == "L6-Syt6"] <- '#666666'

temp=temp %>% mutate(palette = factor(palette, 
                      levels = c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                                 '#666666')))
temp=temp %>% mutate(cluster = factor(cluster, 
                                      levels = c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6',
                                                 'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6')))


diprojection_names=c('DMS+AI','DMS+BLA','DMS+MD','DMS+LH','MD+LH')

plot_list = list()

for (i in 1:5) {
  df= temp %>% select(binary,cluster) %>% 
    filter(binary==diprojection_names[[i]]) %>% 
    group_by(cluster) %>% summarise(., count = n()) %>% 
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


pg=cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=5, align="h")
cowplot::save_plot("diprojection_piechart_plotlist.pdf", pg,base_width=15,base_height = 10)









