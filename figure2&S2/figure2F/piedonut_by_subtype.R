

setwd('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue')

library(moonBook)
library(tidyverse)
library(webr)
library(ggforce)

#exn_meta=read.csv('exn_meta_valid.csv',row.names = 1)

sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
barcode_names=c('AI','DMS','MD','BLA','LH')
sub_clusters=c('L2/3-Calb1','L2/3-Rorb','L5-Bcl6','L5-Htr2c','L5-S100b','L6-Npy','L6-Syt6')

df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
  filter(cluster=='L2/3-Calb1')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_Calb1_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()

df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
    filter(cluster=='L2/3-Rorb')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_Rorb_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()


df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
  filter(cluster=='L5-Bcl6')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_Bcl6_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()

df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
  filter(cluster=='L5-S100b')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_S100b_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()

df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
  filter(cluster=='L6-Npy')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_Npy_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()

df=exn_meta %>% as_tibble() %>% 
  select(all_of(valid_cells_names),cluster) %>% 
  filter(cluster=='L6-Syt6')
df_ai=as.data.frame(table(df$AI_valid))
df_dms=as.data.frame(table(df$DMS_valid))
df_md=as.data.frame(table(df$MD_valid))
df_bla=as.data.frame(table(df$BLA_valid))
df_lh=as.data.frame(table(df$LH_valid))
df_all=rbind(df_ai[1,],df_dms[1,],df_md[1,],df_bla[1,],df_lh[1,])
pdf('piedonut_Syt6_distribution.pdf')
PieDonut(df_all,aes(Var1,count=Freq),r0=0.7,explode=c(1,2,3,4,5),
         showRatioThreshold = 0.01,start=3*pi/2)
dev.off()













  
