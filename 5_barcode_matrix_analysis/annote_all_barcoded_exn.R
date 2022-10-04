

load('sub_trim_seurat_pipe.Robj')
exn.sub=mda
bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)
AI_valid=rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid=rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid=rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
unique_valid_cells=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))


valid_cells_list=vector('list',5)
valid_cells_list[[1]]=AI_valid
valid_cells_list[[2]]=DMS_valid
valid_cells_list[[3]]=MD_valid
valid_cells_list[[4]]=BLA_valid
valid_cells_list[[5]]=LH_valid
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
barcode_names=c('AI','DMS','MD','BLA','LH')
plot_list = list()
exn_meta=read.csv('sub_exn_trim_meta.csv',row.names = 1)

for (i in 1:5) {
  mda <- SetIdent(mda, cells = colnames(mda), value = 'Others')
  mda <- SetIdent(mda, cells = valid_cells_list[[i]], value =valid_cells_names[i])
  temp=as.data.frame(Idents(mda))
  colnames(temp)=valid_cells_names[i]
  exn_meta=cbind(exn_meta,temp)
}
write.csv(exn_meta,file='exn_meta_valid.csv')


sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')

source('PieDonut_modified_v1.R')
for (i in 1:5) {
  temp=exn_meta %>% as_tibble() %>% 
    select(valid_cells_names[i],cluster) %>% 
    filter(.data[[valid_cells_names[[i]]]]==valid_cells_names[[i]])
  df = as.data.frame(table(temp$cluster))
  df$Var1=factor(df$Var1,levels=c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6', 
                                  'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6'))
  p=PieDonut_modified_v1(df,aes(Var1,count=Freq),center_label=barcode_names[i],which_target=i,
                         r0=0.7,explode = c(1,2,3,4,5,6), explodeDonut=TRUE,start=3*pi/2,
                         showRatioThreshold=0.05)
  ggsave(sprintf('piedonut_layer_distribution_%s_cluster_v1.pdf',barcode_names[i]))
}












