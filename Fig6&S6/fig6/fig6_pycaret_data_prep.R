

setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")
library(Seurat)
library(dplyr)


load('pfc.bar.sub.Robj')

seurat_filter_genes = function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
pfc.bar.sub=seurat_filter_genes(pfc.bar.sub,min.cells=3)

exn_meta_valid=read.csv('exn_meta_valid.csv',row.names=1)
meta=exn_meta_valid[colnames(pfc.bar.sub),]
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')

### test DMS ###
Idents(pfc.bar.sub)=meta$DMS_valid
table(Idents(pfc.bar.sub))
dms_markers=FindAllMarkers(pfc.bar.sub,logfc.threshold = 0.25,
                                      min.pct = 0.1,
                                      verbose = TRUE, only.pos = T)
dms_markers_filter= dms_markers %>% filter(p_val_adj<0.000001 & avg_log2FC>0.25)

####  try hvg 150 ############
num_vars=c(50,100,150,300,500,800,1000)
set.seed(200)
for (n in num_vars) {
  pfc.bar.sub=FindVariableFeatures(object = pfc.bar.sub, selection.method = "vst", nfeatures = n)
  var.genes=VariableFeatures(pfc.bar.sub)
  print(table(var.genes=='Syt6'))
  print(table(var.genes=='Foxp2'))
  
  dm=as.matrix(pfc.bar.sub@assays$RNA@data)[var.genes,]
  for (i in valid_cells_names ) {
    df=t(dm)
    df=cbind(df,as.data.frame(meta[,i]))
    colnames(df)[ncol(df)]=('binary')
    df$binary <- factor(df$binary, levels = c("Others",i))
    write.csv(df,file=sprintf('var_%s_%s_binary_normdata.csv',n,i))
  }
  
}

####  try random genes  ############
for (n in num_vars) {
  set.seed(n)
  random.genes <- sample(rownames(pfc.bar.sub),size=n)
  dm=as.matrix(pfc.bar.sub@assays$RNA@data)[random.genes,]
  df=t(dm)
  for (i in valid_cells_names ) {
    df=t(dm)
    df=cbind(df,as.data.frame(meta[,i]))
    colnames(df)[ncol(df)]=('binary')
    write.csv(df,file=sprintf('random_%s_%s_binary_normdata.csv',n,i))
  }
}



