

setwd('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue')

library(Seurat)
library(dplyr)
library(data.table)


load('pfc.bar.sub.Robj')
load('sub_trim_seurat_pipe.Robj')
exn.sub=mda
rm(mda)
seurat_filter_genes = function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
pfc.bar.sub=seurat_filter_genes(pfc.bar.sub,min.cells=3)


barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
exn.sub=subset(exn.sub,features = setdiff(rownames(exn.sub),barcode_genes))
exn.sub=seurat_filter_genes(exn.sub,min.cells=3)
exn_meta=read.csv('exn_meta_valid.csv',row.names=1)


valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')

####  try hvg number of genes ############

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation")


for (n in 1:100) {
  set.seed(n)
  random.cells <- sample(colnames(exn.sub),size=3000)
  temp=subset(exn.sub,cells = random.cells)
  temp=FindVariableFeatures(object = temp, selection.method = "vst", nfeatures = 10000)
  var.genes.temp=VariableFeatures(temp)
  dm=as.matrix(temp@assays$RNA@data)[var.genes.temp,] ## 
  exn_meta_temp=exn_meta[colnames(temp),]
  for (i in valid_cells_names ) {
    df=t(dm)
    df=cbind(df,as.data.frame(exn_meta_temp[,i]))
    colnames(df)[ncol(df)]=('binary')
    df$binary <- factor(df$binary, levels = c("Others",i))
    df=mutate(df, binary_clus = ifelse(binary == i , 1, 0))
    df$binary_clus <- as.factor(df$binary_clus)
    df <- subset(df, select = -c(binary))
    colnames(df)[10001]='binary'
    fwrite(df,file=sprintf('var1500_%s_binary_normdata_3000cells_%s.csv',i,n))
  }
}

####  try random genes  ############

for (n in 1:100) {
  set.seed(n)
  random.cells <- sample(colnames(exn.sub),size=3000)
  random.genes <- sample(rownames(exn.sub),size=50)
  temp=subset(exn.sub,cells = random.cells,features = random.genes)
  dm=as.matrix(temp@assays$RNA@data)[random.genes,] ## 
  exn_meta_temp=exn_meta[colnames(temp),]
  for (i in valid_cells_names ) {
    df=t(dm)
    df=cbind(df,as.data.frame(exn_meta_temp[,i]))
    colnames(df)[ncol(df)]=('binary')
    df$binary <- factor(df$binary, levels = c("Others",i))
    df=mutate(df, binary_clus = ifelse(binary == i , 1, 0))
    df$binary_clus <- as.factor(df$binary_clus)
    df <- subset(df, select = -c(binary))
    colnames(df)[51]='binary'
    fwrite(df,file=sprintf('ran50_%s_binary_normdata_3000cells_%s.csv',i,n))
  }
}

#### all cells ####

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation")

set.seed(123)

for (i in valid_cells_names ) {

  exn.sub=FindVariableFeatures(object = exn.sub, selection.method = "vst", nfeatures = 50)
  var.genes=VariableFeatures(exn.sub)
  dm=as.matrix(exn.sub@assays$RNA@data)[var.genes,]
  df=t(dm)
  df=cbind(df,as.data.frame(exn_meta[,i]))
  colnames(df)[ncol(df)]=('binary')
  df$binary <- factor(df$binary, levels = c("Others",i))
  df=mutate(df, binary_clus = ifelse(binary == i , 1, 0))
  df$binary_clus <- as.factor(df$binary_clus)
  df <- subset(df, select = -c(binary))
  colnames(df)[51]='binary'
  fwrite(df,file=sprintf('var50_%s_binary_normdata_allcells.csv',i))
}

test=read.csv('var50_AI_valid_binary_normdata_allcells.csv',b)







