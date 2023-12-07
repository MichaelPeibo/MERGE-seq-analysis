

library(Seurat)
library(dplyr)
library(data.table)

load('exn.sub_with_valid_projection.Robj')
exn_meta <- read.csv("exn_meta_with_valid_projection.csv", row.names = 1)
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')

####  try hvg number of genes ############

for (n in 1:100) {
  set.seed(n)
  random.cells <- sample(colnames(exn.sub),size=1000)
  temp=subset(exn.sub,cells = random.cells)
  temp=FindVariableFeatures(object = temp, selection.method = "vst", nfeatures = 5000)
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
    colnames(df)[5001]='binary'
    fwrite(df,file=sprintf('./sampling_results/var5000_%s_binary_normdata_1000cells_%s.csv',i,n))
  }
}
### random sampling 1000 cells of all cells for 100 times ###

####  try random genes  ############

for (n in 1:100) {
  set.seed(n)
  random.cells <- sample(colnames(exn.sub),size=1000)
  random.genes <- sample(rownames(exn.sub),size=100)
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
    colnames(df)[101]='binary'
    fwrite(df,file=sprintf('./sampling_results/ran100_%s_binary_normdata_1000cells_%s.csv',i,n))
  }
}