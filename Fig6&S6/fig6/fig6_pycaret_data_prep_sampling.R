

setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

library(Seurat)
library(dplyr)
library(data.table)

load('./data_file/sub_trim_seurat_pipe.Robj')
exn.sub <- mda; rm(mda)
exn.sub <- UpdateSeuratObject(exn.sub)
barcode_genes <- c('barcode0','barcode1','barcode2','barcode3','barcode4')
genes.use <- setdiff(rownames(exn.sub),barcode_genes)
exn.sub <- subset(exn.sub, features = genes.use)
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
exn.sub@meta.data <- exn_meta
bar.mat <- read.csv('./results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
AI_valid <- rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid <- rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid <- rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid <- rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid <- rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
valid_cells_list <- list(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid);
valid_cells_names <- c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
for (i in 1:5) {
  exn.sub <- SetIdent(exn.sub, cells = colnames(exn.sub), value = 'Others')
  exn.sub <- SetIdent(exn.sub, cells = valid_cells_list[[i]], value = valid_cells_names[i])
  temp <- as.data.frame(Idents(exn.sub))
  colnames(temp) <- valid_cells_names[i]
  exn_meta <- cbind(exn_meta,temp)
}
exn.sub@meta.data <- exn_meta

### filter exn.sub by bar.mat cells ###
exn.sub <- subset(exn.sub, cells = rownames(bar.mat))

seurat_filter_genes = function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
exn.sub <- seurat_filter_genes(exn.sub, min.cells=3)
write.csv(exn_meta, file='./results/exn_meta_with_valid_projection.csv')
save(exn.sub,file = './results/exn.sub_with_valid_projection.Robj')

valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
####  try hvg number of genes ############

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
    fwrite(df,file=sprintf('var10000_%s_binary_normdata_3000cells_%s.csv',i,n))
  }
}



#### all cells ####
load('./results/exn.sub_with_valid_projection.Robj')
exn_meta <- read.csv('./results/exn_meta_with_valid_projection.csv',row.names = 1)
exn_meta_temp <- exn_meta[colnames(exn.sub),]
set.seed(123)

for (i in valid_cells_names ) {
  exn.sub=FindVariableFeatures(object = exn.sub, selection.method = "vst", nfeatures = 100)
  var.genes=VariableFeatures(exn.sub)
  dm=as.matrix(exn.sub@assays$RNA@data)[var.genes,]
  df=t(dm)
  df=cbind(df,as.data.frame(exn_meta_temp[,i]))
  colnames(df)[ncol(df)]=('binary')
  df$binary <- factor(df$binary, levels = c("Others",i))
  df=mutate(df, binary_clus = ifelse(binary == i , 1, 0))
  df$binary_clus <- as.factor(df$binary_clus)
  df <- subset(df, select = -c(binary))
  colnames(df)[5001]='binary'
  fwrite(df,file=sprintf('./results/fig6_var100_%s_binary_normdata_allcells.csv',i))
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








