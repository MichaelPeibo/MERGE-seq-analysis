library(Seurat)
library(dplyr)


load('./5_barcode_matrix_analysis/pfc.bar.sub.Robj')


seurat_filter_genes = function(seurat.object,min.cells=3){
  num.cells <- Matrix::rowSums(seurat.object@assays$RNA@counts > 0)
  genes.use <- names(num.cells[which(num.cells >= min.cells)])
  seurat.object=subset(seurat.object,features=genes.use)
  return(seurat.object)
}
pfc.bar.sub=seurat_filter_genes(pfc.bar.sub,min.cells=3)



####  try hvg 150 ############
set.seed(200)

pfc.bar.sub=FindVariableFeatures(object = pfc.bar.sub, selection.method = "vst", nfeatures = 150)
var.genes=VariableFeatures(pfc.bar.sub)
table(var.genes=='Syt6')
table(var.genes=='Foxp2')

dm=as.matrix(pfc.bar.sub@assays$RNA@data)[var.genes,] ## 
df=t(dm)

exn_meta_valid=read.csv('./5_barcode_matrix_analysis/exn_meta_valid.csv',row.names=1)
meta=exn_meta_valid[colnames(pfc.bar.sub),]
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')

for (i in valid_cells_names ) {
	df=t(dm)
	df=cbind(df,as.data.frame(meta[,i]))
	colnames(df)[ncol(df)]=('binary')
	df$binary <- factor(df$binary, levels = c("Others",i))
	write.csv(df,file=sprintf('var150_%s_binary_normdata.csv',i))
}

####  try random genes  ############

set.seed(150)
random.genes <- sample(rownames(pfc.bar.sub),size=150)

dm=as.matrix(pfc.bar.sub@assays$RNA@data)[random.genes,] ## 
df=t(dm)
for (i in valid_cells_names ) {
	df=t(dm)
	df=cbind(df,as.data.frame(meta[,i]))
	colnames(df)[ncol(df)]=('binary')
	write.csv(df,file=sprintf('random_150_%s_binary_normdata.csv',i))
}