### step 3 ###
### Analysis on Excitatory neuron  ###
### do cells filtering and genes filtering as in scanpy ###
sub_meta=read.csv('./3_scanpy_processing_exn/sub_exn_meta.csv',row.names=1,header=T)
sub_vars=read.csv('./3_scanpy_processing_exn/sub_exn_vars.csv',row.names=1,header=T)

load(file='pfc_all.seuratmerge.raw.bar5.Robj')
pfc.sub=subset(mda,cells=rownames(sub_meta),features=rownames(sub_vars))

seuratv3_sub_pipe <- function(mda) {
  mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 50000)
  mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2000)
  cc.genes <- readLines(con = "/home/xupb/regev_lab_cell_cycle_genes.txt")
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]
  mda <- CellCycleScoring(object = mda, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

  seurat_var_genes=VariableFeatures(object = mda)
  barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
  seurat_var_genes=setdiff(seurat_var_genes,barcode_genes)
  length(seurat_var_genes)
  write.csv(seurat_var_genes,file='sub_exn_seurat_var_genes.csv')
  VariableFeatures(object = mda)=seurat_var_genes

  mda <- ScaleData(object = mda, 
                  vars.to.regress = c("nCount_RNA",'nFeature_RNA','facs','percent.mt'),
                  verbose = TRUE)


  mda <- RunPCA(object =  mda, features=seurat_var_genes,npcs = 50, verbose = FALSE)

  print('Saving file in Robj format......')
  mda <- RunUMAP(object = mda, dims = 1:30)
  mda <- RunTSNE(object = mda, dims = 1:30)
  mda <- FindNeighbors(object = mda, dims = 1:30)
  filename = "sub_seurat_pipe.Robj"
  save(mda,file=filename)
  return(mda)
}

pfc.sub=seuratv3_sub_pipe(pfc.sub)
## export pca embeddings
pca_embed=Embeddings(pfc.sub, reduction = "pca")
write.csv(pca_embed,file='sub_exn_seurat_pca_embed.csv')

## convert to loom
head(pfc.sub@meta.data)
meta_info=pfc.sub@meta.data
table(is.na(meta_info))
pfc.sub@graphs <- list()
pfc.sub.loom =as.loom(x = pfc.sub,filename='pfc_exn_sub.loom')
pfc.sub.loom$close_all()
