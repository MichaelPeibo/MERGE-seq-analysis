

library(Seurat)
### on cluster convert to loom file
load('./5_barcode_matrix_analysis/pfc.bar.sub.Robj')
pfc.bar.sub@meta.data$binary=Idents(pfc.bar.sub)
head(pfc.bar.sub@meta.data)
meta_info=pfc.bar.sub@meta.data
table(is.na(meta_info))
pfc.bar.sub@graphs <- list()
pfc.bar.sub.loom=as.loom(x = pfc.bar.sub,filename='pfc.bar.sub.loom')
pfc.bar.sub.loom$close_all()


## go on this analysis on Scanpy
## run file Supp_binary_projection.ipynb ###