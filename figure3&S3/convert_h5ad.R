
library(SeuratDisk)

### on cluster convert to loom file
pfc.bar.sub@meta.data$binary=Idents(pfc.bar.sub)
head(pfc.bar.sub@meta.data)
meta_info=pfc.bar.sub@meta.data
table(is.na(meta_info))

SaveH5Seurat(pfc.bar.sub, filename = "pfc.bar.sub.h5Seurat")
Convert("pfc.bar.sub.h5Seurat", dest = "h5ad")


SaveH5Seurat(exn.sub, filename = "exn.sub.h5Seurat")
Convert("exn.sub.h5Seurat", dest = "h5ad")


pfc.bar.sub@meta.data$binary_cluster=Idents(pfc.bar.sub)
write.csv(pfc.bar.sub@meta.data,file='meta_binary.csv')

DefaultAssay(object = pfc.bar.sub) <- "bar"
SaveH5Seurat(pfc.bar.sub, filename = "pfc.bar.sub.barcodegenes.h5Seurat")
Convert("pfc.bar.sub.barcodegenes.h5Seurat", dest = "h5ad")




pfc.bar.sub@graphs <- list()
pfc.bar.sub.loom=as.loom(x = pfc.bar.sub,filename='pfc.bar.sub.loom')
pfc.bar.sub.loom$close_all()

DefaultAssay(object = pfc.bar.sub) <- "bar"
pfc.bar.sub <- NormalizeData(object = pfc.bar.sub, normalization.method = "CLR")
pfc.bar.sub <- FindVariableFeatures(object = pfc.bar.sub, selection.method = "vst", nfeatures = 5)

pfc.bar.sub@meta.data$binary=Idents(pfc.bar.sub)
head(pfc.bar.sub@meta.data)
meta_info=pfc.bar.sub@meta.data
table(is.na(meta_info))
pfc.bar.sub@graphs <- list()
pfc.bar.sub.loom=as.loom(x = pfc.bar.sub,filename='pfc.bar.sub.bargenes.loom')
pfc.bar.sub.loom$close_all()
