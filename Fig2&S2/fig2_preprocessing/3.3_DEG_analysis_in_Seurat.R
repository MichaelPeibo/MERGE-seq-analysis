##############################
### perform DEG analysis in Seurat ###
##############################

# use seurat to calculate DEGs
## all cells
load('./1_seurat_processing_all_cells/pfc_all.seuratmerge.raw.bar5.Robj')
mda[["percent.mt"]] <- PercentageFeatureSet(object = mda, pattern = "^mt-")
mda[["percent.rp"]] <- PercentageFeatureSet(object = mda, pattern = "^Rps|^Rpl")

mda <- subset(x = mda, subset=nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA>1000 &
	nCount_RNA<60000 & percent.mt < 20)## pt.mt<10 will get 18035 samples

mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 50000)
meta=read.csv('/home/xupb/scRNA_data/mouse_pfc/manuscripts/scanpy/pfc_seurat_merge_meta.csv')
mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2005)
seurat_var_genes=VariableFeatures(object = mda)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
seurat_var_genes=setdiff(seurat_var_genes,barcode_genes)
length(seurat_var_genes)
Idents(mda)=meta$leiden_coarse
### all cluster markers ###
pfc_all_markers_cluster<- FindAllMarkers(object = mda, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(pfc_all_markers_cluster,file = 'pfc_all_markers_cluster.csv')

### excitatory neurons only ###
load('sub_trim_seurat_pipe.Robj')
meta=read.csv('sub_exn_trim_meta.csv',header=T,row.names=1)
Idents(mda)=meta$cluster

### all cluster markers ###
pfc_sub_exn_trim_markers_cluster<- FindAllMarkers(object = mda, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(pfc_sub_exn_trim_markers_cluster,file = 'pfc_sub_exn_trim_markers_cluster.csv')

### all leiden cluster markers ###

Idents(mda)=meta$leiden

pfc_sub_exn_trim_markers_res0.3_pcs30<- FindAllMarkers(object = mda, only.pos = TRUE, min.pct = 0.25,
                                           logfc.threshold = 0.25)
write.csv(pfc_sub_exn_trim_markers_res0.3_pcs30,file = 'pfc_sub_exn_trim_markers_res0.3_pcs30.csv')

