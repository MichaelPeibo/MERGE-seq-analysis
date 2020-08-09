##############################
### perform DEG analysis in Seurat ###
##############################

# use seurat to calculate DEGs

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

