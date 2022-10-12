


library(dplyr)
library(Seurat)
library(Matrix)
### conda activate urd ###
library(MetaNeighbor) ### packageVersion('MetaNeighbor') = 1.6.0
library(SummarizedExperiment)

library(cowplot)
library(matrixStats)
library(ggplot2)
library(ggforce)
library(wesanderson)

load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/pfc_all.seuratmerge.raw.bar5.Robj')
#load('./1_seurat_processing_all_cells/pfc_all.seuratmerge.raw.bar5.Robj')
mda[["percent.mt"]] <- PercentageFeatureSet(object = mda, pattern = "^mt-")
mda[["percent.rp"]] <- PercentageFeatureSet(object = mda, pattern = "^Rps|^Rpl")

mda <- subset(x = mda, subset=nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA>1000 &
	nCount_RNA<60000 & percent.mt < 20)## pt.mt<10 will get 18035 samples

mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 50000)
#meta=read.csv('./2_scanpy_processing_all_cells/pfc_seurat_merge_meta.csv')
meta=read.csv('/home/xupb/scRNA_data/mouse_pfc/manuscripts/scanpy/pfc_seurat_merge_meta.csv')

mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2005)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')


load(file = "/home/xupb/scRNA_data/mouse_pfc/nc_2019/all_cell_inDD_sobj.RData")
pfc=all_cell_inDD_sobj
rm(all_cell_inDD_sobj)
table(pfc@meta.data$CellType)
Idents(pfc)=pfc@meta.data$CellType
pfc@meta.data$celltype=Idents(pfc)

ob.list <- list(mda, pfc)
library(future)
options(future.globals.maxSize = 8000 * 1024^2)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)
mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)
mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)
DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)

save(mda.integrated,file='integrated_all_nc2019.Robj')
setwd('/home/xupb/scRNA_data/mouse_pfc/nc_2019')
load('integrated_all_nc2019.Robj')



var.gene=VariableFeatures(object = mda.integrated)

#combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
combined_mat=as.matrix(mda.integrated@assays$integrated@data)

mda@meta.data$celltype=Idents(mda)=meta$leiden_coarse
pfc@meta.data$celltype=Idents(pfc)

Study_ID = rep(c('1', '2'), c(ncol(mda), ncol(pfc)))
Celltype = c(as.character(mda@meta.data$celltype),as.character(pfc@meta.data$celltype))

dat=SummarizedExperiment(assays=list(counts=combined_mat))
celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)


library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)

pdf('heatmap_pfc_all_integrate_normdata_metaneighbor_nc2019.pdf')
gplots::heatmap.2(celltype_NV,
margins=c(8,8),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 8), rep('deeppink',8)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)

par(lend = 1)      # square line ends for the color legend
legend("top",      # location of the legend on the heatmap plot
    legend = c("1-Bhattacherjee el al.2019", "2-This project"), # category labels
    col = c("darkgreen", "deeppink"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()

library(pheatmap)
ann_row=data.frame(dataset=c(rep("Bhattacherjee", 8), rep('This_study',8)))
rownames(ann_row)=rownames(celltype_NV)

pdf('pheatmap_pfc_all_integrate_normdata_metaneighbor_nc2019.pdf',width=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=4,cutree_cols=4,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Bhattacherjee = wes_palette("Rushmore1")[3], 
            This_study = 'deeppink'))
         )
dev.off()

### ONLY exciatatory neurons ###

###---------------------------------------------------
load(file = "/home/xupb/scRNA_data/mouse_pfc/nc_2019/PFC_excitatory_saline_sobj.RData")
pfc_nc=PFC_excitatory_saline_sobj
rm(PFC_excitatory_saline_sobj)

#pfc_nc_meta=read.csv('./7_metaneighbor_analysis/PFC_exc_saline_meta.csv',row.names=1,header=T)
pfc_nc_meta=read.csv('/home/xupb/scRNA_data/mouse_pfc/nc_2019/PFC_exc_saline_meta.csv',row.names=1,header=T)

Idents(pfc_nc)=pfc_nc_meta$L1_cluster

#load('./3_scanpy_processing_exn/sub_trim_seurat_pipe.Robj')
load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/sub_trim_seurat_pipe.Robj')

pfc=mda
rm(mda)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')

#pfc_meta=read.csv('./3_scanpy_processing_exn/sub_exn_trim_meta.csv',row.names = 1,header = T)
pfc_meta=read.csv('/home/xupb/scRNA_data/mouse_pfc/manuscripts/scanpy/sub_exn_trim_meta.csv',row.names = 1,header = T)

Idents(pfc)=pfc_meta$cluster_v1

ob.list <- list(pfc_nc, pfc)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)

mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)

DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)


save(mda.integrated,file='integrated_exn_nc2019.Robj')
load('integrated_exn_nc2019.Robj')


var.gene=VariableFeatures(object = mda.integrated)

#combined_mat=as.matrix(mda.integrated@assays$integrated@scale.data)
combined_mat=as.matrix(mda.integrated@assays$integrated@data)

pfc_nc@meta.data$celltype=Idents(pfc_nc)
pfc@meta.data$celltype=Idents(pfc)

Study_ID = rep(c('1', '2'), c(ncol(pfc_nc), ncol(pfc)))
Celltype = c(as.character(pfc_nc@meta.data$celltype),as.character(pfc@meta.data$celltype))

dat=SummarizedExperiment(assays=list(counts=combined_mat))
celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)

library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


pdf('heatmap_pfc_exn_integrate_normdata_metaneighbor_nc2019.pdf')
gplots::heatmap.2(celltype_NV,
margins=c(8,8),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 13), rep('deeppink',7)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)

par(lend = 1)      # square line ends for the color legend
legend("top",      # location of the legend on the heatmap plot
    legend = c("1-Bhattacherjee el al.2019", "2-This project"), # category labels
    col = c("darkgreen", "deeppink"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()



library(pheatmap)
ann_row=data.frame(dataset=c(rep("Bhattacherjee", 13), rep('This_study',7)))
rownames(ann_row)=rownames(celltype_NV)

pdf('pheatmap_pfc_exn_integrate_normdata_metaneighbor_nc2019.pdf',width=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=5,cutree_cols=5,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Bhattacherjee = wes_palette("Rushmore1")[3], 
            This_study = 'deeppink'))
         )
dev.off()
