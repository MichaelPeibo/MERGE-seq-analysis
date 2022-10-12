


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


setwd('/home/xyz2020/xpb/pfc_2021')
###---------------------------------------------------

load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/GSE161936_Lui_Nguyen_Seuratobjects.RData')
pfc_cell=Fig1Rbp4
table(Idents(pfc_cell))

newids<- c('Cd44','Figf','Otof','Pld5',
           'Cxcr7','Npr3','Tshz2')
names(newids) <- levels(pfc_cell)
pfc_cell<- RenameIdents(pfc_cell, newids)


load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/GSE161936_Lui_Nguyen_Seuratobjects.RData')

pfc_retro=Fig3Rbp4retro
head(pfc_retro@meta.data)
pfc_retro@meta.data$Cluster=Idents(pfc_retro)

table(pfc_retro@meta.data$Cluster,pfc_retro@meta.data$Cre.line)

table(pfc_retro@meta.data$Batch)
table(pfc_retro@meta.data$Group)
table(pfc_retro@meta.data$Cre.line)

Idents(pfc_retro)=pfc_retro@meta.data$Cre.line
table(Idents(pfc_retro))

pdf('Fig3Rbp4retro umap.pdf')
FeaturePlot(pfc_retro, features = 'Nrip3',reduction = "pca")
dev.off()

pdf('Fig3Rbp4retro vln nrip3.pdf')
VlnPlot(object = pfc_retro, features = 'Nrip3')
dev.off()







load('/home/xupb/scRNA_data/mouse_pfc/manuscripts/sub_trim_seurat_pipe.Robj')
pfc=mda
rm(mda)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
pfc_meta=read.csv('/home/xupb/scRNA_data/mouse_pfc/manuscripts/scanpy/sub_exn_trim_meta.csv',row.names = 1,header = T)
Idents(pfc)=pfc_meta$cluster

ob.list <- list(pfc_cell, pfc)
mda.features <- SelectIntegrationFeatures(object.list = ob.list, nfeatures = 2000)

mda.anchors <- FindIntegrationAnchors(object.list = ob.list, 
    anchor.features = mda.features, verbose = FALSE)

mda.integrated <- IntegrateData(anchorset = mda.anchors,
    verbose = FALSE)

DefaultAssay(mda.integrated) <- "integrated"
mda.integrated <- ScaleData(mda.integrated, verbose = FALSE)

setwd('/home/xyz2020/xpb/pfc_2021')
save(mda.integrated,file='integrated_exn_nc2021.Robj')
load('integrated_exn_nc2021.Robj')


var.gene=VariableFeatures(object = mda.integrated)
combined_mat=as.matrix(mda.integrated@assays$integrated@data)

pfc_cell@meta.data$celltype=Idents(pfc_cell)
pfc@meta.data$celltype=Idents(pfc)

Study_ID = rep(c('1', '2'), c(ncol(pfc_cell), ncol(pfc)))
Celltype = c(as.character(pfc_cell@meta.data$celltype),as.character(pfc@meta.data$celltype))

dat=SummarizedExperiment(assays=list(counts=combined_mat))
celltype_NV=MetaNeighbor::MetaNeighborUS(var_genes = var.gene,
dat = dat,
study_id = Study_ID,
cell_type = Celltype,fast_version=T)

library(gplots)
library(RColorBrewer)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)


pdf('heatmap_pfc_exn_integrate_normdata_metaneighbor_nc2021.pdf')
gplots::heatmap.2(celltype_NV,
margins=c(8,8),
keysize=1,
key.xlab="AUROC",
key.title="NULL",
trace = "none",
density.info =c("none"),
col = cols,
RowSideColors=c(rep("darkgreen", 7), rep('deeppink',7)),
breaks = breaks,
dendrogram="row",     # only draw a row dendrogram
offsetRow=0.1,
offsetCol=0.1,
cexRow = 1.5,
cexCol = 1.5)

par(lend = 1)      # square line ends for the color legend
legend("top",      # location of the legend on the heatmap plot
    legend = c("1-Jan el al.2021", "2-This study"), # category labels
    col = c("darkgreen", "deeppink"),  # color key
    lty= 1,             # line style
    lwd = 10            # line width
)
dev.off()



library(pheatmap)
ann_row=data.frame(dataset=c(rep("Jan", 7), rep('This',7)))
rownames(ann_row)=rownames(celltype_NV)

pdf('pheatmap_pfc_exn_integrate_normdata_metaneighbor_nc2021.pdf',width=10)
pheatmap(celltype_NV, 
         color = cols, cutree_rows=4,cutree_cols=4,
         breaks = breaks, cellwidth = 12, cellheight = 12,
         cluster_rows = TRUE,cluster_cols = TRUE, treeheight_col = 0,
         clustering_method = "ward.D2",
         annotation_row = ann_row,
         annotation_colors = list(dataset = c(Jan = wes_palette("Rushmore1")[3], 
            This = 'deeppink'))
         )
dev.off()
