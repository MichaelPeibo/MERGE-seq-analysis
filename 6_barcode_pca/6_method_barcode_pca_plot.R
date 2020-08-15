

library(tidyverse)
library(Seurat)
library(wesanderson)
library(ggsci)
library(ComplexHeatmap)
library(dendextend)
library(dendsort)

load('./3_scanpy_processing_exn/sub_trim_seurat_pipe.Robj')
exn_meta_valid=read.csv('./5_barcode_matrix_analysis/exn_meta_valid.csv',row.names=1)
bar.mat=read.csv('./5_barcode_matrix_analysis/bar.mat.drop.csv',row.names = 1)

AI_valid=rownames(bar.mat[bar.mat$AI>0,])
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])
MD_valid=rownames(bar.mat[bar.mat$MD>0,])
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])
LH_valid=rownames(bar.mat[bar.mat$LH>0,])
valid_cells_union=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))
mda=subset(mda,cells=valid_cells_union)
meta=exn_meta_valid[colnames(mda),]

bar.mat=bar.mat[rownames(meta),]
rm(mda)

### Seurat pca

temp.bar=bar.mat[colnames(mda),]
colnames(temp.bar)=c("AI", "DMS", "MD",'BLA','LH')
# Transform table
# Add barcode data as a new assay independent from RNA
mda[["bar"]] <- CreateAssayObject(counts =as.data.frame(t(temp.bar)))
DefaultAssay(object = mda) <- "bar"

mda <- NormalizeData(object = mda, normalization.method = "CLR")
mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 5)
mda = ScaleData(object = mda, vars.to.regress = c("nCount_RNA"))
mda <- RunPCA(object = mda,nfeatures.print = 10,npcs = 4)


save(mda,file='binary_seurat_exn.Robj')

load('binary_seurat_exn.Robj')
mda=subset(mda,idents = 'Other',invert=T)

pca_embed=Embeddings(object = mda[["pca"]])
write.csv(pca_embed,file='pca_embed.csv')

col_pal=c("#FAD510","#ED0000FF","#42B540FF","#0099B4FF","#9986A5",
          "#FDAF91FF","#CCBA72","#ADB6B6FF","#5050FFFF",'#ABDDDE')
mda@meta.data$cluster=meta[colnames(mda),]$cluster

pdf('pca plot of barcode by binary cluster.pdf',width = 8)
DimPlot(mda, reduction = "pca",
        pt.size = 2,group.by = 'ident',label = T,label.size = 8,repel = T
) + scale_colour_manual(values=col_pal)
dev.off()

DimPlot(mda, reduction = "pca",
        pt.size = 2,group.by = 'cluster',label = T,label.size = 8
)

plot_list <-vector("list", 6) 
for(i in 1:5) {
  plot_list[[i]] <- FeaturePlot(mda, features = rownames(mda)[i],
                                pt.size = 1,
                                min.cutoff = "q20", max.cutoff = "q80",
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1))
  plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}
plot_list[[6]]=DimPlot(mda, reduction = "pca",
                       pt.size = 2,group.by = 'ident',label = T,label.size = 15,repel = T) + 
  scale_colour_manual(values=col_pal)+ labs(title = "Projection Cluster") + NoLegend()
plot_list[[6]]<- AugmentPlot(plot_list[[6]],dpi = 400)
pl=cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=3, align="h")
cowplot::save_plot("barcoded pc.pdf", pl, base_width=15,base_height = 10)




DefaultAssay(object = mda) <- "RNA"
pca_embed=Embeddings(object = mda[["pca"]])

genes.to.plot=c("Nptxr",'Nrn1','Cck','AI',
                'Nrip3','Ldhb','Crym','LH')
plot_list <-vector("list", 8) 

for(i in 1:8) {
  plot_list[[i]] <- FeaturePlot(mda, features = genes.to.plot[i],
                                pt.size = 1,
                                min.cutoff = "q20", max.cutoff = "q80",reduction = 'pca',
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}
pl=cowplot::plot_grid(plotlist=plot_list, labels = letters[6:13],
                      ncol=4, align="h")
cowplot::save_plot("selected features AI LH.pdf", pl, base_width=20,base_height = 10)


### complexheatmap binary projection ###

df=as.matrix(mda@assays$bar@data)
rownames(df)=c("AI",'DMS','MD','BLA','LH')
pheno = data.frame(Idents(mda),mda@meta.data$cluster)
colnames(pheno)=c('projection_cluster','transcription_cluster')
rownames(pheno)=colnames(mda)
sub_exn_palette=c('#7fc97f', '#beaed4', '#fdc086', '#386cb0', '#f0027f', '#bf5b17',
                  '#666666')
ann_colors = list(
  'projection_cluster'=col_pal,
  'transcription_cluster'=sub_exn_palette
)
names(ann_colors$'projection_cluster')=levels(pheno$'projection_cluster')
names(ann_colors$'transcription_cluster')=levels(pheno$'transcription_cluster')

dend = hclust(dist(df))
dend = color_branches(dend, k = 2)
col<- circlize::colorRamp2(c(0,1,2,3,4),viridis::magma(5))

anno = HeatmapAnnotation(df=pheno,col=ann_colors,
                         annotation_legend_param = list(
                           'projection_cluster' = list(title='Projection',
                                                       legend_direction='horizontal',
                                                       nrow = 3,
                                                       title_gp = gpar(fontsize = 14)),
                           'transcription_cluster' = list(title='Transcription',
                                                          legend_direction='horizontal',
                                                          nrow = 3,
                                                          title_gp = gpar(fontsize = 14))))

p=Heatmap(df,name = "Normalized Counts", 
          show_column_names=F,
          cluster_rows =dend,row_split=2,row_title = NULL,
          #cluster_columns =hclust(dist(t(df))),
          column_split=2,
          column_title = NULL,
          heatmap_legend_param = list(legend_direction = "horizontal"),
          top_annotation = anno,
          #right_annotation=rowAnnotation(df=pheno_row,col=ann_row_colors),
          col = col
)
pdf('complexheatmap binary projection magma.pdf',
    width=15,height=6)
draw(p, heatmap_legend_side = "bottom", merge_legend = TRUE,
     annotation_legend_side = "bottom")
dev.off()











