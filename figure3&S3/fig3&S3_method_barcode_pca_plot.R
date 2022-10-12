

library(tidyverse)
library(Seurat)
library(wesanderson)
library(ggsci)
library(ComplexHeatmap)
library(dendextend)
library(dendsort)


DefaultAssay(object = exn.sub) <- "bar"
exn.sub <- NormalizeData(object = exn.sub, normalization.method = "CLR")
exn.sub <- FindVariableFeatures(object = exn.sub, selection.method = "vst", nfeatures = 5)
exn.sub = ScaleData(object = exn.sub, vars.to.regress = c("nCount_RNA"))
exn.sub <- RunPCA(object = exn.sub,nfeatures.print = 10,npcs = 4)

# save(mda,file='binary_seurat_exn.Robj')
# load('binary_seurat_exn.Robj')

pfc.bar=subset(exn.sub,idents = 'Other',invert=T)

pca_embed=Embeddings(object = pfc.bar[["pca"]])
write.csv(pca_embed,file='pca_embed.csv')

col_pal=c("#FAD510","#ED0000FF","#42B540FF","#0099B4FF","#9986A5",
          "#FDAF91FF","#CCBA72","#ADB6B6FF","#5050FFFF",'#ABDDDE','grey20')

pdf('pca plot of barcode by binary cluster only top10.pdf',width = 8)
DimPlot(pfc.bar, reduction = "pca",
        pt.size = 2,group.by = 'ident',label = T,label.size = 8,repel = T
) + scale_colour_manual(values=col_pal)
dev.off()

plot_list <-vector("list", 6)
for(i in 1:5) {
  plot_list[[i]] <- FeaturePlot(pfc.bar, features = rownames(pfc.bar)[i],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) & 
    theme(plot.title = element_text(size = 20,vjust = 1),legend.position = c(0.1,0.2))
  #plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}
plot_list[[6]]=DimPlot(pfc.bar, reduction = "pca",
                       pt.size = 1,group.by = 'ident',label = T,label.size = 5,repel = T) + 
  scale_colour_manual(values=col_pal)+ labs(title = "Projection Cluster") & NoLegend()

#plot_list[[6]]<- AugmentPlot(plot_list[[6]],dpi = 400)
pl=cowplot::plot_grid(plotlist=plot_list, labels = "auto",
                      ncol=3, align="h")
cowplot::save_plot("barcoded pc black only top10.pdf", pl, base_width=15,base_height = 10)




DefaultAssay(object = exn.sub) <- "RNA"
pca_embed=Embeddings(object = exn.sub[["pca"]])

genes.to.plot=c("Nptxr",'C1ql3','Cck','DMS',
                'Rprm','Npas4','Cyr61','MD')
plot_list <-vector("list", 8) 

for(i in 1:3) {
  plot_list[[i]] <- FeaturePlot(exn.sub, features = genes.to.plot[i],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                reduction = 'pca',
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  #plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}

for(i in 5:7) {
  plot_list[[i]] <- FeaturePlot(exn.sub, features = genes.to.plot[i],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                reduction = 'pca',
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))
  #plot_list[[i]] <- AugmentPlot(plot_list[[i]],dpi = 400)
}
DefaultAssay(object = exn.sub) <- "bar"

plot_list[[4]] <- FeaturePlot(exn.sub, features = 'bDMS',
                              pt.size = 1,
                              #min.cutoff = "q20", max.cutoff = "q80",
                              reduction = 'pca',
                              cols = c("#CECECE", "#CBDAC2", 
                                       RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
  theme(plot.title = element_text(size = 20,vjust = 1),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
#plot_list[[4]] <- AugmentPlot(plot_list[[4]],dpi = 400)

plot_list[[8]] <- FeaturePlot(exn.sub, features = 'bMD',
                              pt.size = 1,
                              #min.cutoff = "q20", max.cutoff = "q80",
                              reduction = 'pca',
                              cols = c("#CECECE", "#CBDAC2", 
                                       RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
  theme(plot.title = element_text(size = 20,vjust = 1),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
#plot_list[[8]] <- AugmentPlot(plot_list[[8]],dpi = 400)


pl=cowplot::plot_grid(plotlist=plot_list, labels = letters[6:13],
                      ncol=4, align="h")
cowplot::save_plot("selected features DMS MD 9368 cells.pdf", pl, base_width=20,base_height = 10)


### complexheatmap binary projection ###
mda=pfc.bar
temp=exn_meta[colnames(pfc.bar.sub),]
temp$binary=Idents(pfc.bar.sub)
head(temp)
df=as.matrix(mda@assays$bar@data)
rownames(df)=c("AI",'DMS','MD','BLA','LH')
mda@meta.data$cluster=as.factor(temp$cluster)
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
          cluster_rows =F,row_title = NULL,
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











