
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

library(tidyverse)
library(Seurat)
library(wesanderson)
library(ggsci)


load('./results/binary_seurat_exn.Robj')

DefaultAssay(object = pfc.bar) <- "RNA"
genes.to.plot=c("Nptxr",'DMS',
                'Rprm','MD')
plot_list <-vector("list", 4) 

plot_list[[1]] <- FeaturePlot(pfc.bar, features = genes.to.plot[1],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                reduction = 'pca',
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))



plot_list[[3]] <- FeaturePlot(pfc.bar, features = genes.to.plot[3],
                                pt.size = 1,
                                #min.cutoff = "q20", max.cutoff = "q80",
                                reduction = 'pca',
                                cols = c("#CECECE", "#CBDAC2", 
                                         RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
    theme(plot.title = element_text(size = 20,vjust = 1),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20))


DefaultAssay(object = pfc.bar) <- "bar"
plot_list[[2]] <- FeaturePlot(pfc.bar, features = 'bDMS',
                              pt.size = 1,
                              #min.cutoff = "q20", max.cutoff = "q80",
                              reduction = 'pca',
                              cols = c("#CECECE", "#CBDAC2", 
                                       RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
  theme(plot.title = element_text(size = 20,vjust = 1),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))

plot_list[[4]] <- FeaturePlot(pfc.bar, features = 'bMD',
                              pt.size = 1,
                              #min.cutoff = "q20", max.cutoff = "q80",
                              reduction = 'pca',
                              cols = c("#CECECE", "#CBDAC2", 
                                       RColorBrewer::brewer.pal(9, "YlGnBu")[3:6])) + 
  theme(plot.title = element_text(size = 20,vjust = 1),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))


pl=cowplot::plot_grid(plotlist=plot_list, labels = NULL,
                      ncol=2, align="h")
cowplot::save_plot("./results/fig6_selected features DMS MD.pdf", pl, base_width=10,base_height = 10)












