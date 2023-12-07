

setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

library(tidyverse)
library(Seurat)
library(EnhancedVolcano)
options(ggrepel.max.overlaps=Inf)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")

load('./data_file/sub_trim_seurat_pipe.Robj')
exn.sub <- mda; rm(mda)
exn.sub <- UpdateSeuratObject(exn.sub)
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
exn.sub@meta.data <- exn_meta

bar.mat <- read.csv('./results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
AI_valid <- rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid <- rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid <- rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid <- rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid <- rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
valid_cells_list <- list(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid);
valid_cells_names <- c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
for (i in 1:5) {
  exn.sub <- SetIdent(exn.sub, cells = colnames(exn.sub), value = 'Others')
  exn.sub <- SetIdent(exn.sub, cells = valid_cells_list[[i]], value = valid_cells_names[i])
  temp <- as.data.frame(Idents(exn.sub))
  colnames(temp) <- valid_cells_names[i]
  exn_meta <- cbind(exn_meta,temp)
}
exn.sub@meta.data <- exn_meta

barcode_genes <- c('barcode0','barcode1','barcode2','barcode3','barcode4')
genes.use <- setdiff(rownames(exn.sub),barcode_genes)
exn.sub <- subset(exn.sub, features = genes.use)

# sub_meta <- exn_meta %>% filter(DMS_valid=='Others' & AI_valid=='Others' & LH_valid=='Others' &
#                   BLA_valid=='Others' & MD_valid=='Others')
#exn.sub <- subset(exn.sub, cells = rownames(sub_meta), invert = TRUE)


### calculate DEGs without barcodes genes###
Idents(exn.sub) <- exn_meta$cluster
pfc_exn_markers_wobarcodes <- FindAllMarkers(exn.sub,logfc.threshold = 0.25,
               test.use = "MAST", min.pct = 0.1,
               verbose = TRUE, only.pos = T)
write.csv(pfc_exn_markers_wobarcodes,file='pfc_exn_markers_wobarcodes.csv')

pfc_exn_markers_wobarcodes=read.csv('../pfc_exn_markers_wobarcodes.csv',row.names = 1)
pfc_exn_markers_wobarcodes %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#### DMS ####
Idents(exn.sub) <- exn_meta$DMS_valid
table(Idents(exn.sub))
dms_markers <- FindMarkers(exn.sub, ident.1 = "DMS_valid", ident.2 = "Others",
                        logfc.threshold = 0.25, min.pct = 0.1,
                        test.use = "MAST", 
                         verbose = TRUE, only.pos = F)
head(dms_markers,n=20)
write.csv(dms_markers,file='./results/fig4_dms_markers.csv')
dms_markers <- read.csv('./results/fig4_dms_markers.csv',row.names = 1)

pdf('./results/fig4_DMS-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(dms_markers,
                     lab = rownames(dms_markers),
                     selectLab=c('C1ql3',"Nptx3","Srgap3",
                      'Syt6','Foxp2','Rprm',
                                'Slc24a3','Rab3c',
                                  'Rorb',
                                 'Ajap1','Cck'),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2,2),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "DMS-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 1e-10,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()


#### AI ####
Idents(exn.sub)=exn_meta$AI_valid
table(Idents(exn.sub))
ai_markers <- FindMarkers(exn.sub,ident.1="AI_valid",ident.2="Others",
                        logfc.threshold = 0.25, min.pct = 0.1,
                        test.use = "MAST", 
                        verbose = TRUE, only.pos = F)
head(ai_markers,n=20)

write.csv(ai_markers,file='./results/fig4_ai_markers.csv')
ai_markers <- read.csv('./results/fig4_ai_markers.csv',row.names = 1)

pdf('./results/fig4_AI-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(ai_markers,
                     lab = rownames(ai_markers),
                     selectLab=c('Nrip3','Foxp2','Sema5a',
                                 'Nptxr','Cck','C1ql3','Slc24a3'
                                 ),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "AI-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 1e-10,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()

#### MD ####
Idents(exn.sub)=exn_meta$MD_valid
table(Idents(exn.sub))
md_markers <- FindMarkers(exn.sub,ident.1="MD_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(md_markers,n=20)
#md_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)
write.csv(md_markers,file='./results/fig4_md_markers.csv')

pdf('./results/fig4_MD-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(md_markers,
                     lab = rownames(md_markers),
                     selectLab=c('Syt6','Foxp2','Crym','Rprm','Cyr61',
                                 'Hs3st4','Sla','Syt6',
                                 'Nptxr','Cck','C1ql3','Npy2r'),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "MD-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 1e-10,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()

#### LH ####
Idents(exn.sub)=exn_meta$LH_valid
table(Idents(exn.sub))

lh_markers <- FindMarkers(exn.sub,ident.1="LH_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(lh_markers,n=20)
#lh_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)
write.csv(lh_markers,file='./results/fig4_lh_markers.csv')
lh_markers <- read.csv('./results/fig4_lh_markers.csv',row.names = 1)

pdf('./results/fig4_LH-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(lh_markers,
                     lab = rownames(lh_markers),
                     selectLab=c('Camk2d','Nrn1','Cck','Pcdh15',
                                 'Crym','Pou3f1','Bcl11b','Nrip3',
                                 "Cyr61"),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "LH-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 1e-10,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()


#### BLA ####
Idents(exn.sub)=exn_meta$BLA_valid
table(Idents(exn.sub))

bla_markers <- FindMarkers(exn.sub,ident.1="BLA_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(bla_markers,n=20)

write.csv(bla_markers,file='./results/fig4_bla_markers.csv')
bla_markers=read.csv('./results/fig4_bla_markers.csv',row.names = 1)

pdf('./results/fig4_BLA-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(bla_markers,
                     lab = rownames(bla_markers),
                     selectLab=c('Foxp2','Rprm','Tle4',
                                 'C1ql3','Nptxr','Cnr1','Syt17',
                                 'Slc24a3'),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "BLA-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 1e-10,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()
