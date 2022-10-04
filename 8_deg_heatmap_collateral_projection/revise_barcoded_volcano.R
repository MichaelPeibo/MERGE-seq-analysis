

setwd('/Users/peiboxu/Desktop/merge-seq analysis/')

library(tidyverse)
library(Seurat)
library(EnhancedVolcano)
options(ggrepel.max.overlaps=Inf)

load('sub_trim_seurat_pipe.Robj')
meta=read.csv('exn_meta_valid.csv',row.names = 1)
mda@meta.data=meta

barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
genes.use=setdiff(rownames(exn.sub),barcode_genes)
exn.sub=subset(exn.sub,features=genes.use)

### calculate DEGs without barcodes genes###
Idents(exn.sub)=meta$cluster
pfc_exn_markers_wobarcodes=FindAllMarkers(mda,logfc.threshold = 0.25,
               test.use = "MAST", min.pct = 0.1,
               verbose = TRUE, only.pos = T)
write.csv(pfc_exn_markers_wobarcodes,file='pfc_exn_markers_wobarcodes.csv')

pfc_exn_markers_wobarcodes=read.csv('pfc_exn_markers_wobarcodes.csv',row.names = 1)
pfc_exn_markers_wobarcodes %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#### DMS ####
Idents(exn.sub)=exn_meta$DMS_valid
table(Idents(exn.sub))
dms_markers=FindMarkers(exn.sub,ident.1="DMS_valid",ident.2="Others",
                        logfc.threshold = 0.25, min.pct = 0.1,
                        test.use = "MAST", 
                         verbose = TRUE, only.pos = F)
head(dms_markers,n=20)
write.csv(dms_markers,file='dms_markers.csv')

dms_markers=read.csv('dms_markers.csv',row.names = 1)
dms_markers_top = dms_markers %>% filter(avg_log2FC>0.5 & p_val_adj<10e-20)
common.genes=intersect(top20$gene,rownames(dms_markers_top))
common.genes
pfc_exn_markers_wobarcodes %>%
  filter(cluster=='L2/3-Rorb' | cluster=='L2/3-Rorb') %>%
  top_n(n = 20, wt = avg_log2FC) -> top20_layer23
common.genes.l23=intersect(top20_layer23$gene,rownames(dms_markers_top))
common.genes.l23

pdf('DMS-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(dms_markers,
                     lab = rownames(dms_markers),
                     selectLab=c('Syt6','Foxp2','Rprm',
                                 'Syt4','Cux2','Slc24a3','Marcksl1','Rab3c',
                                  'Synpr','Tnnc1','Etv1','S100b','Rorb',
                                 'Ajap1','C1ql3','Cck'),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5,2.5),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "DMS-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 10e-20,
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
ai_markers=FindMarkers(exn.sub,ident.1="AI_valid",ident.2="Others",
                        logfc.threshold = 0.25, min.pct = 0.1,
                        test.use = "MAST", 
                        verbose = TRUE, only.pos = F)
head(ai_markers,n=20)
ai_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)
write.csv(ai_markers,file='ai_markers_withbar5.csv')

ai_markers=read.csv('ai_markers.csv',row.names = 1)
ai_markers_top = ai_markers %>% filter(avg_log2FC>0.5 & p_val_adj<10e-20)
common.genes=intersect(top20$gene,rownames(ai_markers_top))
common.genes

pfc_exn_markers_wobarcodes %>%
  filter(cluster=='L6-Npy') %>%
  top_n(n = 20, wt = avg_log2FC) -> top20_layer23
common.genes.l23=intersect(top20_layer23$gene,rownames(ai_markers_top))
common.genes.l23

pdf('AI-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(ai_markers,
                     lab = rownames(ai_markers),
                     selectLab=c('Sparcl1','Nrip3','Foxp2','Sema5a',
                                 'Nptxr','Cck','C1ql3','S100b','Slc24a3',"Etv1"
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
                     pCutoff = 10e-20,
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
md_markers=FindMarkers(exn.sub,ident.1="MD_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(md_markers,n=20)
md_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)
write.csv(md_markers,file='md_markers.csv')

md_markers=read.csv('md_markers.csv',row.names = 1)
md_markers_top = md_markers %>% filter(avg_log2FC>0.5 & p_val_adj<10e-20)
common.genes=intersect(top20$gene,rownames(md_markers_top))
common.genes

pdf('MD-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(md_markers,
                     lab = rownames(md_markers),
                     selectLab=c('Syt6','Foxp2','Crym','Rprm','Cyr61',
                                 'Hs3st4','Sla',
                                 'Nptxr','Cck','C1ql3','Calb1','Npy2r'),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "MD-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 10e-20,
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
lh_markers=FindMarkers(exn.sub,ident.1="LH_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(lh_markers,n=20)
lh_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)
write.csv(lh_markers,file='lh_markers.csv')

lh_markers=read.csv('lh_markers.csv',row.names = 1)
lh_markers_top = lh_markers %>% filter(avg_log2FC>0.5 & p_val_adj<10e-20)
common.genes=intersect(top20$gene,rownames(lh_markers_top))
common.genes

pdf('LH-projecting DEGs volcano.pdf', height = 8,width = 10)
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
                     pCutoff = 10e-20,
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
bla_markers=FindMarkers(exn.sub,ident.1="BLA_valid",ident.2="Others",
                       logfc.threshold = 0.25, min.pct = 0.1,
                       test.use = "MAST", 
                       verbose = TRUE, only.pos = F)
head(bla_markers,n=20)
bla_markers %>% filter(avg_log2FC<(-0.5) & p_val_adj<10e-20)

write.csv(bla_markers,file='bla_markers.csv')
bla_markers=read.csv('bla_markers.csv',row.names = 1)
bla_markers_top = bla_markers %>% filter(avg_log2FC>0.5 & p_val_adj<10e-20)
common.genes=intersect(top20$gene,rownames(bla_markers_top))
common.genes


pdf('BLA-projecting DEGs volcano.pdf', height = 8,width = 10)
plot(EnhancedVolcano(bla_markers,
                     lab = rownames(bla_markers),
                     selectLab=c('Foxp2','Rprm','Tle4',
                                 'C1ql3','Nptxr','Cnr1','Syt17',
                                 'Rorb','Slc24a3',''),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-3,3),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = "BLA-projecting DEGs",
                     xlab ="",
                     subtitle = NULL,
                     pCutoff = 10e-20,
                     FCcutoff = 0.5,
                     pointSize = 1,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     gridlines.major = FALSE, gridlines.minor = FALSE,
                     labSize = 6)) + theme(legend.position = 'none')
dev.off()


lh_valid_cells=rownames(meta[meta$LH_valid=='LH_valid',])
md_valid_cells=rownames(meta[meta$MD_valid=='MD_valid',])

common.cells=intersect(lh_valid_cells,md_valid_cells)





