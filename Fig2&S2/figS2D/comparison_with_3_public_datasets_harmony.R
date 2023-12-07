
library(Seurat)
library(tidyverse)
library(wesanderson)
library(harmony)
setwd('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/')

## function adapted from https://github.com/YiZhang-lab/PFC-MERFISH/blob/main/integrated_merFISH_scRNA.R
scRNA2merFISH_cluster = function(integrated_obj, cluster){

  integrated_scRNA = subset(integrated_obj, subset = tech == "scrna")
  integrated_merFISH = subset(integrated_obj, subset = tech != "scrna")
  
  map_scRNA_merFISH = function(obj_a, obj_b) {
    scRNA2merFISH = data.frame(matrix(nrow = 0, ncol = 3))
    for(i in unique(obj_a[[cluster]][,1])){
      scRNA_cell = colnames(obj_a)[obj_a[[cluster]] == i]
      scRNA_cell_Neighbor = TopNeighbors(integrated_obj@neighbors$RNA.nn, scRNA_cell, n=30)
      
      merFISH_cell = colnames(obj_b)
      scRNA_cell_Neighbor = scRNA_cell_Neighbor[scRNA_cell_Neighbor %in% merFISH_cell]
      
      merFISH_cell_cluster = obj_b[[cluster]][scRNA_cell_Neighbor, ]
      merFISH_cell_cluster =  as.data.frame(table(merFISH_cell_cluster))
      merFISH_cell_cluster$Freq = merFISH_cell_cluster$Freq / sum(merFISH_cell_cluster$Freq)
      merFISH_cell_cluster$scRNA_cluster = i
      scRNA2merFISH = rbind(scRNA2merFISH, merFISH_cell_cluster)
    }
    scRNA2merFISH = reshape2::dcast(scRNA2merFISH, merFISH_cell_cluster~scRNA_cluster, value.var = "Freq")
    rownames(scRNA2merFISH) = scRNA2merFISH$merFISH_cell_cluster
    scRNA2merFISH = scRNA2merFISH[, -1]
    scRNA2merFISH[is.na(scRNA2merFISH)] <- 0
    scRNA2merFISH
  }
  scRNA2merFISH = map_scRNA_merFISH(integrated_scRNA, integrated_merFISH)
  merFISH2scRNA = map_scRNA_merFISH(integrated_merFISH, integrated_scRNA)
  Correspondence_scRNA_merFISH = 0.5 * (merFISH2scRNA + t(scRNA2merFISH)[rownames(merFISH2scRNA), colnames(merFISH2scRNA)])
  #Correspondence_scRNA_merFISH = sqrt(merFISH2scRNA * t(scRNA2merFISH)[rownames(merFISH2scRNA), colnames(merFISH2scRNA)])
  Correspondence_scRNA_merFISH

} 

load("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/sub_trim_seurat_pipe.Robj")
pfc=mda;rm(mda)
pfc_meta=read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/sub_exn_trim_meta.csv',row.names = 1,header = T)
Idents(pfc)=pfc_meta$cluster
table(Idents(pfc))
pfc@meta.data$tech <- "scrna"
pfc@meta.data$L3_cluster <- Idents(pfc)

### dataset 1 ###
#### preprocess subset_pfc_ab_ss_pmacs.R in pmacs cluster####
load("pfc_ab_ss.Robj")
head(pfc_ab@meta.data)
table(pfc_ab@meta.data$region_label)
table(pfc_ab@meta.data$subclass_label)
Idents(pfc_ab) <- pfc_ab@meta.data$subclass_label
table(pfc_ab@meta.data$subclass_label)


pfc_ab <- subset(pfc_ab, idents = c("L2 IT ENTl","L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX","L5 PT CTX","L5/6 IT TPE-ENT", "L5/6 NP CTX", "L6 CT CTX", "L6 IT CTX", "L6b CTX"))
pfc_ab <- NormalizeData(object = pfc_ab, normalization.method = "LogNormalize", scale.factor = 50000)
pfc_ab <- FindVariableFeatures(object = pfc_ab, selection.method = "vst", nfeatures = 2000)
pfc_ab@meta.data$tech <- "scrna_ab"

genes_select <- intersect(rownames(pfc_merfish_sub), rownames(pfc))

integrated <- merge(pfc_ab[genes_select, ], pfc[genes_select,])
integrated <- FindVariableFeatures(integrated) 
integrated <- ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunHarmony(object = integrated, group.by.vars = 'tech',theta = 3, lambda = 0.05,plot_convergence = TRUE)
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)

Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "subclass_label")
saveRDS(Correspondence_scRNA_merFISH, file = "Correspondence_scRNA_with_ab_ss.RDS")
Correspondence_scRNA_merFISH <- readRDS("Correspondence_scRNA_with_ab_ss.RDS")
pdf("Correspondence_scRNA_with_ab_ss.pdf")
pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()
### dataset 2 ###
load('../data_file/GSE161936_Lui_Nguyen_Seuratobjects.RData')
pfc_cell <- Fig1Rbp4;rm(Fig1Rbp4)
pfc_cell <- UpdateSeuratObject(pfc_cell)
table(Idents(pfc_cell))
newids <- c('Cd44','Figf','Otof','Pld5',
           'Cxcr7','Npr3','Tshz2')
names(newids) <- levels(pfc_cell)
pfc_cell <- RenameIdents(pfc_cell, newids)
DefaultAssay(pfc_cell) <- 'RNA'
pfc_cell@meta.data$tech <- "cell2020"
pfc_cell@meta.data$L3_cluster <- Idents(pfc_cell)

genes_select <- intersect(rownames(pfc_cell), rownames(pfc))
integrated <- merge(pfc_cell[genes_select, ], pfc[genes_select,])
integrated <- FindVariableFeatures(integrated) 
integrated <- ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunHarmony(object = integrated, group.by.vars = 'tech',theta = 3, lambda = 0.05,plot_convergence = TRUE)
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)

Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "L3_cluster")
saveRDS(Correspondence_scRNA_merFISH, file = "Correspondence_scRNA_cell2020.RDS")
Correspondence_scRNA_merFISH <- readRDS("Correspondence_scRNA_cell2020.RDS")
pdf("Correspondence_scRNA_cell2020.pdf")
pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()

### dataset 3 ###
load(file = "../data_file/PFC_excitatory_saline_sobj.RData")
pfc_nc=PFC_excitatory_saline_sobj
rm(PFC_excitatory_saline_sobj)
pfc_nc_meta=read.csv('../data_file/PFC_exc_saline_meta.csv',row.names=1,header=T)
Idents(pfc_nc)=pfc_nc_meta$L1_cluster

pfc_nc@meta.data$tech <- "nc2020"
pfc_nc@meta.data$L3_cluster <- Idents(pfc_nc)

genes_select <- intersect(rownames(pfc_nc), rownames(pfc))
integrated <- merge(pfc_nc[genes_select, ], pfc[genes_select,])
integrated <- FindVariableFeatures(integrated) 
integrated <- ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunHarmony(object = integrated, group.by.vars = 'tech',theta = 3, lambda = 0.05,plot_convergence = TRUE)
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)

Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "L3_cluster")
saveRDS(Correspondence_scRNA_merFISH, file = "Correspondence_scRNA_nc2020.RDS")
Correspondence_scRNA_merFISH <- readRDS("Correspondence_scRNA_nc2020.RDS")
pdf("Correspondence_scRNA_nc2020.pdf")
pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()



