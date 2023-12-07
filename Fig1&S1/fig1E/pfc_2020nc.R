library(Seurat)
library(Matrix)
library(tidyverse)
library(SeuratDisk)
library(wesanderson)
library(harmony)
setwd('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/')

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

load('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/pfc_all.seuratmerge.raw.bar5.Robj')
mda[["percent.mt"]] <- PercentageFeatureSet(object = mda, pattern = "^mt-")
mda[["percent.rp"]] <- PercentageFeatureSet(object = mda, pattern = "^Rps|^Rpl")

mda <- subset(x = mda, subset=nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA>1000 &
	nCount_RNA<60000 & percent.mt < 20)## pt.mt<10 will get 18035 samples

mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 50000)
meta=read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/pfc_seurat_merge_meta.csv')

mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2005)
barcode_genes <- c('barcode0','barcode1','barcode2','barcode3','barcode4')

load(file = "/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/all_cell_inDD_sobj.RData")
pfc <- all_cell_inDD_sobj
rm(all_cell_inDD_sobj)
table(pfc@meta.data$CellType)
Idents(pfc)=pfc@meta.data$CellType
pfc@meta.data$celltype=Idents(pfc)

pfc@meta.data$tech <- "nc2020"
mda@meta.data$tech <- "scrna"
mda@meta.data$celltype=Idents(mda)=meta$leiden_coarse

genes_select <- intersect(rownames(pfc), rownames(mda))
integrated <- merge(pfc[genes_select, ], mda[genes_select,])
integrated <- FindVariableFeatures(integrated) 
integrated <- ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunHarmony(object = integrated, group.by.vars = 'tech',theta = 3, lambda = 0.05,plot_convergence = TRUE)
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)
save(integrated,file='integrated_all_nc2020_major_celltype.Robj')
Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "celltype")
saveRDS(Correspondence_scRNA_merFISH, file = "Correspondence_scRNA_with_nc2020_majorcelltype.RDS")
Correspondence_scRNA_merFISH = readRDS(file = "/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/merfish/Correspondence_scRNA_with_nc2020_majorcelltype.RDS")
## change the column order Correspondence_scRNA_merFISH
Correspondence_scRNA_merFISH <- Correspondence_scRNA_merFISH[c(1, 2, 3,4,6,7,5,8),]

pdf("Correspondence_scRNA_nc2020_major_celltype.pdf")
pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()    


