
library(Seurat)
library(MetaNeighbor)
library(SummarizedExperiment)
library(wesanderson)
library(harmony)
library(scCustomize)
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/merfish")

pfc_merfish <- readRDS("PFC.MERFISH.rds")
pfc_merfish
head(pfc_merfish@meta.data)

pfc_merfish_sub <- subset(pfc_merfish, idents=c("L2/3 IT","L4/5 IT","L5 IT","L5 ET","L5/6 NP","L6 CT","L6 IT"))
table(pfc_merfish@meta.data$L3_cluster)
Idents(pfc_merfish_sub) <- pfc_merfish_sub@meta.data$L3_cluster
table(Idents(pfc_merfish_sub))
DimPlot_scCustom(seurat_object = pfc_merfish, figure_plot = TRUE, label.size = 6)


load("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/sub_trim_seurat_pipe.Robj")
pfc=mda;rm(mda)
pfc_meta=read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/data_file/sub_exn_trim_meta.csv',row.names = 1,header = T)
Idents(pfc)=pfc_meta$cluster
table(Idents(pfc))

head(pfc@meta.data)
head(pfc_merfish_sub@meta.data)
pfc_merfish_sub@meta.data$tech <- "merfish"
pfc@meta.data$tech <- "scrna"

pfc@meta.data$L3_cluster <- Idents(pfc)

genes_select <- intersect(rownames(pfc_merfish_sub), rownames(pfc))
integrated <- merge(pfc_merfish_sub[genes_select, ], pfc[genes_select,])
integrated <- FindVariableFeatures(integrated) 
integrated <- ScaleData(integrated, vars.to.regress = c( "nCount_RNA"), verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunHarmony(object = integrated, group.by.vars = 'tech',theta = 3, lambda = 0.05,plot_convergence = TRUE)
integrated <- RunUMAP(integrated, dims = 1:10, reduction = 'harmony')
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30,
                            k.param = 30, return.neighbor = T)

Correspondence_scRNA_merFISH = scRNA2merFISH_cluster(integrated, "L3_cluster")
saveRDS(Correspondence_scRNA_merFISH, file = "Correspondence_scRNA_merFISH.RDS")
Correspondence_scRNA_merFISH <- readRDS("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/code_file/merfish/Correspondence_scRNA_merFISH.RDS")
head(Correspondence_scRNA_merFISH)
## change the column order Correspondence_scRNA_merFISH
Correspondence_scRNA_merFISH <- Correspondence_scRNA_merFISH[, c(1, 2, 3, 5,4, 7, 6)]

pdf("Correspondence_scRNA_merFISH.pdf")
pheatmap::pheatmap(Correspondence_scRNA_merFISH, cluster_rows = F, cluster_cols = F,
                   color = colorRampPalette(c("white", "gray97", "red", "red3", "darkred"))(30), border_color = NA)
dev.off()    


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
