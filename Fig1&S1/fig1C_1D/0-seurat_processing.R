### step 1 ###
# seurat processing of 3 unsort and 1 sort sample datasets

getwd()
setwd("seurat_processing")

library(dplyr)
library(Seurat)
library(cowplot)

# Load the dataset
pfc_1 <- Read10X(data.dir = "pfc_1")
pfc_1 <- CreateSeuratObject(counts = pfc_1, project = "pfc_1", min.cells = 3, min.features = 200)
pfc_1
pfc_1@meta.data$sample <- "pfc_1"
pfc_1@meta.data$facs <- "unsort"
pfc_1 <- RenameCells(pfc_1, add.cell.id = "pfc_1")

pfc_2 <- Read10X(data.dir = "pfc_2")
pfc_2 <- CreateSeuratObject(counts = pfc_2, project = "pfc_2", min.cells = 3, min.features = 200)
pfc_2
pfc_2@meta.data$sample <- "pfc_2"
pfc_2@meta.data$facs <- "unsort"
pfc_2 <- RenameCells(pfc_2, add.cell.id = "pfc_2")

pfc_3 <- Read10X(data.dir = "pfc_3")
pfc_3 <- CreateSeuratObject(counts = pfc_3, project = "pfc_3", min.cells = 3, min.features = 200)
pfc_3
pfc_3@meta.data$sample <- "pfc_3"
pfc_3@meta.data$facs <- "unsort"
pfc_3 <- RenameCells(pfc_3, add.cell.id = "pfc_3")

pfc_4 <- Read10X(data.dir = "pfc_4")
pfc_4 <- CreateSeuratObject(counts = pfc_4, project = "pfc_4", min.cells = 3, min.features = 200)
pfc_4
pfc_4@meta.data$sample <- "pfc_4"
pfc_4@meta.data$facs <- "sort"
pfc_4 <- RenameCells(pfc_4, add.cell.id = "pfc_4")

mda <- merge(x = pfc_1, y = pfc_2, project = "pfc_all")
mda <- merge(x = mda, y = pfc_3, project = "pfc_all")
mda <- merge(x = mda, y = pfc_4, project = "pfc_all")

save(mda,file='pfc_all.seuratmerge.raw.bar5.Robj')

mda[["percent.mt"]] <- PercentageFeatureSet(object = mda, pattern = "^mt-")
mda[["percent.rp"]] <- PercentageFeatureSet(object = mda, pattern = "^Rps|^Rpl")
head(x = mda@meta.data, 5)


mda <- subset(x = mda, subset=nFeature_RNA > 500 & nFeature_RNA < 8000 & nCount_RNA>1000 &
	nCount_RNA<60000 & percent.mt < 20)
mda <- NormalizeData(object = mda, normalization.method = "LogNormalize", scale.factor = 50000)
mda <- FindVariableFeatures(object = mda, selection.method = "vst", nfeatures = 2000)

seurat_var_genes=VariableFeatures(object = mda)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
seurat_var_genes=setdiff(seurat_var_genes,barcode_genes)
length(seurat_var_genes)
write.csv(seurat_var_genes,file='seurat_var_genes.csv')

mda <- ScaleData(object = mda, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt",'facs'))
mda <- RunPCA(object = mda, features = seurat_var_genes,nfeatures.print = 10,npcs = 100)
VariableFeatures(object = mda)=seurat_var_genes

## convert to loom
head(mda@meta.data)
meta_info=mda@meta.data
table(is.na(meta_info))
mda@graphs <- list()
mda.loom=as.loom(x = mda,filename='pfc_seu_merge_loom_bar5.loom')
mda.loom$close_all()
##################