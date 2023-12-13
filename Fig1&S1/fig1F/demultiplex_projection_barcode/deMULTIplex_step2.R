
library(Seurat)
library(ggplot2)
barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
bar.ref=c('CTGCACCGACGCATT','GAAGGCACAGACTTT','GTTGGCTGCAATCCA','AAGACGCCGTCGCAA','TATTCGGAGGACGAC')

pfc_4 <- Read10X(data.dir = "./pfc_4")
pfc_4 <- CreateSeuratObject(counts = pfc_4, project = "pfc_4", min.cells = 3, min.features = 200)
pfc_4
pfc_4_cells=colnames(pfc_4)
readTable_4=read.csv('readTable_pfc_4.csv',row.names=1)
pfc_4_cells=gsub("-1", "", pfc_4_cells)

pfc_1 <- Read10X(data.dir = "./pfc_1")
pfc_1 <- CreateSeuratObject(counts = pfc_1, project = "pfc_1", min.cells = 3, min.features = 200)
pfc_1
pfc_1_cells=colnames(pfc_1)
pfc_1_cells=gsub("-1", "", pfc_1_cells)
readTable_1=read.csv('readTable_pfc_1.csv',row.names=1)

pfc_2 <- Read10X(data.dir = "./pfc_2")
pfc_2 <- CreateSeuratObject(counts = pfc_2, project = "pfc_2", min.cells = 3, min.features = 200)
pfc_2
pfc_2_cells=colnames(pfc_2)
readTable_2=read.csv('readTable_pfc_2.csv',row.names=1)
pfc_2_cells=gsub("-1", "", pfc_2_cells)


pfc_3 <- Read10X(data.dir = "./pfc_3")
pfc_3 <- CreateSeuratObject(counts = pfc_3, project = "pfc_3", min.cells = 3, min.features = 200)
pfc_3
pfc_3_cells=colnames(pfc_3)
readTable_3=read.csv('readTable_pfc_3.csv',row.names=1)
pfc_3_cells=gsub("-1", "", pfc_3_cells)


source('/project/jcreminslab/peibo_projects/readtable/MULTIseq.Align.Suite.v2.R')
bar.table <- MULTIseq.align.v2(readTable_1, pfc_1_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_1_v2_hd2.csv")
bar.table <- MULTIseq.align.v2(readTable_4, pfc_4_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_4_v2_hd2.csv")
bar.table <- MULTIseq.align.v2(readTable_2, pfc_2_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_2_v2_hd2.csv")
bar.table <- MULTIseq.align.v2(readTable_3, pfc_3_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_3_v2_hd2.csv")

