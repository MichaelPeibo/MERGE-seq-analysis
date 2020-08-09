

library(Seurat)
library(deMULTIplex)
library(ggplot2)

barcode_genes=c('barcode0','barcode1','barcode2','barcode3','barcode4')
bar.ref=c('CTGCACCGACGCATT','GAAGGCACAGACTTT','GTTGGCTGCAATCCA','AAGACGCCGTCGCAA','TATTCGGAGGACGAC')


pfc_1 <- Read10X(data.dir = "./data/pfc_1")
pfc_1 <- CreateSeuratObject(counts = pfc_1, project = "pfc_1", min.cells = 3, min.features = 200)
pfc_1
pfc_1_cells=colnames(pfc_1)
readTable <- MULTIseq.preProcess(R1 = "./rawdata/PFC-1_BKDL202572450-1a_1.clean.fq.gz",
                                 R2 = "./rawdata/PFC-1_BKDL202572450-1a_2.clean.fq.gz",
                                 cellIDs = pfc_1_cells,
                                 cell=c(1,16),
                                 umi=c(17,28),
                                 tag=c(31,45))
write.csv(readTable, file="readTable_pfc_1.csv")
readTable=read.csv('readTable_pfc_1.csv',row.names=1)
bar.table <- MULTIseq.align(readTable, pfc_1_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_1.csv")



pfc_2 <- Read10X(data.dir = "./data/pfc_2")
pfc_2 <- CreateSeuratObject(counts = pfc_2, project = "pfc_2", min.cells = 3, min.features = 200)
pfc_2
pfc_2_cells=colnames(pfc_2)

readTable <- MULTIseq.preProcess(R1 = "./rawdata/PFC-2_BKDL202572451-1a_1.clean.fq.gz",
                                 R2 = "./rawdata/PFC-2_BKDL202572451-1a_2.clean.fq.gz",
                                 cellIDs = pfc_2_cells,
                                 cell=c(1,16),
                                 umi=c(17,28),
                                 tag=c(31,45))
write.csv(readTable, file="readTable_pfc_2.csv")
readTable=read.csv('readTable_pfc_2.csv',row.names=1)
bar.table <- MULTIseq.align(readTable, pfc_2_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_2.csv")

pfc_3 <- Read10X(data.dir = "./data/pfc_3")
pfc_3 <- CreateSeuratObject(counts = pfc_3, project = "pfc_3", min.cells = 3, min.features = 200)
pfc_3
pfc_3_cells=colnames(pfc_3)

readTable <- MULTIseq.preProcess(R1 = "./rawdata/PFC-3_BKDL202572452-1a_1.clean.fq.gz",
                                 R2 = "./rawdata/PFC-3_BKDL202572452-1a_2.clean.fq.gz",
                                 cellIDs = pfc_3_cells,
                                 cell=c(1,16),
                                 umi=c(17,28),
                                 tag=c(31,45))
write.csv(readTable, file="readTable_pfc_3.csv")
bar.table <- MULTIseq.align(readTable, pfc_3_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_3.csv")

pfc_4 <- Read10X(data.dir = "./data/pfc_4")
pfc_4 <- CreateSeuratObject(counts = pfc_4, project = "pfc_4", min.cells = 3, min.features = 200)
pfc_4
pfc_4_cells=colnames(pfc_4)
readTable <- MULTIseq.preProcess(R1 = "./rawdata/PFC-4_BKDL202572453-1a_1.clean.fq.gz",
                                 R2 = "./rawdata/PFC-4_BKDL202572453-1a_2.clean.fq.gz",
                                 cellIDs = pfc_4_cells,
                                 cell=c(1,16),
                                 umi=c(17,28),
                                 tag=c(31,45))
write.csv(readTable, file="readTable_pfc_4.csv")
readTable=read.csv('readTable_pfc_4.csv',row.names=1)
bar.table <- MULTIseq.align(readTable, pfc_4_cells, bar.ref)
write.csv(bar.table, file="bar.table_pfc_4.csv")



