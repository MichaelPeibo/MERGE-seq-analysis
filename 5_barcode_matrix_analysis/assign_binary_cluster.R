

### ###
### binary cluster in seurat object
library(Seurat)
library(SeuratDisk)
exn.sub=LoadH5Seurat("exn.sub.h5Seurat")
bar.mat.drop=read.csv('bar.mat.drop.csv',row.names = 1)
barcode_genes=c('AI','DMS','MD','BLA','LH')
temp.bar=bar.mat.drop[colnames(exn.sub),]
colnames(temp.bar)=c("bAI", "bDMS", "bMD",'bBLA','bLH')

# Transform table
# Add barcode data as a new assay independent from RNA
exn.sub[["bar"]] <- CreateAssayObject(counts =as.data.frame(t(temp.bar)))
DefaultAssay(object = exn.sub) <- "bar"

exn.sub

AI=WhichCells(exn.sub, expression = bAI >0 & bDMS == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
BLA=WhichCells(exn.sub, expression = bBLA > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bLH == 0,
               slot='counts')
DMS=WhichCells(exn.sub, expression = bDMS > 0 & bAI == 0 & bMD == 0 & bBLA == 0 & bLH == 0,
               slot='counts')
MD=WhichCells(exn.sub, expression = bMD > 0 & bAI == 0 & bDMS == 0 & bBLA == 0 & bLH == 0,
              slot='counts')
LH=WhichCells(exn.sub, expression = bLH > 0 & bAI == 0 & bDMS == 0 & bMD == 0 & bBLA == 0,
              slot='counts')

AIDMS=WhichCells(exn.sub, expression = bAI > 0 & bDMS > 0 & bMD == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')
MDLH=WhichCells(exn.sub, expression = bMD > 0 & bLH > 0 & bAI == 0 & bDMS == 0 & bBLA == 0,
                slot='counts')
DMSLH=WhichCells(exn.sub, expression = bDMS > 0 & bLH > 0 & bAI == 0 & bMD == 0 & bBLA == 0,
                 slot='counts')
DMSBLA=WhichCells(exn.sub, expression = bDMS > 0 & bBLA > 0 & bAI == 0 & bMD == 0 & bLH == 0,
                  slot='counts')
DMSMD=WhichCells(exn.sub, expression = bDMS > 0 & bMD > 0 & bAI == 0 & bBLA == 0 & bLH == 0,
                 slot='counts')


exn.sub = SetIdent(object = exn.sub, cells = colnames(exn.sub), value = 'Other')
binary_cluster=c('AI','DMS','MD','BLA','LH','DMS+AI','MD+LH','DMS+LH','DMS+BLA','DMS+MD')
binary_cluster_cell_list=list(AI,DMS,MD,BLA,LH,AIDMS,MDLH,DMSLH,DMSBLA,DMSMD)

for (i in 1:10) {
  exn.sub <- SetIdent(object = exn.sub, cells = binary_cluster_cell_list[[i]], value = binary_cluster[i])
}
table(Idents(exn.sub))

exn.sub.sub=subset(exn.sub,idents = 'Other',invert=T)
DefaultAssay(object = exn.sub.sub) <- "RNA"
table(Idents(exn.sub.sub))
save(pfc.bar.sub,file='pfc.bar.sub.Robj')



