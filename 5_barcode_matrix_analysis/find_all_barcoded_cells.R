


### barcoded cell distribution in all celltypes
load('pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all=mda
rm(mda)
pfc_all_meta=read.csv('pfc_seurat_merge_meta.csv',header=T,row.names=1)
pfc_all=subset(pfc_all,cells=rownames(pfc_all_meta))


load('sub_trim_seurat_pipe.Robj')
exn.sub=mda
bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)
AI_valid=rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid=rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid=rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
unique_valid_cells=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))

pfc_all <- SetIdent(pfc_all, cells = colnames(pfc_all), value = 'Unbarcoded')
pfc_all <- SetIdent(pfc_all, cells = unique_valid_cells, value = 'Barcoded')
pfc_all_meta$barcoded=Idents(pfc_all)
write.csv(pfc_all_meta,file='pfc_all_meta_barcoded.csv')
df = pfc_all_meta %>% select(facs,leiden_coarse,barcoded)
df=df %>% 
  group_by(leiden_coarse,barcoded) %>% dplyr::count(name = "CellCount") %>% ungroup() %>%
  mutate(leiden_coarse = factor(leiden_coarse, 
                                levels = c('Excitatory','Microglia', 'Endo',
                                           'OPC','Oligo','Inhibitory',
                                           'Astro','Act.Microglia')))
df
# Create stacked bar graphs with labels
p=ggbarplot(df, x = "leiden_coarse", y = "CellCount",
            color = "barcoded", fill = "barcoded",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",repel=T,lab.size = 7,
            xlab ="Cell type", 
            ylab = 'Cell Count', ggtheme=theme_pubclean()) + 
  rotate_x_text(45)+
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p
ggsave('barcoded cell distribution.pdf',height = 10,width = 10)

pdf('barcoded cell distribution.pdf',width = 10)
p
dev.off()




