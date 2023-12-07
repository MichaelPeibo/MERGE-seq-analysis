


### barcoded cell distribution in all celltypes
load('./data_file/pfc_all.seuratmerge.raw.bar5.Robj')
pfc_all <- mda; rm(mda)
pfc_all_meta <- read.csv('./data_file/pfc_seurat_merge_meta.csv',header = T,row.names = 1)
pfc_all <- subset(pfc_all, cells = rownames(pfc_all_meta))

bar.mat <- read.csv('./results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
## filter pfc_all by bar.mat cells ##
pfc_all <- subset(pfc_all, cells = rownames(bar.mat))
## filter pfc_all_meta by pfc_all cells ##
pfc_all_meta <- pfc_all_meta[colnames(pfc_all),]

AI_valid <- rownames(bar.mat[bar.mat$AI>0,])
DMS_valid <- rownames(bar.mat[bar.mat$DMS>0,])
MD_valid <- rownames(bar.mat[bar.mat$MD>0,])
BLA_valid <- rownames(bar.mat[bar.mat$BLA>0,])
LH_valid <- rownames(bar.mat[bar.mat$LH>0,])
unique_valid_cells <- unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid)); length(unique_valid_cells)

pfc_all <- SetIdent(pfc_all, cells = colnames(pfc_all), value = 'Non-barcoded')
pfc_all <- SetIdent(pfc_all, cells = unique_valid_cells, value = 'Barcoded')
pfc_all_meta$barcoded <- Idents(pfc_all)
table(Idents(pfc_all))
write.csv(pfc_all_meta,file = './results/pfc_all_meta_barcoded.csv')
pfc_all_meta <- read.csv('./results/pfc_all_meta_barcoded.csv',header = T,row.names = 1)
df <- pfc_all_meta %>% dplyr::select(facs,leiden_coarse,barcoded)
df <- df %>% 
  group_by(leiden_coarse,barcoded) %>% dplyr::count(name = "CellCount") %>% ungroup() %>%
  mutate(leiden_coarse = factor(leiden_coarse, 
                                levels = c('Excitatory','Microglia', 'Endo',
                                           'OPC','Oligo','Inhibitory',
                                           'Astro','Act.Microglia')))


# Create stacked bar graphs with labels
p <- ggbarplot(df, x = "leiden_coarse", y = "CellCount",
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
ggsave('./results/fig1_barcoded cell distribution.pdf',height = 10,width = 10)




