

###--------------
### subtype distribution ###

pfc_all_meta_filter = pfc_all_meta[rownames(pfc_all_meta) %in% rownames(exn_meta),]

exn_meta$barcoded=pfc_all_meta_filter$barcoded
df = exn_meta %>% as_tibble() %>% select(cluster,barcoded)

df=df %>% 
  group_by(cluster,barcoded) %>% count(name = "CellCount") %>% ungroup()
df
# Create stacked bar graphs with labels
p=ggbarplot(df, x = "cluster", y = "CellCount",
            color = "barcoded", fill = "barcoded",
            palette = c("#0073C2FF", "grey"),
            label = TRUE, lab.pos = "in", lab.col = "black",lab.size = 7,
            xlab ="Layer",ylab = 'Cell Count',ggtheme=theme_pubclean()) + 
  rotate_x_text(45)+
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p

pdf('barcoded exn distribution.pdf')
p
dev.off()
