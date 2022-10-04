


pfc_all_meta=read.csv('pfc_all_meta_barcoded.csv',row.names = 1,header = T)
#exn_meta=read.csv('sub_exn_trim_meta.csv',row.names = 1,header = T)
exn_barcode_meta=pfc_all_meta[rownames(exn_meta),]
exn_barcode_meta=cbind(exn_barcode_meta,exn_meta$layer_v1,exn_meta$cluster_v1)
colnames(exn_barcode_meta)[15:16]=c('layer','cluster')
load('pfc.bar.sub.Robj')
temp=exn_meta[,c('cluster','layer_v1')]
meta = temp[colnames(pfc.bar.sub),]
meta = meta %>% add_column(Binary_Projection=Idents(pfc.bar.sub)) %>% 
  rownames_to_column() %>% as_tibble() %>%
  select(rowname,layer_v1,cluster,Binary_Projection) %>% mutate_if(is.integer, as.factor) %>%
  group_by(layer_v1,cluster,Binary_Projection) %>% srvyr::summarise(., count = n())
meta
meta$cluster=factor(meta$cluster,levels=c('L2/3-Calb1', 'L2/3-Rorb', 'L5-Bcl6', 
                                          'L5-Htr2c', 'L5-S100b', 'L6-Npy', 'L6-Syt6'))

ggplot(as.data.frame(meta),
       aes(axis1 = cluster, axis2 = Binary_Projection,
           y= count)) +
  scale_x_discrete(limits = c("Layer", "Projection"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = cluster),alpha=0.8,width =  0.3) +
  scale_fill_manual(values=c('#7fc97f', '#beaed4', '#fdc086', 
                             '#386cb0', '#f0027f', '#bf5b17', '#666666'))+
  #geom_flow()+
  geom_stratum(width = 0.3) + geom_text(stat = "stratum", infer.label = TRUE) +
  theme_void()+
  theme(legend.title=element_text(face='bold',size=10),
        legend.text=element_text(face = 'bold',size=15))+
  theme(legend.position = "bottom") +
  ggtitle("Binary projection pattern distribution in layer/subtype")
ggsave('alluvial plot cluster v1.pdf',height = 8,width = 6)

