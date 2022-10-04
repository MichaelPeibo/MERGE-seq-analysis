


DefaultAssay(object = pfc.bar.sub) <- "bar"
valid_cells_names=c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
bar_genes=c("bAI", "bDMS", "bMD",'bBLA','bLH')
n=1
top_deg_list_names <- c('AI','DMS','MD','BLA','LH')
top_deg_list <- vector("list", length(top_deg_list_names))
names(top_deg_list) <- top_deg_list_names

for (i in bar_genes){
  DefaultAssay(object = pfc.bar.sub) <- "bar"
  pfc.bar.sub = SetIdent(object = pfc.bar.sub, cells = colnames(pfc.bar.sub), value = 'Other')
  temp_expr <- FetchData(object = pfc.bar.sub, vars = i,slot='counts')
  temp_cells=colnames(pfc.bar.sub[, which(x = temp_expr > 0)])
  pfc.bar.sub <- SetIdent(object = pfc.bar.sub, cells = temp_cells, value = valid_cells_names[n])
  table(Idents(pfc.bar.sub))
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  temp_markers=FindMarkers(pfc.bar.sub,ident.1=valid_cells_names[n],
                           logfc.threshold = 0.25,
                           test.use = "MAST", min.pct = 0.1,
                           verbose = TRUE, only.pos = T)
  write.csv(temp_markers,file=paste(valid_cells_names[n],'.csv',sep = ''))
  top100 = rownames(temp_markers[1:100,])
  top_deg_list[[n]]=top100
  n=n+1
}


top_deg_list_names <- c('AI','DMS','MD','BLA','LH')
diprojection_names <- c('DMS+AI','MD+LH','DMS+LH','DMS+BLA','DMS+MD')
#volcano_deg_list <- vector("list", 10)
for (i in 1:5){
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  for (n in 1:5) {
    
    temp_markers=FindMarkers(pfc.bar.sub,ident.1=top_deg_list_names[i],ident.2=diprojection_names[n],
                             test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1,
                             verbose = TRUE, only.pos = F)
    write.csv(temp_markers,file=paste(top_deg_list_names[i],'vs',diprojection_names[n],'.csv',sep = ''))
    png(sprintf('EnhancedVolcano %s vs %s .png',top_deg_list_names[i],diprojection_names[n]),
        units="px", width=1600, height=1600, res=300)
    plot(EnhancedVolcano(temp_markers,
                         lab = rownames(temp_markers),
                         x = 'avg_logFC',
                         y = 'p_val_adj',
                         #xlim = c(-3, 3),
                         #ylim = c(0, -log10(10e-250)),
                         drawConnectors = TRUE,
                         widthConnectors = 0.5, 
                         colConnectors = 'grey30',
                         title = sprintf('%s vs %s projection',top_deg_list_names[i],diprojection_names[n]),
                         xlab = bquote(~Log[2]~ .(paste('FC(',top_deg_list_names[[i]],'/',diprojection_names[[n]],')',sep = ''))) ,
                         subtitle = NULL,
                         pCutoff = 10e-20,
                         FCcutoff = 1,
                         pointSize = 2,
                         col=c('black', 'black', 'black', '#cc163a'),
                         colAlpha=1,
                         labSize = 4.0)+ theme(legend.position = 'none')
    )
    dev.off()
  }
}

top_deg_list_names <- c('AI','DMS','MD','BLA','LH')

for (i in 1:5){
  DefaultAssay(object = pfc.bar.sub) <- "RNA"
  for (n in 1:5) {
    if(i==n){
      next
    }
    temp_markers=FindMarkers(pfc.bar.sub,ident.1=top_deg_list_names[i],ident.2=top_deg_list_names[n],
                             test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1,
                             verbose = TRUE, only.pos = F)
    write.csv(temp_markers,file=paste(top_deg_list_names[i],'vs',top_deg_list_names[n],'.csv',sep = ''))
    png(sprintf('EnhancedVolcano %s vs %s .png',top_deg_list_names[i],top_deg_list_names[n]),
        units="px", width=1600, height=1600, res=300)
    plot(EnhancedVolcano(temp_markers,
                         lab = rownames(temp_markers),
                         x = 'avg_logFC',
                         y = 'p_val_adj',
                         #xlim = c(-3, 3),
                         #ylim = c(0, -log10(10e-250)),
                         drawConnectors = TRUE,
                         widthConnectors = 0.5, 
                         colConnectors = 'grey30',
                         title = sprintf('%s vs %s dedicated projection',top_deg_list_names[i],top_deg_list_names[n]),
                         xlab = bquote(~Log[2]~ .(paste('FC(',top_deg_list_names[[i]],'/',top_deg_list_names[[n]],')',sep = ''))) ,
                         subtitle = NULL,
                         pCutoff = 10e-20,
                         FCcutoff = 1,
                         pointSize = 2,
                         col=c('black', 'black', 'black', '#cc163a'),
                         colAlpha=1,
                         labSize = 4.0)+ theme(legend.position = 'none')
    )
    dev.off()
  }
}
UpSetR::upset(UpSetR::fromList(top_deg_list), order.by = "freq")
unique_degs=unique(c(top_deg_list[[1]],top_deg_list[[2]],top_deg_list[[3]],
                     top_deg_list[[4]],top_deg_list[[5]]))

### ###
### heatmap of degs
bar.sub.mat=as.matrix(pfc.bar.sub@assays$RNA@scale.data)[unique_degs,]
col<- circlize::colorRamp2(c(-1,0,1,2,3),viridis::inferno(5))
bar.sub.mat=t(bar.sub.mat)

Heatmap(bar.sub.mat,name = "Scaled Expression", 
        #show_column_names=F,
        show_row_names=F,
        cluster_rows =hclust(dist(bar.sub.mat)),
        #row_order=as.character(unique(deg_genes)[drop=T]), 
        cluster_columns = hclust(dist(t(bar.sub.mat))),
        #top_annotation = HeatmapAnnotation(df=pheno,col=ann_colors),
        col = col
        #viridis::inferno(20)
        #col = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
)

pfc.bar.sub=subset(pfc.bar,idents = 'Other',invert=T)

DefaultAssay(object = pfc.bar.sub) <- "RNA"
table(Idents(pfc.bar.sub))
pfc.bar.sub <- ScaleData(object = pfc.bar.sub, 
                         #features=unique_degs,
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt",'facs'))

binary_cluster_markers=FindAllMarkers(pfc.bar.sub,logfc.threshold = 0.25,
                                      test.use = "MAST", min.pct = 0.1,
                                      verbose = TRUE, only.pos = T)
write.csv(binary_cluster_markers,file='binary_cluster_markers.csv')
# Re-level object@ident
my_levels <- c('AI','DMS','BLA','DMS+AI','DMS+BLA','DMS+MD','DMS+LH','MD','LH','MD+LH')

pfc.bar.sub@active.ident <- factor(x = pfc.bar.sub@active.ident, levels = my_levels)
top10 <- binary_cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10$gene,file='top10_degs.txt', sep = "",row.names=F,col.names=F,quote = FALSE)


bar.sub.mat=as.matrix(pfc.bar.sub@assays$RNA@scale.data)[unique(top10$gene),]
col<- circlize::colorRamp2(c(-1,0,1,2,3),viridis::inferno(5))

p=Heatmap(bar.sub.mat,name = "Scaled Expression", 
          show_column_names=F,
          show_row_names=T,
          cluster_rows =hclust(dist(bar.sub.mat)),
          #row_order=as.character(unique(deg_genes)[drop=T]), 
          cluster_columns = hclust(dist(t(bar.sub.mat))),
          #top_annotation = HeatmapAnnotation(df=pheno,col=ann_colors),
          col = col
          #viridis::inferno(20)
          #col = circlize::colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
)
p = draw(p)
row_order(p)
write.table(unique(top10$gene)[row_order(p)],file='top10_degs.txt', sep = "",row.names=F,col.names=F,quote = FALSE)
