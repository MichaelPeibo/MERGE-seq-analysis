

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
                         x = 'avg_log2FC',
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

temp_markers=read.csv('DMSvsDMS+MD.csv',row.names = 1)
temp_markers=read.csv('DMSvsDMS+LH.csv',row.names = 1)

png(sprintf('EnhancedVolcano %s vs %s .png',top_deg_list_names[i],diprojection_names[n]),
    units="px", width=1600, height=1600, res=300)

pdf('DMS+LH-projecting DEGs volcano.pdf', height = 8,width = 10)
EnhancedVolcano(temp_markers,
                lab = rownames(temp_markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #xlim = c(-3, 3),
                #ylim = c(0, -log10(10e-250)),
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                colConnectors = 'grey30',
                title = sprintf('DMS vs DMS+MD projection'),
                xlab = bquote(~Log[2]~ .(paste('FC(','DMS','/','DMS+MD',')',sep = ''))) ,
                subtitle = NULL,
                pCutoff = 10e-20,
                FCcutoff = 1,
                pointSize = 2,
                col=c('black', 'black', 'black', '#cc163a'),
                colAlpha=1,
                labSize = 4.0)+ theme(legend.position = 'none')
dev.off()



png('DMS+MD-projecting DEGs volcano.png',
    units="px", width=1600, height=1600, res=300)
EnhancedVolcano(temp_markers,
                lab = rownames(temp_markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                #xlim = c(-3, 3),
                #ylim = c(0, -log10(10e-250)),
                selectLab=c('Crym','Rprm','Hs3st4','Hsp90ab1'
                ),
                drawConnectors = TRUE,
                widthConnectors = 0.5, 
                colConnectors = 'grey30',
                title = sprintf('DMS vs DMS+MD projection'),
                xlab = bquote(~Log[2]~ .(paste('FC(','DMS','/','DMS+MD',')',sep = ''))) ,
                subtitle = NULL,
                pCutoff = 10e-20,
                FCcutoff = 1,
                pointSize = 2,
                col=c('black', 'black', 'black', '#cc163a'),
                colAlpha=1,
                labSize = 4.0)+ theme(legend.position = 'none')



pdf('DMS+LH-projecting DEGs volcano.pdf', height = 8,width = 10)
EnhancedVolcano(temp_markers,
                     lab = rownames(temp_markers),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     #xlim = c(-3, 3),
                     #ylim = c(0, -log10(10e-250)),
                selectLab=c('Pou3f1','Bcl11b','Igfbp4','Lypd1'
                ),
                     drawConnectors = TRUE,
                     widthConnectors = 0.5, 
                     colConnectors = 'grey30',
                     title = sprintf('DMS vs DMS+LH projection'),
                     xlab = bquote(~Log[2]~ .(paste('FC(','DMS','/','DMS+LH',')',sep = ''))) ,
                     subtitle = NULL,
                     pCutoff = 10e-20,
                     FCcutoff = 1,
                     pointSize = 2,
                     col=c('black', 'black', 'black', '#cc163a'),
                     colAlpha=1,
                     labSize = 4.0)+ theme(legend.position = 'none')
dev.off()


for (n in 1) {
  
  png('test.png',
      units="px", width=1600, height=1600, res=300)
  plot(EnhancedVolcano(temp_markers,
                       lab = rownames(temp_markers),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       #xlim = c(-3, 3),
                       #ylim = c(0, -log10(10e-250)),
                       selectLab=c('Crym','Rprm','Hs3st4','Hsp90ab1'
                       ),
                       drawConnectors = TRUE,
                       widthConnectors = 0.5, 
                       colConnectors = 'grey30',
                       title = sprintf('DMS vs DMS+MD projection'),
                       xlab = bquote(~Log[2]~ .(paste('FC(','DMS','/','DMS+MD',')',sep = ''))) ,
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





