

### upset plot of binary cluster

#barcode_names=c('AI','DMS','MD','BLA','LH')
upset_list=list(AI=AI_valid,DMS=DMS_valid,MD=MD_valid,BLA=BLA_valid,LH=LH_valid)

pdf('upset valid cell.pdf')
UpSetR::upset(fromList(upset_list), 
              order.by = "freq", sets.bar.color ='deep sky blue',
              point.size = 3.5, line.size = 2, 
              nintersects=10,
              mainbar.y.label = "Cell Count",
              text.scale = c(3, 1.3, 1, 1, 2, 1.5),
              queries = list(list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('DMS'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('LH'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('BLA'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI','DMS'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','BLA'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','DMS'), color='orange',active = T))
              
)
dev.off()
pdf('upset valid cell all groups.pdf')
UpSetR::upset(fromList(upset_list), 
              order.by = "freq", sets.bar.color ='deep sky blue',
              point.size = 3.5, line.size = 2, 
              #nintersects=10,
              mainbar.y.label = "Cell Count",
              text.scale = c(3, 1.3, 1, 1, 2, 1.5),
              queries = list(list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('DMS'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('MD'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('LH'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('BLA'), color=wes_palette("Darjeeling1")[2],active = T),
                             list(query = intersects,
                                  params = list('AI','DMS'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','MD'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('DMS','BLA'), color='orange',active = T),
                             list(query = intersects,
                                  params = list('LH','DMS'), color='orange',active = T))
              
)
dev.off()
##