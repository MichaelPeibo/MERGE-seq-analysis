


### drop negative barcodes counts to zero
colnames(bar.mat)=c('AI','DMS','MD','BLA','LH')
temp=bar.mat
temp = temp %>% mutate(AInew = case_when(AI <= 23 ~ 0,TRUE   ~ as.numeric(AI)),
                       DMSnew=case_when(DMS<=71 ~ 0, TRUE ~ as.numeric(DMS)),
                       MDnew=case_when(MD<=134 ~ 0,TRUE ~ as.numeric(MD)),
                       BLAnew=case_when(BLA<=9 ~ 0, TRUE ~ as.numeric(BLA)),
                       LHnew=case_when(LH<=125 ~ 0, TRUE ~ as.numeric(LH)))
head(temp,n=10)
rownames(temp)=rownames(bar.mat)
temp = temp[,6:10]
colnames(temp)=c('AI','DMS','MD','BLA','LH')

bar.mat.drop=temp
write.csv(bar.mat.drop,file='bar.mat.drop.csv')

df=gather(temp,key = 'barcodes',value='UMI')
head(df)

#pdf('test beeswarm violin log10 drop negative cells.pdf')
p=ggviolin(df,x='barcodes',y='UMI',outlier.shape = NA,
           xlab ="Projection targets", 
           ylab = parse(text = paste0(' ~ log[10] ~','Counts'))) + 
  yscale("log10", .format = TRUE) + geom_quasirandom(size=0.1) + 
  font("xlab", size = 20,face = "bold") +
  font("ylab", size = 20,face = "bold") +
  font("xy.text", size = 20,face = "bold") + 
  font("legend.title",size = 20,face = "bold") + 
  font("legend.text",face = "bold",size = 20)
p
ggsave("beeswarm violin log10 drop negative cells.png", 
       width = 8, height = 8, dpi=1000)

