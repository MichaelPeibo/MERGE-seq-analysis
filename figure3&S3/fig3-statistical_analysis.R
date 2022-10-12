

setwd('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue')

library(tidyverse)
library(plyr)
library(UpSetR)
library(ggpubr)
library(EnhancedVolcano)
library(wesanderson)



load('sub_trim_seurat_pipe.Robj')
exn.sub=mda;rm(mda)
bar.mat=read.csv('bar.mat.drop.csv',row.names = 1)
AI_valid=rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid=rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid=rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid=rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid=rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]

unique_valid_cells=unique(c(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid))
upset_list=list(AI=AI_valid,DMS=DMS_valid,MD=MD_valid,BLA=BLA_valid,LH=LH_valid)


overlapGroups <- function (listInput, sort = TRUE) {
#https://github.com/hms-dbmi/UpSetR/issues/85
  listInputmat    <- fromList(listInput) == 1
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = "+")]] <- myelements
    myelements
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}


li <- overlapGroups(upset_list)
df <- ldply (li, data.frame)
colnames(df)=c('pattern','cellid')
temp_df=as.data.frame(table(df$pattern))
rownames(temp_df)=temp_df$Var1
temp_df[,1]=NULL
  
targets=c('AI','DMS','MD','BLA','LH')
df_pattern <- data.frame(matrix(ncol = 1, nrow = 2^length(targets)))
n=1
for (x in 1:5) {
  temp=combn(targets,x)
  print(temp)
  for (i in 1:ncol(temp)) {
    rownames(df_pattern)[n]=str_c(temp[,i],collapse = '+')
    n=n+1
    #print(n)
  }
}
rownames(df_pattern)
colnames(df_pattern)='Freq'
df_pattern=df_pattern %>% rownames_to_column
temp_df=temp_df %>% rownames_to_column
final_df=left_join(df_pattern,temp_df,by = 'rowname')
final_df[is.na(final_df)] <- 0
final_df[,'Freq.x']=NULL
colnames(final_df)[2]='observed'
write.csv(final_df,file='obeserved_counts.csv')

##calculated based on statistical_analysis_figure3.ipynb
res = numeric(5)
probs = numeric(5)
prob_edge=0.205825443266392
n=4653

for (i in 1:5 ) {
  e1=i
  e2=5-i
  p = (prob_edge ^ e1) * (1 - prob_edge) ^ e2
  print(e1,e2,p)
  exp = as.numeric(p) * n
  res[i] = exp
  probs[i] = p
}
res
probs

exp_df=read.csv('D:/Data/Single cell seq/mouse_pfc/barcode expression enriched from cDNA/barcode matrix/luyue/statistical_analysis/obs_exp_counts.csv',row.names = 1)
rownames(exp_df)=exp_df[,1]
exp_df[,1]=NULL

library(rstatix)
z=vector('character',nrow(exp_df))
effec.size=vector('character',nrow(exp_df))

for (i in 1:nrow(exp_df)){
  if (i >=1 & i <=5) {
    z[i]=binom_test(exp_df$observed[i], 4653, p = probs[1])$p
  } else if ( i >5 & i <=15) {
    z[i]=binom_test(exp_df$observed[i], 4653, p = probs[2])$p
  } else if (i >15 & i <=25) {
    z[i]=binom_test(exp_df$observed[i], 4653, p = probs[3])$p
  } else if (i >25 & i <=30) {
    z[i]=binom_test(exp_df$observed[i], 4653, p = probs[4])$p
  } else {
    z[i]=binom_test(exp_df$observed[i], 4653, p = probs[5])$p
  }
}

exp_df[,'pvalue']=z
exp_df[,'padj']=p.adjust(z, method = 'bonferroni')
exp_df[,'effect_size']=log2((exp_df$observed+1)/(exp_df$expected+1))
exp_df$pvalue=as.numeric(exp_df$pvalue)
exp_df$effect_size=as.numeric(exp_df$effect_size)

write.csv(exp_df,file='exp_df_with_pvalue_effectsize.csv')

exp_df_sig=exp_df[which(exp_df$padj<0.001),]
exp_df_sig=exp_df_sig[order(desc(exp_df_sig$effect_size)),]
exp_df_nonsig=exp_df[which(exp_df$padj>0.001),]

exp_df_sig_nonsig=rbind(exp_df_sig,exp_df_nonsig)
exp_df_reformed=exp_df_sig_nonsig[order(desc(exp_df_sig_nonsig$effect_size)),][,1:2]
exp_df_reformed=exp_df_reformed %>% rownames_to_column() %>% gather('group','count',-rowname)

exp_df_reformed=exp_df_reformed %>% mutate(rowname=factor(rowname,levels=rownames(exp_df_sig_nonsig))) %>%
  mutate(group =factor(group ,levels=c('observed','expected')))

p=ggbarplot(exp_df_reformed, "rowname", "count",
          fill = "group", color = "group", palette = c(wes_palette("BottleRocket2")[2],'grey'),
          label = TRUE,
          position = position_dodge(0.9))+ rotate_x_text(45)
p
ggsave('barplot_obs_exp.pdf',width = 18,height = 7)

pdf('upset_all_group.pdf',width = 15)
UpSetR::upset(fromList(upset_list), empty.intersections = "on",
              order.by = "degree", sets.bar.color ='deep sky blue',
              point.size = 3.5, line.size = 2, 
              #nintersects=10,
              mainbar.y.label = "Cell Count",
              text.scale = c(3, 1.3, 1, 1, 2, 1.5)
)
dev.off()


# 
# nn_df=readxl::read_excel('./statistical_analysis/map-seq-pipeline-master/41593_2020_705_MOESM2_ESM.xlsx',
#                          sheet = 2)
# write.csv(nn_df,file='nn_df.csv')



patterns.counter=vector('character',4)
patterns.counter[1]=sum(exp_df$observed[1:5])
patterns.counter[2]=sum(exp_df$observed[6:15])
patterns.counter[3]=sum(exp_df$observed[16:25])
patterns.counter[4]=sum(exp_df$observed[26:31])

patterns.counter=as.numeric(patterns.counter)
names(patterns.counter)=c('1_target','2_targets','3_targets','>3_targets')
df_counter=patterns.counter %>% as_tibble()
df_counter=cbind(df_counter,rep(sum(df_counter$value),4))
colnames(df_counter)=c('counts','all')
rownames(df_counter)=c('1_target','2_targets','3_targets','>3_targets')
df=df_counter %>% mutate(ratio=counts/all)
df$patterns=c('1 target','2 targets','3 targets','>3 targets')
df$pielabels <- paste(df$patterns, paste(round(df$ratio * 100, 2), "%", sep = ""), sep = "\n")

levels(df$pielabels)=levels(df$patterns)
col_palette=c(wes_palette("BottleRocket2")[3],
              wes_palette("Darjeeling1")[2],
              'orange',wes_palette("IsleofDogs1")[1])

p=ggpie(
  df, "ratio", label = "pielabels", lab.pos = "out",legend='none',repel = TRUE,
  lab.font =  c(6, "bold", "black"), fill = "patterns", color = "black", palette = col_palette
)+ theme(text=element_text(face='bold',size=10))
#p=annotate_figure(p,top = text_grob(sprintf('%s',barcode_names[i]), face = "bold", size = 20))
p
ggsave('pie_targets_ratio.pdf')


# create custom key-value pairs for 'overpresented', 'underpresented', expression by effective size

keyvals <- ifelse(
  exp_df$effect_size < 0 & exp_df$padj < 0.001, wes_palette("Zissou1")[1],
  ifelse(exp_df$effect_size > 0 & exp_df$padj < 0.001, wes_palette("Rushmore1")[5],
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == wes_palette("Rushmore1")[5]] <- 'Overpresented'
names(keyvals)[keyvals == 'black'] <- 'nonsig'
names(keyvals)[keyvals == wes_palette("Zissou1")[1]] <- 'Underpresented'

pdf('volcano_under_over_represented_padj_v1.pdf',width = 23,height = 10)
plot(EnhancedVolcano(exp_df,
                     lab =rownames(exp_df),
                     x = "effect_size",
                     selectLab = rownames(exp_df)[which(names(keyvals) %in% c('Overpresented', 'Underpresented'))],
                     y = "padj",
                     title = NULL,
                     colCustom = keyvals,
                     ylab = bquote(~-Log[10]~adjusted~italic(P)),xlab = NULL,
                     pCutoff = 0.001,
                     FCcutoff = 0,
                     labSize = 4,
                     subtitle = NULL,
                     legendPosition = 'none',
                     xlim = c(-5, 5),
                     boxedLabels=F,
                     ylim = c(0, -log10(10e-80)),
                     colAlpha = 1,
                     hline=0.001,
                     legendLabSize = 10,
                     legendIconSize = 3.0,
                     border = "partial",
                     borderWidth = 1,
                     borderColour = "black",
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE,
                     colConnectors = "black"))
dev.off()

















