


library(ggplot2)
library(ggpubr)
library(readxl)
library(wesanderson)
library(ggbeeswarm)
library(tidyverse)


### LH ###
setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/pacaret_new/lh_results")

num_vars=c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial=c(1:100)
acc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(acc_df)=c('Acc','nVar')
auc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(auc_df)=c('AUC','nVar')
F1_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(F1_df)=c('F1','nVar')

acc_df$nVar=rep(num_vars,length(sample_trial))
auc_df$nVar=rep(num_vars,length(sample_trial))
F1_df$nVar=rep(num_vars,length(sample_trial))

plotlist=vector('list',3)

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    acc_var <- read.csv(sprintf("LH_var_model_results_%s_genes_Acc_%s.csv",n,i),row.names = 1)
    acc_df[j,1]=acc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    auc_var <- read.csv(sprintf("LH_var_model_results_%s_genes_AUC_%s.csv",n,i),row.names = 1)
    auc_df[j,1]=auc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    F1_var <- read.csv(sprintf("LH_var_model_results_%s_genes_F1_%s.csv",n,i),row.names = 1)
    F1_df[j,1]=F1_var[1,1]
    j=j+1
  }
}

p1=ggline(acc_df, x = "nVar", y = "Acc", color= wes_palette("Rushmore1")[4],add = c("mean_sd"),ylim=c(0,1))
p1
plotlist[[1]]=p1

p2=ggline(auc_df, x = "nVar", y = "AUC", color= wes_palette("Rushmore1")[3],add = c("mean_sd"),ylim=c(0,1))
p2
plotlist[[2]]=p2

p3=ggline(F1_df, x = "nVar", y = "F1", color= wes_palette("Rushmore1")[5],add = c("mean_sd"))
p3
plotlist[[3]]=p3

pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=1, align="v")
cowplot::save_plot("lh_metrics.pdf", pg,base_width=7,base_height = 10)

### random genes model metrics ###

sample_trial=c(1:100)
acc_ran_df= data.frame(matrix(ncol = 2,nrow = length(sample_trial)))
auc_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
F1_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
colnames(acc_ran_df)=c('Acc','group')
colnames(auc_ran_df)=c('AUC','group')
colnames(F1_ran_df)=c('F1','group')

acc_ran_df$group=rep('ran',length(sample_trial))
auc_ran_df$group=rep('ran',length(sample_trial))
F1_ran_df$group=rep('ran',length(sample_trial))

acc_var_df=acc_df %>% filter(nVar=='20')
acc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(acc_var_df)[2]='group'
auc_var_df=auc_df %>% filter(nVar=='20')
auc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(auc_var_df)[2]='group'
F1_var_df=F1_df %>% filter(nVar=='20')
F1_var_df$nVar=rep('hvg',length(sample_trial))
colnames(F1_var_df)[2]='group'

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/lh_ran_results")

for (i in sample_trial) {
    temp <- read.csv(sprintf("LH_ran_model_results_20_genes_Acc_%s.csv",i),row.names = 1)
    acc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("LH_ran_model_results_20_genes_AUC_%s.csv",i),row.names = 1)
  auc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("LH_ran_model_results_20_genes_F1_%s.csv",i),row.names = 1)
  F1_ran_df[i,1]=temp[1,1]
}

combined_lh_ran_acc_df=rbind(acc_ran_df,acc_var_df)
ggbarplot(combined_lh_ran_acc_df, x = "group", y = "Acc", color = 'group', palette =c("nrc"),
       add = c("mean_sd"))

combined_lh_ran_auc_df=rbind(auc_ran_df,auc_var_df)
ggbarplot(combined_lh_ran_auc_df, x = "group", y = "AUC", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_lh_ran_f1_df=rbind(F1_ran_df,F1_var_df)
ggbarplot(combined_lh_ran_f1_df, x = "group", y = "F1", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))


### DMS #####
setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/pacaret_new/dms_results")

num_vars=c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial=c(1:100)
acc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(acc_df)=c('Acc','nVar')
auc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(auc_df)=c('AUC','nVar')
F1_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(F1_df)=c('F1','nVar')

acc_df$nVar=rep(num_vars,length(sample_trial))
auc_df$nVar=rep(num_vars,length(sample_trial))
F1_df$nVar=rep(num_vars,length(sample_trial))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    acc_var <- read.csv(sprintf("DMS_var_model_results_%s_genes_Acc_%s.csv",n,i),row.names = 1)
    acc_df[j,1]=acc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    auc_var <- read.csv(sprintf("DMS_var_model_results_%s_genes_AUC_%s.csv",n,i),row.names = 1)
    auc_df[j,1]=auc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    F1_var <- read.csv(sprintf("DMS_var_model_results_%s_genes_F1_%s.csv",n,i),row.names = 1)
    F1_df[j,1]=F1_var[1,1]
    j=j+1
  }
}

p1=ggline(acc_df, x = "nVar", y = "Acc", color= wes_palette("Rushmore1")[4],add = c("mean_sd"),ylim=c(0,1))
plotlist[[1]]=p1
p1

p2=ggline(auc_df, x = "nVar", y = "AUC", color= wes_palette("Rushmore1")[3],add = c("mean_sd"),ylim=c(0,1))
p2
plotlist[[2]]=p2

p3=ggline(F1_df, x = "nVar", y = "F1", color= wes_palette("Rushmore1")[5],add = c("mean_sd"))
p3
plotlist[[3]]=p3

pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=1, align="v")
cowplot::save_plot("dms_metrics.pdf", pg,base_width=7,base_height = 10)

### random genes ###
sample_trial=c(1:100)
acc_ran_df= data.frame(matrix(ncol = 2,nrow = length(sample_trial)))
auc_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
F1_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
colnames(acc_ran_df)=c('Acc','group')
colnames(auc_ran_df)=c('AUC','group')
colnames(F1_ran_df)=c('F1','group')

acc_ran_df$group=rep('ran',length(sample_trial))
auc_ran_df$group=rep('ran',length(sample_trial))
F1_ran_df$group=rep('ran',length(sample_trial))

acc_var_df=acc_df %>% filter(nVar=='50')
acc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(acc_var_df)[2]='group'
auc_var_df=auc_df %>% filter(nVar=='50')
auc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(auc_var_df)[2]='group'
F1_var_df=F1_df %>% filter(nVar=='50')
F1_var_df$nVar=rep('hvg',length(sample_trial))
colnames(F1_var_df)[2]='group'

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/dms_ran_results")

for (i in sample_trial) {
  temp <- read.csv(sprintf("DMS_ran_model_results_50_genes_Acc_%s.csv",i),row.names = 1)
  acc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("DMS_ran_model_results_50_genes_AUC_%s.csv",i),row.names = 1)
  auc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("DMS_ran_model_results_50_genes_F1_%s.csv",i),row.names = 1)
  F1_ran_df[i,1]=temp[1,1]
}

combined_dms_ran_acc_df=rbind(acc_ran_df,acc_var_df)
ggbarplot(combined_lh_ran_acc_df, x = "group", y = "Acc", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_dms_ran_auc_df=rbind(auc_ran_df,auc_var_df)
ggbarplot(combined_lh_ran_auc_df, x = "group", y = "AUC", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_dms_ran_f1_df=rbind(F1_ran_df,F1_var_df)
ggbarplot(combined_lh_ran_f1_df, x = "group", y = "F1", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))


### AI ###

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/pacaret_new/ai_results")
num_vars=c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial=c(1:100)
acc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(acc_df)=c('Acc','nVar')
auc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(auc_df)=c('AUC','nVar')
F1_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(F1_df)=c('F1','nVar')

acc_df$nVar=rep(num_vars,length(sample_trial))
auc_df$nVar=rep(num_vars,length(sample_trial))
F1_df$nVar=rep(num_vars,length(sample_trial))


j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    acc_var <- read.csv(sprintf("AI_var_model_results_%s_genes_Acc_%s.csv",n,i),row.names = 1)
    acc_df[j,1]=acc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    auc_var <- read.csv(sprintf("AI_var_model_results_%s_genes_AUC_%s.csv",n,i),row.names = 1)
    auc_df[j,1]=auc_var[1,1]
    j=j+1
  }
}

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    F1_var <- read.csv(sprintf("AI_var_model_results_%s_genes_F1_%s.csv",n,i),row.names = 1)
    F1_df[j,1]=F1_var[1,1]
    j=j+1
  }
}

p1=ggline(acc_df, x = "nVar", y = "Acc", color= wes_palette("Rushmore1")[4],add = c("mean_sd"),ylim=c(0,1))
plotlist[[1]]=p1
p1

p2=ggline(auc_df, x = "nVar", y = "AUC", color= wes_palette("Rushmore1")[3],add = c("mean_sd"),ylim=c(0,1))
p2
plotlist[[2]]=p2

p3=ggline(F1_df, x = "nVar", y = "F1", color= wes_palette("Rushmore1")[5],add = c("mean_sd"))
p3
plotlist[[3]]=p3

pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=1, align="v")
cowplot::save_plot("ai_metrics.pdf", pg,base_width=7,base_height = 10)


### random genes model metrics ###

sample_trial=c(1:100)
acc_ran_df= data.frame(matrix(ncol = 2,nrow = length(sample_trial)))
auc_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
F1_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
colnames(acc_ran_df)=c('Acc','group')
colnames(auc_ran_df)=c('AUC','group')
colnames(F1_ran_df)=c('F1','group')

acc_ran_df$group=rep('ran',length(sample_trial))
auc_ran_df$group=rep('ran',length(sample_trial))
F1_ran_df$group=rep('ran',length(sample_trial))

acc_var_df=acc_df %>% filter(nVar=='20')
acc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(acc_var_df)[2]='group'
auc_var_df=auc_df %>% filter(nVar=='20')
auc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(auc_var_df)[2]='group'
F1_var_df=F1_df %>% filter(nVar=='20')
F1_var_df$nVar=rep('hvg',length(sample_trial))
colnames(F1_var_df)[2]='group'

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/ai_ran_results")

for (i in sample_trial) {
  temp <- read.csv(sprintf("AI_ran_model_results_20_genes_Acc_%s.csv",i),row.names = 1)
  acc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("AI_ran_model_results_20_genes_AUC_%s.csv",i),row.names = 1)
  auc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("AI_ran_model_results_20_genes_F1_%s.csv",i),row.names = 1)
  F1_ran_df[i,1]=temp[1,1]
}

combined_ai_ran_acc_df=rbind(acc_ran_df,acc_var_df)
ggbarplot(combined_ai_ran_acc_df, x = "group", y = "Acc", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_ai_ran_auc_df=rbind(auc_ran_df,auc_var_df)
ggbarplot(combined_ai_ran_auc_df, x = "group", y = "AUC", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_ai_ran_f1_df=rbind(F1_ran_df,F1_var_df)
ggbarplot(combined_ai_ran_f1_df, x = "group", y = "F1", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))


### MD ###
setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/pacaret_new/md_results")
num_vars=c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial=c(1:100)
acc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(acc_df)=c('Acc','nVar')
auc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(auc_df)=c('AUC','nVar')
F1_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(F1_df)=c('F1','nVar')

acc_df$nVar=rep(num_vars,length(sample_trial))
auc_df$nVar=rep(num_vars,length(sample_trial))
F1_df$nVar=rep(num_vars,length(sample_trial))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    acc_var <- read.csv(sprintf("MD_var_model_results_%s_genes_Acc_%s.csv",n,i),row.names = 1)
    acc_df[j,1]=acc_var[1,1]
    j=j+1
  }
}
ggline(acc_df, x = "nVar", y = "Acc", add = c("jitter","mean_sd"))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    auc_var <- read.csv(sprintf("MD_var_model_results_%s_genes_AUC_%s.csv",n,i),row.names = 1)
    auc_df[j,1]=auc_var[1,1]
    j=j+1
  }
}
ggline(auc_df, x = "nVar", y = "AUC", add = c("jitter","mean_sd"))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    F1_var <- read.csv(sprintf("MD_var_model_results_%s_genes_F1_%s.csv",n,i),row.names = 1)
    F1_df[j,1]=F1_var[1,1]
    j=j+1
  }
}
ggline(F1_df, x = "nVar", y = "F1", add = c("jitter","mean_sd"))

p1=ggline(acc_df, x = "nVar", y = "Acc", color= wes_palette("Rushmore1")[4],add = c("mean_sd"),ylim=c(0,1))
plotlist[[1]]=p1
p1

p2=ggline(auc_df, x = "nVar", y = "AUC", color= wes_palette("Rushmore1")[3],add = c("mean_sd"),ylim=c(0,1))
p2
plotlist[[2]]=p2

p3=ggline(F1_df, x = "nVar", y = "F1", color= wes_palette("Rushmore1")[5],add = c("mean_sd"))
p3
plotlist[[3]]=p3

pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=1, align="v")
cowplot::save_plot("md_metrics.pdf", pg,base_width=7,base_height = 10)

### random genes model metrics ###

sample_trial=c(1:100)
acc_ran_df= data.frame(matrix(ncol = 2,nrow = length(sample_trial)))
auc_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
F1_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
colnames(acc_ran_df)=c('Acc','group')
colnames(auc_ran_df)=c('AUC','group')
colnames(F1_ran_df)=c('F1','group')

acc_ran_df$group=rep('ran',length(sample_trial))
auc_ran_df$group=rep('ran',length(sample_trial))
F1_ran_df$group=rep('ran',length(sample_trial))

acc_var_df=acc_df %>% filter(nVar=='20')
acc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(acc_var_df)[2]='group'
auc_var_df=auc_df %>% filter(nVar=='20')
auc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(auc_var_df)[2]='group'
F1_var_df=F1_df %>% filter(nVar=='20')
F1_var_df$nVar=rep('hvg',length(sample_trial))
colnames(F1_var_df)[2]='group'

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/md_ran_results")

for (i in sample_trial) {
  temp <- read.csv(sprintf("MD_ran_model_results_20_genes_Acc_%s.csv",i),row.names = 1)
  acc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("MD_ran_model_results_20_genes_AUC_%s.csv",i),row.names = 1)
  auc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("MD_ran_model_results_20_genes_F1_%s.csv",i),row.names = 1)
  F1_ran_df[i,1]=temp[1,1]
}

combined_md_ran_acc_df=rbind(acc_ran_df,acc_var_df)
ggbarplot(combined_md_ran_acc_df, x = "group", y = "Acc", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_md_ran_auc_df=rbind(auc_ran_df,auc_var_df)
ggbarplot(combined_md_ran_auc_df, x = "group", y = "AUC", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_md_ran_f1_df=rbind(F1_ran_df,F1_var_df)
ggbarplot(combined_md_ran_f1_df, x = "group", y = "F1", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))


### BLA ###

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/pacaret_new/bla_results")

num_vars=c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial=c(1:100)
acc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(acc_df)=c('Acc','nVar')
auc_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(auc_df)=c('AUC','nVar')
F1_df= data.frame(matrix(ncol = 2, nrow = length(num_vars)*length(sample_trial)))
colnames(F1_df)=c('F1','nVar')

acc_df$nVar=rep(num_vars,length(sample_trial))
auc_df$nVar=rep(num_vars,length(sample_trial))
F1_df$nVar=rep(num_vars,length(sample_trial))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    acc_var <- read.csv(sprintf("BLA_var_model_results_%s_genes_Acc_%s.csv",n,i),row.names = 1)
    acc_df[j,1]=acc_var[1,1]
    j=j+1
  }
}
ggline(acc_df, x = "nVar", y = "Acc", add = c("jitter","mean_sd"))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    auc_var <- read.csv(sprintf("BLA_var_model_results_%s_genes_AUC_%s.csv",n,i),row.names = 1)
    auc_df[j,1]=auc_var[1,1]
    j=j+1
  }
}
ggline(auc_df, x = "nVar", y = "AUC", add = c("jitter","mean_sd"))

j=1
for (i in sample_trial) {
  for (n in num_vars ) {
    F1_var <- read.csv(sprintf("BLA_var_model_results_%s_genes_F1_%s.csv",n,i),row.names = 1)
    F1_df[j,1]=F1_var[1,1]
    j=j+1
  }
}
ggline(F1_df, x = "nVar", y = "F1", add = c("jitter","mean_sd"))

p1=ggline(acc_df, x = "nVar", y = "Acc", color= wes_palette("Rushmore1")[4],add = c("mean_sd"),ylim=c(0,1))
plotlist[[1]]=p1
p1

p2=ggline(auc_df, x = "nVar", y = "AUC", color= wes_palette("Rushmore1")[3],add = c("mean_sd"),ylim=c(0,1))
p2
plotlist[[2]]=p2

p3=ggline(F1_df, x = "nVar", y = "F1", color= wes_palette("Rushmore1")[5],add = c("mean_sd"))
p3
plotlist[[3]]=p3

pg=cowplot::plot_grid(plotlist=plotlist, labels = "auto",
                      ncol=1, align="v")
cowplot::save_plot("bla_metrics.pdf", pg,base_width=7,base_height = 10)

### random genes model metrics ###

sample_trial=c(1:100)
acc_ran_df= data.frame(matrix(ncol = 2,nrow = length(sample_trial)))
auc_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
F1_ran_df= data.frame(matrix(ncol = 2, nrow = length(sample_trial)))
colnames(acc_ran_df)=c('Acc','group')
colnames(auc_ran_df)=c('AUC','group')
colnames(F1_ran_df)=c('F1','group')

acc_ran_df$group=rep('ran',length(sample_trial))
auc_ran_df$group=rep('ran',length(sample_trial))
F1_ran_df$group=rep('ran',length(sample_trial))

acc_var_df=acc_df %>% filter(nVar=='20')
acc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(acc_var_df)[2]='group'
auc_var_df=auc_df %>% filter(nVar=='20')
auc_var_df$nVar=rep('hvg',length(sample_trial))
colnames(auc_var_df)[2]='group'
F1_var_df=F1_df %>% filter(nVar=='20')
F1_var_df$nVar=rep('hvg',length(sample_trial))
colnames(F1_var_df)[2]='group'

setwd("D:/Data/Single cell seq/mouse_pfc/manuscript/github/mpfc_projectome/10_machine_learning_implementation/bla_ran_results")

for (i in sample_trial) {
  temp <- read.csv(sprintf("BLA_ran_model_results_20_genes_Acc_%s.csv",i),row.names = 1)
  acc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("BLA_ran_model_results_20_genes_AUC_%s.csv",i),row.names = 1)
  auc_ran_df[i,1]=temp[1,1]
}
for (i in sample_trial) {
  temp <- read.csv(sprintf("BLA_ran_model_results_20_genes_F1_%s.csv",i),row.names = 1)
  F1_ran_df[i,1]=temp[1,1]
}

combined_bla_ran_acc_df=rbind(acc_ran_df,acc_var_df)
ggbarplot(combined_bla_ran_acc_df, x = "group", y = "Acc", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_bla_ran_auc_df=rbind(auc_ran_df,auc_var_df)
ggbarplot(combined_bla_ran_auc_df, x = "group", y = "AUC", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_bla_ran_f1_df=rbind(F1_ran_df,F1_var_df)
ggbarplot(combined_bla_ran_f1_df, x = "group", y = "F1", color = 'group', palette =c("nrc"),
          add = c("mean_sd"))

combined_ai_ran_acc_df$target=rep('AI',n=nrow(combined_ai_ran_acc_df))
combined_ai_ran_auc_df$target=rep('AI',n=nrow(combined_ai_ran_auc_df))
combined_ai_ran_f1_df$target=rep('AI',n=nrow(combined_ai_ran_f1_df))

combined_lh_ran_acc_df$target=rep('LH',n=nrow(combined_lh_ran_acc_df))
combined_lh_ran_auc_df$target=rep('LH',n=nrow(combined_lh_ran_auc_df))
combined_lh_ran_f1_df$target=rep('LH',n=nrow(combined_lh_ran_f1_df))

combined_dms_ran_acc_df$target=rep('DMS',n=nrow(combined_dms_ran_acc_df))
combined_dms_ran_auc_df$target=rep('DMS',n=nrow(combined_dms_ran_auc_df))
combined_dms_ran_f1_df$target=rep('DMS',n=nrow(combined_dms_ran_f1_df))

combined_bla_ran_acc_df$target=rep('BLA',n=nrow(combined_bla_ran_acc_df))
combined_bla_ran_auc_df$target=rep('BLA',n=nrow(combined_bla_ran_auc_df))
combined_bla_ran_f1_df$target=rep('BLA',n=nrow(combined_bla_ran_f1_df))

combined_md_ran_acc_df$target=rep('MD',n=nrow(combined_md_ran_acc_df))
combined_md_ran_auc_df$target=rep('MD',n=nrow(combined_md_ran_auc_df))
combined_md_ran_f1_df$target=rep('MD',n=nrow(combined_md_ran_f1_df))

combined_all_acc_df=rbind(combined_ai_ran_acc_df,combined_lh_ran_acc_df,combined_dms_ran_acc_df,
                          combined_bla_ran_acc_df,combined_md_ran_acc_df)
combined_all_auc_df=rbind(combined_ai_ran_auc_df,combined_lh_ran_auc_df,combined_dms_ran_auc_df,
                          combined_bla_ran_auc_df,combined_md_ran_auc_df)
combined_all_f1_df=rbind(combined_ai_ran_f1_df,combined_lh_ran_f1_df,combined_dms_ran_f1_df,
                          combined_bla_ran_f1_df,combined_md_ran_f1_df)

### plot ran vs hvg metrics ###
p_acc <- ggviolin(combined_all_acc_df, x = "target", y = "Acc",
                  add = c("boxplot","mean_sd"),
                  ylab='Accuracy',
                  xlab='',
                  size = 0.1,
                  fill = "group",
                  palette =c("nrc"),
                  ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "wilcox.test",paired = F,size=6)
p_acc
ggsave(p_acc,file='acc_all_violinplot.pdf',height = 5,width = 7)

p_auc <- ggviolin(combined_all_auc_df, x = "target", y = "AUC",
                  add = c("boxplot","mean_sd"),
                  ylab='AUC',
                  xlab='',
                  size = 0.1,
                  fill = "group",
                  palette =c("nrc"),
                  ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "wilcox.test",paired = F,size=6)
p_auc
ggsave(p_auc,file='auc_all_violinplot.pdf',height = 5,width = 7)

p_f1 <- ggviolin(combined_all_f1_df, x = "target", y = "F1",
                  add = c("boxplot","mean_sd"),
                  ylab='F1',
                  xlab='',
                  size = 0.1,
                  fill = "group",
                  palette =c("nrc"),
                  ylim = c(0, 0.8)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 0.8, aes(group = group,label = ..p.signif..),
                     method = "wilcox.test",paired = F,size=6)
p_f1
ggsave(p_f1,file='f1_all_violinplot.pdf',height = 5,width = 7)











