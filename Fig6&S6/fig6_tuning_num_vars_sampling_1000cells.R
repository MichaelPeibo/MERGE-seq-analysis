
## this operation is done in pmacs ##
## bsub -Is bash ##
## cd /project/jcreminslab/peibo_projects/shap ##
## micromamba activate shap ##

num_vars <- c(2,5,10,20,50,100,200,300,400,500,1000,2000,5000)
sample_trial <- c(1:100)
valid_names <- c("AI", "DMS", "LH", "MD", "BLA")

for (j in valid_names) {
  combined_df <- data.frame()
  for (i in sample_trial) {
  for (n in num_vars ) {
    metric_df <- read.csv(sprintf("./%s_results_nb/%s_var_model_results_%s_genes_metrics_%s.csv",j,j,n,i))
    metric_df$num_vars <- n
    combined_df <- rbind(combined_df, metric_df)
  }
}
  write.csv(combined_df, sprintf("combined_df_%s.csv",j), row.names = FALSE)
}

### random sampling results ###
sample_trial <- c(1:100)
valid_names <- c("AI", "DMS", "LH", "MD", "BLA")

for (j in valid_names) {
  combined_df <- data.frame()
  for (i in sample_trial) {
    metric_df <- read.csv(sprintf("./%s_results_ran_1000cells/%s_ran_model_results_genes_metrics_%s.csv",j,j,i))
    combined_df <- rbind(combined_df, metric_df)
}
  write.csv(combined_df, sprintf("combined_df_ran_%s.csv",j), row.names = FALSE)
}

################################################
################################################
############## in local computer ###########################
################################################
################################################

setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")
library(ggplot2)
library(readxl)
library(wesanderson)
library(ggbeeswarm)
library(tidyverse)
library(rstatix)
library(coin)

auc_color <- wes_palette("Rushmore1")[4]
acc_color <- wes_palette("Rushmore1")[3]
f1_color <- wes_palette("Rushmore1")[5]
valid_names <- c("AI", "DMS", "LH", "MD", "BLA")

for (j in valid_names) {
  plotlist=vector('list',3)
  combined_df <- read.csv(sprintf("./results/combined_df_%s.csv",j))
  combined_df$num_vars <- as.factor(combined_df$num_vars)
  # Calculate the mean and standard deviation for AUC
  summary_df <- combined_df %>%
    group_by(num_vars) %>%
    summarise(mean_AUC = mean(AUC),
              sd_AUC = sd(AUC),
              se_AUC = sd_AUC / sqrt(n()))

  # Plot the line plot with mean and standard deviation
  p1 <- ggplot(summary_df, aes(x = num_vars, y = mean_AUC)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_AUC - sd_AUC, ymax = mean_AUC + sd_AUC),color = auc_color) +
    labs(x = "num_vars", y = "AUC") +
    geom_line(group = 1,color = auc_color)+
    ylim(0, 1 )+
    theme_Publication()+ 
    theme(axis.ticks.length = unit(0.5, "cm"))
  ggsave(sprintf("./results/fig6_auc_1000cells_nb_100trials_%s.pdf",j), width = 7, height = 5)
  plotlist[[1]] <- p1

  summary_df <- combined_df %>%
    group_by(num_vars) %>%
    summarise(mean_Acc = mean(Accuracy),
              sd_Acc = sd(Accuracy),
              se_Acc = sd_Acc / sqrt(n()))

  p2 <- ggplot(summary_df, aes(x = num_vars, y = mean_Acc)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_Acc - sd_Acc, ymax = mean_Acc + sd_Acc),color = acc_color) +
    labs(x = "num_vars", y = "Accuracy") +
    geom_line(group = 1,color = acc_color)+
    ylim(0, 1)+
    theme_Publication()+ 
    theme(axis.ticks.length = unit(0.5, "cm"))
  ggsave(sprintf("./results/fig6_acc_1000cells_nb_100trials_%s.pdf",j), width = 7, height = 5)
  plotlist[[2]] <- p2

  summary_df <- combined_df %>%
    group_by(num_vars) %>%
    summarise(mean_F1 = mean(F1),
              sd_F1 = sd(F1),
              se_F1 = sd_F1 / sqrt(n()))

  p3 <- ggplot(summary_df, aes(x = num_vars, y = mean_F1)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_F1 - sd_F1, ymax = mean_F1 + sd_F1),color = f1_color) +
    labs(x = "num_vars", y = "F1") +
    geom_line(group = 1,color = f1_color)+
    ylim(min(summary_df$mean_F1 - summary_df$sd_F1), max(summary_df$mean_F1 + summary_df$sd_F1))+
    theme_Publication()+ 
    theme(axis.ticks.length = unit(0.5, "cm"))
  ggsave(sprintf("./results/fig6_f1_1000cells_nb_100trials_%s.pdf",j), width = 7, height = 5)
  plotlist[[3]]=p3

  pg <- cowplot::plot_grid(plotlist=plotlist,
                      nrow=1, align="h")
  cowplot::save_plot(sprintf("./results/fig6_1000cells_nb_100trials_metrics_%s.pdf",j), pg,base_width=21,base_height = 4)
}

####################################
### vs random sampling results #####
####################################
library(ggplot2)
library(ggpubr)
library(ggsci)
library(readxl)
library(wesanderson)
library(ggbeeswarm)
library(tidyverse)
library(rstatix)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

valid_names <- c("AI", "DMS", "LH", "MD", "BLA")
combined_df <- data.frame()
for (j in valid_names) {
  combined_df_ran<- read.csv(sprintf("./results/fig6_1000cells_nb_100trials_metrics_ALL/combined_df_ran_%s.csv",j))
  combined_df_ran <- combined_df_ran %>% mutate(group = "ran") %>% mutate(target = j)
  combined_df_ori <- read.csv(sprintf("./results/fig6_1000cells_nb_100trials_metrics_ALL/combined_df_%s.csv",j))
  combined_df_ori <- combined_df_ori %>% filter(num_vars == 100) %>% mutate(group = "original") %>% mutate(target = j) %>% select(-num_vars)
  combined_df <- rbind(combined_df, combined_df_ran, combined_df_ori)
}
combined_df <- gather(combined_df, Metric, Value, -group, -target)
write.csv(combined_df, "./results/fig6_1000cells_nb_100trials_metrics_ALL/combined_df_var50_vs_random50.csv", row.names = FALSE)
### AUC ###
combined_df_auc <- combined_df %>% filter(Metric == "AUC")
stat.test <- combined_df_auc %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_auc_var_vs_ran_1000cells_nb_100trials-source data file.csv", row.names = FALSE)

combined_df_auc %>% 
  group_by(target) %>%
  wilcox_effsize(Value ~ group,ci = TRUE)

ggplot(combined_df_auc, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size = 0.1) +
  labs(title = "AUC",x = "",y = "AUC",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 1)+  
  ylim(0, 1 ) +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_auc_var_vs_ran_1000cells_nb_100trials.pdf", width = 7, height = 5) 

### ACC ###
combined_df_acc <- combined_df %>% filter(Metric == "Accuracy")
stat.test <- combined_df_acc %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_acc_var_vs_ran_1000cells_nb_100trials-source data file.csv", row.names = FALSE)

ggplot(combined_df_acc, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size = 0.1) +
  labs(title = "Accuracy",x = "",y = "Accuracy",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 1) + 
  ylim(0, 1 ) +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_acc_var_vs_ran_1000cells_nb_100trials.pdf", width = 7, height = 5)

### F1 ###
combined_df_f1 <- combined_df %>% filter(Metric == "F1")
stat.test <- combined_df_f1 %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_f1_var_vs_ran_1000cells_nb_100trials-source data file.csv", row.names = FALSE)

ggplot(combined_df_f1, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size=0.1) +
  labs(title = "F1",x = "",y = "F1",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 0.5) + 
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_f1_var_vs_ran_1000cells_nb_100trials.pdf", width = 7, height = 5)






### test code with AI ###
### NOT RUNNING ### 
combined_df <- read.csv("./results/combined_df_AI.csv")
combined_df$num_vars <- as.factor(combined_df$num_vars)
# Calculate the mean and standard deviation for AUC
summary_df <- combined_df %>%
  group_by(num_vars) %>%
  summarise(mean_AUC = mean(AUC),
            sd_AUC = sd(AUC),
            se_AUC = sd_AUC / sqrt(n()))

# Plot the line plot with mean and standard deviation
ggplot(summary_df, aes(x = num_vars, y = mean_AUC)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_AUC - sd_AUC, ymax = mean_AUC + sd_AUC),color = auc_color) +
  labs(x = "num_vars", y = "AUC") +
  geom_line(group = 1,color = auc_color)+
  ylim(0, 1 )+
  theme_Publication()

summary_df <- combined_df %>%
  group_by(num_vars) %>%
  summarise(mean_Acc = mean(Accuracy),
            sd_Acc = sd(Accuracy),
            se_Acc = sd_Acc / sqrt(n()))

ggplot(summary_df, aes(x = num_vars, y = mean_Acc)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_Acc - sd_Acc, ymax = mean_Acc + sd_Acc),color = acc_color) +
  labs(x = "num_vars", y = "Accuracy") +
  geom_line(group = 1,color = acc_color)+
  ylim(0, 1 )+
  theme_Publication()

summary_df <- combined_df %>%
  group_by(num_vars) %>%
  summarise(mean_F1 = mean(F1),
            sd_F1 = sd(F1),
            se_F1 = sd_F1 / sqrt(n()))

ggplot(summary_df, aes(x = num_vars, y = mean_F1)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_F1 - sd_F1, ymax = mean_F1 + sd_F1),color = f1_color) +
  labs(x = "num_vars", y = "F1") +
  geom_line(group = 1,color = f1_color)+
  ylim(min(summary_df$mean_F1 - summary_df$sd_F1), max(summary_df$mean_F1 + summary_df$sd_F1)) +
  #ylim(0, 0.5)+
  theme_Publication()













