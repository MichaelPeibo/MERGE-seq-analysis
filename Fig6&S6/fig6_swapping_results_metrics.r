
## this operation is done in pmacs ##
## bsub -Is bash ##
## cd /project/jcreminslab/peibo_projects/shap ##
## micromamba activate shap ##

sample_trial <- c(1:100)
valid_names <- c("AI", "DMS", "LH", "MD", "BLA")

for (j in valid_names) {
  combined_df <- data.frame()
  for (i in sample_trial) {
    metric_df <- read.csv(sprintf("./%s_results_swapping_1000cells/%s_var_model_results_genes_metrics_%s.csv",j,j,i))
    combined_df <- rbind(combined_df, metric_df)
  }
  write.csv(combined_df, sprintf("combined_df_swapping_%s.csv",j), row.names = FALSE)
}

################################################
################################################
############## in local computer ###########################
################################################
################################################

library(ggplot2)
library(readxl)
library(wesanderson)
library(ggbeeswarm)
library(tidyverse)
library(ggsci)
library(rstatix)
library(ggpubr)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

valid_names <- c("AI", "DMS", "LH", "MD", "BLA")
combined_df <- data.frame()
for (j in valid_names) {
  combined_df_swapping <- read.csv(sprintf("./results/fig6_swapping_results/combined_df_swapping_%s.csv",j))
  combined_df_swapping <- combined_df_swapping %>% mutate(group = "swap") %>% mutate(target = j)
  combined_df_ori <- read.csv(sprintf("./results/fig6_1000cells_nb_100trials_metrics_ALL/combined_df_%s.csv",j))
  combined_df_ori <- combined_df_ori %>% filter(num_vars == 100) %>% mutate(group = "original") %>% mutate(target = j) %>% select(-num_vars)
  combined_df <- rbind(combined_df, combined_df_swapping, combined_df_ori)
}
combined_df <- gather(combined_df, Metric, Value, -group, -target)
write.csv(combined_df, "./results/fig6_swapping_results/combined_df_swapping_ALL.csv", row.names = FALSE)

### AUC ###
combined_df_auc <- combined_df %>% filter(Metric == "AUC")
stat.test <- combined_df_auc %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_auc_swapping_nb_100trials-source data file.csv", row.names = FALSE)
# combined_df_auc %>% 
#   group_by(target) %>%
#   wilcox_effsize(Value ~ group,ci = TRUE)
ggplot(combined_df_auc, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size = 0.1) +
  labs(title = "AUC",x = "",y = "AUC",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 1)+  
  ylim(0, 1 ) +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_auc_swapping_nb_100trials.pdf", width = 7, height = 5) 

### ACC ###
combined_df_acc <- combined_df %>% filter(Metric == "Accuracy")
stat.test <- combined_df_acc %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_acc_swapping_nb_100trials-source data file.csv", row.names = FALSE)

ggplot(combined_df_acc, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size = 0.1) +
  labs(title = "Accuracy",x = "",y = "Accuracy",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 1) + 
  ylim(0, 1 ) +
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_acc_swapping_nb_100trials.pdf", width = 7, height = 5)
### F1 ###
combined_df_f1 <- combined_df %>% filter(Metric == "F1")
stat.test <- combined_df_f1 %>% 
  group_by(target) %>%
  rstatix::wilcox_test(Value ~ group, detailed = TRUE) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, "./results/fig6_1000cells_nb_100trials_metrics_ALL/fig6_f1_swapping_nb_100trials-source data file.csv", row.names = FALSE)

ggplot(combined_df_f1, aes(x = target, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA,size=0.1) +
  labs(title = "F1",x = "",y = "F1",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 0.0001) +
  scale_color_nejm() + 
  stat_pvalue_manual(stat.test, x = "target", label = "p.adj = {p.adj} ", y.position = 0.6) + 
  theme_Publication() + 
  theme(axis.ticks.length = unit(0.5, "cm"))
ggsave("./results/fig6_f1_swapping_nb_100trials.pdf", width = 7, height = 5)




### test code with AI ###
### NOT RUNNING ### 
combined_df_swapping <- read.csv("./results/combined_df_swapping_AI.csv")
combined_df_swapping <- combined_df_swapping %>% mutate(group = "swap")
combined_df_ori <- read.csv("./results/combined_df_AI.csv")
combined_df_ori <- combined_df_ori %>% filter(num_vars == 100) %>% mutate(group = "original") %>% select(-num_vars)
combined_df <- rbind(combined_df_swapping, combined_df_ori)
long_df_auc <- combined_df %>%
  pivot_longer(cols = -group, names_to = "Metric", values_to = "Value") %>% filter(Metric == "AUC")

ggplot(long_df_auc, aes(x = group, y = Value, color = group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  labs(title = "AUC",x = "",y = "AUC",) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2),size = 1) +
  scale_color_nejm() + 
  stat_compare_means(method = "wilcox") + 
  theme_Publication() 


# Calculate the mean and standard deviation for AUC
summary_df <- combined_df %>%
  group_by(group) %>%
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


combined_df <- read.csv("./results/combined_df_swapping_AI.csv")
head(combined_df)
# Calculate the mean and standard deviation for AUC
summary_df <- combined_df %>%
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



acc_df <- read.csv("./swapping_results/acc_df.csv",row.names = 1)
head(acc_df)
auc_df <- read.csv("./swapping_results/auc_df.csv",row.names = 1)
head(auc_df)
### calculate mean acc_df$Acc ###
acc_df %>% 
  summarise(mean = mean(Acc, na.rm = TRUE), sd = sd(Acc, na.rm = TRUE))
### calculate mean auc_df$AUC ###
auc_df %>% 
  summarise(mean = mean(AUC, na.rm = TRUE), sd = sd(AUC, na.rm = TRUE))  

ggplot(acc_df, aes(x = Subclass, y = mean, fill = Subclass)) +
  geom_bar(stat="identity",fill=getPalette(colourCount)) +
  labs(x = "Dissection Region", y = "Input Reads") +
  ggtitle("Frequency by Variable") + 
  scale_fill_manual(values = getPalette(colourCount)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                 position=position_dodge(.9)) +
  theme_Publication()