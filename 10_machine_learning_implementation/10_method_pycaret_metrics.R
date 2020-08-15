

library(ggplot2)
library(ggpubr)
library(readxl)
library(wesanderson)
library(ggbeeswarm)

metrics <- read_excel("model_metrics_testdata.xlsx")


p1 <- ggviolin(metrics, x = "projection", y = "acc",
               add = c("mean_sd"),#,'mean_range''jitter',
               ylab='Prediction accuracy',
               xlab='',
               width=0,
               size = 0.5,
               color='group',
               position = position_dodge(0.5),
               #add.params = list(shape = "Batch"),
               palette =wes_palette("Darjeeling1")[4:5],
               ylim = c(0.7, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p1
ggsave('violin prediction acc.pdf',height = 10,width = 10)

p2 <- ggviolin(metrics, x = "projection", y = "f1",
               add = c("mean_sd"),
               ylab='F1 score',
               xlab='',
               width=0,
               size = 0.5,
               color='group',
               #position = position_dodge(0.8),
               palette =wes_palette("Darjeeling1")[4:5],
               ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p2
ggsave('violin f1 score.pdf',height = 10,width = 10)

p3 <- ggviolin(metrics, x = "projection", y = "auc",
                add = c("mean_sd"),
                ylab='AUC',
                xlab='',
                width=0,
                size = 0.5,
                color='group',
                position = position_dodge(0.8),
                palette =wes_palette("Darjeeling1")[4:5],
                ylim = c(0.6, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p3
ggsave('violin auc.pdf',height = 10,width = 10)

p4 <- ggviolin(metrics, x = "projection", y = "recall",
                add = c("mean_sd"),
                ylab='Recall',
                xlab='',
                width=0,
                size = 0.5,
                color='group',
                position = position_dodge(0.8),
                palette =wes_palette("Darjeeling1")[4:5],
                ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p4
ggsave('violin recall.pdf',height = 10,width = 10)


p5 <- ggviolin(metrics, x = "projection", y = "precision",
                add = c("mean_sd"),
                ylab='Precision',
                xlab='',
                width=0,
                size = 0.5,
                color='group',
                position = position_dodge(0.8),
                palette =wes_palette("Darjeeling1")[4:5],
                ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p5
ggsave('violin precision.pdf',height = 10,width = 10)

p6 <- ggviolin(metrics, x = "projection", y = "kappa",
                add = c("mean_sd"),
                ylab='Kappa',
                xlab='',
                width=0,
                size = 0.5,
                color='group',
                position = position_dodge(0.8),
                palette =wes_palette("Darjeeling1")[4:5],
                ylim = c(0, 1)
) + theme_pubr() + labs_pubr() + 
  stat_compare_means(label.y = 1, aes(group = group,label = ..p.signif..),
                     method = "t.test",size=6)
p6
ggsave('violin kappa.pdf',height = 10,width = 10)





