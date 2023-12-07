
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")

## y=0.999, umi egfp ==0, nonneuron and sort-neg
## barcode 0, ai, 28; 26
## barcode 1, dms, 83; 101
## barcode 2, md, 114; 61
## barcode 3, bla, 35; 25
## barcode 4, lh, 103; 72


### drop negative barcodes counts to zero
colnames(bar.mat) <- c('AI','DMS','MD','BLA','LH')
temp <- bar.mat
temp <- temp %>% mutate(AInew = case_when(AI <= 28 ~ 0, TRUE   ~ as.numeric(AI)),
                       DMSnew=case_when(DMS <= 101 ~ 0, TRUE ~ as.numeric(DMS)),
                       MDnew=case_when(MD <= 114 ~ 0, TRUE ~ as.numeric(MD)),
                       BLAnew=case_when(BLA <= 35 ~ 0, TRUE ~ as.numeric(BLA)),
                       LHnew=case_when(LH <= 103 ~ 0, TRUE ~ as.numeric(LH)))
head(temp, n = 10)
rownames(temp) <- rownames(bar.mat)
temp <- temp[, 6:10]
colnames(temp) <- c('AI','DMS','MD','BLA','LH')

bar.mat.drop <- temp

write.csv(bar.mat.drop, file =  './results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv')

bar.mat.drop <- read.csv('/Users/peiboxu/Desktop/merge-seq analysis/elife_revision/results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
df <- gather(bar.mat.drop, key = 'barcodes',value = 'UMI')
head(df)

ggplot(df, aes(x = barcodes, y = UMI)) +
  geom_violin(outlier.shape = NA) +  # Create violin plot without outliers
  geom_jitter(size = 0.1, width = 0.2) +  # Add quasirandom points
  scale_y_continuous(trans = 'log10', labels = scales::label_log())+
  labs(x = "Projection targets", y = "Counts") +
  theme_Publication(base_size = 12) +
  theme(axis.title.x = element_text(size = 20, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20, face = "bold"),  # Bold axis text
    legend.title = element_text(size = 20, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 20, face = "bold"))

ggsave("beeswarm violin log10 drop negative cells.png", 
       width = 8, height = 8, dpi=1000)

