
setwd("/Users/peiboxu/Desktop/merge-seq analysis/elife_revision")

library(tidyverse)
library(Seurat)
library(EnhancedVolcano)
library(scCustomize)
options(ggrepel.max.overlaps=Inf)
source("~/Desktop/Memory/analysis/ggplot_theme_Publication.R")

### calculate DEGs with of dediprojection vs collateral projection###
load('./data_file/sub_trim_seurat_pipe.Robj')
exn.sub <- mda; rm(mda)
exn.sub <- UpdateSeuratObject(exn.sub)
barcode_genes <- c('barcode0','barcode1','barcode2','barcode3','barcode4')
genes.use <- setdiff(rownames(exn.sub),barcode_genes)
exn.sub <- subset(exn.sub, features = genes.use)
exn_meta <- read.csv('./results/sub_exn_trim_meta_barcoded.csv',row.names = 1)
exn.sub@meta.data <- exn_meta
bar.mat <- read.csv('./results/bar.mat.drop_ecdf.0.999.filter_top.05.quantile.csv', row.names = 1)
AI_valid <- rownames(bar.mat[bar.mat$AI>0,])[rownames(bar.mat[bar.mat$AI>0,]) %in% colnames(exn.sub)]
DMS_valid <- rownames(bar.mat[bar.mat$DMS>0,])[rownames(bar.mat[bar.mat$DMS>0,]) %in% colnames(exn.sub) ]
MD_valid <- rownames(bar.mat[bar.mat$MD>0,])[rownames(bar.mat[bar.mat$MD>0,]) %in% colnames(exn.sub)]
BLA_valid <- rownames(bar.mat[bar.mat$BLA>0,])[rownames(bar.mat[bar.mat$BLA>0,]) %in% colnames(exn.sub)]
LH_valid <- rownames(bar.mat[bar.mat$LH>0,])[rownames(bar.mat[bar.mat$LH>0,]) %in% colnames(exn.sub)]
valid_cells_list <- list(AI_valid,DMS_valid,MD_valid,BLA_valid,LH_valid);
valid_cells_names <- c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')
for (i in 1:5) {
  exn.sub <- SetIdent(exn.sub, cells = colnames(exn.sub), value = 'Others')
  exn.sub <- SetIdent(exn.sub, cells = valid_cells_list[[i]], value = valid_cells_names[i])
  temp <- as.data.frame(Idents(exn.sub))
  colnames(temp) <- valid_cells_names[i]
  exn_meta <- cbind(exn_meta,temp)
}
exn.sub@meta.data <- exn_meta

### filter exn.sub by bar.mat cells ###
exn.sub <- subset(exn.sub, cells = rownames(bar.mat))


### dms ###
### filter rows of bar.mat based on the colnames of exn.sub ###
bar.mat.exb <- bar.mat[colnames(exn.sub),]

### extract exn.sub@meta.data$DMS_valid and bar.mat.exb$DMS and combine into a dataframe ###
dms_df <- data.frame(
    DMS_valid = exn.sub@meta.data$DMS_valid,
    DMS = bar.mat.exb$DMS, row.names = rownames(exn.sub@meta.data))
head(dms_df)

correlation <- cor(dms_df$DMS_valid, dms_df$DMS)

library(ggplot2)

# Convert "DMS_valid" to a factor
# Create a violin plot
plot <- ggplot(dms_df, aes(x = DMS_valid, y = DMS)) +
  geom_violin(trim = FALSE, scale = "width", fill = "#66C2A5") +
  geom_boxplot(width = 0.1, fill = "white", color = "#1F78B4") +
  theme_bw() +
  labs(x = "DMS_valid", y = "DMS", title = "Correlation between DMS_valid and DMS")

# Adjust plot appearance and labels for publication-level output
plot + theme(plot.title = element_text(size = 12, face = "bold"),
             axis.text = element_text(size = 10),
             axis.title = element_text(size = 11, face = "bold"))

# Assuming your dataframe is named "df" and the column names are "column1" and "column2"
library(ggplot2)

# Line plot
ggplot(dms_df, aes(x = DMS_valid, y = DMS)) +
  geom_line() +
  labs(x = "Column 1", y = "Column 2") +
  ggtitle("Line Plot")

# Dot plot
## transform dms_df to long format ##
dms_df_long <- gather(dms_df, key = "DMS_valid", value = "DMS")

ggplot(dms_df, aes(x = DMS_valid, y = DMS, color = DMS_valid)) +
  geom_point() +
  labs(x = "Column 1", y = "Column 2")

