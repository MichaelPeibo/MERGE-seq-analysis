

## this operation is done in pmacs ##
## bsub -Is bash ##
## cd /project/jcreminslab/peibo_projects/shap ##
## micromamba activate shap ##

library(Seurat)
library(dplyr)
library(data.table)

load('./results/exn.sub_with_valid_projection.Robj')
head(exn.sub@meta.data)
exn_meta <- read.csv('./results/exn_meta_with_valid_projection.csv',row.names = 1)
exn_meta <- exn_meta[colnames(exn.sub),]

#### 1000 cells ####
valid_cells_names <- c('AI_valid','DMS_valid','MD_valid','BLA_valid','LH_valid')

for (i in valid_cells_names) {
  for (n in 1:100) {
    # Set the seed for reproducibility (optional)
    set.seed(n)
    random.cells <- sample(colnames(exn.sub),size=1000)
    temp=subset(exn.sub,cells = random.cells)
    temp=FindVariableFeatures(object = temp, selection.method = "vst", nfeatures = 100)
    var.genes=VariableFeatures(temp)
    exn_meta_temp <- exn_meta[colnames(temp),]
    # Get the indices of "AI_valid" and "Others" cells
    valid_indices <- which(exn_meta_temp[,i] == i)
    others_indices <- which(exn_meta_temp[,i] == "Others")
    # Randomly shuffle the indices within each category
    shuffled_valid <- sample(valid_indices)
    shuffled_others <- sample(others_indices)
    # Randomly swap cells between "AI_valid" and "others"
    num_swaps <- min(length(valid_indices), length(others_indices))
    for (j in 1:num_swaps) {
      exn_meta_temp[,i][shuffled_valid[j]] <- "Others"
      exn_meta_temp[,i][shuffled_others[j]] <- i
    }
    df=t(as.matrix(temp@assays$RNA@data)[var.genes,])
    df=cbind(df,as.data.frame(exn_meta_temp[,i]))
    colnames(df)[ncol(df)]=('binary')
    df$binary <- factor(df$binary, levels = c("Others",i))
    df=mutate(df, binary_clus = ifelse(binary == i , 1, 0))
    df$binary_clus <- as.factor(df$binary_clus)
    df <- subset(df, select = -c(binary))
    colnames(df)[101]='binary'
    fwrite(df,file=sprintf('./swapping_results/fig6_var100_%s_swapping_%s_1000cells.csv',i,n))
  }
}
