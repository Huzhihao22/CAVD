library(bestNormalize)
library(stringr)
library(parallel)
library(tidyverse)
library(dplyr)
orderNorm_gene <- function(gene.name, expr.mat)
{
  expr_norm <- expr.mat[,gene.name] %>% orderNorm()
  df <- data.frame(gene = expr_norm$x.t)
  colnames(df) <- gene.name  
  return(df)
}
psulist_name <- list.files("/path/pseudobulk/") 
psulist_file_list <- list()
for(i in psulist_name){
  j = str_split_fixed(i,"\\.",2)[,1]
  psulist_file <- read.csv(paste0("/path/pseudobulk/",i),row.names = 1)
  psulist_file_list[[j]] <- psulist_file
}
use_celltype = c()
for(i in names(psulist_file_list)){
  if(dim(psulist_file_list[[i]])[1] > 20) 
  use_celltype = c(use_celltype,i)
  
}

write.csv(use_celltype,file = '/path/use_celltype.csv',quote = F)

for (celltype in use_celltype) {
  
  pseudobulk_matrix_path1 <- paste0('/path/pseudobulk/', celltype, '.csv')
  
  pseudobulk_matrix_1 <- read.csv(pseudobulk_matrix_path1, row.names = 1)
  
  norm_list1 <- mclapply(colnames(pseudobulk_matrix_1), orderNorm_gene, expr.mat = pseudobulk_matrix_1, mc.cores = 30)
  
  expr1 <- bind_cols(norm_list1)
  row.names(expr1) <- row.names(pseudobulk_matrix_1)
  output_path1 <- paste0('/path/normal_dis/', celltype, '.csv')
  write.csv(expr1, file = output_path1)
  
  pseudobulk_matrix_path2 <- paste0('/path/top2000_pseudobulk/', celltype, '.csv')
  pseudobulk_matrix_2 <- read.csv(pseudobulk_matrix_path2, row.names = 1)
  norm_list2 <- mclapply(colnames(pseudobulk_matrix_2), orderNorm_gene, expr.mat = pseudobulk_matrix_2, mc.cores = 30)
  expr2 <- bind_cols(norm_list2)
  row.names(expr2) <- row.names(pseudobulk_matrix_2)
  output_path2 <- paste0('/path/top2000_normal_dis/', celltype, '.csv')
  write.csv(expr2, file = output_path2)
  cat("Processed and saved for celltype:", celltype, "\n")
}
