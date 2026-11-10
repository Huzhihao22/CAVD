library(qvalue)
celltypes = list.files("/path/eQTL_result")
lead_perm_list <- list()
for(i in celltypes){
  perm_temp <- read.csv((paste0("/path/eQTL_result/",i,"/all_lead_perm.csv")))
  qvalue_obj <- qvalue(perm_temp$pval_beta)
  perm_temp$celltype_qvalue <- qvalue_obj$qvalues
  lead_perm_list[[i]] <- perm_temp
}
lead_perm_merge <- do.call(rbind,lead_perm_list)
lead_perm_merge$study_qvalue <- qvalue(lead_perm_merge$pval_beta)$qvalues
