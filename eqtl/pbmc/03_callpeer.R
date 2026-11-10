library(peer)
mapvaluesss <- function (x, from, to, warn_missing = TRUE) 
{
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x) && !is.null(x)) {
    stop("`x` must be an atomic vector or NULL.")
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}
usecell <- c("B_intermediate","B_memory","B_naive","CD14_Mono","CD16_Mono","CD4_CTL","CD4_Naive","CD4_Proliferating",
             "CD4_TCM","CD4_TEM","CD8_Naive","CD8_TCM","CD8_TEM","cDC1","cDC2","dnT","gdT","HSPC","ILC","MAIT","NK_CD56bright",
             "NK_Proliferating","NK","pDC","Plasmablast","Treg")

blood_with_snp <- read.csv("/path/blood_with_snp.csv",row.names = 1)
blood_with_snp$name <- as.character(blood_with_snp$name)
blood_with_snp$snp_file <- as.character(blood_with_snp$snp_file)
snp_pc_mat <- as.data.frame(read.table("/path/snp.eigenvec"))
rownames(snp_pc_mat) <- snp_pc_mat$V2
snp_pc_mat <- snp_pc_mat[,-c(1,2)]
colnames(snp_pc_mat) <- paste0("PC_",seq(1,20,1))
snp_pc_mat_blood <- snp_pc_mat[rownames(snp_pc_mat) %in% blood_with_snp$snp_file,]
for(i in usecell){
  print(sprintf('%s STARTTTTTTTTTTTTT', i))
  expr_mat <- as.data.frame(read.csv(paste0("/path/top2000_normal_dis/",i,".csv"),row.names = 1))
  if("HUMwmfrN_C_B_3326543_2" %in% rownames(expr_mat)){
    expr_mat <- expr_mat[rownames(expr_mat) != "HUMwmfrN_C_B_3326543_2",]
  }
  cova_mat <- as.data.frame(read.csv("/path/blood_with_snp.csv",row.names = 1))
  rownames(cova_mat) <- cova_mat$name
  cosample <- intersect(rownames(expr_mat),rownames(cova_mat))
  expr_mat <- expr_mat[cosample,]
  cova_mat <- cova_mat[cosample,c("sex","age")]
  cova_mat$sex_adj <- 1
  cova_mat[cova_mat$sex == "female","sex_adj"] <- 0
  cova_mat = cova_mat[,c("sex_adj","age")]
  print(table(rownames(expr_mat) == rownames(cova_mat)))
  rownames(expr_mat) <- mapvaluesss(rownames(expr_mat),
                                    from = blood_with_snp$name,
                                    to = blood_with_snp$snp_file)
  rownames(cova_mat) <- mapvaluesss(rownames(cova_mat),
                                    from = blood_with_snp$name,
                                    to = blood_with_snp$snp_file)
  snp_pc_mat_use <- snp_pc_mat_blood[rownames(cova_mat),]
  print(table(rownames(snp_pc_mat_use) == rownames(cova_mat)))
  cova_mat <- cbind(cova_mat,snp_pc_mat_use[,c(1,2)])
  
  
  ##### PEER-------------------------------------------------------------------------
  # Set PEER paramaters based on the instructions from PEER package website
  model = PEER()
  PEER_setPhenoMean(model, as.matrix(expr_mat))
  PEER_setCovariates(model, as.matrix(cova_mat))
  dim(PEER_getPhenoMean(model))
  # Set to generate 10 PEER factors: K=15
  K <- 15
  PEER_setNk(model,K)
  PEER_getNk(model)
  PEER_setNmax_iterations(model, 5000)
  # Calculate
  PEER_update(model)
  ##### Observing output--------------------------------------------------------------
  pf.names <- paste('pf',1:K, sep = '')
  # factors: NxK matrix
  f <- PEER_getX(model)
  f_df <- data.frame(f)
  colnames(f_df) <- c(colnames(cova_mat),pf.names)
  f_df$sample <- rownames(expr_mat)
  f_df <- f_df[,c(ncol(f_df),1:(ncol(f_df)-1))]
  write.csv(f_df, sprintf('/path/peer/factor/%s.csv',i))
  
  # weights: GxK matrix
  w <- PEER_getW(model)
  w_df <- data.frame(w)
  colnames(w_df) <- c(colnames(cova_mat),pf.names)
  w_df$geneid <- colnames(expr_mat)
  w_df <- w_df[,c(ncol(w_df),1:(ncol(w_df)-1))]
  write.csv(w_df, sprintf('/path/peer/weights/%s.csv',i))
  
  # precision: Kx1 matrix
  p <- PEER_getAlpha(model)
  p_df <- data.frame(p)
  p_df$peer_factor <- c(colnames(cova_mat),pf.names)
  p_df <- p_df[,c(2,1)]
  write.csv(p_df, sprintf('/path/peer/precision/%s.csv',i))
  
  # residuals: NxG matrix
  r <- PEER_getResiduals(model)
  r_df <- data.frame(r)
  colnames(r_df) <- colnames(expr_mat)
  r_df <- cbind(sample = rownames(expr_mat), r_df)
  write.csv(r_df, sprintf('/path/peer/residuals/%s.csv',i))
  print(sprintf('%s DONEEEEEEEEEEEEEE', i))
}
