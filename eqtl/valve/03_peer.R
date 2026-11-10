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
usecell <- c("Endo_C1_NETA1","Endo_C2_HAPLN1","Endo_C4_PI16","M_DC_CD1C","M_Macro_C3","M_Macro_HSPA1A",
             "M_Macro_IL1B","M_Macro_MRC1","M_Macro_NEAT1",
             "M_Macro_SPP1","M_Macro_TXNIP","M_Mast_TPSB2","NK_NCAM1",
             "T_CD4_Memory_GPR183","T_CD4_TH17_CCL20","T_CD8_EM_CD69","T_CD8_EM_DUSP4",
             "T_STR_HSPA1A","T_Trans_TXNIP","VIC_C0_POSTN","VIC_C1_MYH11",
             "VIC_C2_NEAT1","VIC_C3_RPS11","VIC_C4_C7","VIC_C5_ISG15",
             "VIC_C6_PRELP","VIC_C7_LIF","VIC_C9_S100B")
valve_snp_sample <- read.csv("/path/valve_with_snp.csv",row.names = 1)
valve_snp_sample$name <- as.character(valve_snp_sample$name)
valve_snp_sample$snp_file <- as.character(valve_snp_sample$snp_file)
snp_pc_mat <- as.data.frame(read.table("/path/snp.eigenvec"))
rownames(snp_pc_mat) <- snp_pc_mat$V2
snp_pc_mat <- snp_pc_mat[,-c(1,2)]
colnames(snp_pc_mat) <- paste0("PC_",seq(1,20,1))
snp_pc_mat_valve <- snp_pc_mat[rownames(snp_pc_mat) %in% valve_snp_sample$snp_file,]
for(i in usecell){
  print(sprintf('%s START', i))
  expr_mat <- as.data.frame(read.csv(paste0("/path/top2000_normal_dis/",i,".csv"),row.names = 1))
  cova_mat <- as.data.frame(read.csv("/path/valve_with_snp.csv",row.names = 1))
  rownames(cova_mat) <- cova_mat$name
  cosample <- intersect(rownames(expr_mat),rownames(cova_mat))
  expr_mat <- expr_mat[cosample,]
  cova_mat <- cova_mat[cosample,c("sex","age")]
  cova_mat$sex_adj <- 1
  cova_mat[cova_mat$sex == "female","sex_adj"] <- 0
  cova_mat = cova_mat[,c("sex_adj","age")]
  print(table(rownames(expr_mat) == rownames(cova_mat)))
  rownames(expr_mat) <- mapvaluesss(rownames(expr_mat),
                                    from = valve_snp_sample$name,
                                    to = valve_snp_sample$snp_file)
  rownames(cova_mat) <- mapvaluesss(rownames(cova_mat),
                                    from = valve_snp_sample$name,
                                    to = valve_snp_sample$snp_file)
  snp_pc_mat_use <- snp_pc_mat_valve[rownames(cova_mat),]
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
  print(sprintf('%s DONE', i))
}
