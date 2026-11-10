# VIC GOenrich ------------------------------------------------------------
mk_vic <- read.csv("./mk_vic.csv")
mk_vic_sig <- mk_vic[mk_vic$pvals_adj < 0.05 & mk_vic$pct_nz_group > 0.25,]
mk_vic_sig <- mk_vic_sig[order(mk_vic_sig$group),]
mk_vic_sig$group <- factor(mk_vic_sig$group,levels = unique(as.character(mk_vic_sig$group)))
run_go_enrichment <- function(cluster_id) {
  mk_temp <- mk_vic_sig[mk_vic_sig$group == cluster_id,]
  mk_temp <- mk_temp[order(mk_temp$logfoldchanges,decreasing = T),]
  cluster_genes <- mk_temp[1:50,"names"]
  entrez_ids <- bitr(cluster_genes, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)$ENTREZID
  go_res <- enrichGO(gene = entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP",  # 可选 "CC", "MF"
                     pAdjustMethod = "BH",
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
  return(go_res)
}
all_clusters <- as.character(unique(mk_vic_sig$group))
go_results <- lapply(all_clusters, run_go_enrichment)
names(go_results) <- all_clusters
all_go_results <- bind_rows(
  lapply(names(go_results), function(x) {
    go_results[[x]]@result %>% 
      mutate(source = x)
  })
)
top_go_terms <- lapply(names(go_results), function(cluster) {
  res <- go_results[[cluster]]@result %>%
    arrange(p.adjust) %>% 
    head(5) %>% 
    mutate(Cluster = cluster,
           Term = Description,
           `-log10(p.adjust)` = -log10(p.adjust)) %>% 
    select(Cluster, Term, `-log10(p.adjust)`)
}) %>% bind_rows()
use_term <- unique(top_go_terms$Term)
all_go_results_select <- all_go_results[all_go_results$Description %in% use_term,c("Description","pvalue","source")]
all_go_results_select$pvalue <- (-log10(all_go_results_select$pvalue))
all_go_results_select_wide <- all_go_results_select %>%
  pivot_wider(
    names_from = source,
    values_from = pvalue
  )
all_go_results_select_wide <- as.data.frame(all_go_results_select_wide)
all_go_results_select_wide[is.na(all_go_results_select_wide)] = 0
rownames(all_go_results_select_wide) <- all_go_results_select_wide$Description
all_go_results_select_wide <- all_go_results_select_wide[,-1]
all_go_results_select_wide_adj <- all_go_results_select_wide
all_go_results_select_wide_adj[all_go_results_select_wide_adj >6] = 6
pheatmap::pheatmap(all_go_results_select_wide_adj[use_term,],
                   color = colorRampPalette(rev(brewer.pal(name = "RdBu",n = 11)[2:8]))(100),
                   border_color = NA,cluster_rows = F,cluster_cols = F)

# monocle3 ----------------------------------------------------------------
library(Seurat)
library(monocle3)
new_cell_data_set <- function (expression_data, cell_metadata = NULL, gene_metadata = NULL) 
{
  assertthat::assert_that(methods::is(expression_data, "matrix") || 
                            is_sparse_matrix(expression_data), msg = paste("Argument expression_data must be a", 
                                                                           "matrix - either sparse from the", "Matrix package or dense"))
  if (!is.null(cell_metadata)) {
    assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data), 
                            msg = paste("cell_metadata must be NULL or have", 
                                        "the same number of rows as columns", "in expression_data"))
    assertthat::assert_that(!is.null(row.names(cell_metadata)) & 
                              all(row.names(cell_metadata) == colnames(expression_data)), 
                            msg = paste("row.names of cell_metadata must be equal to colnames of", 
                                        "expression_data"))
  }
  if (!is.null(gene_metadata)) {
    assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data), 
                            msg = paste("gene_metadata must be NULL or have", 
                                        "the same number of rows as rows", "in expression_data"))
    assertthat::assert_that(!is.null(row.names(gene_metadata)) & 
                              all(row.names(gene_metadata) == row.names(expression_data)), 
                            msg = paste("row.names of gene_metadata must be equal to row.names of", 
                                        "expression_data"))
  }
  if (is.null(cell_metadata)) {
    cell_metadata <- data.frame(cell = colnames(expression_data), 
                                row.names = colnames(expression_data))
  }
  if (!("gene_short_name" %in% colnames(gene_metadata))) {
    warning("gene_metadata must contain a column verbatim ", 
            "named 'gene_short_name' for certain functions.")
  }
  sce <- SingleCellExperiment(list(counts = methods::as(expression_data, 
                                                        "dgCMatrix")), rowData = gene_metadata, colData = cell_metadata)
  cds <- methods::new("cell_data_set", assays = SummarizedExperiment::Assays(list(counts = methods::as(expression_data, 
                                                                                                       "dgCMatrix"))), colData = colData(sce), int_elementMetadata = SingleCellExperiment::int_elementMetadata(sce), 
                      int_colData = SingleCellExperiment::int_colData(sce), 
                      int_metadata = SingleCellExperiment::int_metadata(sce), 
                      metadata = S4Vectors::metadata(sce), NAMES = NULL, elementMetadata = elementMetadata(sce)[, 
                                                                                                                0], rowRanges = rowRanges(sce))
  S4Vectors::metadata(cds)$cds_version <- Biobase::package.version("monocle3")
  clusters <- stats::setNames(S4Vectors::SimpleList(), character(0))
  cds <- estimate_size_factors(cds)
  cds
}
seurat_metadata <- vic_seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(seurat_data))
rownames(gene_annotation) <- rownames(seurat_data)
gene_annotation <- gene_annotation[1:100,,drop=F]

cds <- new_cell_data_set(expression_data = seurat_data,
                         cell_metadata = seurat_metadata,
                         gene_metadata = gene_annotation)
seurat_embed <- Embeddings(vic_seurat, reduction = "Xumap_")
table(colnames(cds) == rownames(seurat_embed))
cds@int_colData$reducedDims$UMAP <- seurat_embed
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
vic_monocle3 <- order_cells(cds)
vic_ps_time <- as.data.frame(pseudotime(vic_monocle3))
vic_seurat$monocle_time = vic_ps_time$`pseudotime(vic_monocle3)`
FeaturePlot(vic,features = "monocle_time")
vic_gene_time = FetchData(vic,vars = c("monocle_time","IL6","LIF"))
vic_gene_time_long <- vic_gene_time %>%
  pivot_longer(
    cols = -monocle_time, 
    names_to = "gene",
    values_to = "score"
  )
ggplot(vic_gene_time_long,aes(x = monocle_time,y = score,color=gene))+
  geom_smooth()+
  theme_classic2()+
  ylim(0, NA)

vic_hmk <- read.csv("./vic_hmk_genescore.csv",row.names = 1)
table(rownames(vic_hmk) == colnames(vic_seurat))
vic_seurat@meta.data <- cbind(vic_seurat@meta.data,vic_hmk[,34:83])
vic_path_time = FetchData(vic,vars = c("monocle_time","HALLMARK_IL6_JAK_STAT3_SIGNALING_score","HALLMARK_INFLAMMATORY_RESPONSE_score"))
vic_path_time_long <- vic_path_time %>%
  pivot_longer(
    cols = -monocle_time, 
    names_to = "path",
    values_to = "score"
  )
ggplot(vic_path_time_long,aes(x = monocle_time,y = score,color=path))+
  geom_smooth()+
  theme_classic2()

# miloR -------------------------------------------------------------------
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scuttle)
library(plyr)
library(ggbeeswarm)
library(ggpubr)
library(readr)
plotDAbeeswarm_self = function (da.res, group.by = NULL, alpha = 0.1, subset.nhoods = NULL) 
{
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", 
           group.by, ")?")
    }
    if (is.numeric(da.res[, group.by])) {
    }
    da.res <- mutate(da.res, group_by = da.res[, group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  
  if (!is.factor(da.res[, "group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, 
                                               levels = unique(group_by)))
  }
  
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods, ]
  }
  
  beeswarm_pos <- ggplot_build(
    da.res %>% 
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>% 
      arrange(group_by) %>% 
      ggplot(aes(group_by, logFC)) + 
      geom_quasirandom()
  )
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  n_groups <- unique(da.res$group_by) %>% length()
  
  da.res <- da.res %>% 
    mutate(signif_type = case_when(
      SpatialFDR < alpha & logFC > 0 ~ 1,   
      SpatialFDR < alpha & logFC < 0 ~ 2,   
      TRUE ~ 3                             
    )) 
  
  da.res %>% 
    arrange(group_by) %>% 
    mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>% 
    mutate(pos_x = pos_x, pos_y = pos_y) %>% 
    ggplot(aes(pos_x, pos_y)) + 
    geom_point(aes(color = factor(signif_type)), size = 1) + 
    scale_color_manual(values = c("2" = "#440154FF", "1" =  "#FDE725FF", "3" = "lightgrey")) +
    guides(color = "none", shape = "none") +  
    xlab("") + 
    ylab("Log Fold Change") + 
    scale_x_continuous(
      breaks = seq(1, n_groups), 
      labels = setNames(levels(da.res$group_by), seq(1, n_groups))
    ) + 
    coord_flip() + 
    theme_classic2(base_size = 22) + 
    theme(strip.text.y = element_text(angle = 0))+
    theme(axis.text.y = element_blank())
}
#VIC
vic_sce <- as.SingleCellExperiment(vic_seurat)
vic.milo <- Milo(vic_sce)
reducedDim(vic.milo, "PCA") <- reducedDim(vic_sce, "XPCA_")
reducedDim(vic.milo, "UMAP") <- reducedDim(vic_sce, "XUMAP_")
vic.milo <- buildGraph(vic.milo, k = 10, d = 30)
vic.milo <- makeNhoods(vic.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
vic.milo <- countCells(vic.milo, meta.data = data.frame(colData(vic.milo)), samples="sample_id") # make the nhood X sample counts matrix
plotNhoodSizeHist(vic.milo)
thy.design <- data.frame(colData(vic.milo))[,c("sample_id","disease")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$sample_id
thy.design <- thy.design[colnames(nhoodCounts(vic.milo)), , drop=FALSE]
vic.milo <- calcNhoodDistance(vic.milo, d=30)
contrast.1 <- c("diseasedisease - diseasenormal")
da_results1 <- testNhoods(vic.milo, design = ~ 0 + disease, design.df = thy.design, 
                          model.contrasts = contrast.1,
                          fdr.weighting="graph-overlap", norm.method="TMM")
vic.milo_nood_anno <- annotateNhoods(vic.milo,da.res = da_results1,coldata_col = "vic_celltype")
vic.milo_nood_anno2 <- vic.milo_nood_anno
vic.milo_nood_anno2[vic.milo_nood_anno2$subtype_fraction < 0.8,"vic_celltype"] <- "mix"
plotDAbeeswarm_self(vic.milo_nood_anno2,group.by = "vic_celltype")

#Mye
mye_sce <- as.SingleCellExperiment(mye_seurat)
mye.milo <- Milo(mye_sce)
reducedDim(mye.milo, "PCA") <- reducedDim(mye_sce, "XPCA_")
reducedDim(mye.milo, "UMAP") <- reducedDim(mye_sce, "XUMAP_")
mye.milo <- buildGraph(mye.milo, k = 10, d = 30)
mye.milo <- makeNhoods(mye.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
mye.milo <- countCells(mye.milo, meta.data = data.frame(colData(mye.milo)), samples="sample_id") # make the nhood X sample counts matrix
plotNhoodSizeHist(mye.milo)
thy.design <- data.frame(colData(mye.milo))[,c("sample_id","disease")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$sample_id
thy.design <- thy.design[colnames(nhoodCounts(mye.milo)), , drop=FALSE]
mye.milo <- calcNhoodDistance(mye.milo, d=30)
contrast.1 <- c("diseasedisease - diseasenormal")
da_results1 <- testNhoods(mye.milo, design = ~ 0 + disease, design.df = thy.design, 
                          model.contrasts = contrast.1,
                          fdr.weighting="graph-overlap", norm.method="TMM")
mye.milo <- buildNhoodGraph(mye.milo)
mye.milo_nood_anno <- annotateNhoods(mye.milo,da.res = da_results1,coldata_col = "mye_celltype")
mye.milo_nood_anno2 <- mye.milo_nood_anno
mye.milo_nood_anno2[mye.milo_nood_anno2$subtype_fraction < 0.8,"mye_celltype"] <- "mix"
plotDAbeeswarm_self(mye.milo_nood_anno2,group.by = "mye_celltype")

#tnk
tnk_sce <- as.SingleCellExperiment(tnk_seurat)
tnk.milo <- Milo(tnk_sce)
reducedDim(tnk.milo, "PCA") <- reducedDim(tnk_sce, "XPCA_")
reducedDim(tnk.milo, "UMAP") <- reducedDim(tnk_sce, "XUMAP_")
tnk.milo <- buildGraph(tnk.milo, k = 10, d = 30)
tnk.milo <- makeNhoods(tnk.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
tnk.milo <- countCells(tnk.milo, meta.data = data.frame(colData(tnk.milo)), samples="sample_id") # make the nhood X sample counts matrix
plotNhoodSizeHist(tnk.milo)
thy.design <- data.frame(colData(tnk.milo))[,c("sample_id","disease")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$sample_id
thy.design <- thy.design[colnames(nhoodCounts(tnk.milo)), , drop=FALSE]
tnk.milo <- calcNhoodDistance(tnk.milo, d=30)
contrast.1 <- c("diseasedisease - diseasenormal")
da_results1 <- testNhoods(tnk.milo, design = ~ 0 + disease, design.df = thy.design, 
                          model.contrasts = contrast.1,
                          fdr.weighting="graph-overlap", norm.method="TMM")
tnk.milo <- buildNhoodGraph(tnk.milo)
tnk.milo_nood_anno <- annotateNhoods(tnk.milo,da.res = da_results1,coldata_col = "t_celltype")
tnk.milo_nood_anno2 <- tnk.milo_nood_anno
tnk.milo_nood_anno2[tnk.milo_nood_anno2$subtype_fraction < 0.8,"t_celltype"] <- "mix"
plotDAbeeswarm_self(tnk.milo_nood_anno2,group.by = "t_celltype")

#endo
endo_sce <- as.SingleCellExperiment(endo_seurat)
endo.milo <- Milo(endo_sce)
reducedDim(endo.milo, "PCA") <- reducedDim(endo_sce, "XPCA_")
reducedDim(endo.milo, "UMAP") <- reducedDim(endo_sce, "XUMAP_")
endo.milo <- buildGraph(endo.milo, k = 10, d = 30)
endo.milo <- makeNhoods(endo.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
endo.milo <- countCells(endo.milo, meta.data = data.frame(colData(endo.milo)), samples="sample_id") # make the nhood X sample counts matrix
plotNhoodSizeHist(endo.milo)
thy.design <- data.frame(colData(endo.milo))[,c("sample_id","disease")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$sample_id
thy.design <- thy.design[colnames(nhoodCounts(endo.milo)), , drop=FALSE]
endo.milo <- calcNhoodDistance(endo.milo, d=30)
contrast.1 <- c("diseasedisease - diseasenormal")
da_results1 <- testNhoods(endo.milo, design = ~ 0 + disease, design.df = thy.design, 
                          model.contrasts = contrast.1,
                          fdr.weighting="graph-overlap", norm.method="TMM")
endo.milo <- buildNhoodGraph(endo.milo)
endo.milo_nood_anno <- annotateNhoods(endo.milo,da.res = da_results1,coldata_col = "endo_celltype")
endo.milo_nood_anno2 <- endo.milo_nood_anno
endo.milo_nood_anno2[endo.milo_nood_anno2$subtype_fraction < 0.8,"endo_celltype"] <- "mix"
plotDAbeeswarm_self(endo.milo_nood_anno2,group.by = "endo_celltype")

#pbmc
pbmc_all_sce <- as.SingleCellExperiment(pbmc_all_seurat)
pbmc_all.milo <- Milo(pbmc_all_sce)
reducedDim(pbmc_all.milo, "PCA") <- reducedDim(pbmc_all_sce, "XPCA_")
reducedDim(pbmc_all.milo, "UMAP") <- reducedDim(pbmc_all_sce, "XUMAP_")
plotUMAP(pbmc_all.milo,colour_by="predicted.celltype.l2") + plotUMAP(pbmc_all.milo, colour_by="disease")
pbmc_all.milo <- buildGraph(pbmc_all.milo, k = 10, d = 30)
pbmc_all.milo <- makeNhoods(pbmc_all.milo, prop = 0.1, k = 10, d=30, refined = TRUE, refinement_scheme="graph") # make nhoods using graph-only as this is faster
pbmc_all.milo <- countCells(pbmc_all.milo, meta.data = data.frame(colData(pbmc_all.milo)), samples="sample_id") # make the nhood X sample counts matrix
plotNhoodSizeHist(pbmc_all.milo)

thy.design <- data.frame(colData(pbmc_all.milo))[,c("sample_id","group")]
thy.design <- distinct(thy.design)
rownames(thy.design) <- thy.design$sample_id
thy.design <- thy.design[colnames(nhoodCounts(pbmc_all.milo)), , drop=FALSE]

pbmc_all.milo <- calcNhoodDistance(pbmc_all.milo, d=30)
contrast.1 <- c("groupdisease - groupnormal")
da_results1 <- testNhoods(pbmc_all.milo, design = ~ 0 + group, design.df = thy.design, 
                          model.contrasts = contrast.1,
                          fdr.weighting="graph-overlap", norm.method="TMM")
nood_anno <- annotateNhoods(pbmc_all.milo,da.res = da_results1,coldata_col = "predicted.celltype.l2")
nood_anno_2 <- nood_anno
nood_anno_2[nood_anno_2$predicted.celltype.l2_fraction < 0.8,"predicted.celltype.l2"] <- "mix"
nood_anno_2_sub = nood_anno_2[nood_anno_2$predicted.celltype.l2 != "mix",]
nood_anno_2_sub$predicted.celltype.l2 <- factor(nood_anno_2_sub$predicted.celltype.l2,
                                                levels = rev(c("CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL","CD4 Proliferating","Treg","CD8 Naive","CD8 TCM",
                                                               "CD8 TEM","CD8 Proliferating","NK","NK_CD56bright","NK Proliferating","ILC","MAIT","dnT",
                                                               "gdT","B naive","B memory","B intermediate","Plasmablast","CD14 Mono","CD16 Mono","cDC1",
                                                               "cDC2","pDC","ASDC","HSPC")))
plotDAbeeswarm_self(nood_anno_2_sub,group.by = "predicted.celltype.l2")


# cellchat ----------------------------------------------------------------
#稍微降采样至每个细胞类型最多100,000细胞，避免无法计算
Idents(sc_valve_disease) <- sc_valve_disease$subtype
sc_valve_disease <- subset(sc_valve_disease,downsample=100000)
cellChat_disease <- createCellChat(object = sc_valve_disease, group.by = "subtype", assay = "RNA")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)
cellChat_disease@DB <- CellChatDB.use
cellChat_disease <- subsetData(cellChat_disease) # This step is necessary even if using the whole database
cellChat_disease <- identifyOverExpressedGenes(cellChat_disease)
cellChat_disease <- identifyOverExpressedInteractions(cellChat_disease)
cellChat_disease <- computeCommunProb(cellChat_disease, type = "triMean")
cellChat_disease <- filterCommunication(cellChat_disease, min.cells = 5)
cellChat_disease <- computeCommunProbPathway(cellChat_disease)
cellChat_disease <- aggregateNet(cellChat_disease)
cellChat_disease <- netAnalysis_computeCentrality(cellChat_disease, slot.name = "netP")

cellChat_normal <- createCellChat(object = sc_valve_normal, group.by = "subtype", assay = "RNA")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB)
cellChat_normal@DB <- CellChatDB.use
cellChat_normal <- subsetData(cellChat_normal) 
cellChat_normal <- identifyOverExpressedGenes(cellChat_normal)
cellChat_normal <- identifyOverExpressedInteractions(cellChat_normal)
cellChat_normal <- computeCommunProb(cellChat_normal, type = "triMean")
cellChat_normal <- filterCommunication(cellChat_normal, min.cells = 5)
cellChat_normal <- computeCommunProbPathway(cellChat_normal)
cellChat_normal <- aggregateNet(cellChat_normal)
cellChat_normal <- netAnalysis_computeCentrality(cellChat_normal, slot.name = "netP")
object.list <- list(normal = cellChat_normal, disease = cellChat_disease)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#绘制Role
netAnalysis_signalingRole_scatter(cellChat_normal)+
  netAnalysis_signalingRole_scatter(cellChat_disease)

#绘制Interactions数
compareInteractions(cellchat, show.legend = F, group = c(1,2))+
  compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

#绘制top差异通路
p = rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, 
            stacked = T, do.stat = TRUE,color.use = c("#BAD9E9","#FCD7C2"))
p_data <- ggplot_build(p)$plot$data
p_data$name <- as.character(p_data$name)
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, 
        stacked = T, do.stat = TRUE,color.use = c("#BAD9E9","#FCD7C2"),signaling = p_data[c(1:30,197:216),1])

#IL6 heatmap
pathways.show <- c("IL6") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#IL6 dotplot
p = netVisual_bubble(cellchat,sources.use = "VIC_C7_LIF", signaling = "IL6",comparison = c(1, 2), angle.x = 45)
data_plt <- ggplot_build(p)$plot$data
ggplot(data_plt, aes(x=source.target, y=interaction_name_2, fill=prob, size=pval)) +
  geom_point(pch=21,stroke=0.1) +
  scale_radius(range = c(0, 20)) +
  scale_fill_gradientn(name='Commun.Prob', colors=rev(brewer.pal(11, 'RdBu'))) +
  theme_bw(base_size=14) +
  geom_vline(xintercept=c(2.5,4.5,6.5,7.5,9.5,11.5,13.5,15.5,17.5), linetype="dashed", color = "grey20",size=0.1)





# cellratio nmf -----------------------------------------------------------
library(NMF)
valve_meta_VIC <- valve_obs[valve_obs$main_celltype == "VIC",]
valve_meta_Myeloid <- valve_obs[valve_obs$main_celltype == "Myeloid",]
valve_meta_TNK <- valve_obs[valve_obs$main_celltype == "TNK",]
valve_meta_EC <- valve_obs[valve_obs$main_celltype == "EC",]
cellratio_VIC <- as.data.frame.matrix(prop.table(table(valve_meta_VIC$sample_id,valve_meta_VIC$subtype),margin = 1))
cellratio_Myeloid <- as.data.frame.matrix(prop.table(table(valve_meta_Myeloid$sample_id,valve_meta_Myeloid$subtype),margin = 1))
cellratio_TNK <- as.data.frame.matrix(prop.table(table(valve_meta_TNK$sample_id,valve_meta_TNK$subtype),margin = 1))
cellratio_EC <- as.data.frame.matrix(prop.table(table(valve_meta_EC$sample_id,valve_meta_EC$subtype),margin = 1))
cellratio_merge = cbind(cellratio_VIC,cellratio_Myeloid,cellratio_TNK,cellratio_EC)


nmf_input = cellratio_merge
nmf.res_raw = nmf(t(nmf_input), rank=3:10, nrun=20, seed=123456)
plot(nmf.res_raw)
nrank = 3
nmf.res = nmf(t(nmf_input), rank=nrank, nrun=20, seed=123456)
w = basis(nmf.res)
colnames(w) = paste0('NMF', seq(1,nrank))
w_ = w %>% t %>% scale
Heatmap(w_, width = 11, height = 4, column_km = nrank,
        cluster_rows = F, clustering_method_columns = 'single',
        col =c(colorRampPalette(c("#98A6C8","white"))(25),colorRampPalette(c("white","#940E26"))(25))) 
h = coef(nmf.res)
h_t <- t(h)
h_t <- as.data.frame(h_t)
colnames(h_t) <- c("NMF1","NMF2","NMF3")
h_t$group <- rownames(h_t)
h_t$group <- mapvalues(h_t$group,from = name_type$name,to = name_type$disease)
h_t_melt <- reshape2::melt(h_t)
ggplot(h_t_melt, aes(x = group, y = value, color = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 1, width = 0.2,size=3) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif" 
  ) +
  scale_color_manual(values = c("#98A6C8", "#940E26")) +
  theme(legend.position = "none") +
  theme_classic2(base_size = 20)
