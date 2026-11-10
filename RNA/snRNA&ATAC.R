# snRNA -------------------------------------------------------------------
library(Seurat)
library(harmony)
atac_rna_files <- list.files("/path_to_snRNA")
atac_rna_obj <- list()
for(i in atac_rna_files){
  data_temp <- Read10X(paste0("/path_to_snRNA/",i,"/02.count/filter_matrix/"),gene.column = 1)
  data_temp_seurat <- CreateSeuratObject(data_temp,min.cells = 3,min.features = 3)
  data_temp_seurat$sample_id = i
  atac_rna_obj[[i]] = data_temp_seurat
}
atac_rna_all <- merge(atac_rna_obj[[1]],y = atac_rna_obj[-1])
atac_rna_all[["percent.mt"]] <- PercentageFeatureSet(atac_rna_all, pattern = "^MT-")
atac_rna_all_filter <- subset(atac_rna_all,nFeature_RNA >200 & nFeature_RNA < 7000 & percent.mt < 15)
atac_rna_all_filter <- NormalizeData(atac_rna_all_filter,normalization.method = "LogNormalize", scale.factor = 1e4)
atac_rna_all_filter <- FindVariableFeatures(atac_rna_all_filter, selection.method = 'vst', nfeatures = 2000)
atac_rna_all_filter <- ScaleData(atac_rna_all_filter)
atac_rna_all_filter <- RunPCA(atac_rna_all_filter, features = VariableFeatures(object = atac_rna_all_filter))
atac_rna_all_filter <- RunHarmony(atac_rna_all_filter,group.by.vars = "sample_id")
atac_rna_all_filter <- FindNeighbors(atac_rna_all_filter, dims = 1:30,reduction = "harmony")
atac_rna_all_filter <- FindClusters(atac_rna_all_filter,resolution = 0.4)
atac_rna_all_filter <- RunUMAP(atac_rna_all_filter, dims = 1:30,reduction = "harmony")
mk_round1 <- FindAllMarkers(atac_rna_all_filter)


# doublet remove ----------------------------------------------------------
atac_rna_all_filter_sub <- subset(atac_rna_all_filter,seurat_clusters != 6)
atac_rna_all_filter_sub <- DietSeurat(atac_rna_all_filter_sub)
atac_rna_all_filter_sub <- NormalizeData(atac_rna_all_filter_sub,normalization.method = "LogNormalize", scale.factor = 1e4)
atac_rna_all_filter_sub <- FindVariableFeatures(atac_rna_all_filter_sub, selection.method = 'vst', nfeatures = 2000)
atac_rna_all_filter_sub <- ScaleData(atac_rna_all_filter_sub)
atac_rna_all_filter_sub <- RunPCA(atac_rna_all_filter_sub, features = VariableFeatures(object = atac_rna_all_filter_rm_eys))
atac_rna_all_filter_sub <- RunHarmony(atac_rna_all_filter_sub,group.by.vars = "sample_id")
atac_rna_all_filter_sub <- FindNeighbors(atac_rna_all_filter_sub, dims = 1:30,reduction = "harmony")
atac_rna_all_filter_sub <- FindClusters(atac_rna_all_filter_sub,resolution = 0.4)
atac_rna_all_filter_sub <- RunUMAP(atac_rna_all_filter_sub, dims = 1:30,reduction = "harmony")
DimPlot(atac_rna_all_filter_sub,label = T)
mk_round2 <- FindAllMarkers(atac_rna_all_filter_sub)

atac_rna_all_filter_sub$celltype <- atac_rna_all_filter_sub$seurat_clusters
atac_rna_all_filter_sub$celltype <- atac_rna_all_filter_sub(snrna_all$celltype,
                                from = seq(0,13,1),
                                to = c("VIC","VEC","VEC","VIC","Myeloid","VIC",
                                       "VEC","TNK","Myeloid","VEC","VEC","VEC","VIC","Myeloid"))



# snATAC -------------------------------------------------------------------
library(ArchR)
library(RColorBrewer)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(threads = 8) 
pathToMacs2 <- "/PATH/macs2"

ArrowFiles <- createArrowFiles(
  inputFiles = inputfile_list,
  sampleNames = atac_files_all,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)
proj_all <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archRProject_dir",
  copyArrows = TRUE
)
proj_all_rm_db <- filterDoublets(proj_all)
proj_all_rm_db <- addIterativeLSI(
  ArchRProj = proj_all_rm_db,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = T
)
proj_all_rm_db <- addHarmony(
  ArchRProj = proj_all_rm_db,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = T
)
proj_all_rm_db <- addClusters(
  input = proj_all_rm_db,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = T
)
proj_all_rm_db <- addUMAP(
  ArchRProj = proj_all_rm_db, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 40, 
  minDist = 0.4, 
  metric = "cosine",
  force = T
)
markersGS <- getMarkerFeatures(
  ArchRProj = proj_all_rm_db, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
proj_all_rm_db$celtype = mapvalues(proj_all_rm_db$celtype,
                                   from=c("C1","C2","C3","C4","C5","C6","C7","C8","C9"),
                                   to=c("VIC","VIC","VIC","VIC","VIC","VEC","VEC","Myeloid","TNK"))
proj_all_rm_db <- addReproduciblePeakSet(
  ArchRProj = proj_all_rm_db, 
  groupBy = "celtype", 
  pathToMacs2 = pathToMacs2)
PeakSet <- getPeakSet(proj_all_rm_db)
proj_all_rm_db <- addPeakMatrix(proj_all_rm_db)
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj_all_rm_db, 
  useMatrix = "PeakMatrix", 
  groupBy = "celtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

proj_all_rm_db <- addArchRAnnotations(ArchRProj = proj_all_rm_db, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(
  seMarker = markerPeaks,
  ArchRProj = proj_all_rm_db,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 10, transpose = TRUE)
peak_heatmap_data = heatmapEncode@matrix
peak_heatmap_data = peak_heatmap_data[c(3,1,2,4),c(16:25,1:15,26:35)]
pheatmap::pheatmap(peak_heatmap_data,color = colorRampPalette(rev(brewer.pal(name = "RdBu",n = 11)[2:8]))(100),
                   border_color = "lightgrey",treeheight_row = 0,treeheight_col = 0,
                   cluster_rows = F,cluster_cols = F)
p_track <- plotBrowserTrack(
  ArchRProj = proj_all_rm_db, 
  groupBy = "celtype", 
  geneSymbol = c("PECAM1","MRC1","IL7R","COL1A1"),
  features =  getPeakSet(proj_all_rm_db),
  upstream = 50000,
  downstream = 50000,
  plotSummary = c("bulkTrack","geneTrack")
)
grid::grid.draw(p_track$PECAM1)
grid::grid.draw(p_track$MRC1)
grid::grid.draw(p_track$IL7R)
grid::grid.draw(p_track$COL1A1)

