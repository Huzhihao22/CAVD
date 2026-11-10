###### run in R #####
###### RunAzimuth #####
seuobj_list <- list()
options(future.globals.maxSize = 10*1024*1024*1024)
j=1
for(i in blood_sample){
  data <- Read10X(data.dir = "path_of_i",gene.column = 1)
  seuobj <- CreateSeuratObject(counts = data,min.cells = 3,min.features = 3)
  seuobj[["percent.mt"]] <- PercentageFeatureSet(seuobj, pattern = "^MT-")
  seuobj$sample_id = i
  seuobj = subset(seuobj,percent.mt <= 15 & nFeature_RNA >= 200 & nFeature_RNA <= 6000)
  seuobj <- RunAzimuth(seuobj, reference = azimuth_pbmc,verbose = F)
  seuobj_list[[i]] <- seuobj
  print(paste0("第",j,"个完成"))
  j = j+1
}
blood_metadata <- seuobj_list[[1]]@meta.data
blood_metadata$cellid_raw <- rownames(blood_metadata)
for(i in 2:length(seuobj_list)){
  meta_i = seuobj_list[[i]]@meta.data
  meta_i$cellid_raw = rownames(meta_i)
  blood_metadata = rbind(blood_metadata,meta_i)
}  

###### run in PYTHON #####
import os
import scanpy as sc
import pandas as pd
from collections import defaultdict
import anndata
import scanpy.external as sce
import numpy as np

# 自定义read10x
def self_read10x(data_dir):
    adata_raw = anndata.read_mtx(data_dir + '/matrix.mtx.gz').T
    features = pd.read_csv(data_dir + '/features.tsv.gz', header=None, sep='\t')
    adata_raw.var_names = anndata.utils.make_index_unique(pd.Index(features[0].values))
    adata_raw.obs_names = pd.read_csv(data_dir + '/barcodes.tsv.gz', header=None)[0].values
    return adata_raw

# 读取临床信息
sample_info = pd.read_csv("sample.csv")
print(sample_info.head())
print(len(sample_info))

#筛选血液样本
blood_sample_name  = sample_info.loc[sample_info["type"] == 'blood', "name"]
print(blood_sample_name.head())
print(len(blood_sample_name))

# azimuth注释信息
azimuth_anno = pd.read_csv("./azimuth_anno_metadata.csv")

# 读取样本
adata_dict = defaultdict(sc.AnnData)
for sample in blood_sample_name:
  data_dir = os.path.join("path", sample, "output","attachment", "splice_matrix")
  adata = self_read10x(data_dir)
  sc.pp.filter_cells(adata, min_genes=3)
  sc.pp.filter_genes(adata, min_cells=3)
  adata.obs['sample_id'] = sample
  adata.var["mt"] = adata.var_names.str.startswith("MT-")
  adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
  adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
  sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
  sc.pp.filter_cells(adata, min_genes=200)
  sc.pp.filter_cells(adata, max_genes=6000)
  adata = adata[adata.obs.pct_counts_mt < 15, :].copy()
  azimuth_anno_use = azimuth_anno[azimuth_anno['sample_id'] == sample].copy()
  azimuth_anno_use = azimuth_anno_use.set_index('cellid_raw')
  adata.obs = adata.obs.join(azimuth_anno_use[['predicted.celltype.l1.score', 'predicted.celltype.l1','predicted.celltype.l2.score','predicted.celltype.l2','predicted.celltype.l3.score','predicted.celltype.l3',"group"]], how='left')
  adata_dict[sample] = adata
    
print(len(adata_dict))

# 合并所有样本
adatas = list(adata_dict.values())
merged_adata = sc.concat(adatas,join = 'outer')
merged_adata.obs_names_make_unique()


sc.pp.scrublet(merged_adata,batch_key="sample_id")
adata_rm_db = merged_adata[merged_adata.obs.predicted_doublet == False,:].copy()
adata_rm_db.layers["counts"] = adata_rm_db.X.copy()

adata_rm_db.X = adata_rm_db.layers["counts"].copy()
sc.pp.normalize_total(adata_rm_db, target_sum=1e4)
sc.pp.log1p(adata_rm_db)
sc.pp.highly_variable_genes(adata_rm_db, n_top_genes=2500,batch_key="sample_id")
adata_rm_db.raw = adata_rm_db.copy()
sc.pp.scale(adata_rm_db, max_value=10)
sc.tl.pca(adata_rm_db, svd_solver="arpack")
sce.pp.harmony_integrate(adata_rm_db,"sample_id")
sc.pp.neighbors(adata_rm_db,n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_rm_db)
sc.tl.leiden(adata_rm_db,flavor="igraph", n_iterations=2,resolution=1,key_added="res_1")
sc.pl.umap(adata_rm_db, color=["res_1"],legend_loc="on data")

#Eryth,Platelet and doublet remove
adata_rm_db_sub = adata_rm_db[
    (~adata_rm_db.obs['res_1'].isin(["11", "16","20", "22"])) & 
    (~adata_rm_db.obs['predicted.celltype.l2'].isin(["Eryth", "Platelet"]))
].copy()
adata_rm_db_sub.X = adata_rm_db_sub.layers["counts"].copy()
sc.pp.normalize_total(adata_rm_db_sub, target_sum=1e4)
sc.pp.log1p(adata_rm_db_sub)
sc.pp.highly_variable_genes(adata_rm_db_sub, n_top_genes=3000)
adata_rm_db_sub.raw = adata_rm_db_sub.copy()
sc.pp.scale(adata_rm_db_sub, max_value=10)
sc.tl.pca(adata_rm_db_sub, svd_solver="arpack")
sce.pp.harmony_integrate(adata_rm_db_sub,"sample_id")
sc.pp.neighbors(adata_rm_db_sub,n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(adata_rm_db_sub)
sc.pl.umap(adata_rm_db_sub, color=["predicted.celltype.l2"],legend_loc="on data")
