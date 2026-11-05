import os
import scanpy as sc
import pandas as pd
from collections import defaultdict
import anndata
import scanpy.external as sce
import numpy as np

# 读取临床信息
sample_info = pd.read_csv("sample.csv")
print(sample_info.head())
print(len(sample_info))

#筛选瓣膜样本
tissue_sample_name  = sample_info.loc[sample_info["type"] == 'valve', "name"]

# 自定义read10x
def self_read10x(data_dir):
    adata_raw = anndata.read_mtx(data_dir + '/matrix.mtx.gz').T
    features = pd.read_csv(data_dir + '/features.tsv.gz', header=None, sep='\t')
    adata_raw.var_names = anndata.utils.make_index_unique(pd.Index(features[0].values))
    adata_raw.obs_names = pd.read_csv(data_dir + '/barcodes.tsv.gz', header=None)[0].values
    return adata_raw

# 读取样本
adata_dict = defaultdict(sc.AnnData)
for sample in tissue_sample_name:
    data_dir = os.path.join("path", sample, "output", "filter_matrix")
    adata = self_read10x(data_dir)
    sc.pp.filter_cells(adata, min_genes=3)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs['sample_id'] = sample
    adata_dict[sample] = adata
print(len(adata_dict))

# 合并所有样本
adatas = list(adata_dict.values())
merged_adata = sc.concat(adatas,join = 'outer')
merged_adata.obs_names_make_unique()

merged_adata

# 质控
merged_adata.var["mt"] = merged_adata.var_names.str.startswith("MT-")
merged_adata.var["ribo"] = merged_adata.var_names.str.startswith(("RPS", "RPL"))
merged_adata.var["hb"] = merged_adata.var_names.str.contains("^HB[^(P)]")
sc.pp.calculate_qc_metrics(
    merged_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)
sc.pl.violin(
    merged_adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0,
    multi_panel=True,
)

sc.pp.filter_cells(merged_adata, min_genes=200)
sc.pp.filter_cells(merged_adata, max_genes=7000)
merged_adata = merged_adata[merged_adata.obs.pct_counts_mt < 15, :].copy()
sc.pp.scrublet(merged_adata,batch_key="sample_id")

merged_adata.layers["counts"] = merged_adata.X.copy()


disease_dict = sample_info.set_index("name")["disease"].to_dict()
merged_adata.obs["disease"] = merged_adata.obs["sample_id"].map(disease_dict)


adata_rm_db = merged_adata[merged_adata.obs.predicted_doublet == False,:].copy()
sc.pp.normalize_total(adata_rm_db, target_sum=1e4)
sc.pp.log1p(adata_rm_db)
sc.pp.highly_variable_genes(adata_rm_db, n_top_genes=2000)
adata_rm_db.raw = adata_rm_db.copy()
sc.pp.scale(adata_rm_db, max_value=10)
sc.tl.pca(adata_rm_db, svd_solver="arpack")
sce.pp.harmony_integrate(adata_rm_db,"sample_id")
sc.pp.neighbors(adata_rm_db,n_pcs=30, use_rep='X_pca_harmony')
sc.tl.umap(adata_rm_db)
sc.tl.leiden(adata_rm_db,flavor="igraph", n_iterations=2,resolution=1,key_added="res_1")
sc.pl.umap(adata_rm_db, color=["res_1"],legend_loc="on data")

##### VIC anno
VIC_cell = adata_rm_db[adata_rm_db.obs['res_1'].isin(["0","1","2","3","4","5","6","7","8","9","10"])].copy()
VIC_cell.X = VIC_cell.layers["counts"].copy()
sc.pp.normalize_total(VIC_cell, target_sum=1e4)
sc.pp.log1p(VIC_cell)
sc.pp.highly_variable_genes(VIC_cell, n_top_genes=2000)
VIC_cell.raw = VIC_cell.copy()
sc.pp.scale(VIC_cell, max_value=10)
sc.tl.pca(VIC_cell, svd_solver="arpack")
sce.pp.harmony_integrate(VIC_cell,"sample_id")
sc.pp.neighbors(VIC_cell,n_pcs=20,use_rep='X_pca_harmony')
sc.tl.umap(VIC_cell)
sc.tl.leiden(VIC_cell,flavor="igraph", n_iterations=2,resolution=0.7,key_added="vic_res_07")
sc.pl.umap(VIC_cell, color=["vic_res_07"],legend_loc="on data")
sc.tl.rank_genes_groups(VIC_cell, "vic_res_07", method="wilcoxon",pts=True)
VIC_cell.obs["vic_celltype"] = VIC_cell.obs["vic_res_07"].replace({"0":"VIC_C0_POSTN",
                                                                  "1":"VIC_C1_MYH11",
                                                                  "2":"VIC_C2_NEAT1",
                                                                  "3":"VIC_C3_RPS11",
                                                                  "4":"VIC_C4_C7",
                                                                  "5":"VIC_C5_ISG15",
                                                                  "6":"VIC_C6_PRELP",
                                                                  "7":"VIC_C7_LIF",
                                                                  "8":"VIC_C7_LIF",
                                                                  "9":"VIC_C8_PLCG2",
                                                                  "10":"VIC_C9_S100B",
                                                                  "11":"VIC_C2_NEAT1"})

##### Myeloid anno
mye = adata_rm_db[adata_rm_db.obs['res_1'].isin(["12","13"])].copy()
sc.pp.normalize_total(mye, target_sum=1e4)
sc.pp.log1p(mye)
sc.pp.highly_variable_genes(mye, n_top_genes=3000)
mye.raw = mye.copy()
sc.pp.scale(mye, max_value=10)
sc.tl.pca(mye, svd_solver="arpack")
sce.pp.harmony_integrate(mye,"sample_id")
sc.pp.neighbors(mye,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(mye)
sc.tl.leiden(mye,flavor="igraph", n_iterations=2,resolution=0.6,key_added="mye_res_06")
sc.pl.umap(mye, color=["mye_res_06"],legend_loc="on data")
sc.tl.rank_genes_groups(mye, "mye_res_06", method="wilcoxon",pts=True)

# doublet
mye_sub = mye[~mye.obs['mye_res_1'].isin(["10","16","17"])].copy()
mye_sub.X = mye_sub.layers["counts"].copy()
sc.pp.normalize_total(mye_sub, target_sum=1e4)
sc.pp.log1p(mye_sub)
sc.pp.highly_variable_genes(mye_sub, n_top_genes=3000)
mye_sub.raw = mye_sub.copy()
sc.pp.scale(mye_sub, max_value=10)
sc.tl.pca(mye_sub, svd_solver="arpack")
sce.pp.harmony_integrate(mye_sub,"sample_id")
sc.pp.neighbors(mye_sub,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(mye_sub)
sc.tl.leiden(mye_sub,flavor="igraph", n_iterations=2,resolution=0.6,key_added="mye_sub_res_06")
sc.pl.umap(mye_sub, color=["mye_sub_res_06"],legend_loc="on data")
sc.tl.rank_genes_groups(mye_sub, "mye_sub_res_06", method="wilcoxon",pts=True)

mye_sub_sub = mye_sub[~mye_sub.obs['mye_sub_res_06'].isin(["11"])].copy()
mye_sub_sub.X = mye_sub_sub.layers["counts"].copy()
sc.pp.normalize_total(mye_sub_sub, target_sum=1e4)
sc.pp.log1p(mye_sub_sub)
sc.pp.highly_variable_genes(mye_sub_sub, n_top_genes=3000)
mye_sub_sub.raw = mye_sub_sub.copy()
sc.pp.scale(mye_sub_sub, max_value=10)
sc.tl.pca(mye_sub_sub, svd_solver="arpack")
sce.pp.harmony_integrate(mye_sub_sub,"sample_id")
sc.pp.neighbors(mye_sub_sub,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(mye_sub_sub)
sc.tl.leiden(mye_sub_sub,flavor="igraph", n_iterations=2,resolution=0.6,key_added="mye_sub_sub_res_06")
sc.pl.umap(mye_sub_sub, color=["mye_sub_sub_res_06"],legend_loc="on data")
sc.tl.rank_genes_groups(mye_sub_sub, "mye_sub_sub_res_06", method="wilcoxon",pts=True)
mye_sub_sub.obs["mye_celltype"] = mye_sub_sub.obs["mye_sub_sub_res_06"].replace({"0":"M_Macro_TXNIP",
                                                             "1":"M_Macro_C3",
                                                             "2":"M_Macro_IL1B",
                                                             "3":"M_Macro_MRC1",
                                                             "4":"M_Macro_NEAT1",
                                                             "5":"M_Macro_HSPA1A",
                                                             "6":"M_Cycling_MKI67",
                                                             "7":"M_DC_CD1C",
                                                             "8":"M_Macro_SPP1",
                                                             "9":"M_Mast_TPSB2",
                                                             "10":"M_Macro_S100A8"})

##### T cell anno
tcell = adata_rm_db[adata_rm_db.obs['res_1'].isin(["14"])].copy()
tcell.X = tcell.layers["counts"].copy()
sc.pp.normalize_total(tcell, target_sum=1e4)
sc.pp.log1p(tcell)
sc.pp.highly_variable_genes(tcell, n_top_genes=3000)
tcell.raw = tcell.copy()
sc.pp.scale(tcell, max_value=10)
sc.tl.pca(tcell, svd_solver="arpack")
sce.pp.harmony_integrate(tcell,"sample_id",max_iter_harmony=50)
sc.pp.neighbors(tcell,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(tcell)
sc.tl.leiden(tcell,flavor="igraph", n_iterations=2,resolution=1,key_added="tcell_res_1")
sc.pl.umap(tcell, color=["tcell_res_1"],legend_loc="on data")
sc.tl.rank_genes_groups(tcell, "tcell_res_1", method="wilcoxon",pts=True)

#doublet remove
tcell_sub = tcell[~tcell.obs['tcell_res_1'].isin(["3","4","13","16"])].copy()
tcell_sub.X = tcell_sub.layers["counts"].copy()
sc.pp.normalize_total(tcell_sub, target_sum=1e4)
sc.pp.log1p(tcell_sub)
sc.pp.highly_variable_genes(tcell_sub, n_top_genes=3000)
tcell_sub.raw = tcell_sub.copy()
sc.pp.scale(tcell_sub, max_value=10)
sc.tl.pca(tcell_sub, svd_solver="arpack")
sce.pp.harmony_integrate(tcell_sub,"sample_id",max_iter_harmony=50)
sc.pp.neighbors(tcell_sub,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(tcell_sub)
sc.tl.leiden(tcell_sub,flavor="igraph", n_iterations=2,resolution=1,key_added="tcell_sub_res_1")
sc.pl.umap(tcell_sub, color=["tcell_sub_res_1"],legend_loc="on data")
sc.tl.rank_genes_groups(tcell_sub, "tcell_sub_res_1", method="wilcoxon",pts=True)

tcell_sub_sub = tcell_sub[~tcell_sub.obs['tcell_sub_res_1'].isin(["13"])].copy()
tcell_sub_sub.X = tcell_sub_sub.layers["counts"].copy()
tcell_sub_sub.X = tcell_sub_sub.layers["counts"].copy()
sc.pp.normalize_total(tcell_sub_sub, target_sum=1e4)
sc.pp.log1p(tcell_sub_sub)
sc.pp.highly_variable_genes(tcell_sub_sub, n_top_genes=2500)
tcell_sub_sub.raw = tcell_sub_sub.copy()
sc.pp.scale(tcell_sub_sub, max_value=10)
sc.tl.pca(tcell_sub_sub, svd_solver="arpack")
sce.pp.harmony_integrate(tcell_sub_sub,"sample_id",max_iter_harmony=50)
sc.pp.neighbors(tcell_sub_sub,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(tcell_sub_sub)
sc.tl.leiden(tcell_sub_sub,flavor="igraph", n_iterations=2,resolution=0.7,key_added="tcell_sub_sub_res_07")
sc.pl.umap(tcell_sub_sub, color=["tcell_sub_sub_res_07"],legend_loc="on data")
sc.tl.rank_genes_groups(tcell_sub_sub, "tcell_sub_sub_res_07", method="wilcoxon",pts=True)
tcell_sub_sub.obs["t_celltype"] = tcell_sub_sub.obs["tcell_sub_sub_res_07"].replace({"0":"T_CD4_IFR_ISG15",
                                                                                     "1":"T_CD8_EM_DUSP4",
                                                                                     "2":"T_MT2A",
                                                                                     "3":"T_STR_HSPA1A",
                                                                                     "4":"T_CD4_Memory_GPR183",
                                                                                     "5":"T_CD8_EM_CD69",
                                                                                     "6":"T_Trans_TXNIP",
                                                                                     "7":"T_CD4_Memory_GPR183",
                                                                                     "8":"T_CD4_TH17_CCL20",
                                                                                     "9":"NK_NCAM1",
                                                                                     "10":"ILC_KIT",
                                                                                     "11":"NK_FCGR3A",
                                                                                    "12":"T_Cycling_MKI67"})



##### endo cell anno
endo = adata_rm_db[adata_rm_db.obs['res_1'].isin(["11"])].copy()
endo.X = endo.layers["counts"].copy()
sc.pp.normalize_total(endo, target_sum=1e4)
sc.pp.log1p(endo)
sc.pp.highly_variable_genes(endo, n_top_genes=3000)
endo.raw = endo.copy()
sc.pp.scale(endo, max_value=10)
sc.tl.pca(endo, svd_solver="arpack")
sce.pp.harmony_integrate(endo,"sample_id")
sc.pp.neighbors(endo,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(endo)
sc.tl.leiden(endo,flavor="igraph", n_iterations=2,resolution=0.6,key_added="endo_res_06")
sc.pl.umap(endo, color=["endo_res_06"],legend_loc="on data")
sc.tl.rank_genes_groups(endo, "endo_res_06", method="wilcoxon",pts=True)

#doublet remove
endo_sub = endo[~endo.obs['endo_res_06'].isin(["9","10","13"])].copy()
endo_sub.X = endo_sub.layers["counts"].copy()
sc.pp.normalize_total(endo_sub, target_sum=1e4)
sc.pp.log1p(endo_sub)
sc.pp.highly_variable_genes(endo_sub, n_top_genes=3000)
endo_sub.raw = endo_sub.copy()
sc.pp.scale(endo_sub, max_value=10)
sc.tl.pca(endo_sub, svd_solver="arpack")
sce.pp.harmony_integrate(endo_sub,"sample_id")
sc.pp.neighbors(endo_sub,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(endo_sub)
sc.tl.leiden(endo_sub,flavor="igraph", n_iterations=2,resolution=0.4,key_added="endo_sub_res_04")
sc.pl.umap(endo_sub, color=["endo_sub_res_04"],legend_loc="on data")
sc.tl.rank_genes_groups(endo_sub, "endo_sub_res_04", method="wilcoxon",pts=True)
endo_sub.obs["endo_celltype"] = endo_sub.obs["endo_sub_res_04"].replace({"0":"Endo_C0_SELP",
                                                                        "1":"Endo_C1_NETA1",
                                                                        "2":"Endo_C2_HAPLN1",
                                                                        "3":"Endo_C3_HSPA1A",
                                                                        "4":"Endo_C4_PI16",
                                                                        "5":"Endo_C5_CLDN5",
                                                                        "6":"Endo_C1_NETA1",
                                                                        "7":"Endo_C6_CD74",
                                                                        "8":"Endo_C7_MKI67",
                                                                        "9":"Endo_C8_COL4A1"})
