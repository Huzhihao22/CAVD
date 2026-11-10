import scanpy as sc
import pandas as pd
import decoupler as dc
adata = sc.read_h5ad("/path/valve.h5ad")
adata.X = adata.layers["counts"].copy()
valve_snp = pd.read_csv("/path/valve_with_snp.csv")
valve_with_snp = valve_snp["name"]
tss_bed = pd.read_csv('/path/tss.bed',sep = ' ')
adata = adata[:,adata.var_names.isin(tss_bed['gene_id'])]
adata = adata[adata.obs['sample_id'].isin(valve_with_snp),:]
def preprocess_RNA_data(adata_use):
    for celltype in adata_use.obs['subtype'].unique():
        print(f"正在处理 {celltype}", flush=True)
        #分小类
        adata_use_celltype = adata_use[adata_use.obs['subtype'] == celltype]
        #算伪Bulk
        pdata = dc.pp.pseudobulk(
            adata_use_celltype,
            sample_col='sample_id',
            groups_col=None,
            mode='mean')
        dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000)
        
        #生成伪bulk矩阵
        pseudo_matrix = pd.DataFrame(pdata.X)
        pseudo_matrix.columns = pdata.var_names
        pseudo_matrix.index = pdata.obs.index

        #只选取在90%的样本中都表达的特征
        non_zero_ratio = (pseudo_matrix != 0).mean()
        columns_to_keep = non_zero_ratio[non_zero_ratio >= 0.9].index
        pseudo_matrix = pseudo_matrix[columns_to_keep]

        # 计算每列的均值和标准差
        means = pseudo_matrix.mean()
        stds = pseudo_matrix.std()
        # 计算变异系数
        cv = (stds / means).abs() * 100
        # 按变异系数从大到小排序
        sorted_columns = cv.sort_values(ascending=False).index
        # 选取变异系数最高的前2000列
        top_2000_columns = sorted_columns[:min(2000,pseudo_matrix.shape[1])]
        if pseudo_matrix.shape[1] > 0:
            pseudo_matrix.to_csv(f'pseudobulk/{celltype}.csv')
            pseudo_matrix[top_2000_columns].to_csv(f'top2000_pseudobulk/{celltype}.csv')
        print(f"处理完成 {celltype}", flush=True)
preprocess_RNA_data(adata)
