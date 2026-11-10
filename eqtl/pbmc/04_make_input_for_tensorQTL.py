import os
import pandas as pd
import numpy as np
tss = pd.read_csv('/path/tss.bed',sep = ' ')
tss['gene_id_map'] = tss['gene_id'].str.replace('-', '.', regex=False)
gene_map_dict = pd.Series(tss.gene_id.values, index=tss.gene_id_map).to_dict()
celltype = ["B_intermediate","B_memory","B_naive","CD14_Mono","CD16_Mono","CD4_CTL","CD4_Naive","CD4_Proliferating",
             "CD4_TCM","CD4_TEM","CD8_Naive","CD8_TCM","CD8_TEM","cDC1","cDC2","dnT","gdT","HSPC","ILC","MAIT","NK_CD56bright",
             "NK_Proliferating","NK","pDC","Plasmablast","Treg"]
for cell in celltype:
    print(f'processing_{cell}')
    #make_bed_file_eQTL
    ## Prepare phenotype data
    norm_expr = pd.read_csv(f'/path/normal_dis_rename/{cell}.csv', index_col=0)
    filtered_map_dict = {key: gene_map_dict[key] for key in norm_expr.columns}
    norm_expr.rename(columns=filtered_map_dict, inplace=True)

    norm_expr = norm_expr.T
    norm_expr = norm_expr.reset_index()  

    pheno_file =  pd.merge(tss, norm_expr, right_on='index', left_on='gene_id', how='inner')
    pheno_file = pheno_file.drop(['index'], axis=1).rename(columns={'chr':'#chr'})
    pheno_file['#chr'] = pheno_file['#chr'].str.replace('chr','')   

    if norm_expr.shape[0] != pheno_file.shape[0]:
        raise ValueError("The number of rows in the two dataframes are not equal.")
    else:
        pheno_file.to_csv(f'/path/bed_file/{cell}.bed', sep='\t', index=False)
