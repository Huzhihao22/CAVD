import os
import pandas as pd
import tensorqtl
import numpy as np
from tensorqtl import genotypeio, cis, trans
pr = genotypeio.PlinkReader("/path_to_genome")
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
genotype_df.to_parquet('/path/genome_pqt/sample_genotype.parquet')
variant_df.to_parquet('/path/genome_pqt/sample_variant_df.parquet')
tss = pd.read_csv('/path/tss.bed',sep = ' ')
tss['gene_id_map'] = tss['gene_id'].str.replace('-', '.', regex=False)
gene_map_dict = pd.Series(tss.gene_id.values, index=tss.gene_id_map).to_dict()
celltype = ["Endo_C1_NETA1","Endo_C2_HAPLN1","Endo_C4_PI16","M_DC_CD1C","M_Macro_C3","M_Macro_HSPA1A",
             "M_Macro_IL1B","M_Macro_MRC1","M_Macro_NEAT1",
             "M_Macro_SPP1","M_Macro_TXNIP","M_Mast_TPSB2","NK_NCAM1",
             "T_CD4_Memory_GPR183","T_CD4_TH17_CCL20","T_CD8_EM_CD69","T_CD8_EM_DUSP4",
             "T_STR_HSPA1A","T_Trans_TXNIP","VIC_C0_POSTN","VIC_C1_MYH11",
             "VIC_C2_NEAT1","VIC_C3_RPS11","VIC_C4_C7","VIC_C5_ISG15",
             "VIC_C6_PRELP","VIC_C7_LIF","VIC_C9_S100B"]
for cell in celltype:
    print(f'processing_{cell}')
    norm_expr = pd.read_csv(f'/path/{cell}.csv', index_col=0)
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
