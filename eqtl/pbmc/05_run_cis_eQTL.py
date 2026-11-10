import os
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post,genotypeio
import argparse
parser = argparse.ArgumentParser(description='Define arguments')
parser.add_argument('--cell', type=str, required=True)
args = parser.parse_args()
cell = args.cell
print(f'processing_{cell}')

device = torch.device("cpu")
## Use tensorqtl
nominal_out_dir = f'/path/eQTL_result/{cell}/'
if not os.path.exists(nominal_out_dir):
    os.makedirs(nominal_out_dir)
    print('dir created')
#read_genotype_df and vatiant_df
pr = genotypeio.PlinkReader("/path_to_genome")
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
#read_expression_file and cov
expression_bed = f'/path/bed_file/{cell}.bed'
covariates_file = f'/path/peer/factor/{cell}.csv'
# load phenotypes and covariates
# align phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
del phenotype_df[phenotype_df.columns[0]]
covariates_df = pd.read_csv(covariates_file, sep=',', index_col=1)
del covariates_df[covariates_df.columns[0]] 
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, maf_threshold=0.05, seed=1234)
cis_df['celltype'] = cell
cis_df.to_csv(f'{nominal_out_dir}all_lead_perm.csv')
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, 'all', covariates_df, maf_threshold=0.05, output_dir=nominal_out_dir)
print(f'{cell}_done')
