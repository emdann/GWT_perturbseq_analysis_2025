import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import pertpy
from tqdm import tqdm

adata_state = sc.read_h5ad('pbulk_test.h5ad')
design_formula = '~ log10_n_cells + donor_id + target'
n_cpus = 3
st = adata_state.obs['culture_condition'].unique()[0]

model = pertpy.tl.PyDESeq2(adata_state, design=design_formula)
model.fit(n_cpus = n_cpus, quiet=True)

all_targets = adata_state.obs['target'].unique().tolist()
all_targets.remove('NTC')

all_res_df = pd.DataFrame()

# # for t in tqdm(all_targets, desc="Testing targets"):
# all_contrasts = {t:(model.cond(target = t) - model.cond(target = 'NTC')) for t in all_targets}
# res_df = model.test_contrasts(all_contrasts[all_targets[3]], n_cpus=n_cpus)
# res_df['culture_condition'] = st
# res_df['contrast'] = t
# all_res_df = pd.concat([all_res_df, res_df])

# all_res_df = all_res_df.reset_index().drop('index', axis=1)

# all_res_df.to_csv('results.csv')
