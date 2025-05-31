'''
Script to correct failed demultiplexing when the wrong config.csv was passed to cellranger
'''
import scanpy as sc
import numpy as np
import pandas as pd
import os,sys
import json
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Correct failed demultiplexing from cellranger')
    parser.add_argument('input_dir', type=str, help='Path to input directory containing cellranger output')
    parser.add_argument('--datadir', type=str, default='/mnt/oak/users/emma/data/GWT/',
                        help='Base data directory (default: %(default)s)')
    parser.add_argument('--experiment', type=str, default='CD4iR1_Psomagen',
                        help='Experiment name (default: %(default)s)')
    return parser.parse_args()

SAMPLE_BARCODE_MAP = {
    'CD4i_R1_D1_Rest': ['BC001', 'BC002', 'BC003', 'BC004'],
    'CD4i_R1_D2_Rest': ['BC005', 'BC006', 'BC007', 'BC008'],
    'CD4i_R1_D1_Stim8hr': ['BC009', 'BC010', 'BC011', 'BC012'],
    'CD4i_R1_D2_Stim8hr': ['BC013', 'BC014', 'BC015', 'BC016']
}

def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    lane_name = input_dir.name

    # Read cells per tag mapping
    with open(input_dir / 'multi/cells_per_tag.json', 'r') as f:
        cells_per_tag = json.load(f)

    # Read cellranger output
    ad = sc.read_10x_h5(input_dir / f'per_sample_outs/{lane_name}/count/sample_filtered_feature_bc_matrix.h5', 
                        gex_only=False)

    # Assign tags to cells based on cells_per_tag dictionary
    ad.obs['tag'] = np.nan
    for tag, barcodes in cells_per_tag.items():
        ad.obs.loc[barcodes, 'tag'] = tag
    assert ad.obs['tag'].isna().sum() == 0, f"Found {ad.obs['tag'].isna().sum()} cells with missing tag assignments"

    # Assign sample_id based on tag and sample_barcode_map
    ad.obs['sample_id'] = np.nan
    for sample_id, barcodes in SAMPLE_BARCODE_MAP.items():
        mask = ad.obs['tag'].isin(barcodes)
        ad.obs.loc[mask, 'sample_id'] = sample_id
    assert ad.obs['sample_id'].isna().sum() == 0, f"Found {ad.obs['sample_id'].isna().sum()} cells with missing sample assignments"

    # Create output directory
    output_dir = Path(args.datadir) / args.experiment / 'cellranger_outs' / lane_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create separate AnnData objects for each sample and save as h5 files
    for sample_id in SAMPLE_BARCODE_MAP.keys():
        mask = ad.obs['sample_id'] == sample_id
        sample_adata = ad[mask].copy()
        sample_adata.obs = pd.DataFrame(index=sample_adata.obs.index)
        
        # Save as h5 file
        output_file = f"{sample_id}_sample_filtered_feature_bc_matrix.h5"
        sample_adata.write_h5ad(output_dir / output_file)
        print(f"Saved {output_file}")

if __name__ == '__main__':
    main()