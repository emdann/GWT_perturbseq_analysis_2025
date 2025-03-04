'''
Script to download sample metadata from Google Drive in machine-readable format.
'''

import pandas as pd
import yaml
import subprocess
import os

# Parse command line arguments
import argparse
parser = argparse.ArgumentParser(description='Process sample metadata')
parser.add_argument('--config', type=str, default='../../metadata/experiments_config.yaml',
                   help='Path to experiments config YAML file')
parser.add_argument('--experiment_name', type=str, default=None,
                   help='Name of specific experiment to process. If not specified, processes all experiments')
parser.add_argument('--force_metadata_download', action='store_true',
                   help='Force download of metadata file even if it exists locally')
parser.add_argument('--datadir', type=str, default='/oak/stanford/groups/pritch/users/emma/data/GWT/',
                   help='Path to data directory')
args = parser.parse_args()

# Define the output path for the xlsx file
datadir = args.datadir
xlsx_path = os.path.join(datadir, 'GWT_sample_metadata.xlsx')

# Create metadata directory
metadata_dir = os.path.join(datadir, 'sample_metadata/')
os.makedirs(metadata_dir, exist_ok=True)

if args.force_metadata_download:
    # Force download the file using rclone
    subprocess.run(['rclone', 'copy', 'gdrive:GWT_perturbseq_analysis/metadata/GWT_sample_metadata.xlsx', datadir], check=True)

try:
    # Read the metadata description sheet from the downloaded Excel file 
    df = pd.read_excel(xlsx_path, sheet_name='metadata_description')
except FileNotFoundError:
    # Download the file using rclone if not found
    subprocess.run(['rclone', 'copy', 'gdrive:GWT_perturbseq_analysis/metadata/GWT_sample_metadata.xlsx', datadir], check=True)
    df = pd.read_excel(xlsx_path, sheet_name='metadata_description')
# Read config file for all experiments
with open(args.config, 'r') as f:
    experiment_config = yaml.safe_load(f)

# Select columns to read in csv
read_cols = df['COLUMN_NAME'][df.READ_H5AD].tolist()

# Get list of experiments to process
experiments = [args.experiment_name] if args.experiment_name else experiment_config.keys()

for experiment in experiments:
    try:
        # Read each experiment sheet from the Excel file
        metadata_df = pd.read_excel(xlsx_path, sheet_name=experiment)[read_cols]
        metadata_df.to_csv(f'{metadata_dir}/GWT_sample_metadata.{experiment}.csv')
    except Exception as e:
        print(f"Error processing {experiment}: {e}")