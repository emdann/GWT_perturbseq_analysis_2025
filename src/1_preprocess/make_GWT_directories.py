import os
import yaml
import argparse

from preprocess import _convert_oak_path

def create_directories(config_path, experiment_name):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    exp_config = config[experiment_name]

    base_dir = _convert_oak_path(exp_config['datadir'])
    
    # Create main directories
    exp_dir = os.path.dirname(base_dir) if base_dir.endswith('/') else base_dir
    cellranger_dir = os.path.join(exp_dir, "cellranger_outs")
    os.makedirs(cellranger_dir, exist_ok=True)
    
    # Create lane subdirectories
    for lane_id in exp_config.get('lane_ids', []):
        os.makedirs(os.path.join(cellranger_dir, lane_id), exist_ok=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create directory structure for scRNA-seq processing')
    parser.add_argument('experiment', help='Experiment name')
    parser.add_argument('--config', default='../../metadata/experiments_config.yaml', help='Path to experiments_config.yaml')
    
    args = parser.parse_args()
    create_directories(args.config, args.experiment)