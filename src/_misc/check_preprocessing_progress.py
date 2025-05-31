#!/usr/bin/env python3

import os
import glob
import argparse
from collections import defaultdict
import sys

# ANSI color codes
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    YELLOW = '\033[93m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'  # End color

def colored_text(text, color):
    """Return colored text for terminal output"""
    return f"{color}{text}{Colors.END}"

def find_samples_and_lanes(datadir, experiment_name):
    """
    Find all unique sample names and lane IDs by scanning the cellranger_outs directory
    for input .h5 files that should be processed by step 1.
    """
    pattern = os.path.join(datadir, experiment_name, "cellranger_outs", "*", "*.h5")
    files = glob.glob(pattern)
    
    samples_lanes = set()
    for file_path in files:
        basename = os.path.basename(file_path)
        lane_id = os.path.basename(os.path.dirname(file_path))
        
        # Extract sample name from the h5 filename
        # Assuming format like: {sample_name}_sample_filtered_feature_bc_matrix.h5
        if basename.endswith('_sample_filtered_feature_bc_matrix.h5'):
            sample_name = basename.replace('_sample_filtered_feature_bc_matrix.h5', '')
            samples_lanes.add((sample_name, lane_id))
    
    return sorted(samples_lanes)

def check_file_exists(filepath):
    """Check if a file exists and return colored status symbol"""
    if os.path.exists(filepath):
        return colored_text("✓", Colors.GREEN)
    else:
        return colored_text("✗", Colors.RED)

def check_pipeline_status(datadir, experiment_name):
    """
    Check the status of all pipeline outputs for each sample-lane combination.
    """
    samples_lanes = find_samples_and_lanes(datadir, experiment_name)
    
    if not samples_lanes:
        print(f"No samples found for experiment {experiment_name}")
        print(f"Searched in: {os.path.join(datadir, experiment_name, 'tmp', '*_CD4i_R1_Ultima.*.scRNA.h5ad')}")
        return
    
    print(colored_text(f"Pipeline Status for Experiment: {experiment_name}", Colors.BOLD + Colors.BLUE))
    print(colored_text(f"Data Directory: {datadir}", Colors.CYAN))
    print(colored_text("=" * 80, Colors.BOLD))
    
    # Header with colors
    header = (f"{colored_text('Sample Name', Colors.BOLD):<30} "
              f"{colored_text('Lane', Colors.BOLD):<18} "
              f"{colored_text('Step1-scRNA', Colors.BOLD):<22} "
              f"{colored_text('Step1-sgRNA', Colors.BOLD):<22} "
              f"{colored_text('Step2-assign', Colors.BOLD):<22} "
              f"{colored_text('Step3-postQC', Colors.BOLD):<22}")
    print(header)
    print(colored_text("-" * 80, Colors.BOLD))
    
    # Track missing files
    missing_step1_scrna = []
    missing_step1_sgrna = []
    missing_step2_assign = []
    missing_step3_postqc = []
    
    for sample_name, lane_id in samples_lanes:
        # Define file paths for each step
        step1_scrna = os.path.join(datadir, experiment_name, "tmp", 
                                  f"{sample_name}_CD4i_R1_Ultima.{lane_id}.scRNA.h5ad")
        
        step1_sgrna = os.path.join(datadir, experiment_name, 
                                  f"{sample_name}_CD4i_R1_Ultima.{lane_id}.sgRNA.h5ad")
        
        step2_assign = os.path.join(datadir, experiment_name, 
                                   f"{sample_name}_CD4i_R1_Ultima.{lane_id}.sgrna_assignment.csv")
        
        step3_postqc = os.path.join(datadir, experiment_name, "tmp", 
                                   f"{sample_name}_CD4i_R1_Ultima.{lane_id}.scRNA.postQC.h5ad")
        
        # Check file existence
        status1_scrna = check_file_exists(step1_scrna)
        status1_sgrna = check_file_exists(step1_sgrna)
        status2_assign = check_file_exists(step2_assign)
        status3_postqc = check_file_exists(step3_postqc)
        
        # Print status row with colors
        row = (f"{colored_text(sample_name, Colors.WHITE):<30} "
               f"{colored_text(lane_id, Colors.YELLOW):<18} "
               f"{status1_scrna:<22} "
               f"{status1_sgrna:<22} "
               f"{status2_assign:<22} "
               f"{status3_postqc:<22}")
        print(row)
        
        # Track missing files (check the actual file existence, not the colored status)
        sample_lane = f"{sample_name}.{lane_id}"
        if not os.path.exists(step1_scrna):
            missing_step1_scrna.append(sample_lane)
        if not os.path.exists(step1_sgrna):
            missing_step1_sgrna.append(sample_lane)
        if not os.path.exists(step2_assign):
            missing_step2_assign.append(sample_lane)
        if not os.path.exists(step3_postqc):
            missing_step3_postqc.append(sample_lane)
    
    # Summary with colors
    print("\n" + colored_text("=" * 80, Colors.BOLD))
    print(colored_text("SUMMARY", Colors.BOLD + Colors.MAGENTA))
    print(colored_text("=" * 80, Colors.BOLD))
    
    total_samples = len(samples_lanes)
    print(f"Total sample-lane combinations: {colored_text(str(total_samples), Colors.BOLD)}")
    
    # Color-coded completion stats
    step1_scrna_complete = total_samples - len(missing_step1_scrna)
    step1_sgrna_complete = total_samples - len(missing_step1_sgrna)
    step2_complete = total_samples - len(missing_step2_assign)
    step3_complete = total_samples - len(missing_step3_postqc)
    
    def format_completion(complete, total):
        if complete == total:
            return colored_text(f"{complete}/{total}", Colors.GREEN + Colors.BOLD)
        elif complete == 0:
            return colored_text(f"{complete}/{total}", Colors.RED + Colors.BOLD)
        else:
            return colored_text(f"{complete}/{total}", Colors.YELLOW + Colors.BOLD)
    
    print(f"Step 1 (scRNA) complete: {format_completion(step1_scrna_complete, total_samples)}")
    print(f"Step 1 (sgRNA) complete: {format_completion(step1_sgrna_complete, total_samples)}")
    print(f"Step 2 (assignment) complete: {format_completion(step2_complete, total_samples)}")
    print(f"Step 3 (post-QC) complete: {format_completion(step3_complete, total_samples)}")
    
    # Print missing files with colors
    if missing_step1_scrna:
        print(f"\n{colored_text(f'Missing Step 1 (scRNA) outputs ({len(missing_step1_scrna)}):', Colors.RED + Colors.BOLD)}")
        for sample in missing_step1_scrna:
            print(f"  {colored_text('- ' + sample, Colors.RED)}")
    
    if missing_step1_sgrna:
        print(f"\n{colored_text(f'Missing Step 1 (sgRNA) outputs ({len(missing_step1_sgrna)}):', Colors.RED + Colors.BOLD)}")
        for sample in missing_step1_sgrna:
            print(f"  {colored_text('- ' + sample, Colors.RED)}")
    
    if missing_step2_assign:
        print(f"\n{colored_text(f'Missing Step 2 (assignment) outputs ({len(missing_step2_assign)}):', Colors.RED + Colors.BOLD)}")
        for sample in missing_step2_assign:
            print(f"  {colored_text('- ' + sample, Colors.RED)}")
    
    if missing_step3_postqc:
        print(f"\n{colored_text(f'Missing Step 3 (post-QC) outputs ({len(missing_step3_postqc)}):', Colors.RED + Colors.BOLD)}")
        for sample in missing_step3_postqc:
            print(f"  {colored_text('- ' + sample, Colors.RED)}")
    
    if not any([missing_step1_scrna, missing_step1_sgrna, missing_step2_assign, missing_step3_postqc]):
        print(f"\n{colored_text('🎉 All pipeline steps completed for all samples!', Colors.GREEN + Colors.BOLD)}")

def main():
    parser = argparse.ArgumentParser(description="Check pipeline status for scRNA-seq data processing")
    parser.add_argument("experiment_name", help="Name of the experiment")
    parser.add_argument("--datadir", default="/oak/stanford/groups/pritch/users/emma/data/GWT/", 
                       help="Base data directory (default: /oak/stanford/groups/pritch/users/emma/data/GWT/)")
    
    args = parser.parse_args()
    
    # Check if base directory exists
    if not os.path.exists(args.datadir):
        print(f"Error: Data directory does not exist: {args.datadir}")
        sys.exit(1)
    
    experiment_dir = os.path.join(args.datadir, args.experiment_name)
    if not os.path.exists(experiment_dir):
        print(f"Error: Experiment directory does not exist: {experiment_dir}")
        sys.exit(1)
    
    check_pipeline_status(args.datadir, args.experiment_name)

if __name__ == "__main__":
    main()