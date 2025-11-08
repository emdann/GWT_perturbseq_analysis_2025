#!/bin/bash

# =============================================================================
# GWT PerturbSeq Analysis - Plot Sync Script
# Syncs plot files from src/ subfolders to Google Drive using rclone
# =============================================================================

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
SRC_DIR="${REPO_ROOT}/src"

# Google Drive configuration
GDRIVE_REMOTE="gdrive"  # Change this to your rclone remote name
GDRIVE_BASE_PATH="GWT_perturbseq_analysis/plots"

# Logging
LOG_DIR="${SCRIPT_DIR}/logs"
LOG_FILE="${LOG_DIR}/plot_sync_$(date +%Y%m%d_%H%M%S).log"
SUMMARY_LOG="${LOG_DIR}/sync_summary.log"

# rclone options
RCLONE_OPTIONS="--checksum --ignore-times --progress --log-level INFO"
DRY_RUN=false
BANDWIDTH_LIMIT=""

# =============================================================================
# Functions
# =============================================================================

log() {
    local level="$1"
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $*" | tee -a "$LOG_FILE"
}

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Sync plot files from GWT_perturbseq_analysis/src/ to Google Drive

OPTIONS:
    -r, --remote NAME       Google Drive remote name (default: gdrive)
    -p, --path PATH         Base path on Google Drive (default: GWT_perturbseq_analysis/plots)
    -n, --dry-run          Show what would be synced without actually syncing
    -b, --bandwidth LIMIT  Bandwidth limit (e.g., 10M, 100k)
    -v, --verbose          Verbose output
    -h, --help             Show this help message

EXAMPLES:
    $0                                    # Basic sync
    $0 --dry-run                         # Preview what would be synced
    $0 --remote my_gdrive --verbose      # Use custom remote with verbose output
    $0 --bandwidth 50M                   # Limit bandwidth to 50MB/s

REQUIREMENTS:
    - rclone must be installed and configured with a Google Drive remote
    - Run 'rclone config' to set up your Google Drive connection first

EOF
}

check_requirements() {
    log "INFO" "Checking requirements..."
    
    # Check if rclone is installed
    if ! command -v rclone &> /dev/null; then
        log "ERROR" "rclone is not installed. Please install it first."
        exit 1
    fi
    
    # Check if the remote exists
    if ! rclone listremotes | grep -q "^${GDRIVE_REMOTE}:$"; then
        log "ERROR" "rclone remote '${GDRIVE_REMOTE}' not found."
        log "INFO" "Available remotes:"
        rclone listremotes | sed 's/^/  /'
        log "INFO" "Run 'rclone config' to set up a Google Drive remote."
        exit 1
    fi
    
    # Check if src directory exists
    if [[ ! -d "$SRC_DIR" ]]; then
        log "ERROR" "Source directory not found: $SRC_DIR"
        exit 1
    fi
    
    log "INFO" "Requirements check passed"
}

create_log_dir() {
    mkdir -p "$LOG_DIR"
    log "INFO" "Logging to: $LOG_FILE"
}

count_plots() {
    local dir="$1"
    find "$dir" -type f \( -name "*.pdf" -o -name "*.png" \) 2>/dev/null | wc -l
}

sync_analysis_folder() {
    local analysis_name="$1"
    local src_path="$2"
    local dest_path="${GDRIVE_REMOTE}:${GDRIVE_BASE_PATH}/${analysis_name}"
    
    local plot_count
    plot_count=$(count_plots "$src_path")
    
    if [[ $plot_count -eq 0 ]]; then
        log "INFO" "No plots found in $analysis_name, skipping"
        return 0
    fi
    
    log "INFO" "Syncing $analysis_name ($plot_count plots): $src_path -> $dest_path"
    
    local sync_cmd="rclone sync \"$src_path\" \"$dest_path\" --include \"*.pdf\" --include \"*.png\" $RCLONE_OPTIONS"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        sync_cmd="$sync_cmd --dry-run"
        log "INFO" "[DRY RUN] Would execute: $sync_cmd"
    fi
    
    if [[ -n "$BANDWIDTH_LIMIT" ]]; then
        sync_cmd="$sync_cmd --bwlimit $BANDWIDTH_LIMIT"
    fi
    
    # Execute the sync command
    if eval "$sync_cmd" 2>&1 | tee -a "$LOG_FILE"; then
        log "INFO" "Successfully synced $analysis_name"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - $analysis_name: $plot_count plots synced" >> "$SUMMARY_LOG"
        return 0
    else
        log "ERROR" "Failed to sync $analysis_name"
        return 1
    fi
}

main_sync() {
    log "INFO" "Starting plot sync from $SRC_DIR"
    log "INFO" "Destination: ${GDRIVE_REMOTE}:${GDRIVE_BASE_PATH}"
    
    local failed_syncs=0
    local total_plots=0
    
    # Define analysis folders and their destinations
    declare -A analysis_folders=(
        ["DE_analysis"]="3_DE_analysis"
        ["polarization_signatures"]="4_polarization_signatures"
        ["functional_interaction"]="6_functional_interaction"
        ["1k1k_analysis"]="7_1k1k_analysis"
        ["perturbation_prediction_LM"]="5_perturbation_prediction_LM"
        ["preprocess_results"]="1_preprocess/results"
    )
    
    # Sync each analysis folder
    for dest_name in "${!analysis_folders[@]}"; do
        local src_folder="${analysis_folders[$dest_name]}"
        local full_src_path="${SRC_DIR}/${src_folder}"
        
        if [[ -d "$full_src_path" ]]; then
            local folder_plots
            folder_plots=$(count_plots "$full_src_path")
            total_plots=$((total_plots + folder_plots))
            
            if ! sync_analysis_folder "$dest_name" "$full_src_path" "$dest_name"; then
                failed_syncs=$((failed_syncs + 1))
            fi
        else
            log "WARN" "Analysis folder not found: $full_src_path"
        fi
    done
    
    # Summary
    log "INFO" "Sync completed!"
    log "INFO" "Total plots found: $total_plots"
    log "INFO" "Failed syncs: $failed_syncs"
    
    if [[ $failed_syncs -gt 0 ]]; then
        log "ERROR" "Some syncs failed. Check the log for details."
        return 1
    else
        log "INFO" "All syncs completed successfully!"
        return 0
    fi
}

# =============================================================================
# Parse command line arguments
# =============================================================================

while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--remote)
            GDRIVE_REMOTE="$2"
            shift 2
            ;;
        -p|--path)
            GDRIVE_BASE_PATH="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -b|--bandwidth)
            BANDWIDTH_LIMIT="$2"
            shift 2
            ;;
        -v|--verbose)
            RCLONE_OPTIONS="$RCLONE_OPTIONS --log-level DEBUG"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log "ERROR" "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# =============================================================================
# Main execution
# =============================================================================

main() {
    create_log_dir
    log "INFO" "=== GWT PerturbSeq Plot Sync Started ==="
    log "INFO" "Repository root: $REPO_ROOT"
    log "INFO" "Source directory: $SRC_DIR"
    log "INFO" "Remote: $GDRIVE_REMOTE"
    log "INFO" "Destination path: $GDRIVE_BASE_PATH"
    log "INFO" "Dry run: $DRY_RUN"
    
    check_requirements
    
    if main_sync; then
        log "INFO" "=== Sync completed successfully ==="
        exit 0
    else
        log "ERROR" "=== Sync completed with errors ==="
        exit 1
    fi
}

# Run main function
main "$@"