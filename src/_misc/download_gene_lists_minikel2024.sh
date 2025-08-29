#!/bin/bash

# Base URL for raw files from the GitHub repository
BASE_URL="https://raw.githubusercontent.com/ericminikel/genetic_support/main/data/gene_lists"

# List of TSV files to download (based on the repository structure)
TSV_FILES=(
   "enzymes.tsv"
   "catalytic_receptors.tsv"
   "ion_channels.tsv"
   "kinases.tsv"
   "nuclear_receptors.tsv"
   "rhodop_gpcr.tsv"
   "transporters.tsv"
)

# Create directory to store downloaded files
for file in "${TSV_FILES[@]}"; do
    echo "Downloading $file..."
    curl -L -o "metadata/gene_lists/$file" "$BASE_URL/$file"
    
    # Check if download was successful
    if [ $? -eq 0 ]; then
        echo "✓ Successfully downloaded $file"
    else
        echo "✗ Failed to download $file"
    fi
done
