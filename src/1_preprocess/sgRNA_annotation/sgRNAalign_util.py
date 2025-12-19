
import pandas as pd
from Bio.Seq import Seq
import numpy as np
import csv
import pysam

def convert_csv_to_fasta(input_csv_file, output_fasta_file):
    """
    Reads a CSV file with 'name' and 'sgRNA' columns and converts it
    into a FASTA file, taking only the last 19 characters of the sgRNA sequence, since hCRISPRi-v2 only uses 19bp protospacer.

    Args:
        input_csv_file (str): The path to the input CSV file.
        output_fasta_file (str): The path for the output FASTA file.
    """
    print(f"Reading from '{input_csv_file}'...")
    
    try:
        with open(input_csv_file, mode='r', encoding='utf-8') as csv_file:
            # Use DictReader to easily access columns by their header name
            csv_reader = csv.DictReader(csv_file)
            
            # Open the output file for writing
            with open(output_fasta_file, mode='w', encoding='utf-8') as fasta_file:
                record_count = 0
                for row in csv_reader:
                    # Get the name and sequence from the current row
                    # .strip() removes any accidental leading/trailing whitespace
                    name = row.get('name', '').strip()
                    # Get the full sequence, strip whitespace, and then take the last 19 characters
                    sequence = row.get('sgRNA', '').strip()[1:]

                    # Ensure both name and sequence exist before writing
                    if name and sequence:
                        # Write the header line, starting with '>'
                        fasta_file.write(f">{name}\n")
                        # Write the truncated sequence on the next line
                        fasta_file.write(f"{sequence}\n")
                        record_count += 1
                
                print(f"Successfully converted {record_count} records.")
                print(f"FASTA file saved as '{output_fasta_file}'")

    except FileNotFoundError:
        print(f"Error: The file '{input_csv_file}' was not found.")
        print("Please make sure the file exists and is in the correct directory.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def sam_to_dataframe(sam_file_path):
    """
    Parses a SAM file and extracts alignment information into a pandas DataFrame.

    Args:
        sam_file_path (str): The path to the SAM file.

    Returns:
        pandas.DataFrame: A DataFrame containing the alignment information.
    """
    alignments = []
    # Open the SAM file for reading ('r' mode for SAM)
    with pysam.AlignmentFile(sam_file_path, "r") as samfile:
        # Optionally, you can access the header information
        header = samfile.header
        # print("Header:", header)

        # Iterate over each alignment in the SAM file
        for read in samfile.fetch():
            alignments.append({
                'sgRNA': read.query_name,
                'chromosome': read.reference_name,
                'pos': read.reference_start,  # 0-based
                'seq': str(Seq(read.query_sequence).reverse_complement()) if read.is_reverse else read.query_sequence,
                'strand': '-' if read.is_reverse else '+',
            })

    # Create a pandas DataFrame from the list of dictionaries
    df = pd.DataFrame(alignments)
    return df

def calculate_tss(row):
    """Calculates the Transcription Start Site (TSS) list."""
    is_in_ts = isinstance(row['ts_starts'], (np.ndarray, list)) or isinstance(row['ts_starts'], (np.ndarray, list))
    if is_in_ts:
        return row['ts_starts'] if row['strand'] == '+' else row['ts_ends']
    else:
        return [row['start']] if row['strand'] == '+' else [row['end']]

def calculate_cds(row):
    """Calculates the Coding Sequence (CDS) start sites list."""
    has_cds = isinstance(row['cds_starts'], (np.ndarray, list)) or isinstance(row['cds_starts'], (np.ndarray, list))
    if has_cds:
        return row['cds_starts'] if row['strand'] == '+' else row['cds_ends']
    else:
        return []

def is_gene_nearby(sgrna_row, genes_df, distance_threshold=30000):
    """
    Checks if there is a transcription start site (TSS) within a certain distance of an sgRNA.

    Args:
        sgrna_row (pd.Series): A row from the sgRNA dataframe.
        genes_df (pd.DataFrame): The dataframe containing gene information.
        distance_threshold (int): The maximum distance to consider a gene "nearby".

    Returns:
        bool: True if a nearby gene is found, False otherwise.
    """
    # Filter genes to only those on the same chromosome
    genes_on_same_chrom = genes_df[genes_df['chromosome'] == sgrna_row['chromosome']]
    
    # Check each gene on the same chromosome
    for _, gene_row in genes_on_same_chrom.iterrows():
        # First, check for proximity to any TSS for the current gene
        for tss in gene_row['tss']:
            if abs(sgrna_row['pos'] - tss) <= distance_threshold:
                return True # Found a nearby TSS

        # If no nearby TSS was found for this gene, check if the sgRNA is within the gene body
        if gene_row['gene_start'] <= sgrna_row['pos'] <= gene_row['gene_end']:
            return True # Found sgRNA within a gene body

    # If the loop completes without finding any match, return False
    return False

def find_closest_target_info(row, distance_threshold=2000):
    sgrna_pos = row['pos']
    tss_list = row['tss']
    gene_id = row['gene_id']
    gene_name = row['gene_name']
    if not isinstance(tss_list, (np.ndarray, list)):
        return pd.Series([np.nan, np.nan, np.nan])
    distances = [abs(sgrna_pos - tss_pos) for tss_pos in tss_list]
    min_distance = min(distances)
    # Check if sgRNA is within certain distance from of tss
    if min_distance <= distance_threshold:
        return pd.Series([gene_id, gene_name, min_distance])
    else:
        return pd.Series([np.nan, np.nan, np.nan])

def find_nearby_genes(row, all_genes_df, distance_threshold=30000):
    """
    Find nearby genes within certain distance on the same chromosome.
    """
    sgrna_chrom = row['chromosome_norm']
    sgrna_pos = row['pos']

    # Filter for genes on the same chromosome
    same_chrom_genes = all_genes_df[all_genes_df['chromosome_norm'] == sgrna_chrom].copy()

    # Calculate minimum distance for each gene
    same_chrom_genes['distance'] = same_chrom_genes['tss'].apply(
        lambda tss_list: min([abs(sgrna_pos - tss) for tss in tss_list])
    )

    # Filter by distance
    nearby_genes = same_chrom_genes[same_chrom_genes['distance'] <= distance_threshold]
    
    return list(nearby_genes['gene_id'])

def find_nearest_genes(row, all_genes_df, distance_threshold=30000):
    """
    Find nearest genes within certain distance on the same chromosome.
    """
    sgrna_chrom = row['chromosome_norm']
    sgrna_pos = row['pos']
    target_gene_id = row['designed_target_gene_id']
    target_gene_name = row['designed_target_gene_name']

    # Filter for genes on the same chromosome
    same_chrom_genes = all_genes_df[all_genes_df['chromosome_norm'] == sgrna_chrom].copy()

    # Calculate minimum distance for each gene
    same_chrom_genes['distance'] = same_chrom_genes['tss'].apply(
        lambda tss_list: min([abs(sgrna_pos - tss) for tss in tss_list])
    )

    # Filter by distance
    nearby_genes = same_chrom_genes[same_chrom_genes['distance'] <= distance_threshold]

    if nearby_genes.empty:
        nearest_gene_id = np.nan
        nearest_gene_name = np.nan
        nearest_gene_dist = np.nan
    else:
        nearest_gene_id = nearby_genes.loc[nearby_genes['distance'].idxmin(), 'gene_id']
        nearest_gene_name = nearby_genes.loc[nearby_genes['distance'].idxmin(), 'gene_name']
        nearest_gene_dist = nearby_genes.loc[nearby_genes['distance'].idxmin(), 'distance']

    filtered = nearby_genes[(nearby_genes['gene_id'] != target_gene_id) & (nearby_genes['gene_name'] != target_gene_name)]
    if filtered.empty:
        nearest_nontarget_gene_id = np.nan
        nearest_nontarget_gene_name = np.nan
        nearest_nontarget_gene_dist = np.nan
    else:
        nearest_nontarget_gene_id = filtered.loc[filtered['distance'].idxmin(), 'gene_id']
        nearest_nontarget_gene_name = filtered.loc[filtered['distance'].idxmin(), 'gene_name']
        nearest_nontarget_gene_dist = filtered.loc[filtered['distance'].idxmin(), 'distance']

    # Filter by distance
    nearby_genes = same_chrom_genes[same_chrom_genes['distance'] <= distance_threshold]
    
    return pd.Series([nearest_gene_id, nearest_gene_name, nearest_gene_dist,
        nearest_nontarget_gene_id, nearest_nontarget_gene_name, nearest_nontarget_gene_dist])

def is_near_cds_start(pos, cds_starts_list, distance_threshold=2000):
    """Checks if a position is within 2000bp of any CDS start in the list."""
    try:
        for start_pos in cds_starts_list:
            if abs(pos - start_pos) <= distance_threshold:
                return True
    except TypeError:
        return False
    return False

def find_nearest_nontarget(sgrna_row, genes_df):
    """
    Finds the nearest non-target gene for an sgRNA by checking its primary
    and other alignment positions against a gene database.
    
    Args:
        sgrna_row (pd.Series): A single row from the sgrna_df_final DataFrame.
        genes_df (pd.DataFrame): The genes_df_subset DataFrame for searching.
        
    Returns:
        pd.Series: A series containing the ID, name, and distance of the 
                   nearest non-target gene.
    """
    # Get the ID of the gene this sgRNA was designed to target
    designed_target_id = sgrna_row['designed_target_gene_id']
    
    # --- 1. Collect all genomic positions to search ---
    search_locations = []
    
    # Add the primary sgRNA alignment location
    search_locations.append((sgrna_row['chromosome'], sgrna_row['pos']))
    
    # Add other alignment locations if they exist and are valid lists
    other_pos = sgrna_row['other_alignment_pos']
    other_chrom = sgrna_row['other_alignment_chromosome']
    
    if isinstance(other_pos, (np.ndarray, list)) and isinstance(other_chrom, (np.ndarray, list)):
        search_locations.extend(list(zip(other_chrom, other_pos)))

    # --- 2. Search for the closest non-target gene across all locations ---
    overall_best_dist = float('inf')
    overall_best_gene_id = np.nan
    overall_best_gene_name = np.nan

    for chrom, pos in search_locations:
        # Filter genes to the correct chromosome and exclude the designed target
        candidate_genes = genes_df[
            (genes_df['chromosome'] == chrom) &
            (genes_df['gene_id'] != designed_target_id)
        ].copy() # Use .copy() to prevent SettingWithCopyWarning
        
        # If no potential non-target genes exist on this chromosome, skip
        if candidate_genes.empty:
            continue
            
        # For each candidate gene, find the minimum distance from the sgRNA
        # position to any of its transcription start sites (TSS)
        # The lambda function efficiently handles the list of TSSs for each gene
        distances = candidate_genes['tss'].apply(
            lambda tss_list: min([abs(pos - tss) for tss in tss_list])
        )
        
        # Find the closest gene at this specific location
        if not distances.empty:
            local_min_dist = distances.min()
            
            # If this location's best is better than our overall best, update it
            if local_min_dist < overall_best_dist:
                overall_best_dist = local_min_dist
                # Get the details of the new best gene
                best_gene_series = candidate_genes.loc[distances.idxmin()]
                overall_best_gene_id = best_gene_series['gene_id']
                overall_best_gene_name = best_gene_series['gene_name']

    # Convert infinity back to NaN for cleaner output if no gene was found
    if overall_best_dist == float('inf'):
        overall_best_dist = np.nan

    return pd.Series([overall_best_gene_id, overall_best_gene_name, overall_best_dist])