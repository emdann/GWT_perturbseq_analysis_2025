#!/usr/bin/env python3
"""
OpenTargets Autoimmune Gene Analysis
Query OpenTargets for genes with genetic evidence for autoimmune disorders
and test for enrichment in gene clusters.
"""

import pandas as pd
import numpy as np
import requests
import json
import time
from typing import List, Dict, Set, Tuple
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

def download_google_sheet(sheet_url: str, save_path: str = None) -> pd.DataFrame:
    """Download data from Google Sheets URL and return as DataFrame"""
    # Extract sheet ID and convert to CSV download URL
    if '/spreadsheets/d/' in sheet_url:
        # Extract sheet ID from URL
        sheet_id = sheet_url.split('/spreadsheets/d/')[1].split('/')[0]
        
        # Convert to CSV download URL
        csv_url = f"https://docs.google.com/spreadsheets/d/{sheet_id}/export?format=csv&gid=0"
        
        try:
            print(f"Downloading data from Google Sheets: {sheet_id}")
            df = pd.read_csv(csv_url)
            print(f"Downloaded {df.shape[0]} rows and {df.shape[1]} columns")
            
            # Save to file if path provided
            if save_path:
                df.to_csv(save_path, index=False)
                print(f"Saved downloaded data to: {save_path}")
            
            return df
        except Exception as e:
            print(f"Error downloading Google Sheet: {e}")
            raise
    else:
        raise ValueError("Invalid Google Sheets URL format")

class OpenTargetsQuerier:
    """Class to query OpenTargets GraphQL API for genetic evidence"""
    
    def __init__(self):
        self.base_url = "https://api.platform.opentargets.org/api/v4/graphql"
        self.autoimmune_diseases = [
            "EFO_0005140",  # autoimmune disease
            "EFO_0000685",  # rheumatoid arthritis (3822 associations)
            "MONDO_0007915",  # systemic lupus erythematosus
            "EFO_0003767",  # inflammatory bowel disease (8551 associations)
            "MONDO_0005301",  # multiple sclerosis
            "MONDO_0005147",  # type 1 diabetes
            "EFO_0000676",  # psoriasis
            "EFO_0003898",  # ankylosing spondylitis (1449 associations)
            "MONDO_0004979",  # asthma
            "EFO_0003779", # Hashimotos
            "EFO_0000384",  # Crohn's disease (6497 associations) - part of IBD
            "EFO_0000729",  # ulcerative colitis (9068 associations) - part of IBD
            "EFO_0001060",  # celiac disease (1439 associations)
        ]
    
    def query_genetic_evidence(self, disease_id: str, page_size: int = 1000) -> Dict:
        """Query genetic evidence for a specific disease using pagination"""
        
        query = """
        query GeneticEvidence($diseaseId: String!, $pageIndex: Int!, $pageSize: Int!) {
          disease(efoId: $diseaseId) {
            id
            name
            associatedTargets(page: {index: $pageIndex, size: $pageSize}) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """
        
        all_results = []
        page_index = 0
        
        while True:
            variables = {
                "diseaseId": disease_id,
                "pageIndex": page_index,
                "pageSize": page_size
            }
            
            response = requests.post(
                self.base_url,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code != 200:
                print(f"HTTP error querying {disease_id}: {response.status_code}")
                try:
                    error_detail = response.json()
                    print(f"Error details: {error_detail}")
                except:
                    print(f"Response text: {response.text}")
                return None
            
            result = response.json()
            if 'errors' in result:
                print(f"GraphQL errors for {disease_id}: {result['errors']}")
                return None
                
            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                break
                
            associated_targets = disease_data['associatedTargets']
            current_rows = associated_targets['rows']
            
            if not current_rows:
                break
                
            all_results.extend(current_rows)
            
            # If we got fewer rows than requested, we're done
            if len(current_rows) < page_size:
                break
                
            page_index += 1
            
            # Add small delay between requests
            time.sleep(0.1)
        
        # Return result in same format as original
        if all_results:
            return {
                'data': {
                    'disease': {
                        'id': disease_data['id'],
                        'name': disease_data['name'],
                        'associatedTargets': {
                            'rows': all_results,
                            'count': len(all_results)
                        }
                    }
                }
            }
        else:
            return None
    
    def get_autoimmune_genes(self, min_score: float = 0.1) -> Tuple[Set[str], Dict[str, float], pd.DataFrame, Dict[str, Set[str]]]:
        """Get all genes with genetic evidence for autoimmune disorders
        
        Only includes genes with evidence from:
        - genetic_association (GWAS)
        - gene_burden (Gene Burden)
        - somatic (ClinVar and similar)
        """
        
        # Define the genetic evidence types we want
        genetic_evidence_types = {
            'genetic_association',  # GWAS associations
            'gene_burden',          # Gene burden studies  
            'somatic'               # Includes ClinVar and somatic mutations
        }
        
        all_genes = set()
        gene_scores = {}
        detailed_results = []
        disease_gene_sets = {}
        
        print("Querying OpenTargets for autoimmune genetic evidence...")
        print(f"Filtering for evidence types: {genetic_evidence_types}")
        
        for disease_id in self.autoimmune_diseases:
            print(f"  Querying {disease_id}...")
            
            result = self.query_genetic_evidence(disease_id)
            if not result or 'data' not in result:
                continue
                
            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                continue
                
            print(f"    Found {len(disease_data['associatedTargets']['rows'])} associations")
            
            disease_name = disease_data['name']
            disease_genes_above_threshold = set()
            genetic_evidence_count = 0
            
            for row in disease_data['associatedTargets']['rows']:
                target = row['target']
                overall_score = row['score']
                gene_symbol = target['approvedSymbol']
                gene_id = target['id']
                gene_name = target['approvedName']
                
                # Check if gene has any of our desired genetic evidence types
                has_genetic_evidence = False
                best_genetic_score = 0
                evidence_types_found = []
                
                for datatype_score in row['datatypeScores']:
                    evidence_type = datatype_score['id']
                    evidence_score = datatype_score['score']
                    
                    if evidence_type in genetic_evidence_types and evidence_score >= min_score:
                        has_genetic_evidence = True
                        evidence_types_found.append(evidence_type)
                        if evidence_score > best_genetic_score:
                            best_genetic_score = evidence_score
                
                # Only include genes with genetic evidence
                if has_genetic_evidence:
                    genetic_evidence_count += 1
                    
                    # Add to detailed results
                    detailed_results.append({
                        'disease_efo': disease_id,
                        'disease_name': disease_name,
                        'gene_symbol': gene_symbol,
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'association_score': overall_score,
                        'genetic_evidence_score': best_genetic_score,
                        'genetic_evidence_types': ','.join(evidence_types_found)
                    })
                    
                    # Add to gene sets
                    all_genes.add(gene_symbol)
                    disease_genes_above_threshold.add(gene_symbol)
                    
                    # Keep track of best genetic score for each gene
                    if gene_symbol not in gene_scores or best_genetic_score > gene_scores[gene_symbol]:
                        gene_scores[gene_symbol] = best_genetic_score
            
            # Store disease-specific gene set
            disease_gene_sets[disease_name] = disease_genes_above_threshold
            
            print(f"    Genes with genetic evidence: {genetic_evidence_count}")
            
            # Add small delay to be respectful to API
            time.sleep(0.5)
        
        # Create detailed DataFrame
        detailed_df = pd.DataFrame(detailed_results).sort_values(
            ['disease_name', 'genetic_evidence_score'], 
            ascending=[True, False]
        )
        
        print(f"Found {len(all_genes)} unique genes with genetic evidence (score >= {min_score})")
        print(f"Total disease-gene associations with genetic evidence: {len(detailed_df)}")
        print(f"Disease-specific gene sets created for {len(disease_gene_sets)} diseases")
        
        return all_genes, gene_scores, detailed_df, disease_gene_sets


class EnrichmentAnalyzer:
    """Class to perform enrichment analysis of autoimmune genes in clusters"""
    
    def __init__(self, clustering_data_source: str):
        """Load clustering data from file path or Google Sheets URL"""
        if clustering_data_source.startswith('https://docs.google.com/spreadsheets/'):
            # Download from Google Sheets and save locally
            self.df = download_google_sheet(clustering_data_source, save_path='cluster_data_downloaded.csv')
        elif clustering_data_source.endswith('.parquet'):
            # Load from local parquet file
            self.df = pd.read_parquet(clustering_data_source)
        elif clustering_data_source.endswith('.csv'):
            # Load from local CSV file
            self.df = pd.read_csv(clustering_data_source)
        else:
            raise ValueError("Unsupported data source format. Use .parquet, .csv file path, or Google Sheets URL")
        
        print(f"Loaded clustering data: {self.df.shape[0]} clusters")
        
        # Extract all unique genes across clusters
        self.all_cluster_genes = set()
        for cluster_genes in self.df['cluster_member']:
            # Handle string representation of lists if needed
            if isinstance(cluster_genes, str):
                # If stored as string representation of list, eval it (be careful with this)
                try:
                    cluster_genes = eval(cluster_genes)
                except:
                    # If that fails, try splitting by comma
                    cluster_genes = [gene.strip() for gene in cluster_genes.split(',')]
            self.all_cluster_genes.update(cluster_genes)
        
        print(f"Total unique genes in clusters: {len(self.all_cluster_genes)}")
    
    def test_enrichment(self, autoimmune_genes: Set[str]) -> pd.DataFrame:
        """Test each cluster for enrichment of autoimmune genes"""
        
        results = []
        
        # Background: all genes used for clustering (universe of all possible genes in clustering)
        background_genes = self.all_cluster_genes
        background_autoimmune = background_genes.intersection(autoimmune_genes)
        background_size = len(background_genes)
        background_autoimmune_count = len(background_autoimmune)
        
        print(f"Background (clustering universe): {background_size} total genes, {background_autoimmune_count} autoimmune")
        
        for idx, row in self.df.iterrows():
            cluster_id = row['cluster']
            cluster_genes_raw = row['cluster_member']
            
            # Handle string representation of lists if needed
            if isinstance(cluster_genes_raw, str):
                try:
                    cluster_genes = set(eval(cluster_genes_raw))
                except:
                    cluster_genes = set([gene.strip() for gene in cluster_genes_raw.split(',')])
            else:
                cluster_genes = set(cluster_genes_raw)
            
            cluster_size = len(cluster_genes)
            
            # Find autoimmune genes in this cluster
            cluster_autoimmune = cluster_genes.intersection(autoimmune_genes)
            cluster_autoimmune_count = len(cluster_autoimmune)
            
            # Fisher's exact test
            # Contingency table:
            # [cluster_autoimmune, cluster_not_autoimmune]
            # [background_autoimmune - cluster_autoimmune, background_not_autoimmune - cluster_not_autoimmune]
            
            a = cluster_autoimmune_count  # autoimmune in cluster
            b = cluster_size - cluster_autoimmune_count  # not autoimmune in cluster
            c = background_autoimmune_count - cluster_autoimmune_count  # autoimmune not in cluster
            d = (background_size - background_autoimmune_count) - b  # not autoimmune not in cluster
            
            if a > 0:  # Only test clusters with at least one autoimmune gene
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
                
                # Calculate confidence intervals for odds ratio
                # Using standard error of log odds ratio
                se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
                ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se_log_or)
                ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se_log_or)
                
                results.append({
                    'cluster': cluster_id,
                    'cluster_size': cluster_size,
                    'autoimmune_count': cluster_autoimmune_count,
                    'autoimmune_genes': list(cluster_autoimmune),
                    'odds_ratio': odds_ratio,
                    'odds_ratio_ci_lower': ci_lower,
                    'odds_ratio_ci_upper': ci_upper,
                    'p_value': p_value,
                    'enrichment_fold': (a / cluster_size) / (background_autoimmune_count / background_size) if background_autoimmune_count > 0 else 0
                })
        if not results:
            print("No clusters with autoimmune genes found!")
            return pd.DataFrame()
        
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        _, results_df['p_adj_fdr'], _, _ = multipletests(
            results_df['p_value'], 
            method='fdr_bh'
        )
        
        # Sort by p-value
        results_df = results_df.sort_values('p_value')
        
        return results_df
    
    def test_disease_specific_enrichment(self, disease_gene_sets: Dict[str, Set[str]]) -> Dict[str, pd.DataFrame]:
        """Test each cluster for enrichment of genes from each disease separately"""
        
        disease_results = {}
        
        for disease_name, disease_genes in disease_gene_sets.items():
            if len(disease_genes) == 0:
                continue
                
            print(f"\nTesting enrichment for {disease_name} ({len(disease_genes)} genes)...")
            
            results = []
            
            # Background: all genes used for clustering (universe of all possible genes in clustering)
            background_genes = self.all_cluster_genes
            background_disease = background_genes.intersection(disease_genes)
            background_size = len(background_genes)
            background_disease_count = len(background_disease)
            
            if background_disease_count == 0:
                print(f"  No {disease_name} genes found in clustering universe, skipping...")
                continue
                
            print(f"  Background (clustering universe): {background_size} total genes, {background_disease_count} {disease_name} genes")
            
            for idx, row in self.df.iterrows():
                cluster_id = row['cluster']
                cluster_genes_raw = row['cluster_member']
                
                # Handle string representation of lists if needed
                if isinstance(cluster_genes_raw, str):
                    try:
                        cluster_genes = set(eval(cluster_genes_raw))
                    except:
                        cluster_genes = set([gene.strip() for gene in cluster_genes_raw.split(',')])
                else:
                    cluster_genes = set(cluster_genes_raw)
                
                cluster_size = len(cluster_genes)
                
                # Find disease genes in this cluster
                cluster_disease = cluster_genes.intersection(disease_genes)
                cluster_disease_count = len(cluster_disease)
                
                # Only test clusters with at least one disease gene
                if cluster_disease_count > 0:
                    # Fisher's exact test
                    a = cluster_disease_count  # disease genes in cluster
                    b = cluster_size - cluster_disease_count  # not disease genes in cluster
                    c = background_disease_count - cluster_disease_count  # disease genes not in cluster
                    d = (background_size - background_disease_count) - b  # not disease genes not in cluster
                    
                    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
                    
                    # Calculate confidence intervals for odds ratio
                    # Using standard error of log odds ratio
                    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
                    ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se_log_or)
                    ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se_log_or)
                    
                    results.append({
                        'disease': disease_name,
                        'cluster': cluster_id,
                        'cluster_size': cluster_size,
                        'disease_gene_count': cluster_disease_count,
                        'disease_genes': list(cluster_disease),
                        'odds_ratio': odds_ratio,
                        'odds_ratio_ci_lower': ci_lower,
                        'odds_ratio_ci_upper': ci_upper,
                        'p_value': p_value,
                        'enrichment_fold': (a / cluster_size) / (background_disease_count / background_size) if background_disease_count > 0 else 0
                    })
            
            if not results:
                print(f"  No clusters with {disease_name} genes found!")
                continue
            
            results_df = pd.DataFrame(results)
            
            # Multiple testing correction
            _, results_df['p_adj_fdr'], _, _ = multipletests(
                results_df['p_value'], 
                method='fdr_bh'
            )
            
            # Sort by p-value
            results_df = results_df.sort_values('p_value')
            
            disease_results[disease_name] = results_df
            
            # Print summary
            significant_clusters = results_df[results_df['p_adj_fdr'] < 0.1]
            nominal_sig = results_df[results_df['p_value'] < 0.05]
            
            print(f"  Clusters tested: {len(results_df)}")
            print(f"  Significantly enriched (FDR < 0.1): {len(significant_clusters)}")
            print(f"  Nominally significant (p < 0.05): {len(nominal_sig)}")
        
        return disease_results


def main():
    """Main analysis pipeline"""
    
    # Initialize
    querier = OpenTargetsQuerier()
    
    # Use Google Sheets URL for clustering data
    clustering_url = "https://docs.google.com/spreadsheets/d/1qNfDOBRoNc6U0DVOlD3E5I0iFmxXRo-Ma6UXUOtCglI/edit?pli=1&gid=0#gid=0"
    analyzer = EnrichmentAnalyzer(clustering_url)
    
    # Load existing autoimmune gene data instead of querying OpenTargets
    try:
        print("Loading existing autoimmune gene data...")
        detailed_disease_gene_df = pd.read_csv('disease_gene_associations_detailed.csv')
        detailed_disease_gene_df = detailed_disease_gene_df[detailed_disease_gene_df['genetic_evidence_score'] > 0.2]
        autoimmune_df = pd.read_csv('autoimmune_genes_opentargets.csv')
        
        autoimmune_genes = set(autoimmune_df['gene'].tolist())
        gene_scores = dict(zip(autoimmune_df['gene'], autoimmune_df['max_score']))
        
        # Reconstruct disease gene sets
        disease_gene_sets = {}
        for disease in detailed_disease_gene_df['disease_name'].unique():
            disease_genes = detailed_disease_gene_df[
                detailed_disease_gene_df['disease_name'] == disease
            ]['gene_symbol'].unique()
            disease_gene_sets[disease] = set(disease_genes)
        
        print(f"Loaded {len(autoimmune_genes)} unique autoimmune genes")
        print(f"Loaded {len(detailed_disease_gene_df)} disease-gene associations") 
        print(f"Loaded {len(disease_gene_sets)} disease-specific gene sets")
        
    except FileNotFoundError as e:
        print(f"Saved data not found ({e}), querying OpenTargets...")
        # Get autoimmune genes from OpenTargets
        autoimmune_genes, gene_scores, detailed_disease_gene_df, disease_gene_sets = querier.get_autoimmune_genes(min_score=0.1)
        detailed_disease_gene_df = detailed_disease_gene_df[detailed_disease_gene_df['genetic_evidence_score'] > 0.2]
    
    # Save detailed disease-gene table
    detailed_disease_gene_df.to_csv('disease_gene_associations_detailed.csv', index=False)
    print(f"Saved {len(detailed_disease_gene_df)} disease-gene associations to disease_gene_associations_detailed.csv")
    
    # Perform overall enrichment analysis (union of all autoimmune genes)
    print("\n" + "="*80)
    print("OVERALL AUTOIMMUNE ENRICHMENT ANALYSIS (Union of all diseases)")
    print("="*80)
    results_df = analyzer.test_enrichment(autoimmune_genes)
    
    if len(results_df) > 0:
        # Save results
        results_df.to_csv('cluster_autoimmune_enrichment_results.csv', index=False)
        print(f"Saved overall enrichment results to cluster_autoimmune_enrichment_results.csv")
        
        significant_clusters = results_df[results_df['p_adj_fdr'] < 0.1]
        nominal_sig = results_df[results_df['p_value'] < 0.05]
        
        print(f"Clusters tested: {len(results_df)}")
        print(f"Significantly enriched clusters (FDR < 0.1): {len(significant_clusters)}")
        print(f"Nominally significant (p < 0.05): {len(nominal_sig)}")
        
        if len(nominal_sig) > 0:
            print(f"\nTop results (p < 0.05):")
            print(nominal_sig.head(8)[['cluster', 'cluster_size', 'autoimmune_count', 
                                     'enrichment_fold', 'p_value', 'p_adj_fdr']].round(4))
    else:
        print("No overall enrichment results to report")
    
    # Perform disease-specific enrichment analysis
    print("\n" + "="*80)
    print("DISEASE-SPECIFIC ENRICHMENT ANALYSIS")
    print("="*80)
    
    disease_results = analyzer.test_disease_specific_enrichment(disease_gene_sets)
    
    # Save disease-specific results
    all_disease_results = []
    for disease_name, disease_df in disease_results.items():
        disease_df_copy = disease_df.copy()
        disease_df_copy['disease'] = disease_name
        all_disease_results.append(disease_df_copy)
    
    if all_disease_results:
        combined_disease_results = pd.concat(all_disease_results, ignore_index=True)
        combined_disease_results = combined_disease_results[['disease', 'cluster', 'cluster_size', 
                                                           'disease_gene_count', 'disease_genes', 
                                                           'odds_ratio', 'odds_ratio_ci_lower', 'odds_ratio_ci_upper', 'p_value', 'enrichment_fold', 'p_adj_fdr']]
        combined_disease_results.to_csv('cluster_disease_specific_enrichment_results.csv', index=False)
        print(f"Saved disease-specific enrichment results to cluster_disease_specific_enrichment_results.csv")
    else:
        print("No disease-specific enrichment results to report")


if __name__ == "__main__":
    main()