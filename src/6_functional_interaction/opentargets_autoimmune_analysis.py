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
from typing import List, Dict, Set, Tuple, Optional
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

@dataclass
class AnalysisConfig:
    """Configuration parameters for the analysis"""
    # Data sources
    clustering_data_path: str = "results/clustering_nde30ntotal75_reg_and_downstream.csv"

    # Column names in the input file
    cluster_id_col: str = "cluster"
    cluster_member_col: str = "cluster_member_reg"
    downstream_member_col: str = "cluster_member_downstream"

    # Filtering parameters
    min_genetic_evidence_score: float = 0.1  # Initial OpenTargets query threshold
    genetic_evidence_score_filter: float = 0.2  # Post-filtering threshold for analysis

    # Statistical parameters
    fdr_threshold: float = 0.1
    nominal_p_threshold: float = 0.05

    # Output files
    disease_gene_associations_file: str = 'disease_gene_associations_detailed.csv'
    cluster_enrichment_file: str = 'cluster_autoimmune_enrichment_results.csv'
    cluster_disease_enrichment_file: str = 'cluster_disease_specific_enrichment_results.csv'
    downstream_enrichment_file: str = 'downstream_autoimmune_enrichment_results.csv'
    downstream_disease_enrichment_file: str = 'downstream_disease_specific_enrichment_results.csv'

    # Verbosity
    verbose: bool = True


# Genetic evidence types to include
GENETIC_EVIDENCE_TYPES = {
    'genetic_association',  # GWAS associations
    'gene_burden',          # Gene burden studies
    'somatic'               # Includes ClinVar and somatic mutations
}

# Autoimmune disease IDs from OpenTargets
AUTOIMMUNE_DISEASES = [
    "EFO_0005140",  # autoimmune disease
    "EFO_0000685",  # rheumatoid arthritis
    "MONDO_0007915",  # systemic lupus erythematosus
    "EFO_0003767",  # inflammatory bowel disease
    "MONDO_0005301",  # multiple sclerosis
    "MONDO_0005147",  # type 1 diabetes
    "EFO_0000676",  # psoriasis
    "EFO_0003898",  # ankylosing spondylitis
    "MONDO_0004979",  # asthma
    "EFO_0003779",  # Hashimotos
    "EFO_0000384",  # Crohn's disease
    "EFO_0000729",  # ulcerative colitis
    "EFO_0001060",  # celiac disease
]


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def perform_fisher_test(a: int, b: int, c: int, d: int) -> Tuple[float, float, float, float]:
    """
    Perform Fisher's exact test and calculate confidence intervals

    Args:
        a: Target genes in cluster
        b: Non-target genes in cluster
        c: Target genes not in cluster
        d: Non-target genes not in cluster

    Returns:
        odds_ratio, ci_lower, ci_upper, p_value
    """
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

    # Calculate confidence intervals for odds ratio
    # Add 0.5 to all cells if any are zero (Haldane-Anscombe correction)
    if 0 in [a, b, c, d]:
        a_adj, b_adj, c_adj, d_adj = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        se_log_or = np.sqrt(1/a_adj + 1/b_adj + 1/c_adj + 1/d_adj)
    else:
        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)

    ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se_log_or)
    ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se_log_or)

    return odds_ratio, ci_lower, ci_upper, p_value


def apply_fdr_correction(df: pd.DataFrame, p_col: str = 'p_value', fdr_col: str = 'p_adj_fdr') -> pd.DataFrame:
    """Apply FDR correction to p-values in DataFrame"""
    if len(df) == 0:
        return df

    _, df[fdr_col], _, _ = multipletests(df[p_col], method='fdr_bh')
    return df.sort_values(p_col)


def load_or_parse_gene_list(gene_list_raw) -> Set[str]:
    """Parse gene list from various formats (string, list, etc.)"""
    if isinstance(gene_list_raw, str):
        try:
            return set(eval(gene_list_raw))
        except:
            return set([gene.strip() for gene in gene_list_raw.split(',')])
    else:
        return set(gene_list_raw)


# ============================================================================
# OPENTARGETS QUERIER
# ============================================================================

class OpenTargetsQuerier:
    """Class to query OpenTargets GraphQL API for genetic evidence"""

    def __init__(self, diseases: List[str] = AUTOIMMUNE_DISEASES):
        self.base_url = "https://api.platform.opentargets.org/api/v4/graphql"
        self.diseases = diseases

    def query_genetic_evidence(self, disease_id: str, page_size: int = 1000) -> Optional[Dict]:
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
            variables = {"diseaseId": disease_id, "pageIndex": page_index, "pageSize": page_size}
            response = requests.post(self.base_url, json={"query": query, "variables": variables},
                                   headers={"Content-Type": "application/json"})

            if response.status_code != 200:
                return None

            result = response.json()
            if 'errors' in result:
                return None

            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                break

            current_rows = disease_data['associatedTargets']['rows']
            if not current_rows:
                break

            all_results.extend(current_rows)

            if len(current_rows) < page_size:
                break

            page_index += 1
            time.sleep(0.1)

        if all_results:
            return {
                'data': {
                    'disease': {
                        'id': disease_data['id'],
                        'name': disease_data['name'],
                        'associatedTargets': {'rows': all_results, 'count': len(all_results)}
                    }
                }
            }
        return None

    def get_autoimmune_genes(self, min_score: float = 0.1,
                            evidence_types: Set[str] = GENETIC_EVIDENCE_TYPES,
                            verbose: bool = True) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
        """Get all genes with genetic evidence for autoimmune disorders"""
        detailed_results = []
        disease_gene_sets = {}

        if verbose:
            print(f"Querying OpenTargets for genetic evidence (types: {evidence_types})...")

        for disease_id in self.diseases:
            if verbose:
                print(f"  Querying {disease_id}...")

            result = self.query_genetic_evidence(disease_id)
            if not result or 'data' not in result:
                continue

            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                continue

            disease_name = disease_data['name']
            disease_genes = set()

            for row in disease_data['associatedTargets']['rows']:
                target = row['target']

                # Check for genetic evidence
                best_genetic_score = 0
                evidence_types_found = []

                for datatype_score in row['datatypeScores']:
                    evidence_type = datatype_score['id']
                    evidence_score = datatype_score['score']

                    if evidence_type in evidence_types and evidence_score >= min_score:
                        evidence_types_found.append(evidence_type)
                        best_genetic_score = max(best_genetic_score, evidence_score)

                # Only include genes with genetic evidence
                if evidence_types_found:
                    detailed_results.append({
                        'disease_efo': disease_id,
                        'disease_name': disease_name,
                        'gene_symbol': target['approvedSymbol'],
                        'gene_id': target['id'],
                        'gene_name': target['approvedName'],
                        'association_score': row['score'],
                        'genetic_evidence_score': best_genetic_score,
                        'genetic_evidence_types': ','.join(evidence_types_found)
                    })
                    disease_genes.add(target['approvedSymbol'])

            disease_gene_sets[disease_name] = disease_genes
            time.sleep(0.5)

        detailed_df = pd.DataFrame(detailed_results).sort_values(
            ['disease_name', 'genetic_evidence_score'], ascending=[True, False]
        )

        if verbose:
            print(f"Found {len(detailed_df['gene_symbol'].unique())} unique genes, "
                  f"{len(detailed_df)} associations across {len(disease_gene_sets)} diseases")

        return detailed_df, disease_gene_sets


# ============================================================================
# ENRICHMENT ANALYZER
# ============================================================================

class EnrichmentAnalyzer:
    """Class to perform enrichment analysis of disease genes in clusters"""

    def __init__(self, config: AnalysisConfig):
        self.config = config
        self._load_clustering_data()
        self._load_downstream_genes()

    def _load_clustering_data(self):
        """Load clustering data from file"""
        data_path = self.config.clustering_data_path

        if data_path.endswith('.parquet'):
            self.df = pd.read_parquet(data_path)
        elif data_path.endswith('.csv'):
            self.df = pd.read_csv(data_path)
        else:
            raise ValueError("Unsupported data source format. Use .csv or .parquet")

        # Extract all unique cluster member genes (regulators)
        self.all_cluster_genes = set()
        for cluster_genes in self.df[self.config.cluster_member_col]:
            if pd.notna(cluster_genes):
                self.all_cluster_genes.update(load_or_parse_gene_list(cluster_genes))

        if self.config.verbose:
            print(f"Loaded clustering data: {self.df.shape[0]} clusters, {len(self.all_cluster_genes)} unique regulator genes")

    def _load_downstream_genes(self):
        """Load downstream genes from the same file"""
        self.downstream_genes_by_cluster = {}
        self.all_downstream_genes = set()

        # Extract downstream genes for each cluster
        for _, row in self.df.iterrows():
            cluster_id = row[self.config.cluster_id_col]
            downstream_genes_raw = row[self.config.downstream_member_col]

            if pd.notna(downstream_genes_raw):
                downstream_genes = load_or_parse_gene_list(downstream_genes_raw)
                self.downstream_genes_by_cluster[cluster_id] = downstream_genes
                self.all_downstream_genes.update(downstream_genes)

        if self.config.verbose:
            print(f"Loaded downstream genes: {len(self.all_downstream_genes)} unique genes across {len(self.downstream_genes_by_cluster)} clusters")

    def _test_enrichment_core(self, target_genes: Set[str], background_genes: Set[str],
                             cluster_genes_dict: Dict[int, Set[str]],
                             label: str = 'target', apply_fdr: bool = True) -> pd.DataFrame:
        """Core enrichment testing logic"""
        background_target = background_genes.intersection(target_genes)
        background_size = len(background_genes)
        background_target_count = len(background_target)

        if background_target_count == 0:
            return pd.DataFrame()

        results = []
        for cluster_id, cluster_genes in cluster_genes_dict.items():
            cluster_size = len(cluster_genes)
            cluster_target = cluster_genes.intersection(target_genes)
            cluster_target_count = len(cluster_target)

            # Skip clusters with no overlap with target genes
            if cluster_target_count == 0:
                continue

            a = cluster_target_count
            b = cluster_size - cluster_target_count
            c = background_target_count - cluster_target_count
            d = (background_size - background_target_count) - b

            odds_ratio, ci_lower, ci_upper, p_value = perform_fisher_test(a, b, c, d)

            results.append({
                'cluster': cluster_id,
                'cluster_size': cluster_size,
                f'{label}_count': cluster_target_count,
                f'{label}_genes': list(cluster_target),
                'odds_ratio': odds_ratio,
                'odds_ratio_ci_lower': ci_lower,
                'odds_ratio_ci_upper': ci_upper,
                'p_value': p_value,
                'enrichment_fold': (a / cluster_size) / (background_target_count / background_size)
            })

        if not results:
            return pd.DataFrame()

        results_df = pd.DataFrame(results)

        # Only apply FDR correction if requested
        if apply_fdr:
            return apply_fdr_correction(results_df)
        else:
            return results_df

    def test_enrichment(self, target_genes: Set[str]) -> pd.DataFrame:
        """Test cluster enrichment for target genes (regulators)"""
        cluster_genes_dict = {
            row[self.config.cluster_id_col]: load_or_parse_gene_list(row[self.config.cluster_member_col])
            for _, row in self.df.iterrows()
            if pd.notna(row[self.config.cluster_member_col])
        }
        return self._test_enrichment_core(target_genes, self.all_cluster_genes,
                                         cluster_genes_dict, 'autoimmune')

    def test_downstream_enrichment(self, target_genes: Set[str]) -> pd.DataFrame:
        """Test downstream gene enrichment"""
        if not self.downstream_genes_by_cluster:
            return pd.DataFrame()

        results_df = self._test_enrichment_core(target_genes, self.all_downstream_genes,
                                               self.downstream_genes_by_cluster, 'autoimmune')

        # Rename cluster_size to downstream_size for clarity
        if 'cluster_size' in results_df.columns:
            results_df.rename(columns={'cluster_size': 'downstream_size'}, inplace=True)

        return results_df

    def test_disease_specific_enrichment(self, disease_gene_sets: Dict[str, Set[str]]) -> pd.DataFrame:
        """Test disease-specific enrichment for each disease (regulators)

        Returns combined DataFrame with FDR correction applied across all diseases
        """
        cluster_genes_dict = {
            row[self.config.cluster_id_col]: load_or_parse_gene_list(row[self.config.cluster_member_col])
            for _, row in self.df.iterrows()
            if pd.notna(row[self.config.cluster_member_col])
        }

        all_results = []
        for disease_name, disease_genes in disease_gene_sets.items():
            if len(disease_genes) == 0:
                continue

            # Get results without FDR correction
            results_df = self._test_enrichment_core(disease_genes, self.all_cluster_genes,
                                                   cluster_genes_dict, 'disease_gene', apply_fdr=False)

            if len(results_df) > 0:
                results_df['disease'] = disease_name
                all_results.append(results_df)

        if not all_results:
            return pd.DataFrame()

        # Combine all results and apply FDR correction once across all tests
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df = apply_fdr_correction(combined_df)

        return combined_df

    def test_disease_specific_downstream_enrichment(self, disease_gene_sets: Dict[str, Set[str]]) -> pd.DataFrame:
        """Test disease-specific downstream enrichment

        Returns combined DataFrame with FDR correction applied across all diseases
        """
        if not self.downstream_genes_by_cluster:
            return pd.DataFrame()

        all_results = []
        for disease_name, disease_genes in disease_gene_sets.items():
            if len(disease_genes) == 0:
                continue

            # Get results without FDR correction
            results_df = self._test_enrichment_core(disease_genes, self.all_downstream_genes,
                                                   self.downstream_genes_by_cluster, 'disease_gene', apply_fdr=False)

            if len(results_df) > 0:
                results_df['disease'] = disease_name
                # Rename cluster_size to downstream_size
                if 'cluster_size' in results_df.columns:
                    results_df.rename(columns={'cluster_size': 'downstream_size'}, inplace=True)
                all_results.append(results_df)

        if not all_results:
            return pd.DataFrame()

        # Combine all results and apply FDR correction once across all tests
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df = apply_fdr_correction(combined_df)

        return combined_df


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def load_or_query_disease_genes(config: AnalysisConfig, querier: OpenTargetsQuerier) -> Tuple[Set[str], pd.DataFrame, Dict[str, Set[str]]]:
    """Load disease genes from file or query OpenTargets"""
    try:
        if config.verbose:
            print(f"Loading existing data from {config.disease_gene_associations_file}...")

        detailed_disease_gene_df = pd.read_csv(config.disease_gene_associations_file)
        detailed_disease_gene_df = detailed_disease_gene_df[
            detailed_disease_gene_df['genetic_evidence_score'] > config.genetic_evidence_score_filter
        ]

        autoimmune_genes = set(detailed_disease_gene_df['gene_symbol'].unique())

        disease_gene_sets = {
            disease: set(detailed_disease_gene_df[detailed_disease_gene_df['disease_name'] == disease]['gene_symbol'].unique())
            for disease in detailed_disease_gene_df['disease_name'].unique()
        }

        if config.verbose:
            print(f"Loaded {len(autoimmune_genes)} genes, {len(detailed_disease_gene_df)} associations")

    except FileNotFoundError:
        if config.verbose:
            print("File not found, querying OpenTargets...")

        detailed_df_raw, disease_gene_sets_raw = querier.get_autoimmune_genes(
            min_score=config.min_genetic_evidence_score, verbose=config.verbose
        )

        detailed_disease_gene_df = detailed_df_raw[
            detailed_df_raw['genetic_evidence_score'] > config.genetic_evidence_score_filter
        ]

        autoimmune_genes = set(detailed_disease_gene_df['gene_symbol'].unique())

        disease_gene_sets = {
            disease: set(detailed_disease_gene_df[detailed_disease_gene_df['disease_name'] == disease]['gene_symbol'].unique())
            for disease in detailed_disease_gene_df['disease_name'].unique()
        }

        detailed_df_raw.to_csv(config.disease_gene_associations_file, index=False)
        if config.verbose:
            print(f"Saved to {config.disease_gene_associations_file}")

    return autoimmune_genes, detailed_disease_gene_df, disease_gene_sets


def run_enrichment_analysis(analyzer: EnrichmentAnalyzer, autoimmune_genes: Set[str],
                           disease_gene_sets: Dict[str, Set[str]], config: AnalysisConfig):
    """Run all enrichment analyses and save results"""

    # 1. Overall cluster enrichment (regulators)
    if config.verbose:
        print("\n" + "="*80)
        print("CLUSTER ENRICHMENT ANALYSIS (REGULATORS)")
        print("="*80)

    results_df = analyzer.test_enrichment(autoimmune_genes)
    if len(results_df) > 0:
        results_df.to_csv(config.cluster_enrichment_file, index=False)
        if config.verbose:
            sig = len(results_df[results_df['p_adj_fdr'] < config.fdr_threshold])
            print(f"Tested {len(results_df)} clusters, {sig} significant (FDR < {config.fdr_threshold})")

    # 2. Disease-specific cluster enrichment (regulators)
    disease_results_df = analyzer.test_disease_specific_enrichment(disease_gene_sets)
    if len(disease_results_df) > 0:
        disease_results_df.to_csv(config.cluster_disease_enrichment_file, index=False)
        if config.verbose:
            n_diseases = disease_results_df['disease'].nunique()
            sig = len(disease_results_df[disease_results_df['p_adj_fdr'] < config.fdr_threshold])
            print(f"Saved {len(disease_results_df)} results across {n_diseases} diseases, {sig} significant")

    # 3. Downstream enrichment
    downstream_df = analyzer.test_downstream_enrichment(autoimmune_genes)
    if len(downstream_df) > 0:
        downstream_df.to_csv(config.downstream_enrichment_file, index=False)
        if config.verbose:
            sig = len(downstream_df[downstream_df['p_adj_fdr'] < config.fdr_threshold])
            print(f"Tested {len(downstream_df)} clusters, {sig} significant")

    # 4. Disease-specific downstream enrichment
    downstream_disease_df = analyzer.test_disease_specific_downstream_enrichment(disease_gene_sets)
    if len(downstream_disease_df) > 0:
        downstream_disease_df.to_csv(config.downstream_disease_enrichment_file, index=False)
        if config.verbose:
            n_diseases = downstream_disease_df['disease'].nunique()
            sig = len(downstream_disease_df[downstream_disease_df['p_adj_fdr'] < config.fdr_threshold])
            print(f"Saved {len(downstream_disease_df)} results across {n_diseases} diseases, {sig} significant")


def main():
    """Main analysis pipeline"""
    # Initialize configuration
    config = AnalysisConfig()

    # Initialize querier
    querier = OpenTargetsQuerier()

    # Load clustering data and downstream genes
    analyzer = EnrichmentAnalyzer(config)

    # Load or query disease genes
    autoimmune_genes, detailed_df, disease_gene_sets = load_or_query_disease_genes(config, querier)

    # Run all enrichment analyses
    run_enrichment_analysis(analyzer, autoimmune_genes, disease_gene_sets, config)

    if config.verbose:
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)


if __name__ == "__main__":
    main()
