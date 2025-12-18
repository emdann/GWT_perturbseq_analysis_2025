#!/usr/bin/env python3
"""
Download Disease Gene Associations from OpenTargets
Query OpenTargets API for genes associated with specified diseases
and save detailed association scores.
"""

import pandas as pd
import requests
import time
import argparse
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass
import warnings
warnings.filterwarnings('ignore')


# ============================================================================
# CONFIGURATION
# ============================================================================

@dataclass
class DownloadConfig:
    """Configuration parameters for downloading disease genes"""
    # Filtering parameters
    min_genetic_evidence_score: float = 0.1  # Minimum threshold for genetic evidence

    # Output file
    output_file: str = 'disease_gene_associations_detailed.csv'

    # API parameters
    page_size: int = 1000
    request_delay: float = 0.5  # seconds between requests

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
    "EFO_0000274", # atopic eczema
    # Neg controls
    "EFO_0001645", # CAD
    "EFO_0001365", # macular degeneration
    "EFO_0000249", # alzheimers
    "EFO_0003884" # chronic kidney disease
]


# ============================================================================
# OPENTARGETS QUERIER
# ============================================================================

class OpenTargetsQuerier:
    """Class to query OpenTargets GraphQL API for genetic evidence"""

    def __init__(self, diseases: List[str], config: DownloadConfig = None):
        """
        Initialize the querier

        Args:
            diseases: List of disease EFO/MONDO IDs to query
            config: Configuration object
        """
        self.base_url = "https://api.platform.opentargets.org/api/v4/graphql"
        self.diseases = diseases
        self.config = config or DownloadConfig()

    def query_genetic_evidence(self, disease_id: str) -> Optional[Dict]:
        """
        Query genetic evidence for a specific disease using pagination

        Args:
            disease_id: Disease EFO/MONDO ID

        Returns:
            Dictionary with disease data and all associated targets
        """
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
                "pageSize": self.config.page_size
            }

            response = requests.post(
                self.base_url,
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"}
            )

            if response.status_code != 200:
                if self.config.verbose:
                    print(f"    Error: HTTP {response.status_code}")
                return None

            result = response.json()
            if 'errors' in result:
                if self.config.verbose:
                    print(f"    Error: {result['errors']}")
                return None

            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                break

            current_rows = disease_data['associatedTargets']['rows']
            if not current_rows:
                break

            all_results.extend(current_rows)

            # Check if we've retrieved all results
            if len(current_rows) < self.config.page_size:
                break

            page_index += 1
            time.sleep(0.1)  # Small delay between pages

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
        return None

    def download_disease_genes(self,
                              evidence_types: Set[str] = GENETIC_EVIDENCE_TYPES) -> Tuple[pd.DataFrame, Dict[str, Set[str]]]:
        """
        Download all genes with genetic evidence for specified diseases

        Args:
            evidence_types: Set of evidence types to include

        Returns:
            Tuple of (detailed_df, disease_gene_sets)
            - detailed_df: DataFrame with all associations and scores
            - disease_gene_sets: Dictionary mapping disease name to set of gene symbols
        """
        detailed_results = []
        disease_gene_sets = {}

        if self.config.verbose:
            print(f"Querying OpenTargets for {len(self.diseases)} diseases...")
            print(f"Evidence types: {evidence_types}")
            print(f"Min genetic evidence score: {self.config.min_genetic_evidence_score}")
            print()

        for i, disease_id in enumerate(self.diseases, 1):
            if self.config.verbose:
                print(f"[{i}/{len(self.diseases)}] Querying {disease_id}...")

            result = self.query_genetic_evidence(disease_id)

            if not result or 'data' not in result:
                if self.config.verbose:
                    print(f"    No data returned")
                continue

            disease_data = result['data']['disease']
            if not disease_data or not disease_data['associatedTargets']:
                if self.config.verbose:
                    print(f"    No associated targets")
                continue

            disease_name = disease_data['name']
            disease_genes = set()
            associations_count = 0

            for row in disease_data['associatedTargets']['rows']:
                target = row['target']

                # Check for genetic evidence
                best_genetic_score = 0
                evidence_types_found = []

                for datatype_score in row['datatypeScores']:
                    evidence_type = datatype_score['id']
                    evidence_score = datatype_score['score']

                    if evidence_type in evidence_types and evidence_score >= self.config.min_genetic_evidence_score:
                        evidence_types_found.append(evidence_type)
                        best_genetic_score = max(best_genetic_score, evidence_score)

                # Only include genes with genetic evidence above threshold
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
                    associations_count += 1

            disease_gene_sets[disease_name] = disease_genes

            if self.config.verbose:
                print(f"    Found {len(disease_genes)} genes, {associations_count} associations")

            # Delay between disease queries to avoid rate limiting
            if i < len(self.diseases):
                time.sleep(self.config.request_delay)

        # Create DataFrame and sort by disease and genetic evidence score
        detailed_df = pd.DataFrame(detailed_results).sort_values(
            ['disease_name', 'genetic_evidence_score'],
            ascending=[True, False]
        )

        if self.config.verbose:
            print()
            print("="*80)
            print(f"DOWNLOAD COMPLETE")
            print(f"Total unique genes: {len(detailed_df['gene_symbol'].unique())}")
            print(f"Total associations: {len(detailed_df)}")
            print(f"Diseases with results: {len(disease_gene_sets)}")
            print("="*80)

        return detailed_df, disease_gene_sets


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main function to download disease gene associations"""

    parser = argparse.ArgumentParser(
        description='Download disease-gene associations from OpenTargets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download autoimmune disease genes
  python download_disease_genes.py

  # Download with custom disease list
  python download_disease_genes.py --diseases EFO_0000685 MONDO_0007915

  # Download with custom output file and min score
  python download_disease_genes.py --output my_genes.csv --min-score 0.2
        """
    )

    parser.add_argument(
        '--diseases',
        nargs='+',
        default=None,
        help='Disease EFO/MONDO IDs to query (default: autoimmune diseases)'
    )

    parser.add_argument(
        '--output', '-o',
        default='disease_gene_associations_detailed.csv',
        help='Output CSV file path (default: disease_gene_associations_detailed.csv)'
    )

    parser.add_argument(
        '--min-score',
        type=float,
        default=0.1,
        help='Minimum genetic evidence score threshold (default: 0.1)'
    )

    parser.add_argument(
        '--page-size',
        type=int,
        default=1000,
        help='API page size for pagination (default: 1000)'
    )

    parser.add_argument(
        '--delay',
        type=float,
        default=0.5,
        help='Delay between disease queries in seconds (default: 0.5)'
    )

    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress verbose output'
    )

    args = parser.parse_args()

    # Use default AUTOIMMUNE_DISEASES if no diseases specified
    diseases = args.diseases if args.diseases is not None else AUTOIMMUNE_DISEASES

    # Create configuration
    config = DownloadConfig(
        min_genetic_evidence_score=args.min_score,
        output_file=args.output,
        page_size=args.page_size,
        request_delay=args.delay,
        verbose=not args.quiet
    )

    # Initialize querier
    querier = OpenTargetsQuerier(diseases=diseases, config=config)

    # Download disease genes
    detailed_df, disease_gene_sets = querier.download_disease_genes()

    # Save to file
    detailed_df.to_csv(config.output_file, index=False)

    if config.verbose:
        print(f"\nSaved results to: {config.output_file}")
        print(f"\nSummary by disease:")
        print("-" * 80)
        for disease_name in sorted(disease_gene_sets.keys()):
            n_genes = len(disease_gene_sets[disease_name])
            print(f"  {disease_name}: {n_genes} genes")


if __name__ == "__main__":
    main()
