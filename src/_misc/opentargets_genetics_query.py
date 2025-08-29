#!/usr/bin/env python3

import requests
import pandas as pd
import sys

def get_genetic_associations(disease_id):
    """Get genes with genetic evidence for a disease"""
    
    query = """
    query diseaseGenes($diseaseId: String!, $pageIndex: Int!, $pageSize: Int!) {
      disease(efoId: $diseaseId) {
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
            datasourceScores {
              id
              score
            }
          }
        }
      }
    }
    """
    
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    
    # Get all associations with pagination
    all_associations = []
    disease_name = None
    page_index = 0
    page_size = 3000
    
    while True:
        response = requests.post(url, json={
            "query": query,
            "variables": {"diseaseId": disease_id, "pageIndex": page_index, "pageSize": page_size}
        })
        
        data = response.json()
        
        # Debug: print the response
        if response.status_code != 200:
            print(f"HTTP Error {response.status_code}: {response.text}")
            return None, pd.DataFrame()
        
        if 'errors' in data:
            print(f"GraphQL Errors: {data['errors']}")
            return None, pd.DataFrame()
        
        if 'data' not in data:
            print(f"No 'data' field in response: {data}")
            return None, pd.DataFrame()
        
        disease = data['data']['disease']
        if disease_name is None:
            disease_name = disease['name']
        
        associations = disease['associatedTargets']['rows']
        
        if not associations:
            break
            
        all_associations.extend(associations)
        
        # Check if we have all results
        if len(associations) < page_size:
            break
            
        page_index += 1
    
    # Filter for genetic evidence only
    genetic_genes = []
    for assoc in all_associations:
        print(assoc)
        for dtype in assoc['datasourceScores']:
            print(dtype)
            if dtype['id'] in ['gwas_credible_sets', 'gene_burden'] and dtype['score'] > 0:
                genetic_genes.append({
                    'ensembl_id': assoc['target']['id'],
                    'gene_symbol': assoc['target']['approvedSymbol'],
                    'gene_name': assoc['target']['approvedName'],
                    'overall_score': assoc['score'],
                    'genetic_score': dtype['score']
                })
                break
    
    df = pd.DataFrame(genetic_genes)
    df = df.sort_values('genetic_score', ascending=False)
    
    print(f"Disease: {disease_name}")
    print(f"Total associations: {len(all_associations)}")
    print(f"Genetic associations: {len(genetic_genes)}")
    
    return disease_name, df

def main():
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Query OpenTargets Genetics API for disease associations')
    parser.add_argument('disease_id', nargs='?', default="MONDO_0004784",
                      help='Disease ID (default: MONDO_0004784)')
    parser.add_argument('-o', '--output-dir', default='/mnt/oak/users/emma/bin/GWT_perturbseq_analysis/metadata/',
                      help='Output directory for results (default: current directory)')
    
    args = parser.parse_args()
    
    result = get_genetic_associations(args.disease_id)
    if result[0] is None:  # Error occurred
        return
        
    disease_name, df = result
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save to CSV
    filename = os.path.join(args.output_dir, f"{args.disease_id}_genetic_associations.csv")
    df.to_csv(filename, index=False)
    
    print(f"Saved {len(df)} genes to {filename}")
    print("\nTop 10 genes:")
    print(df.head(10)[['gene_symbol', 'genetic_score']].to_string(index=False))

if __name__ == "__main__":
    main()