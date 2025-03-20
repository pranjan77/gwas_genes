import json
import pandas as pd
import numpy as np
import argparse
import logging


def analyze_snps_and_genes(gene_file, snp_file, distance_threshold=5000, pvalue_threshold=1.0, output_prefix="Populus_trichocarpa", 
                          save_output=True, verbose=True, save_gene_function=False):
    """
    Analyze SNPs and genes to identify associations based on distance and p-value thresholds.
    
    Parameters:
    -----------
    gene_file : str
        Path to the gene JSON file
    snp_file : str
        Path to the SNP JSON file
    distance_threshold : int, optional
        Distance threshold in base pairs (default: 5000)
    pvalue_threshold : float, optional
        P-value threshold (default: 1.0, i.e., no filtering)
    output_prefix : str, optional
        Prefix for output file names (default: "Populus_trichocarpa")
    save_output : bool, optional
        Whether to save output files (default: True)
    verbose : bool, optional
        Whether to print progress information (default: True)
    save_gene_function : bool, optional
        Whether to save the complete gene function file for the entire genome (default: False)
    
    Returns:
    --------
    dict
        Dictionary containing DataFrames and summary information:
        - 'gene_df': DataFrame of gene information
        - 'snp_gene_df': DataFrame of SNP-gene associations
        - 'gene_snp_df': DataFrame of gene-centric view
        - 'summary': Dictionary with summary statistics
        - 'output_files': Dictionary with output file paths
    """
    if verbose:
        print(f"Using distance threshold: {distance_threshold} bp")
        print(f"Using p-value threshold: {pvalue_threshold}")

    # Define helper function to calculate distance
    def calculate_distance(snp_pos, gene):
        """
        Calculate the distance between an SNP and a gene.
        For SNPs outside the gene: minimum distance to either start or end
        For SNPs inside the gene: relative position within the gene
        """
        # If SNP is within the gene
        if gene['is_within_gene']:
            # Return distance as internal position
            return 0
        else:
            # If outside, return minimum distance to either start or end
            return min(
                abs(snp_pos - gene['start']),
                abs(snp_pos - gene['end'])
            )

    # Load the gene data
    with open(gene_file, 'r') as file:
        data = json.load(file)

    # Create a list to store each gene's data
    records = []

    for j in data['features']:
        gene_id = j['id']
        
        location = j['location'][0]
        chr_num = location[0]
        start = location[1]
        orientation = location[2]
        length = location[3]

        if orientation == "-":
            end = start 
            start = start - length
        else:
            start = location[1] - 1
            end = start + length
        
        # Process functions
        try:
            function = " ".join(j['functions'])
        except Exception as e:
            function = None
        
        # Process GO terms
        try:
            go_terms = ", ".join(j['ontology_terms']['GO'].keys())
        except Exception as e:
            go_terms = None

        records.append({
            'gene_id': gene_id,
            'function': function,
            'go_terms': go_terms,
            'location': location,
            "chr": chr_num,
            "orientation": orientation,
            "start": start,
            "end": end
        })

    # Create a DataFrame from the records list
    gene_df = pd.DataFrame(records)

    # Load SNP data from JSON file
    with open(snp_file, 'r') as file:
        snp_data = json.load(file)
    logging.info(f"SNP Data: {snp_data}")

    # Create a DataFrame for SNPs based on the JSON structure
    snp_records = []
    for association in snp_data.get('association_details', []):
        for result in association.get('association_results', []):
            if len(result) >= 5:  # Ensure we have enough elements in each result
                snp_records.append({
                    'chr': result[0],                  # Chromosome
                    'snp_id': result[1],               # SNP ID
                    'pos': result[2],                  # Position
                    'pvalue': result[3],               # P-value
                    'additional_value': result[4]      # Additional value (may represent something specific)
                })

    # Create DataFrame from SNP records
    gwas_df = pd.DataFrame(snp_records)

    logging.info(f"GWAS DataFrame: {gwas_df}")

    # Filter SNPs by p-value if threshold is specified
    orig_count = len(gwas_df)
    if pvalue_threshold < 1.0:
        gwas_df = gwas_df[gwas_df['pvalue'] <= pvalue_threshold]
        filtered_count = len(gwas_df)
        if verbose:
            print(f"Filtered SNPs by p-value threshold {pvalue_threshold}: {orig_count} → {filtered_count} SNPs")
    else:
        filtered_count = orig_count
        if verbose:
            print(f"No p-value filtering applied. Using all {len(gwas_df)} SNPs.")

    # Create a list to store SNP data with nearest gene info
    snp_gene_records = []

    # Process each SNP
    for _, snp in gwas_df.iterrows():
        snp_chr = snp['chr']
        snp_pos = snp['pos']
        pvalue = snp['pvalue'] if 'pvalue' in snp and not pd.isna(snp['pvalue']) else None
        
        # Filter genes on the same chromosome
        genes_on_chr = gene_df[gene_df['chr'] == snp_chr].copy()
        
        if len(genes_on_chr) > 0:
            # Add is_within_gene column first to simplify calculations
            genes_on_chr['is_within_gene'] = genes_on_chr.apply(
                lambda gene: (snp_pos >= gene['start']) and (snp_pos <= gene['end']),
                axis=1
            )
            
            # Calculate distance to each gene
            genes_on_chr['distance'] = genes_on_chr.apply(
                lambda gene: calculate_distance(snp_pos, gene),
                axis=1
            )
            
            # Filter genes within the distance threshold
            nearby_genes = genes_on_chr[genes_on_chr['distance'] <= distance_threshold]
            
            if len(nearby_genes) > 0:
                # Sort by distance to get closest genes first
                nearby_genes = nearby_genes.sort_values('distance')
                
                # Process each nearby gene
                for _, gene in nearby_genes.iterrows():
                    # Get is_within_gene from the calculated column
                    is_within_gene = gene['is_within_gene']
                    
                    # Determine SNP position relative to gene (3', 5', or within)
                    gene_orientation = gene['orientation']
                    
                    if is_within_gene:
                        # Simplify to just "within gene" without percentage
                        snp_position_category = "within gene"
                    else:
                        # For "+" orientation: 5' is upstream (before start), 3' is downstream (after end)
                        # For "-" orientation: 5' is downstream (after end), 3' is upstream (before start)
                        if gene_orientation == "+":
                            if snp_pos < gene['start']:
                                snp_position_category = "5'"  # Upstream of gene
                            else:  # snp_pos > gene['end']
                                snp_position_category = "3'"  # Downstream of gene
                        else:  # gene_orientation == "-"
                            if snp_pos > gene['end']:
                                snp_position_category = "5'"  # Downstream of gene (for - orientation)
                            else:  # snp_pos < gene['start']
                                snp_position_category = "3'"  # Upstream of gene (for - orientation)
                    
                    snp_gene_records.append({
                        'snp_chr': snp_chr,
                        'snp_id': snp['snp_id'],
                        'snp_pos': snp_pos,
                        'pvalue': pvalue,
                        'gene_id': gene['gene_id'],
                        'gene_start': gene['start'],
                        'gene_end': gene['end'],
                        'gene_orientation': gene_orientation,
                        'distance': gene['distance'],
                        'is_within_gene': is_within_gene,
                        'snp_position_category': snp_position_category,
                        'gene_function': gene['function'],
                        'gene_go_terms': gene['go_terms']
                    })
            else:
                # No genes within threshold distance
                snp_gene_records.append({
                    'snp_chr': snp_chr,
                    'snp_id': snp['snp_id'],
                    'snp_pos': snp_pos,
                    'pvalue': pvalue,
                    'gene_id': None,
                    'gene_start': None,
                    'gene_end': None,
                    'gene_orientation': None,
                    'distance': None,
                    'is_within_gene': False,
                    'snp_position_category': None,
                    'gene_function': None,
                    'gene_go_terms': None
                })
        else:
            # No genes on this chromosome
            snp_gene_records.append({
                'snp_chr': snp_chr,
                'snp_id': snp['snp_id'],
                'snp_pos': snp_pos,
                'pvalue': pvalue,
                'gene_id': None,
                'gene_start': None,
                'gene_end': None,
                'gene_orientation': None,
                'distance': None,
                'is_within_gene': False,
                'snp_position_category': None,
                'gene_function': None,
                'gene_go_terms': None
            })

    # Create a DataFrame for SNPs with nearest gene info
    snp_gene_df = pd.DataFrame(snp_gene_records)

    # Create gene-centric view
    # Step 1: Filter out rows with no gene_id
    valid_snp_gene_df = snp_gene_df[snp_gene_df['gene_id'].notna()]

    # Step 2: Create a dictionary to store gene-centric information
    gene_snp_dict = {}

    # Process each valid SNP-gene association
    for _, row in valid_snp_gene_df.iterrows():
        gene_id = row['gene_id']
        snp_info = f"{row['snp_id']}"
        
        if row['pvalue'] is not None:
            snp_info += f" (p={row['pvalue']})"
        
        # Add position info differently based on whether SNP is within gene or not
        if row['is_within_gene']:
            snp_info += f" [within gene]"
        else:
            # Format distance without 'bp' and with the position category
            snp_info += f" [{int(row['distance'])}, {row['snp_position_category']}]"
        
        # Initialize the gene entry if it doesn't exist
        if gene_id not in gene_snp_dict:
            gene_snp_dict[gene_id] = {
                'gene_id': gene_id,
                'chr': row['snp_chr'],
                'gene_start': row['gene_start'],
                'gene_end': row['gene_end'],
                'gene_orientation': row['gene_orientation'],
                'gene_function': row['gene_function'],
                'gene_go_terms': row['gene_go_terms'],
                'associated_snps': [],
                'snp_count': 0,
                'min_pvalue': None
            }
            
            # Initialize min_pvalue only if we have a valid p-value
            if row['pvalue'] is not None:
                gene_snp_dict[gene_id]['min_pvalue'] = row['pvalue']
        
        # Add this SNP to the gene's list
        gene_snp_dict[gene_id]['associated_snps'].append(snp_info)
        gene_snp_dict[gene_id]['snp_count'] += 1
        
        # Update minimum p-value if applicable and if we have a valid p-value
        if row['pvalue'] is not None:
            if gene_snp_dict[gene_id]['min_pvalue'] is None:
                gene_snp_dict[gene_id]['min_pvalue'] = row['pvalue']
            else:
                gene_snp_dict[gene_id]['min_pvalue'] = min(gene_snp_dict[gene_id]['min_pvalue'], row['pvalue'])

    # Convert the dictionary to a list of records
    gene_snp_records = []
    for gene_id, gene_data in gene_snp_dict.items():
        # Join all SNP information into a single string
        gene_data['associated_snps'] = ", ".join(gene_data['associated_snps'])
        
        # Add the gene to the records
        gene_snp_records.append(gene_data)

    # Create DataFrame from records
    gene_snp_df = pd.DataFrame(gene_snp_records)

    # Add thresholds to output filenames
    threshold_info = f"_d{distance_threshold}"
    if pvalue_threshold < 1.0:
        threshold_info += f"_p{pvalue_threshold:.0e}"

    # Generate output filenames with thresholds
    gene_function_file = f"{output_prefix}_gene_function.csv"
    snp_analysis_file = f"{output_prefix}_SNP_GWAS_analysis{threshold_info}.csv"
    gene_centric_file = f"{output_prefix}_gene_centric_GWAS{threshold_info}.csv"

    # Save files if requested
    saved_files = []
    if save_output:
        # Save SNP-centric view (SNP-gene associations)
        snp_gene_df.to_csv(snp_analysis_file, index=False)
        saved_files.append(snp_analysis_file)
        
        # Save gene-centric view
        gene_snp_df.to_csv(gene_centric_file, index=False)
        saved_files.append(gene_centric_file)
        
        # Optionally save the complete gene function file
        if save_gene_function:
            gene_df.to_csv(gene_function_file, index=False)
            saved_files.append(gene_function_file)
        
        if verbose:
            print(f"Results saved to:")
            for file in saved_files:
                print(f"  - {file}")

    # Prepare summary information
    summary = {
        'total_snps': orig_count,
        'filtered_snps': filtered_count,
        'total_genes': len(gene_df),
        'snp_gene_associations': len(snp_gene_df),
        'gene_centric_records': len(gene_snp_df),
        'distance_threshold': distance_threshold,
        'pvalue_threshold': pvalue_threshold
    }
    
    if verbose:
        print(f"Processed {filtered_count} SNPs and {len(gene_df)} genes.")
        print(f"Found {len(valid_snp_gene_df)} SNP-gene associations within {distance_threshold} bp and p-value ≤ {pvalue_threshold}.")
        print(f"Generated {len(gene_snp_df)} gene-centric records with SNP counts and min p-values.")

    # Return results dictionary
    return {
        'gene_df': gene_df,
        'snp_gene_df': snp_gene_df,
        'gene_snp_df': gene_snp_df,
        'summary': summary,
        'output_files': {
            'gene_function_file': gene_function_file if save_gene_function else None,
            'snp_analysis_file': snp_analysis_file,
            'gene_centric_file': gene_centric_file
        }
    }


if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Analyze SNPs and genes to identify associations based on distance and p-value thresholds')
    parser.add_argument('--distance', type=int, default=5000, help='Distance threshold in base pairs (default: 5000)')
    parser.add_argument('--pvalue', type=float, default=1.0, help='P-value threshold (default: 1.0, i.e., no filtering)')
    parser.add_argument('--gene_file', type=str, default="Populus_trichocarpa.json", help='Gene JSON file path')
    parser.add_argument('--snp_file', type=str, default="snps.json", help='SNP JSON file path')
    parser.add_argument('--output_prefix', type=str, default="Populus_trichocarpa", help='Prefix for output file names')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress messages')
    parser.add_argument('--save_gene_function', action='store_true', help='Save the complete gene function file for the entire genome')

    args = parser.parse_args()

    # Call the function with parsed arguments
    analyze_snps_and_genes(
        gene_file=args.gene_file,
        snp_file=args.snp_file,
        distance_threshold=args.distance,
        pvalue_threshold=args.pvalue,
        output_prefix=args.output_prefix,
        verbose=not args.quiet,
        save_gene_function=args.save_gene_function
    )

