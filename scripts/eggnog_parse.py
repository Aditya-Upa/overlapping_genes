import os
import pandas as pd

from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
#import pre-existing functions
from overlapped import parse_gff_for_genes_and_cds, find_overlaps 

# --------------------------
# Enhanced Analysis Pipeline
# --------------------------

def parse_emapper(emapper_file):
    """Parse eggNOG annotations with error handling"""
    try:
        df = pd.read_csv(emapper_file, sep='\t', comment='#', header=None,
                        names=[
                            'query', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs',
                            'max_annot_lvl', 'COG_category', 'Description', 'Preferred_name',
                            'GOs', 'EC', 'KEGG_ko', 'KEGG_Pathway', 'KEGG_Module',
                            'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy',
                            'BiGG_Reaction', 'PFAMs'
                        ])
        
                # Convert '-' to NaN for better pandas handling
        # df.replace('-', pd.NA, inplace=True)

        df['primary_OG'] = df['eggNOG_OGs'].str.split('@').str[0].str.split(',').str[0]
        return df[['query', 'primary_OG', 'COG_category', 'KEGG_Pathway', 'GOs']].dropna()
        
    except Exception as e:
        print(f"Error parsing {emapper_file}: {str(e)}")
        return pd.DataFrame()

def load_strain_data(strain):
    """Integrated data loader for a single strain"""
    base_path = f"{strain}"
    
    # Load GFF data
    gff_path = f"gff_files/{base_path}.gff"
    genes, gene_to_protein = parse_gff_for_genes_and_cds(gff_path)
    #print(gene_to_protein)
    
    # Find 4nt overlaps
    overlaps = find_overlaps(genes, overlap_length=4)
    
    # Load eggNOG annotations
    emapper_path = f"eggnog_files/{base_path}.emapper.annotations"
    emapper_df = parse_emapper(emapper_path)

    # Create annotation mapping (gene_id -> function)
    annotation_map = {}
    for gene_id, protein_id in gene_to_protein.items():
        annotations = emapper_df[emapper_df['query'] == protein_id]
        if not annotations.empty:
            annotation_map[gene_id] = {
                'COG': annotations['COG_category'].values[0],
                'KEGG': annotations['KEGG_Pathway'].values[0],
                'GO': annotations['GOs'].values[0],
                'protein_id': protein_id
            }
   
    return {
        'overlaps': overlaps,
        'annotations': annotation_map,
        'gene_to_protein': gene_to_protein
    }

def analyze_functional_sharing(strain_data):
    """Calculate functional sharing statistics, ignoring '-' and NaN"""
    stats = defaultdict(int)
    cog_counts = defaultdict(int)
    kegg_counts = defaultdict(int)
    
    for gene1, gene2 in strain_data['overlaps']:
        ann1 = strain_data['annotations'].get(gene1)
        ann2 = strain_data['annotations'].get(gene2)
        
        if not ann1 or not ann2:
            stats['pairs_missing_annotations'] += 1
            continue
            
        stats['total_pairs'] += 1
        
        # COG comparison (ignore '-' and NaN)
        cog1 = ann1['COG']
        cog2 = ann2['COG']
        if pd.notna(cog1) and pd.notna(cog2) and cog1 != '-' and cog2 != '-':
            if cog1 == cog2:
                stats['shared_cog'] += 1
                cog_counts[cog1] += 1
        else:
            stats['missing_cog'] += 1
            
        # KEGG comparison (handle '-' as empty)
        kegg1 = str(ann1['KEGG']).strip()
        kegg2 = str(ann2['KEGG']).strip()
        kegg1 = set() if kegg1 in ['-', 'nan'] else set(kegg1.split(','))
        kegg2 = set() if kegg2 in ['-', 'nan'] else set(kegg2.split(','))
        
        common_kegg = kegg1 & kegg2
        if common_kegg:
            stats['shared_kegg'] += 1
            for path in common_kegg:
                kegg_counts[path] += 1
                
        # GO comparison (handle '-' as empty)
        go1 = str(ann1['GO']).strip()
        go2 = str(ann2['GO']).strip()
        go1 = set() if go1 in ['-', 'nan'] else set(go1.split(','))
        go2 = set() if go2 in ['-', 'nan'] else set(go2.split(','))
        
        if go1 & go2:
            stats['shared_go'] += 1
            
    return {
        'stats': dict(stats),
        'cog_counts': dict(cog_counts),
        'kegg_counts': dict(kegg_counts)
    }

# def analyze_functional_sharing(strain_data):
#     """Calculate functional sharing statistics"""
#     stats = defaultdict(int)
#     cog_counts = defaultdict(int)
#     kegg_counts = defaultdict(int)
    
#     for gene1, gene2 in strain_data['overlaps']:
#         ann1 = strain_data['annotations'].get(gene1)
#         ann2 = strain_data['annotations'].get(gene2)
        
#         if not ann1 or not ann2:
#             stats['pairs_missing_annotations'] += 1
#             continue
            
#         stats['total_pairs'] += 1
        
#         # COG comparison
#         if ann1['COG'] == ann2['COG'] and pd.notna(ann1['COG']):
#             stats['shared_cog'] += 1
#             cog_counts[ann1['COG']] += 1
            
#         # KEGG comparison
#         kegg1 = set(str(ann1['KEGG']).split(','))
#         kegg2 = set(str(ann2['KEGG']).split(','))
#         common_kegg = kegg1 & kegg2
#         if common_kegg:
#             stats['shared_kegg'] += 1
#             for path in common_kegg:
#                 kegg_counts[path] += 1
                
#         # GO comparison
#         go1 = set(str(ann1['GO']).split(','))
#         go2 = set(str(ann2['GO']).split(','))
#         if go1 & go2:
#             stats['shared_go'] += 1
            
#     return {
#         'stats': dict(stats),
#         'cog_counts': dict(cog_counts),
#         'kegg_counts': dict(kegg_counts)
#     }

#WRITE THE OVERLAP FILES TO A TEXT FILE IN A FORMAT WHICH WILL BE USED BY R SCRIPT NEXT

def save_overlaps_for_r(overlaps_data, filename="overlap_genes.txt"):
    """Save overlapping protein IDs in R-compatible format.
    
    overlaps_data is expected to be a dictionary with keys as strains and values as a dictionary:
        {
            'overlaps': list of (gene1, gene2) tuples,
            'gene_to_protein': { gene_id: protein_id, ... }
        }
    """
    with open(filename, 'w') as f:
        f.write("# Overlapping protein IDs for clusterProfiler analysis\n")
        f.write("overlap_genes <- c(\n")
        
        unique_proteins = set()
        for strain, data in overlaps_data.items():
            overlaps = data['overlaps']
            gene_to_protein = data['gene_to_protein']
            for gene1, gene2 in overlaps:
                protein1 = gene_to_protein.get(gene1)
                protein2 = gene_to_protein.get(gene2)
                if protein1:
                    unique_proteins.add(f'"{strain}|{protein1}"')
                if protein2:
                    unique_proteins.add(f'"{strain}|{protein2}"')
        
        # Write in chunks for readability
        chunk_size = 10
        proteins_list = sorted(list(unique_proteins))
        for i in range(0, len(proteins_list), chunk_size):
            chunk = proteins_list[i:i+chunk_size]
            f.write("  " + ",\n  ".join(chunk) + ("," if i+chunk_size < len(proteins_list) else "") + "\n")
        
        f.write(")\n")
        
def full_analysis():
    strains = [
        "GCA_000025685.1", "GCA_010692905.1", "GCA_012726105.1",
        "GCA_029225785.1", "GCA_029232225.1"
    ]
    
    all_stats = {}
    overlaps_data = {}
    
    for strain in strains:
        print(f"\nProcessing {strain}...")
        strain_data = load_strain_data(strain)
        
        #Store overlaps with strain prefix
        #overlaps_data[strain] = strain_data['overlaps']
        overlaps_data[strain] = {
    'overlaps': strain_data['overlaps'],
    'gene_to_protein': strain_data['gene_to_protein']
}

        #print(strain_data['overlaps'])
        if not strain_data['overlaps']:
            print(f"No overlaps found for {strain}")
            continue
            
        results = analyze_functional_sharing(strain_data)
        
        # Print strain summary
        print(f"Total 4nt overlaps: {results['stats']['total_pairs']}")
        print(f"Shared COG: {results['stats'].get('shared_cog', 0)}")
        print(f"Shared KEGG: {results['stats'].get('shared_kegg', 0)}")
        print(f"Shared GO: {results['stats'].get('shared_go', 0)}")
        
        all_stats[strain] = results

    # Save overlaps for clusterProfiler analysis using R
    save_overlaps_for_r(overlaps_data)        
 
    # Generate visualizations
    plot_combined_results(all_stats)
    
def plot_combined_results(all_stats):
    """Visualize results across all strains"""
    os.makedirs("results", exist_ok=True)
    
    # Aggregate data
    combined = defaultdict(int)
    cog_dist = defaultdict(int)
    kegg_dist = defaultdict(int)
    
    for strain, data in all_stats.items():
        for k, v in data['stats'].items():
            combined[k] += v
        for cog, count in data['cog_counts'].items():
            cog_dist[cog] += count
        for kegg, count in data['kegg_counts'].items():
            kegg_dist[kegg] += count
            
    # Plotting
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))

    # Functional sharing
    metrics = ['shared_cog', 'shared_kegg', 'shared_go']
    values = [combined[m]/combined['total_pairs'] for m in metrics]
    sns.barplot(x=metrics, y=values, ax=ax[0], palette="Blues_d")
    ax[0].set_title("Functional Category Sharing")
    ax[0].set_ylabel("Proportion of Overlapping Pairs")
    
    #COG distribution
    if kegg_dist:
        cog_df=pd.DataFrame.from_dict(cog_dist, orient='index').reset_index()
        cog_df.columns = ['COG', 'Count']
        cog_df = cog_df.sort_values('Count', ascending=False).head(10)
        sns.barplot(x='Count', y='COG', data=cog_df, ax=ax[1], palette="Accent")
        ax[1].set_title("Top 10 Shared COG Categories")

    plt.tight_layout()
    plt.savefig("results/functional_analysis.png")
    plt.close()

if __name__ == "__main__":
    full_analysis()
