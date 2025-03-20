#check how many of gene pairs share the same orthologous group 

import os
from BCBio import GFF
from collections import defaultdict, Counter
import csv


def create_protein_mapping(gff_file):
    """Extract protein_id to locus_tag mappings from GFF"""
    id_map = {}
    with open(gff_file) as fh:
        for rec in GFF.parse(fh):
            for feat in rec.features:
                if feat.type == "gene":
                    locus_tag = feat.qualifiers['locus_tag']
                    locus_tag = ''.join(locus_tag)  # convert list into string
                    for subfeat in feat.sub_features:
                        if "protein_id" in subfeat.qualifiers:
                            protein_id = subfeat.qualifiers["protein_id"][0]
                            id_map[protein_id] = locus_tag
    return id_map


def parse_genes(gff_file):
    """Parse gene coordinates from GFF file"""
    genes = []
    with open(gff_file) as fh:
        for rec in GFF.parse(fh):
            for feat in rec.features:
                if feat.type == "gene":
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    locus_tag = feat.qualifiers['locus_tag']
                    locus_tag = ''.join(locus_tag)  # convert list to string
                    genes.append((start, end, locus_tag))  # (start, end, locus_tag)
    return genes


def find_overlaps(genes, overlap_length):
    """Find locus tag pairs with a specific length of overlap"""
    sorted_genes = sorted(genes, key=lambda x: x[0])
    overlaps = []

    for i in range(len(sorted_genes)):
        start1, end1, tag1 = sorted_genes[i]
        for j in range(i + 1, len(sorted_genes)):
            start2, end2, tag2 = sorted_genes[j]
            if start2 >= end1:
                break

            overlap_start = max(start1, start2)
            overlap_end = min(end1, end2)
            if overlap_end - overlap_start == overlap_length:
                overlaps.append((tag1, tag2))
    return overlaps


def process_strains(gff_dir, overlap_length):
    """Process all strains' GFF files with names"""
    strain_data = []
    for strain_file in sorted(os.listdir(gff_dir)):  # Sort for consistent order
        if not strain_file.endswith(".gff"):
            continue

        gff_path = os.path.join(gff_dir, strain_file)
        protein_map = create_protein_mapping(gff_path)
        genes = parse_genes(gff_path)
        overlaps = find_overlaps(genes, overlap_length)

        strain_data.append({
            'name': os.path.splitext(strain_file)[0],  # Store filename as strain name
            'protein_map': protein_map,
            'overlaps': overlaps
        })
    return strain_data


def process_orthogroups(orthofile, strain_data):
    """Convert orthogroups with proper strain name matching"""
    ortho_dict = defaultdict(lambda: defaultdict(set))
    
    with open(orthofile) as fh:
        header = next(fh)
        tsv_strain_names = header.strip().split('\t')[1:]
        
        # Create mapping from TSV names to strain_data indices
        strain_index_map = {}
        for idx, sd in enumerate(strain_data):
            if sd['name'] in tsv_strain_names:
                strain_index_map[tsv_strain_names.index(sd['name'])] = idx
        
        for line in fh:
            parts = line.strip().split('\t')
            group_id = parts[0]
            
            for tsv_idx, strain_proteins in enumerate(parts[1:]):
                # Map TSV column index to strain_data index
                strain_idx = strain_index_map.get(tsv_idx, None)
                if strain_idx is None:
                    continue
                
                strain_map = strain_data[strain_idx]['protein_map']
                for prot_id in strain_proteins.split(','):
                    prot_id = prot_id.strip()
                    if prot_id and prot_id in strain_map:
                        ortho_dict[group_id][strain_idx].add(strain_map[prot_id])
    
    #print(ortho_dict)
    return (ortho_dict)

def analyze_conservation(strain_data, ortho_dict):
    """Analyze overlaps with direct orthogroup checking"""
    # Track overlaps per strain
    overlap_tracker = defaultdict(set)
    for strain_idx, data in enumerate(strain_data):
        for gene1, gene2 in data['overlaps']:
            sorted_pair = tuple(sorted((gene1, gene2)))
            overlap_tracker[sorted_pair].add(strain_idx)

    # Analyze conservation
    report = []
    for pair, strain_indices in overlap_tracker.items():
        gene1, gene2 = pair
        conserved_ogs = set()
        
        # Check each orthogroup in relevant strains
        for og_id, og_data in ortho_dict.items():
            for strain_idx in strain_indices:
                if strain_idx in og_data:
                    og_genes = og_data[strain_idx]
                    if gene1 in og_genes and gene2 in og_genes:
                        conserved_ogs.add(og_id)
                        break  # Only need one strain match per OG
        
        report.append({
            'gene_pair': pair,
            'strains': sorted(strain_indices),
            'num_strains': len(strain_indices),
            'orthogroups': sorted(conserved_ogs)
        })
    
    return sorted(report, key=lambda x: (-x['num_strains'], x['gene_pair']))

def verify_mappings(ortho_dict, strain_data):
    for og_id, og_data in ortho_dict.items():
        for strain_idx, loci in og_data.items():
            if not loci:
                print(f"Warning: Empty loci for OG {og_id} in strain {strain_idx}")


def main():
    overlap_length = 4
    strain_data = process_strains("gff_files", overlap_length)
    ortho_dict = process_orthogroups("Orthogroups_old.tsc", strain_data)
    verify_mappings(ortho_dict, strain_data)
    report = analyze_conservation(strain_data, ortho_dict)

    with open("overlap_report.csv", "w", newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(["Gene1", "Gene2", "Strains", "Orthogroups"])
        for entry in report:
            genes = sorted(entry['gene_pair'])
            writer.writerow([
                genes[0], genes[1],
                ";".join([strain_data[i]['name'] for i in entry['strains']]),
                ";".join(entry['orthogroups'])
            ])

    print(f"Generated report with {len(report)} overlap pairs in overlap_report.csv")


if __name__ == "__main__":
    main()


