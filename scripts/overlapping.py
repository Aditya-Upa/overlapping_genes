#this function when run within a directory containing a gff file, finds the gff file inside the directory and parses it.
#it counts the total number of overlapping genes within the gff file as well as the number of genes which are overlapping by exactly 4-nt and returns both values
#USAGE: python overlapping.py
#Make sure to use it inside a directory which contains a gff file. 

import os
from BCBio import GFF
from collections import defaultdict

def find_gff_file(directory):
    """Finds the GFF file inside the given directory."""
    for file in os.listdir(directory):
        if file.endswith(".gff"):
            return os.path.join(directory, file)
    return None

def parse_gff(gff_file):
    """Group genes by chromosome with correct coordinates"""
    chrom_genes = defaultdict(list)
    with open(gff_file) as infile:
        for rec in GFF.parse(infile):
            for feature in rec.features:
                if feature.type == "gene":
                    chrom = rec.id
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    chrom_genes[chrom].append((start, end))
    return chrom_genes

def count_overlaps(chrom_genes):
    """Count all overlaps including exact 4-nt overlaps"""
    total_overlaps = 0
    exact_4nt = 0

    for chrom, genes in chrom_genes.items():
        # Sort genes by start position
        sorted_genes = sorted(genes, key=lambda x: x[0])

        # Compare each gene with subsequent genes
        for i in range(len(sorted_genes)):
            start1, end1 = sorted_genes[i]

            for j in range(i+1, len(sorted_genes)):
                start2, end2 = sorted_genes[j]

                # Stop checking if subsequent genes can't overlap
                if start2 >= end1:
                    break

                # Calculate actual overlap
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                overlap = overlap_end - overlap_start

                if overlap > 0:
                    total_overlaps += 1
                    if overlap == 4:
                        exact_4nt += 1
    return total_overlaps, exact_4nt

# Usage

directory = os.getcwd()
gff_file = find_gff_file(directory)

#gff_file = "GCA_000092025.1_ASM9202v1_genomic.gff"

chrom_genes = parse_gff(gff_file)
overlaps, exact_4 = count_overlaps(chrom_genes)


if gff_file:
    genes = parse_gff(gff_file)
    print(f"Total overlapping genes: {overlaps}")
    print(f"Genes overlapping exactly by 4 nucleotides: {exact_4}")
else:
    print("No GFF file found in the directory.")
