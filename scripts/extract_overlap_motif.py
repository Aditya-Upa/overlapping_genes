
#This script iterates over a list of strains and creates a FASTA file (one per strain) containing the 4‑nt overlap sequences along with ±10 nt flanking regions for later MEME analysis.


import os
import re
from Bio import SeqIO
from collections import defaultdict

# --- Configuration ---
# List of assembly IDs (used for file naming)
STRAINS = [
    "GCA_000025685.1",
    "GCA_010692905.1",
    "GCA_012726105.1",
    "GCA_029225785.1",
    "GCA_029232225.1"
]

GENOME_DIR = "genome_files"  # Directory containing .fna files
GFF_DIR = "gff_files"      # Directory containing .gff files
OUTPUT_DIR = "motif_analysis_4nt" # Output directory for FASTA files

TARGET_OVERLAP_LENGTH = 4  # We want exactly 4-nt overlaps
FLANK_SIZE = 10            # +/- 10 nt flanking region

# --- Functions ---

def parse_gff_genes(gff_path):
    """
    Parses a GFF3 file to extract gene features.

    Args:
        gff_path (str): Path to the GFF file.

    Returns:
        dict: A dictionary where keys are contig IDs and values are lists of
              gene tuples, sorted by start position.
              Gene tuple: (start, end, strand, gene_id)
              Returns None if the file cannot be processed.
    """
    genes_by_contig = defaultdict(list)
    try:
        with open(gff_path, 'r') as gff_file:
            for line in gff_file:
                if line.startswith('#') or line.strip() == '':
                    continue # Skip comments and empty lines

                parts = line.strip().split('\t')
                if len(parts) != 9:
                    print(f"  Warning: Skipping malformed GFF line: {line.strip()}")
                    continue

                feature_type = parts[2]
                # Focus on 'gene' features. Adjust if your GFF uses 'CDS' primarily
                # and lacks 'gene' features covering the full CDS extent.
                if feature_type.lower() == 'gene':
                    contig_id = parts[0]
                    try:
                        # GFF is 1-based, inclusive
                        start = int(parts[3])
                        end = int(parts[4])
                    except ValueError:
                        print(f"  Warning: Skipping GFF line with non-integer coordinates: {line.strip()}")
                        continue

                    strand = parts[6]
                    attributes = parts[8]

                    # Extract a gene ID - look for 'ID=' first, then 'locus_tag='
                    gene_id = None
                    id_match = re.search(r'ID=([^;]+)', attributes)
                    if id_match:
                        gene_id = id_match.group(1)
                    else:
                        locus_match = re.search(r'locus_tag=([^;]+)', attributes)
                        if locus_match:
                            gene_id = locus_match.group(1)

                    if not gene_id:
                        # Fallback: Use coordinates if no ID found
                        gene_id = f"gene_{contig_id}_{start}_{end}"
                        # print(f"  Warning: Could not find 'ID' or 'locus_tag' for gene at {contig_id}:{start}-{end}. Using generated ID.")

                    genes_by_contig[contig_id].append((start, end, strand, gene_id))

        # Sort genes within each contig by start position
        for contig_id in genes_by_contig:
            genes_by_contig[contig_id].sort(key=lambda x: x[0])

        return genes_by_contig

    except FileNotFoundError:
        print(f"  Error: GFF file not found at {gff_path}")
        return None
    except Exception as e:
        print(f"  Error parsing GFF file {gff_path}: {e}")
        return None

def find_overlaps_and_extract(genes_by_contig, genome_sequences, strain_id):
    """
    Finds adjacent genes with TARGET_OVERLAP_LENGTH, extracts sequences.

    Args:
        genes_by_contig (dict): Output from parse_gff_genes.
        genome_sequences (dict): Dictionary of SeqRecords from SeqIO.to_dict.
        strain_id (str): Identifier for the current strain.

    Returns:
        list: A list of tuples, where each tuple contains:
              (fasta_header, extracted_sequence_string)
    """
    extracted_sequences = []
    pair_counter = 0

    for contig_id, genes in genes_by_contig.items():
        if contig_id not in genome_sequences:
            print(f"  Warning: Contig '{contig_id}' found in GFF but not in FASTA for strain {strain_id}. Skipping.")
            continue

        contig_seq = genome_sequences[contig_id].seq
        contig_len = len(contig_seq)

        # Iterate through adjacent gene pairs
        for i in range(len(genes) - 1):
            gene1_start, gene1_end, _, gene1_id = genes[i]
            gene2_start, gene2_end, _, gene2_id = genes[i+1]

            # Calculate overlap coordinates (1-based)
            overlap_start = max(gene1_start, gene2_start)
            overlap_end = min(gene1_end, gene2_end)

            # Check if there is an overlap and if it's the target length
            # overlap_start <= overlap_end ensures they actually overlap
            if overlap_start <= overlap_end:
                overlap_length = overlap_end - overlap_start + 1
                if overlap_length == TARGET_OVERLAP_LENGTH:
                    pair_counter += 1

                    # Define extraction boundaries (1-based), clamped to contig limits
                    extract_start = max(1, overlap_start - FLANK_SIZE)
                    extract_end = min(contig_len, overlap_end + FLANK_SIZE)

                    # Extract sequence (Python slicing is 0-based, exclusive end)
                    # Adjust 1-based coordinates to 0-based index
                    sequence_str = str(contig_seq[extract_start - 1 : extract_end])

                    # Create FASTA header
                    header = f">{strain_id}_pair{pair_counter}|{gene1_id}-{gene2_id}|{contig_id}:{extract_start}-{extract_end}_Overlap:{overlap_start}-{overlap_end}"
                    extracted_sequences.append((header, sequence_str))

    return extracted_sequences

# --- Main Execution ---
def main():
    """
    Main function to process all strains.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Starting analysis. Output directory: {OUTPUT_DIR}")
    print(f"Looking for exactly {TARGET_OVERLAP_LENGTH}nt overlaps with +/- {FLANK_SIZE}nt flanks.")

    total_sequences_written = 0

    for strain in STRAINS:
        print(f"\nProcessing strain: {strain}")
        genome_path = os.path.join(GENOME_DIR, f"{strain}.fna")
        gff_path = os.path.join(GFF_DIR, f"{strain}.gff")

        # 1. Check input files exist
        if not os.path.exists(genome_path):
            print(f"  Error: Genome file not found: {genome_path}. Skipping strain.")
            continue
        if not os.path.exists(gff_path):
            print(f"  Error: GFF file not found: {gff_path}. Skipping strain.")
            continue

        # 2. Load Genome Sequences
        genome_sequences = None
        try:
            # Use to_dict for easy access by contig ID
            genome_sequences = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
            if not genome_sequences:
                print(f"  Error: No sequences found in genome file: {genome_path}. Skipping strain.")
                continue
            print(f"  Loaded {len(genome_sequences)} contigs from {genome_path}")
        except Exception as e:
            print(f"  Error reading genome file {genome_path}: {e}. Skipping strain.")
            continue

        # 3. Parse GFF File
        genes_by_contig = parse_gff_genes(gff_path)
        if genes_by_contig is None: # Error occurred during parsing
             print(f"  Failed to parse GFF for {strain}. Skipping strain.")
             continue
        if not genes_by_contig:
            print(f"  No 'gene' features found or parsed from {gff_path}. Skipping strain.")
            continue
        num_genes = sum(len(glist) for glist in genes_by_contig.values())
        print(f"  Parsed {num_genes} genes from {gff_path} across {len(genes_by_contig)} contigs.")


        # 4. Find overlaps and Extract Sequences
        extracted_sequences = []
        try:
            extracted_sequences = find_overlaps_and_extract(genes_by_contig, genome_sequences, strain)
        except Exception as e:
            print(f"  Error during overlap finding/extraction for {strain}: {e}")
            # Continue to next strain or handle as needed

        # 5. Write Output FASTA File
        if extracted_sequences:
            output_fasta_path = os.path.join(OUTPUT_DIR, f"{strain}_overlaps_{TARGET_OVERLAP_LENGTH}nt.fasta")
            try:
                with open(output_fasta_path, "w") as f_out:
                    for header, seq in extracted_sequences:
                        # Ensure sequence wraps if needed (optional, MEME handles long lines)
                        # Typically not needed for short sequences like these (24 nt)
                        f_out.write(f"{header}\n{seq}\n")
                print(f"  Success: Wrote {len(extracted_sequences)} sequences to {output_fasta_path}")
                total_sequences_written += len(extracted_sequences)
            except Exception as e:
                print(f"  Error writing output FASTA file {output_fasta_path}: {e}")
        else:
            print(f"  No {TARGET_OVERLAP_LENGTH}-nt overlaps found meeting criteria for {strain}.")

    print(f"\nAnalysis complete. Total sequences written across all strains: {total_sequences_written}")

if __name__ == "__main__":
    main()
