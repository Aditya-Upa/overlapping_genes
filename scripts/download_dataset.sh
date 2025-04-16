#download all the information for an NCBI accession id using ncbi datasets

for dir in "GCA_000025685.1" "GCA_010692905.1" "GCA_012726105.1" "GCA_029225785.1" "GCA_029232225.1"; do
echo $dir
datasets download genome accession $dir  --include gff3,protein,genome,seq-report

7z x ncbi_dataset.zip #unzip step

mv ncbi_dataset/data/$dir .
rm -r ncbi_dataset*

#transfer all the genomic fasta files into a directory called genome_files present in the current working directory
cp $dir/*fna genome_files/$dir.fna


done
