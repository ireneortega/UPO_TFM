#!/bin/bash

echo "This script will BLAST the spacers of the CRISPR-Cas system type I-C against the IMGV database"

# Remove sequences that contain N nucleotide from the fasta file with the spacers representing each type of CRISPR system (example: "Pseudomonas aeruginosa I-C System CRISPR-Cas spacers.fasta"). This step is necessary as BLAST cannot handle sequences with Ns
echo "The spacers will N nucleotides will be removed to avoid problems with BLAST"

input_file="Pseudomonas aeruginosa I-C System CRISPR-Cas spacers.fasta"
input_file_base=$(basename "$input_file")
spacers_no_Ns="${input_file_base%%.fasta} no Ns.fasta"

# Loop through the fasta file
while read -r header; do
    read -r sequence # Read the sequence

    # Check if the sequence contains "N"
    if [[ ! "$sequence" =~ N ]]; then
    # if not, write the header and sequence to the output file spacers_no_Ns
    echo "$header" >> "$spacers_no_Ns"
    echo "$sequence" >> "$spacers_no_Ns"
    fi

done < "$input_file"


# Run BLAST for the CRISPR system and filter results based on minimum 95 % identity and 100 % coverage
blastn -task blastn-short -query "$spacers_no_Ns" -db /home/iortsan/tfm/part2/IMGVR/IMGVR_nucl_PA_1L.fasta -outfmt '6 qseqid sseqid pident qcovs qcovhsp length qlen slen evalue qstart qend sstart send' -max_target_seqs 10000 -num_threads 40 | awk '$3 >= 95 && $4 >= 100' > blast_IMGVR_PA_vs_I-C_system.tab

echo "BLAST spacers vs IMGVR was done"
