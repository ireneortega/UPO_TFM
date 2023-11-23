#!/bin/bash

echo "This script will count the frequency of the consensus upstream PAM sequence in the plasmids (for spacers CRISPR subtype I-C)”

# Get subjects from BLAST output
cut -f2 blast_IMGVR_PA_vs_I-C_system.tab > tmp.txt

# Remove duplicates subject from tmp.txt
sort tmp.txt | uniq > tmp_uniq.txt

# From the IMGVR database, remove sequences not aligning to any of the spacers (recognized sequences will be kept)
grep -A 1 -f tmp_uniq.txt IMGVR_nucl_PA_1L.fasta >> IMGVR_nucl_PA_1L_I-C_blast.fasta

# From the IMGVR database, remove sequences aligning the spacers (non-recognized sequences will be kepts)
cp /home/iortsan/tfm/part2/IMGVR_validation/IMGVR_nucl_PA_1L.fasta IMGVR_nucl_PA_1L_I-C_no_blast.fasta

while IFS= read -r line; do 
    sed -i "/$line/,+1d" IMGVR_nucl_PA_1L_I-C_no_blast.fasta
done < tmp_uniq.txt


# Create function to calculate the frequency of PAM sequence
function PAM_frequency {

  # Write the header to the output file
  echo -e "IMGVR_sequence\tIMGVR_sequence_length\tcounts_PAM_forw\tcounts_PAM_rev_comp\ttotal_counts_PAM\ttotal_K-nucleotides_PAM\tfrequency_PAM_seq_%\tGC_counts\tGC_content_%" > "$output_filename"

  # Read the multifasta database file
  while read line; do

    # If the line starts with a ">" character, it's a header line
    if [[ $line =~ ^\> ]]; then
      # Store the header line
      header=$(echo "$line" | sed 's/>//g') # Remove ">" that causes echo split information into two lines after the header

    else
      # Calculate the length of the sequence
      sequence_length=${#line}

      # Calculate the number of occurrences of the consensus upstream PAM sequence (forward) in the line
      count_forw=$(grep -o $PAM_sequence <<< "$line" | wc -l)

      # Calculate the number of occurences of the consensus upstream PAM sequence (reverse complement) in the line
      count_rev_comp=$(grep -o $PAM_sequence_rev_comp <<< "$line" | wc -l)

      # Calculate the total number of occurrences of the consensus upstream PAM sequence (forward and reverse complement) in the line
      total_PAM_counts=$(awk "BEGIN { print $count_forw + $count_rev_comp }")

      # Calculate the number of <di | tri>nucleotides the sequence can have (n-1 for dinucleotides; n-2 for dinucleotides; etc)
      total_K_nucleotides=$((sequence_length - PAM_length + 1)) # $(( )) is required to perform arithmetic operations

      # Calculate the frequency of the consensus upstream PAM sequence in the line
      frequency=$(awk "BEGIN { print $total_PAM_counts/$total_K_nucleotides * 100 }")

      # Calculate the GC content
      gc_count=$(grep -o -i "[GC]" <<< $line | wc -l)
      gc_content=$(awk "BEGIN { print $gc_count/$sequence_length * 100 }")

      # Print the header and frequency in a single line separated by a tab
      echo -e "$header\t$sequence_length\t$count_forw\t$count_rev_comp\t$total_PAM_counts\t$total_K_nucleotides\t$frequency %\t${gc_count}\t${gc_content}" >> "$output_filename"

    fi
  done < "$file"

}

# Calculate the frequency of the corresponding consensus upstream PAM sequence in each genome
PAM_sequence="TTC" # Define the consensus PAM sequence you want to count
PAM_sequence_rev_comp=$(echo "$PAM_sequence" | tr 'ATCG' 'TAGC' | rev)
PAM_length=3 # Indicate the length of the consensus PAM sequence

for file in IMGVR_*_blast.fasta; do
  output_filename="${file/IMGVR_nucl_PA_1L/freq_upPAM_IMGVR}"
  output_filename="${output_filename/.fasta/.tab}"

  #Call the function to calculate the frequency of PAM sequence
  PAM_frequency 

done

echo "Script finished!
