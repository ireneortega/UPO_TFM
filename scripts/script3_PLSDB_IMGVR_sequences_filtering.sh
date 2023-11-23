#!/bin/bash

echo "This script will filter the IMGVR database and extract the viruses related to Pseudomonas aeruginosa"

## PART A: IMG/VR
# Filter the rows from the IMGVR_all_Host_information-high_confidence.tsv file with the sequences metadata that contain "s__Pseudomonas aeruginosa" (s__ means species). 
echo "Processing the IMGVR_all_Host_information-high_confidence.tsv file to filter the rows that correspond to Pseudomonas aeruginosa..."

grep "s__Pseudomonas aeruginosa" IMGVR_all_Host_information-high_confidence.tsv > IMGVR_Host_PA_information.tsv

echo "IMGVR_all_Host_information-high_confidence.tsv file processed!"

# Extract the values of the first field after the tab in the filtered tsv file. These values correspond to the fasta headers of the sequences that need to be filtered extracted from the IMGVR database with the fasta sequences
uvig_values=$(awk '{print $1}' IMGVR_Host_PA_information.tsv)

# Extract the sequences from the IMGVR database (IMGVR_all_Host_information-high_confidence.fna)
echo "Extracting the sequences that correspond to Pseudomonas aeruginosa from the IMGVR database..."

INDEX="1" # to know the progress of the for loop
for uvig_value in $uvig_values
do
    echo ${INDEX}\ $uvig_value

    line=$(grep -n -m 1 "^>${uvig_value}\b" IMGVR_all_nucleotides-high_confidence.fna | cut -d: -f1) # Line number that correspond to the uvig_value in IMGVR_all_nucleotides-high_confidence.fna. The parameter -m 1 means stops searching once the first result is found

    sed "${line}q;d" IMGVR_all_nucleotides-high_confidence.fna >> IMGVR_nucl_PA.fasta # Print header

    ((line++))
    
    awk -v start="$line" 'NR>=start{if (/^>/) exit; else print}' IMGVR_all_nucleotides-high_confidence.fna >> IMGVR_nucl_PA.fasta # Print the next lines (nucleotides) up to a new header sequence starting with “>” is reached 

    let INDEX=${INDEX}+1

done

echo "Sequences corresponding to Pseudomonas aeruginosa from the IMGVR database extracted!"


# Convert IMGVR_nucl_PA.fasta (multifasta multiline format) to oneline format
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' IMGVR_nucl_PA.fasta > IMGVR_nucl_PA_1L.fasta
echo "The file IMGVR_nucl_PA.fasta was converted to oneline format"


# Short FASTA headers (max length for makeblastdb is 50)
sed -i 's/|.*//' IMGVR_nucl_PA_1L.fasta
echo "The headers of the fasta sequences were shortened for BLAST"

# Make BLAST database with the IMGVR database for Pseudomonas aeruginosa
echo "Creating BLAST database with the IMGVR database for Pseudomonas aeruginosa"
makeblastdb -in IMGVR_nucl_PA_1L.fasta -parse_seqids -dbtype nucl



## PART B: PLSDB
# Extract the headers of the sequences that correspond to Pseudomonas aeruginosa in the PLSDB database
echo "Extracting the sequences that correspond to Pseudomonas aeruginosa from the PLSDB..."

grep "Pseudomonas aeruginosa" plsdb.fna > plsdb_Host_PA_information.tsv

PA_headers=$(awk '{print $1}' plsdb_Host_PA_information.tsv)

INDEX="1" # to know the progress of the for loop
for PA_header in $PA_headers
do
    echo ${INDEX}\ $PA_header

    line=$(grep -n -m 1 "^${PA_header}\b" plsdb.fna | cut -d: -f1) # Line number that correspond to the uvig_value in PLSDB database. The parameter -m 1 means stops searching once the first result is found

    sed "${line}q;d" plsdb.fna >> plsdb_Host_PA.fasta # Print header

    ((line++))
    
    awk -v start="$line" 'NR>=start{if (/^>/) exit; else print}' plsdb.fna >> plsdb_Host_PA.fasta # # Print the next lines (nucleotides) up to a new header sequence starting with “>” is reached

    let INDEX=${INDEX}+1

done

echo "Sequences corresponding to Pseudomonas aeruginosa from the PLSDB extracted!"


#Convert plsdb_Host_PA.fasta (multifasta multiline format) to oneline format
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' plsdb_Host_PA.fasta > plsdb_Host_PA_1L.fasta
echo "The file plsdb_Host_PA.fasta was converted to oneline format"


# Short FASTA headers (max length for makeblastdb is 50)
sed -i 's/|.*//' plsdb_Host_PA_1L.fasta
echo "The headers of the fasta sequences were shortened for BLAST"

# Make BLAST database with the PLSDB for Pseudomonas aeruginosa
echo "Creating BLAST database with the PLSDB for Pseudomonas aeruginosa"
makeblastdb -in plsdb_Host_PA_1L.fasta -parse_seqids -dbtype nucl


echo "Script finished!"
