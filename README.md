# PAM prediction in Pseudomonas aeruginosa

Here you will the scripts and data associated with the research article "Search for PAM sequences associated with CRISPR-Cas systems in _Pseudomonas aeruginosa_ and their enrichment in plasmids and phages" [under preparation]. A brief comment on the purpose of each script is provided here:
- **script1_spacers_df.R**
Construction of the _Pseudomonas aeruginosa_ spacers dataframe for the `df2fasta()` function of the Spacer2PAM library. The spacers were collected from the output of CRISPRCasFinder and filtered based on known CRISPR-Cas array orientation and evidence level equal to 4, and known subtype determined by CRISPRCas-Typer.

- **script2_PAM_prediction.R**
After the information regarding each spacer has been collected, the PAM for each CRISPR-Cas subtype will be predicted using Spacer2PAM.

- **script3_PLSDB_IMGVR_sequences_filtering.sh**
From the PLSDB database v2020_06_23_v2 and the IMG/VR v3 high-quality genomes database, the _Pseudomonas aeruginosa_ sequences will be filtered.

- **script4_plasmids_viruses_BLAST.sh**
The spacers representing each _Pseudomonas aeruginosa_ CRISPR-Cas subtype will be blasted against the _P. aeruginosa_ plasmids and viruses from the PLSDB and IMG/VR databases, respectively (an example is provided for _P. aeruginosa_ CRISPR-Cas subtype I-C and IMG/VR database).

- **script5_DNA_logos_plasmids_viruses.R**
For the _P. aeruginosa_ plasmids and viruses from the PLSDB and IMG/VR databases, respectively, that are recognized by each _Pseudomonas aeruginosa_ CRISPR-Cas system subtype, the DNA logo will be constructed (an example is provided for _P. aeruginosa_ CRISPR-Cas subtype I-C and IMG/VR database).

- **script6_PAM_freq_GC.sh**
Determination of the occurrence of the PAM and GC content in the foreign sequences (plasmids and viruses from the PLSDB and IMG/VR databases, respectively) recognized by each _Pseudomonas aeruginosa_ CRISPR-Cas system (an example is provided for _P. aeruginosa_ CRISPR-Cas subtype I-C and IMG/VR database).
