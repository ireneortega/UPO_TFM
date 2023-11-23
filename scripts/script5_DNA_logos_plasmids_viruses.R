library(dplyr)

blast_IMGVR <- read.csv(file = "I-C/blast_IMGVR_PA_vs_I-C_system.tab", sep = "\t", header = FALSE, dec = ".") # Read previous BLAST file
names(blast_IMGVR) <- c("qseqid", "sseqid", "pident", "qcovs", "qcovhsp", "length", "qlen", "slen", "evalue", "qstart", "qend", "sstart", "send") # Add column headers to BLAST file

# PART 1. Add subject sequence (that is the virus sequence) to the BLAST file 
blast_IMGVR$sseqid <- gsub("ref\\||gb\\|", "", blast_IMGVR$sseqid) # Edit subject sequence ID column (sseqid): BLAST file (BLAST changed NZ_CP025052.1 to ref|NZ_CP025052.1|)
blast_IMGVR$sseqid <- gsub("\\|", "", blast_IMGVR$sseqid)

# Apply same procedure as Spacer2PAM: remove hits of the same spacer to the same subject and keep the best alignment (low evalue), that is the first result
blast_IMGVR_unique <- blast_IMGVR %>% 
  dplyr::distinct(qseqid, sseqid, .keep_all = TRUE)

# Read IMGVR database (be careful last line needs to be empty)
IMGVR_list <- seqinr::read.fasta(file = "IMGVR_nucl_PA_1L.fasta", as.string = TRUE)

# Get headers and sequences from IMGVR_list and put into dataframe
sseqid <- c()
for (IMGVR_seq in 1:length(IMGVR_list)){
  sseqid[IMGVR_seq] = attr(IMGVR_list[[IMGVR_seq]], "name")
}

genomeSequence <- c()
for (IMGVR_seq in 1:length(IMGVR_list)){
  genomeSequence[IMGVR_seq] = unlist(IMGVR_list[[IMGVR_seq]])
}

IMGVR_data <- data.frame(sseqid, genomeSequence)
print("The data from the IMGVR database was converted to a dataframe")
# Join BLAST dataframe and IMGVR_data dataframe: add the genome sequence to the BLAST output
blast_withGenome <- blast_IMGVR_unique %>%
  dplyr::left_join(IMGVR_data, by="sseqid")

print("The dataframe with the IMGVR information was added to the BLAST output dataframe")


# PART 2. Extract the flanking nucleotides of the region of the subject than align with the spacer

# Function to determine if alignment is on plus or minus strand of the subject (plus is TRUE, minus is FALSE)
whichStrand = function (start, end) {
  ifelse(end - start > 0, TRUE, FALSE)
}

# Use previous function to add orientation of subject strand
blast_withGenome_stranded <- blast_withGenome %>%
  dplyr::mutate(strand_subject = whichStrand(sstart, send))

print("The orientation of the subject strand was determined for each alignment")


# Extract flankLength nucleotides upstream and downstream of alignment from genome sequence based on direction of alignment
flankLength <- 10 # Default PAM sequence length used by Spacer2PAM (keep same parameters)

# Functions to extract position of flankLength nucleotides depending on orientation of alignment)
upStrandAccount <- function (seq, strand_subject, substart, querystart, length) {
  ifelse(strand_subject == TRUE, substr(seq, substart - querystart + 1 - flankLength, substart - querystart), substr(seq, substart + querystart - length - flankLength, substart + querystart - 1 - length))
}

downStrandAccount <- function (seq, strand_subject, substart, querystart, length) {
  ifelse(strand_subject == TRUE, substr(seq, substart - querystart + length + 1, substart - querystart + length + flankLength), substr(seq, substart + querystart, substart + querystart + flankLength - 1))
}

revcomplement <- function (strand_subject, flank) {
  ifelse(strand_subject == FALSE, sapply(lapply(strsplit(chartr("atgc","tacg",as.character(flank)), NULL), rev), paste, collapse=""), sapply(lapply(strsplit(chartr("atgc","tacg",as.character(flank)), NULL), rev), paste, collapse=""))
}

# Use previous functions to extract flankLength nucleotides
blast_withGenome_stranded_upnDown <- blast_withGenome_stranded %>%
  dplyr::mutate(upstream = upStrandAccount(genomeSequence, strand_subject, sstart, qstart, length))%>%
  dplyr::mutate(downstream = downStrandAccount(genomeSequence, strand_subject, sstart, qstart, length))%>%
  dplyr::mutate(upstream.rev = revcomplement(strand_subject, upstream))%>%
  dplyr::mutate(downstream.rev = revcomplement(strand_subject, downstream))

print("Extracting flanking nucleotides upstream and downstream for each hit")

# Create function to extract corresponding flanking nucleotides
whichFlank = function(orientation_array, strand_subject, flank, revcompflank){
  ifelse((orientation_array == "Forward") == strand_subject, flank, revcompflank)
}

# The orientation of array needs to be added to the blast_withGenome_stranded_upnDown dataframe
spacers_subtype_info <- read.csv(file = "spacers_subtype_info.tab", sep = "\t", header = TRUE) # Read spacers dataframe

# Filter spacers representing this CRISPR array (I-C)
df_spacers_I_C <- subset(spacers_subtype_info, Array.Orientation != "ND")
df_spacers_I_C <- subset(df_spacers_I_C, Evidence_lev == "4")
df_spacers_I_C <- subset(df_spacers_I_C, Subtype == "I-C")

# Add new column qseqid that contains the spacer ID as it appears in the BLAST file (and now in blast_withGenome_stranded_upnDown dataframe). This can be done by merging the string subtype+A+column_array+S+column spacer
df_spacers_I_C$qseqid <- paste("I-CA", df_spacers_I_C$Array, "S", df_spacers_I_C$Spacer, sep ="")

# Add column Array.Orientation from the df_spacers_I_C dataframe to the BLAST file (now blast_withGenome_stranded_upnDown dataframe)
blast_withGenome_stranded_upnDown$Array.Orientation <- NA # Create empty column so far
for (i in 1:nrow(blast_withGenome_stranded_upnDown)) {
  qseqid_value <- blast_withGenome_stranded_upnDown$qseqid[i]
  orientation_value <- df_spacers_I_C$Array.Orientation[df_spacers_I_C$qseqid == qseqid_value]
  blast_withGenome_stranded_upnDown$Array.Orientation[i] <- orientation_value
}
print("The orientation of the array was added")

# Create list of upstream PAM sequences: get corresponding flanking nucleotides
print("Creating list of upstream PAM sequences...")
alignmentUp <- c()
if(nrow(blast_withGenome_stranded_upnDown) >= 1) {
  for (i in 1:nrow(blast_withGenome_stranded_upnDown)) {
    alignmentUp[i] <- as.character(whichFlank(blast_withGenome_stranded_upnDown$Array.Orientation[i],blast_withGenome_stranded_upnDown$strand_subject[i],blast_withGenome_stranded_upnDown$upstream[i],blast_withGenome_stranded_upnDown$downstream.rev[i]))
    
    ID_subject <- append(ID_subject, as.character(blast_withGenome_stranded_upnDown$sseqid[i]))
  }
} else{alignmentUp = c("")}

# Remove those flanking sequences with a length shorter than 10 nt
alignmentUp <- alignmentUp[nchar(alignmentUp) == flankLength]


# PART 3. Generate table of nucleotides upstream (PAM) protospacer for each subject
print("Creating dataframe with PAM nucleotides for each subject...")
## Create dataframe with ID_subject and alignmentUp
length(alignmentUp) == length(ID_subject) # Check that the number of elements in these lists is equal, otherwise the dataframe would be unbalanced (needs to be [TRUE])
ID_alignmentUp_df <- data.frame(ID_subject, alignmentUp)

write.table(ID_alignmentUp_df, file = "I-C/I-C_ID_Upstream_PAM_nucl.tab", sep = "\t", row.names = FALSE) # save results


# PART 4. Plot upstream WebLogo and save
print("Plotting DNA seq logo...")

upstreamLogoOutput <- "I-C/Protospacers I-C Upstream PAM frequence IMGVR.png" # output filename

library(ggplot2)

if(length(alignmentUp)>=1){
  ggplot2::ggplot() + 
    ggseqlogo::geom_logo(as.character(toupper(alignmentUp)), seq_typ="DNA") + 
    ggseqlogo::theme_logo() + 
    ggplot2::theme(axis.text.x = element_text(size=32),axis.text.y = element_text(size=32), axis.title=element_text(size=48) ,axis.ticks = element_line(size = 1.5), axis.ticks.length = unit(20, "pt"))
  
  ggplot2::ggsave(upstreamLogoOutput, width = 10, height = 7, units = "in")
}
