library(taxonomizr)
library(Spacer2PAM)
library("jsonlite")
library(dplyr)

# Create list to store the location of input information: strain ID, JSON file from CRISPRCasFinder and crisprs_all.tab from CRISPRCasTyper 
strains_list <- list.dirs(path = './ccfinder/', full.names = FALSE, recursive = FALSE)
json_list <- list.files(list.dirs(path = './ccfinder/', full.names = TRUE, recursive = FALSE), pattern = "json", full.names = TRUE)
subtype_list <- list.files(list.dirs(path = './cctyper/', full.names = TRUE, recursive = FALSE), pattern = "crisprs_all.tab", full.names = TRUE)

# PART 1: create empty dataframe to store the spacers information required for the df2fasta() function
df_spacers <- data.frame(Strain = character(),
                         Spacers = character(),
                         Array.Orientation = character(),
                         Repeat = character(),
                         Array = character(),
                         Spacer = int(),
                         Evidence_lev = int(),
                         Contig = character(),
                         Start = character(),
                         End = character())

for (i in 1:length(json_list)) { # For each JSON file, that is, for each strain… 
  
  # PART 1A: extract spacers information from CRISPRCasFinder
  
  strain_name <- strains_list[i] # Get strain ID
  
  json_file <- json_list[i] # Get JSON file associated with the strain
  json_file <- fromJSON(txt = json_file) # Read JSON file (jsonlite package required as with rjson package an error was encoutered)
  json_data_frame <- as.data.frame(json_file) # Convert JSON file to a dataframe
  
  all_contigs_lista <- json_data_frame$Sequences.Crisprs # From the JSON dataframe, extract the list of dataframes that contain information about the CRISPR arrays in all contigs of the genome
  contigs_with_CRISPR_lista <- Filter(function(x) dim(x)[1] > 0, all_contigs_lista) # From the previous list, remove empty dataframes that correspond to contigs without CRISPR arrays
  
  if (length(contigs_with_CRISPR_lista) > 0) { # When multiple contigs harbor CRISPR array(s…
    for (contig in 1:length(contigs_with_CRISPR_lista)) { # For each contig with a CRISPR array…
      contig_with_CRISPR <- contigs_with_CRISPR_lista[contig] # Select one of these contigs
      array_s_ID_per_contig <- contig_with_CRISPR[[1]]$Name # Get name of CRISPR array(s) within contig. The length of this variable determines whether the contig contains multiple CRISPR arrays or not
      
      if ((length(array_s_ID_per_contig)) == 1) { # If the contig has just one CRISPR array…
        CRISPR_array <- contig_with_CRISPR # …the CRISPR array within contig is in fact the CRISPR array
        
        array_orientation <- CRISPR_array[[1]]$Potential_Orientation # Get CRISPR array orientation and modify orientation coding (+ means Forwards and – means Reverse)
        if (array_orientation == "+") {
          array_orientation <- "Forward"
        } else if (array_orientation == "-") {
          array_orientation <- "Reverse"
        }
        
        DR_seq <- CRISPR_array[[1]]$DR_Consensus # Get repeat sequence within CRSPR array
        array_ID <- CRISPR_array[[1]]$Name # Get CRISPR array ID
        
        contig_id <- sub("_[0-9]+$", "", array_ID) # Get name of contig that contains the CRISPR array, that is equal to the CRISPR array ID except last underscore followed by a number
        start_position <- CRISPR_array[[1]]$Start # Get start position of CRISPR array
        end_position <- CRISPR_array[[1]]$End # Get end position of CRISPR array
        
        evidence_level <- CRISPR_array[[1]]$Evidence_Level # Get evidence_level of CRISPR array. This will be used to filter out those with a value below 4.
        
        array_regions <- CRISPR_array[[1]]$Regions # Get regions of the CRISPR array (spacers and repeat)
        spacers_seq <- array_regions[[1]][array_regions[[1]]$Type == "Spacer", ]$Sequence # Get spacer sequences
        
        for (i in 1:length(spacers_seq)) { # For each spacer within the CRISPR array, get information of interest and add to a new dataframe
          df_single_spacer <- data.frame(Strain = strain_name,
                                         Spacers = spacers_seq[i], # Get the corresponding spacer sequence
                                         Array.Orientation = array_orientation,
                                         Repeat = DR_seq,
                                         Array = array_ID,
                                         Spacer = i, # It beeds to be a number
                                         Evidence_lev = evidence_level, 
                                         Contig = contig_id,
                                         Start = start_position,
                                         End = end_position)
          
          df_spacers <- rbind(df_spacers, df_single_spacer) # Add the previous dataframe to the initial one
        } # Close for-loop (create dataframe with the information of the spacer)
        
      } else { # …or when more than one CRISPR array are found in a contig…
        for (array_ID in array_s_ID_per_contig) { # For each of these CRISPR arrays in the same contig…
          CRISPR_array <- contig_with_CRISPR[[1]][contig_with_CRISPR[[1]]$Name == array_ID, ]
          
          contig_id <- sub("_[0-9]+$", "", array_ID) # Get name of contig that contains the CRISPR array, that is equal to the CRISPR array ID except last underscore followed by a number
          start_position <- CRISPR_array$Start # Get Get start position of CRISPR array
          end_position <- CRISPR_array$End # Get end start position of CRISPR array
          
          array_orientation <- CRISPR_array$Potential_Orientation # Get CRISPR array orientation and modify orientation coding (+ means Forwards and – means Reverse)      
          if (array_orientation == "+") {
            array_orientation <- "Forward"
          } else if (array_orientation == "-") {
            array_orientation <- "Reverse"
          }
          
          DR_seq <- CRISPR_array$DR_Consensus # Get repeat sequence within CRSPR array
          
          evidence_level <- CRISPR_array$Evidence_Level # Get evidence_level of CRISPR array. This will be used to filter out those with a value below 4
          
          array_regions <- CRISPR_array$Regions # Get regions of the CRISPR array (spacers and repeat)
          spacers_seq <- array_regions[[1]][array_regions[[1]]$Type == "Spacer", ]$Sequence # Get spacer sequences
          
          for (i in 1:length(spacers_seq)) { # For each spacer within the CRISPR array, get information of interest and add to a new dataframe
            df_single_spacer <- data.frame(Strain = strain_name,
                                           Spacers = spacers_seq[i], # Get the corresponding spacer sequence
                                           Array.Orientation = array_orientation,
                                           Repeat = DR_seq,
                                           Array = array_ID,
                                           Spacer = i, # It beeds to be a number
                                           Evidence_lev = evidence_level, 
                                           Contig = contig_id,
                                           Start = start_position,
                                           End = end_position)
            
            df_spacers <- rbind(df_spacers, df_single_spacer) # Add the previous dataframe to the initial one
          } # Close for-loop (create dataframe with the information of the spacer)
          
          
        } # Close for-loop when more than one CRISPR array is found in the same contig 
      } # Close if-statement that considers the number of CRISPR arrays within contig
      
      
    } # Close for-loop for each contig that contains CRISPR array
  } # Close if-statement that considers if the strain has some CRISPR array
  
  print(paste0("Processed spacers for ", strain_name))
  
} # Close for-loop for processing each strain 


print("Adding subtype information for spacers dataframe... Remember that CRISPR arrays both detected by CRISPRCasFinder and CRISPRCasTyper but with position mismatches of 100 or less in Start or End or/and mistmatches in Repeat Sequence are considered equal arrays for subtype assessing")

subtype_data <- data.frame()
for (subtype_file in subtype_list) { # For each CRISPRCasTyper output file
  single_subtype_data <- read.table(file = subtype_file, header = TRUE, sep = "\t") # Read the CRISPRCasTyper output file
  single_subtype_data$Contig <- sub("\\..*", "", single_subtype_data$Contig) # Remove "." and subsequent caracters in the column Contig
  single_subtype_data <- single_subtype_data[c('Contig', 'Start', 'End', 'Consensus_repeat', 'Subtype')] # Keep only columns of interest for the comparison between CRISPRCasFinder and CRISPRCasTyper, and the column with the Subtype (empty)
  colnames(single_subtype_data)[which(names(single_subtype_data) == "Consensus_repeat")] <- "Repeat"
  subtype_data <- rbind(subtype_data, single_subtype_data)
}

df_spacers_subtype_raw <- left_join(df_spacers, subtype_data, by = c("Contig", "Start", "End", "Repeat")) # Add the Subtype when Contig, Start, End and Repeat are equal. WARNING: those spacers within CRISPR arrays whose Start or End is different, or Repeat sequence is not identical, between the two tools are missed


## Add the Subtype for the spacers having been missed: starting point
df_spacers_subtype <- merge(df_spacers_subtype_raw, subtype_data, by = c("Contig"), all.x = TRUE)

df_spacers_subtype <- df_spacers_subtype %>% # Keep the spacers without Subtype defined so far
  filter(is.na(Subtype.x))

df_spacers_subtype$Subtype.x <- ifelse(abs(df_spacers_subtype$Start.x - df_spacers_subtype$Start.y) <= 100 & abs(df_spacers_subtype$End.x - df_spacers_subtype$End.y) <= 100 & !is.na(df_spacers_subtype$Repeat.y) | grepl(df_spacers_subtype$Repeat.x, df_spacers_subtype$Repeat.y) == TRUE | grepl(df_spacers_subtype$Repeat.y, df_spacers_subtype$Repeat.x) == TRUE, 
                                       df_spacers_subtype$Subtype.y, 
                                       df_spacers_subtype$Subtype.x) # Replace NA value in column Subtype.x (spacers without Subtype defined) by the columna Subtype.y from the Subtype column in CRISPRCasTyper when the positions of the array between the two tools differ 100 bp or less, and when the Repeat sequence of one tool is contained in the other tool (this allows certain mistmaches).

df_spacers_subtype <- unique(df_spacers_subtype) # Remove duplicates just created when the previous conditions are not satisfied

# Remove columns from CRISPRCasTyper as their information has been added to the dataframe with the spacers information from CRISPRCasFinder
df_spacers_subtype$Start.y <- NULL
df_spacers_subtype$End.y <- NULL
df_spacers_subtype$Subtype.y <- NULL
df_spacers_subtype$Repeat.y <- NULL

# Rename columns (remove subscripts)
colnames(df_spacers_subtype)[colnames(df_spacers_subtype) == 'Start.x'] <- 'Start' 
colnames(df_spacers_subtype)[colnames(df_spacers_subtype) == 'End.x'] <- 'End'
colnames(df_spacers_subtype)[colnames(df_spacers_subtype) == 'Subtype.x'] <- 'Subtype'
colnames(df_spacers_subtype)[colnames(df_spacers_subtype) == 'Repeat.x'] <- 'Repeat'

## Add the Subtype for the spacers having been missed: ending point



df_spacers_subtype <- left_join(df_spacers_subtype_raw, df_spacers_subtype, by = names(df_spacers_subtype)[1:(ncol(df_spacers_subtype)-1)]) # Add the new spacers within already known CRISPR system to the initial spacers within defined Subtype 

df_spacers_subtype$Subtype.x <- ifelse(is.na(df_spacers_subtype$Subtype.x), df_spacers_subtype$Subtype.y, df_spacers_subtype$Subtype.x) # Merge Subtype columns: add the new Subtype (Subtype.y) to the existing Subtype column (Subtype.x)

df_spacers_subtype$Subtype.y <- NULL # Remove column Subtype.y as its information has already been added
colnames(df_spacers_subtype)[colnames(df_spacers_subtype) == 'Subtype.x'] <- 'Subtype' # Rename columns (remove subscripts)

# Remove rows whose value in all except last column (Subtype) is duplicated and the values of the columns Subtype is not missing (NA), that is the Subtype is defined
df_spacers_subtype <- df_spacers_subtype[order(df_spacers_subtype$Subtype),] # As duplicated() function in the next step only keeps the first duplicated observation, the column Subtype will be reorder so that the NA rows for the Subtype column are move to the end and therefore will be removed
df_spacers_subtype <- df_spacers_subtype[!duplicated(df_spacers_subtype[,colnames(df_spacers_subtype) != "Subtype"]),]

# In the previous reorderin process, the spacers dataframe has been disordere. The dataframe will be ordered by index
df_spacers_subtype$index <- as.numeric(row.names(df_spacers_subtype))
df_spacers_subtype <- df_spacers_subtype[order(df_spacers_subtype$index), ]
df_spacers_subtype$index <- NULL

print("Subtype added to spacers dataframe!")


# Save spacers dataframes
write.table(df_spacers, file = "./spacers_info.tab", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(df_spacers_subtype, file = "./spacers_subtype_info.tab", sep = "\t", row.names = FALSE, col.names = TRUE)

print("Dataframe with spacers information and subtype saved to file spacers_subtype_info.tab")


# FILTER SPACERS. From the dataframe with spacers information and CRISPR array, remove rows that matches array with unknown orientation, evidence level below 4 and unknown subtype. 
df_spacers_filtered <- df_spacers_subtype[df_spacers_subtype$Array.Orientation != "ND" & df_spacers_subtype$Evidence_lev == "4" & !(is.na(df_spacers_subtype$Subtype)), ]

write.table(df_spacers_filtered, file = "./spacers_for_df2fasta.tab", sep = "\t", row.names = FALSE, col.names = TRUE)



# PARTE 2: convert spacers dataframe in fasta format required for the df2fasta() function
print("Spacers with evidence level equal to 4, known direction and subtype will be proccessed")

# Remove unnecessary columns for df2fasta()
df_spacers_filtered_sql <- df_spacers_filtered %>%
  dplyr::select(Spacers, Array.Orientation, Repeat, Array, Spacer, Subtype)

subtype_list_df <- split(df_spacers_filtered_sql, df_spacers_filtered_sql$Subtype) # Create individual list of dataframes for each CRISPR array subtype

for (CRISPR_subtype in names(subtype_list_df)) { # For each CRISPR array subtype
  subtype_list_df[[CRISPR_subtype]]$Subtype <- NULL # Remove subtype column
  
  setCRISPRInfo(genus = "Pseudomonas",
                species = "aeruginosa",
                strain = CRISPR_subtype,
                crisprSystemNumber = "CRISPR-Cas") # It will only be considered for the fasta filename
  
  df2fasta(subtype_list_df[[CRISPR_subtype]]) # Fasta file will be created in the working directory
  
  print(paste0("FASTA file for ", CRISPR_subtype, " CRISPR-Cas array subtype created!"))
  
}
