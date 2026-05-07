# parse orthofinder results to get hte 1-to-1 orthologs between the 4 spp.
# Define file names
single_copy_file <- "../Orthogroups_SingleCopyOrthologues.txt"
orthogroups_file <- "../Orthogroups.tsv"
output_file <- "orthologs.txt"

# Read single copy orthogroup IDs
single_copy_ids <- readLines(single_copy_file)

# Read the Orthogroups.tsv file
orthogroups_data <- read.table(orthogroups_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter the data to keep only rows with IDs in single_copy_ids
filtered_data <- orthogroups_data[orthogroups_data$Orthogroup %in% single_copy_ids, ]

# Write the filtered data to a new TSV file
write.table(filtered_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
