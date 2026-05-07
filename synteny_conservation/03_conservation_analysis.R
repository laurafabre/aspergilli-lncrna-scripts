# Load necessary libraries
library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(gridExtra)
library(svglite)
library(RColorBrewer)
library(amelia)
library(ggridges)
library(dunn.test)
library(ggsignif)
my_palette <- unname(amelia_colors())

setwd("~/orthology_synteny/synteny")

# Function to read BED file and create a data frame of gene coordinates
read_bed_file <- function(bed_file) {
  bed_df <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed_df) <- c("chrom", "start", "end", "gene_id", "strand")
  return(bed_df)
}

# Paths to the BED files
bed_file_afla <- "../bed_files/afla_lncrna_intergenic.bed"
bed_file_afum <- "../bed_files/afum_lncrna_intergenic.bed"
bed_file_anid <- "../bed_files/anid_lncrna_intergenic.bed"
bed_file_anig <- "../bed_files/anig_lncrna_intergenic.bed"

# Read the BED files
afla_bed <- read_bed_file(bed_file_afla)
afum_bed <- read_bed_file(bed_file_afum)
anid_bed <- read_bed_file(bed_file_anid)
anig_bed <- read_bed_file(bed_file_anig)

#################################### ORTHOBASE SYNTHENY #################################

# Read the synteny results file
file_path <- "intergenic/intergenic.fam"  # Modify this path
gene_data <- read_delim(file_path, delim = "\t", col_names = FALSE)
orthologs <- read_delim(file = 'orthologs.txt')

# Step 2: Assign column names for clarity
colnames(gene_data) <- c("family", "gene_id")

# Step 3: Separate gene_id into components (gene ID and species)
# Create a new dataframe with separate species columns
gene_data <- gene_data %>%
  separate(gene_id, into = c("gene_id", "species"), sep = "\\|", extra = "merge") %>%
  mutate(species = sub(".*\\|", "", species))  # Extract only the species part
gene_data <- gene_data %>%
  mutate(gene_id = paste0(gene_id, "_", species)) 

write_csv(gene_data, file = 'families_lnc_intergenic_ID.csv')

# Step 4: Create a dataframe for each species
# Extract unique gene IDs for each species
gene_list <- gene_data %>%
  group_by(species) %>%
  summarise(unique_genes = list(unique(gene_id)))


## Optional: If you need to separate the genes for specific species
#syntenic_afla <- gene_data %>% filter(species == "afla") %>% pull(gene_id) %>% unique()
#syntenic_afum <- gene_data %>% filter(species == "afum") %>% pull(gene_id) %>% unique()
#syntenic_anid <- gene_data %>% filter(species == "anid") %>% pull(gene_id) %>% unique()
#syntenic_anig <- gene_data %>% filter(species == "anig") %>% pull(gene_id) %>% unique()


# Function to extract non-conserved lncRNA for each species
get_non_conserved_lncrna <- function(conserved_genes, bed_file) {
  bed_data <- read_bed_file(bed_file)
  non_conserved_genes <- setdiff(bed_data$gene_id, conserved_genes)
  return(non_conserved_genes)
}

# Get non-conserved lncRNA for each species
non_conserved_afla <- setdiff(afla_bed$gene_id, syntenic_afla)
non_conserved_afum <- setdiff(afum_bed$gene_id, syntenic_afum)
non_conserved_anid <- setdiff(anid_bed$gene_id, syntenic_anid)
non_conserved_anig <- setdiff(anig_bed$gene_id, syntenic_anig)


all_data <- read.csv('../all_data_prediction.csv')
all_data_lncrna <- read.csv('../all_lncrna_prediction.csv', sep = ',')

# filter only the intergenic lncrna
filtered_data <- all_data_lncrna %>%
  filter(Class_code == "u") 

filtered_data <- filtered_data %>%
  mutate(synteny = case_when(
    ID %in% syntenic_afla ~ "conserved",
    ID %in% syntenic_afum ~ "conserved",
    ID %in% syntenic_anid ~ "conserved",
    ID %in% syntenic_anig ~ "conserved",
    TRUE ~ "not conserved" # If not in any of the syntenic lists
  ))

write_csv(filtered_data, file = 'intergenic_with_orthology.csv')


# Statatics and boostrap
filtered_length <- all_data %>%
  mutate(Type_2= case_when(
    Type_2 == 'protein' ~ 'PC genes',
    Type_2 == 'intergenic'  ~ "LincRNA",
    Type_2 == "antisense" ~ "LncRNA AS",
    TRUE ~ Type_2 # If not in any of the syntenic lists
  ))
filtered_length$Type_2<-factor(filtered_length$Type_2,levels = c("PC genes","LincRNA","LncRNA AS"))

# boostrap with replacemnet
# Step 1: Find the minimum group size
# Find the minimum group size for each species and Type_2
min_size_per_species <- filtered_length %>%
  group_by(Spp, Type_2) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Spp) %>%
  summarise(min_size = min(n)) # Get the minimum size per species

# Merge the minimum sizes back with the original dataset
balanced_data_length<- filtered_length %>%
  left_join(min_size_per_species, by = "Spp") %>%
  group_by(Spp, Type_2) %>%
  sample_n(size = first(min_size), replace = TRUE) %>%  # Use the species-specific min size
  ungroup()

# Step 3: Ensure factor levels remain correct
balanced_data_length$Type_2 <- factor(balanced_data_length$Type_2, 
                               levels = c("PC genes", "LincRNA", "LncRNA AS"))

balanced_data_length$Species <- factor(balanced_data_length$Species, 
                                levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))

# Step 1: Perform Dunn's Test for each species
dunn_results <- list()
for(speciess in unique(balanced_data_length$Species)) {
  # Filter data for the current species
  data_species <- balanced_data_length %>% filter(Species == speciess)
  
  # Ensure the filter produced some data
  if(nrow(data_species) > 0) {
    # Perform Dunn's Test using the filtered data
    dunn_test <- dunn.test(
      x = log2(data_species$Length),  # Log-transformed values for the current species
      g = data_species$Type_2,              # Grouping variable for the current species
      method = "bonferroni"
    )
    
    # Store results in a list, using the species name as the key
    dunn_results[[speciess]] <- dunn_test
  } else {
    # Handle cases where there is no data for a species
    cat("No data for species:", speciess, "\n")
  }
}

# Extract pairwise p-values (and format them)
# We will loop over the Dunn test results for each species

# Initialize an empty data frame for the results
pval_labels <- data.frame()

# Loop through each species in the Dunn's test results
for(species in names(dunn_results)) {
  # Get the adjusted p-values and comparisons
  pvals <- dunn_results[[species]]$P.adjusted
  comparisons <- dunn_results[[species]]$comparisons
  
  # Create a data frame for the current species
  temp_df <- data.frame(
    comparison = comparisons,     # Comparison names
    pval = pvals,                 # Adjusted p-values
    Species = species             # Add the species name
  )
  
  # Bind the temporary data frame to the main pval_labels data frame
  pval_labels <- rbind(pval_labels, temp_df)
}

pairwise_results_length <- balanced_data_length %>%
  group_by(Species) %>%
  rstatix::pairwise_wilcox_test(Length ~ Type_2, p.adjust.method = "fdr")

# plot
my_comparisons <- list(
  c("PC genes", "LincRNA"),  
  c("PC genes", "LncRNA AS"),
  c("LincRNA", "LncRNA AS")
)

y_positions <- c(18.5, 17.5, 16.5)

# Create the plot with pairwise comparisons
length_plot_1 <- ggplot(filtered_length, aes(x = Type_2, y = log2(Length), fill = Type_2)) +
  geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 2) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 2, 
              vjust = 0.5, 
              size = 0.5) +
  guides(fill = guide_legend(title = "")) +
  labs(title = "Length", y = "Log2(Length)", x = "") +
  scale_fill_manual(values =  c(my_palette[6], my_palette[2], "#1b9e77"))+  #c(my_palette[6], my_palette[2],my_palette[1], my_palette[8], my_palette[4])
  theme(legend.position = "none")+
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) + # Remove minor gridlines)
  panel_border()

ggsave(length_plot_1, filename = "../../paper/plots/length_plot.svg", width = 8, height = 4, dpi = 600, device = "svg")
ggsave(length_plot_1, filename = "../../paper/plots/length_plot.tiff", width = 8, height = 4, dpi = 600, device = "tiff", compression = "lzw")


# conservation analysis
conserved_ids <- unique(c(syntenic_afla, syntenic_afum, syntenic_anid, syntenic_anig))

filtered_data <- all_data %>%
  mutate(Type_2= case_when(
    Type_2 == 'intergenic' & ID %in% conserved_ids ~ "intergenic conserved",
    Type_2 == "intergenic" & !(ID %in% conserved_ids) ~ "intergenic non-conserved",
    TRUE ~ Type_2 # If not in any of the syntenic lists
  ))
filtered_data$Type_2<-factor(filtered_data$Type_2,levels = c("protein","intergenic conserved","intergenic non-conserved","antisense"))

intergenic_lncrna <- all_data_lncrna %>%
  filter(Class_code == "u")  # Select only antisense (x) class

# Step 1: Filter only antisense lncRNAs from all_data_lncrna
antisense_lncrna <- all_data_lncrna %>%
  filter(Class_code == "x")  # Select only antisense (x) class

# Step 2: Function to classify antisense lncRNAs as conserved or not conserved
classify_conservation <- function(df, orthologs_df) {
  df %>%
    rowwise() %>%  # Process each row individually
    mutate(Conservation_Status = case_when(
      Spp == "afla" & Gene %in% orthologs_df$afla ~ "conserved",
      Spp == "afum" & Gene %in% orthologs_df$afum ~ "conserved",
      Spp == "anid" & Gene %in% orthologs_df$anid ~ "conserved",
      Spp == "anig" & Gene %in% orthologs_df$anig ~ "conserved",
      TRUE ~ "not_conserved"  # Default to not conserved if no match is found
    )) %>%
    ungroup()  # Ungroup after rowwise operation
}

# Step 3: Apply the function to classify conservation status
antisense_lncrna_with_orthologs <- classify_conservation(antisense_lncrna, orthologs)

write_csv(antisense_lncrna_with_orthologs, file = 'antisense_with_orthology.csv')

conserved_antisense <- antisense_lncrna_with_orthologs %>%
  filter(Conservation_Status == "conserved") %>%
  pull(ID)

filtered_data <- all_data %>%
  mutate(Type_2 = case_when(
    # Intergenic lncRNA classification
    Type_2 == 'intergenic' & ID %in% conserved_ids ~ "intergenic conserved",
    Type_2 == "intergenic" & !(ID %in% conserved_ids) ~ "intergenic non-conserved",
    
    # Antisense lncRNA classification
    Type_2 == 'antisense' & ID %in% conserved_antisense ~ "antisense conserved",
    Type_2 == 'antisense' & !(ID %in% conserved_antisense) ~ "antisense non-conserved",
    
    # Default case to keep Type_2 unchanged
    TRUE ~ Type_2
  ))

# Adjust the factor levels to include the new antisense categories
filtered_data$Type_2 <- factor(filtered_data$Type_2, 
                               levels = c("protein", 
                                          "intergenic conserved", 
                                          "intergenic non-conserved", 
                                          "antisense conserved", 
                                          "antisense non-conserved"))


my_palette <- unname(amelia_colors())

# Define comparision
# my_comparisons <- list(
#   c("protein", "intergenic_conserved"), 
#   c("protein", "intergenic_not_conserved"),
#   c("protein", "antisense_conserved"),
#   c("protein", "antisense_not_conserved"),
#   c("intergenic_conserved", "intergenic_not_conserved"),
#   c("intergenic_conserved", "antisense_conserved"),
#   c("intergenic_conserved", "antisense_not_conserved"),
#   c("intergenic_not_conserved", "antisense_conserved"),
#   c("intergenic_not_conserved", "antisense_not_conserved"),
#   c("antisense_conserved", "antisense_not_conserved")
# )

my_comparisons <- list(
  c("protein", "intergenic conserved"),  
  c("protein", "antisense conserved"),
  c("intergenic conserved", "intergenic non-conserved"),
  c("antisense conserved", "antisense non-conserved")
)

# Step 1: Perform Dunn's Test for each species
dunn_results <- list()

for(species in unique(filtered_data$Species)) {
  data_species <- filtered_data %>% filter(Species == species)
  dunn_test <- dunn.test(x = log2(data_species$Length), 
                         g = data_species$Type_2, 
                         method = "bonferroni")
  
  # Store results in a list
  dunn_results[[species]] <- dunn_test
}

# Step 2: Extract pairwise p-values (and format them)
# We will loop over the Dunn test results for each species

# Initialize an empty data frame for the results
pval_labels <- data.frame()

# Loop through each species in the Dunn's test results
for(species in names(dunn_results)) {
  # Get the adjusted p-values and comparisons
  pvals <- dunn_results[[species]]$P.adjusted
  comparisons <- dunn_results[[species]]$comparisons
  
  # Create a data frame for the current species
  temp_df <- data.frame(
    comparison = comparisons,     # Comparison names
    pval = pvals,                 # Adjusted p-values
    Species = species             # Add the species name
  )
  
  # Bind the temporary data frame to the main pval_labels data frame
  pval_labels <- rbind(pval_labels, temp_df)
}

# Adjust these values based on your data range
y_positions <- c(15, 16, 14, 14)

# Create the plot with pairwise comparisons
length_plot <- ggplot(filtered_data, aes(x = Type_2, y = log2(Length), fill = Type_2)) +
  geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 2) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 3, 
              vjust = 0.5, 
              size = 0.5) +
  guides(fill = guide_legend(title = "")) +
  labs(title = "Length", y = "Log2(Length)", x = "") +
  scale_fill_manual(values = c(my_palette[6], my_palette[1], my_palette[2],my_palette[8], my_palette[4])) +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  panel_border()

ggsave(filename = "LENGTH_boxplot_convervation.svg" , plot = length_plot, width = 14, height = 7)

# gc
GC_content <- read.csv('../GC_intergenic.csv', sep = ',')
GC_content$ID <- sub("\\|.*", "", GC_content$gene)  # Remove everything after the first "|"
GC_content$ID <- paste(GC_content$ID, GC_content$spp, sep = "_")  # Append species

GC_filtered_data <- GC_content %>%
  mutate(Class_code_2= case_when(
    Class_code_2 == 'intergenic' & ID %in% conserved_ids ~ "intergenic conserved",
    Class_code_2 == "intergenic" & !(ID %in% conserved_ids) ~ "intergenic non-conserved",
    # Antisense lncRNA classification
    Class_code_2== 'antisense' & ID %in% conserved_antisense ~ "antisense conserved",
    Class_code_2 == 'antisense' & !(ID %in% conserved_antisense) ~ "antisense non-conserved",
    TRUE ~ Class_code_2 # If not in any of the syntenic lists
  ))

GC_filtered_data$Class_code_2<-factor(GC_filtered_data$Class_code_2,
                                      levels = c("protein", 
                                                 "intergenic conserved", 
                                                 "intergenic non-conserved", 
                                                 "antisense conserved", 
                                                 "antisense non-conserved"))

# Step 1: Perform Dunn's Test for each species
dunn_results <- list()

for(species in unique(GC_filtered_data$Species)) {
  data_species <- GC_filtered_data %>% filter(Species == species)
  dunn_test <- dunn.test(x = data_species$GC, 
                         g = data_species$Class_code_2, 
                         method = "bonferroni")
  
  # Store results in a list
  dunn_results[[species]] <- dunn_test
}

# Step 2: Extract pairwise p-values (and format them)
# We will loop over the Dunn test results for each species

# Step 1: Initialize an empty data frame for the results
pval_labels <- data.frame()

# Step 2: Loop through each species in the Dunn's test results
for(species in names(dunn_results)) {
  # Get the adjusted p-values and comparisons
  pvals <- dunn_results[[species]]$P.adjusted
  comparisons <- dunn_results[[species]]$comparisons
  
  # Create a data frame for the current species
  temp_df <- data.frame(
    comparison = comparisons,     # Comparison names
    pval = pvals,                 # Adjusted p-values
    Species = species             # Add the species name
  )
  
  # Bind the temporary data frame to the main pval_labels data frame
  pval_labels <- rbind(pval_labels, temp_df)
}

y_positions <- c(0.95, 1, 0.8, 0.8)
GC_plot <- ggplot(GC_filtered_data, aes(x=Class_code_2, y=GC,fill=Class_code_2)) +
    geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 0.1) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 3, 
              vjust = 0.5, 
              size = 0.5) +
  
  guides(fill = guide_legend(title = "")) +
  labs(title = "GC", 
       y = "GC content", x="") + 
  scale_fill_manual(values = c( my_palette[6], my_palette[1], my_palette[2], my_palette[8],my_palette[4]))+theme(legend.position = "none")+
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))  +
  panel_border()

ggsave(filename = "gc_boxplot_convervation.svg" , plot = GC_plot, width = 14, height = 7)



# number
filtered_data_number <- filtered_data %>%
  group_by(Type_2) %>%          # Group by Type_2
  mutate(n = n()) %>%          # Count occurrences and create a new column "n"
  ungroup()      

filtered_data_number <- filtered_data_number %>% filter(!filtered_data_number$Type == 'pc')
# Count occurrences for each Type_2 and Species
counted_data <- filtered_data_number %>%
  group_by(Type_2, Species) %>%  # Group by Type_2 and Species
  summarise(n = n(), .groups = 'drop')  # Count occurrences and drop grouping

# Reorder Type_2 for consistent plotting
counted_data$Type_2 <- factor(counted_data$Type_2,
                              levels = c("intergenic conserved", 
                                         "intergenic non-conserved", 
                                         "antisense conserved", 
                                         "antisense non-conserved"))

# Create a bar plot using the counted data
number_plot <- ggplot(counted_data, aes(x = Type_2, y = n, fill = Type_2)) + 
  facet_grid(. ~ Species) +  # Facet by Species
  geom_bar(stat = "identity") +  # Use 'identity' to plot precomputed counts
  geom_text(aes(label = n), size = 5, position = position_stack(vjust = 0.5)) +  # Add text labels
  scale_fill_manual(values = c( my_palette[2], my_palette[1],  my_palette[8], my_palette[4])) +
  guides(fill = guide_legend(title = "Class code")) +
  ylab("Number of lncRNAs") + xlab("") +
  scale_x_discrete(position = "bottom") +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", aspect.ratio = 1.5) +
  panel_border()

ggsave(filename = "number_lncrna_conserved.svg", plot = number_plot, width = 12)

# Create broader categories (e.g., intergenic, antisense)
filtered_data_2 <- filtered_data %>%
  mutate(Type_main = case_when(
    grepl("intergenic", Type_2) ~ "intergenic",
    grepl("antisense", Type_2) ~ "antisense",
    TRUE ~ Type_2  # Keep other types unchanged
  ),
  conserved_status = case_when(
    Type_2 == "intergenic conserved" ~ "conserved",
    Type_2 == "antisense conserved" ~ "conserved",
    Type_2 == "intergenic non-conserved" ~ "non-conserved",
    Type_2 == "antisense non-conserved" ~ "non-conserved"
  ))

filtered_data_2 <- filtered_data_2 %>% filter(!filtered_data_2$Type == 'pc')
# Group by the new main category and species
counted_data <- filtered_data_2 %>%
  group_by(Type_main, conserved_status, Species) %>%
  summarise(n = n(), .groups = 'drop')

# Calculate the percentages for each conserved_status within each Type_main and Species group
counted_data <- filtered_data_2 %>%
  group_by(Type_main, Species) %>%
  summarise(total = n(), .groups = 'drop') %>%  # Calculate total for each Type_main-Species
  right_join(
    filtered_data_2 %>%
      group_by(Type_main, conserved_status, Species) %>%
      summarise(n = n(), .groups = 'drop'), 
    by = c("Type_main", "Species")
  ) %>%
  mutate(percentage = (n / total) * 100) 


# Plot the stacked bar plot with percentages
stacked_plot <- ggplot(counted_data, aes(x = Type_main, y = percentage, fill = conserved_status)) + 
  facet_grid(. ~ Species) +  # Facet by Species
  geom_bar(stat = "identity", position = "stack") +  # Stack bars by conserved_status
  geom_text(aes(label = sprintf("%.1f%%", percentage)),  # Format percentages to 1 decimal place
            size = 4, position = position_stack(vjust = 0.5)) +  # Add percentage labels inside the bars
  scale_fill_manual(values = c(my_palette[8], my_palette[1])) +  # Custom colors for conserved_status
  guides(fill = guide_legend(title = "Conserved Status")) +
  ylab("Percentage of lncRNAs") + xlab("") +
  scale_x_discrete(position = "bottom") +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",  # Adjust legend position
        aspect.ratio = 1.5) +
  panel_border()

ggsave(filename = "number_lncrna_stacked_conserved.svg", plot = stacked_plot, width = 12)


# expression
expression <- read.csv('../expression_prediction.csv', sep = ',')
expression <-expression %>%
  separate(gene, into = c("gene", "species"), sep = "\\|", extra = "merge") %>%
  mutate(species = sub(".*\\|", "", Spp))

expression <- expression %>% mutate(gene= paste0(gene, "_", Spp))
colnames(expression)[3] <- 'ID'

# data for plot 1
expression$Type_2<-factor(expression$Type_2,
                                        levels = c("protein",
                                                   "intergenic",
                                                   "antisense"))
expression$Species <- factor(expression$Species, levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))


# data for plot with conservation info
expression_filtered_data <- expression %>%
  mutate(Type_2= case_when(
    Type_2 == 'protein' ~ 'PC genes',
    Type_2 == 'intergenic' & ID %in% conserved_ids ~ "LincRNA conserved",
    Type_2 == "intergenic" & !(ID %in% conserved_ids) ~ "LincRNA non-conserved",
    # # Antisense lncRNA classification
    Type_2 == 'antisense' ~ "LncRNA AS",
    
    # Default case to keep Type_2 unchanged
    TRUE ~ Type_2
  ))
expression_filtered_data$Type_2<-factor(expression_filtered_data$Type_2,
                                        levels = c("PC genes", 
                                                   "LincRNA conserved", 
                                                   "LincRNA non-conserved",
                                                   "LncRNA AS"))
expression_filtered_data$Species <- factor(all_data$Species, levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))

expression <- expression %>%
  mutate(Type_2= case_when(
    Type_2 == 'protein' ~ 'PC genes',
    Type_2 == 'intergenic'  ~ "LincRNA",
    Type_2 == 'antisense' ~ "LncRNA AS",
    
    # Default case to keep Type_2 unchanged
    TRUE ~ Type_2
  ))

# boostrap with replacemnet
# Step 1: Find the minimum group size
min_size_per_species <- expression %>%
  group_by(species, Type_2) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(species) %>%
  summarise(min_size = min(n)) # Get the minimum size per species


# Step 2: Perform sampling with replacement for each group
# Merge the minimum sizes back with the original dataset
balanced_data_expression <- expression %>%
  left_join(min_size_per_species, by = "species") %>%
  group_by(species, Type_2) %>%
  sample_n(size = first(min_size), replace = TRUE) %>%  # Use the species-specific min size
  ungroup()

# Step 3: Ensure factor levels remain correct
balanced_data_expression$Type_2 <- factor(balanced_data_expression$Type_2, 
                                          levels = c("PC genes", "LincRNA", 
                                                    "LncRNA AS"))

balanced_data_expression$Species <- factor(balanced_data_expression$Species, 
                                           levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))

# Step 1: Perform Dunn's Test for each species
dunn_results <- list()
for(speciess in unique( balanced_data_expression$Species)) {
  # Filter data for the current species
  data_species <-  balanced_data_expression %>% filter(Species == speciess)
  
  # Ensure the filter produced some data
  if(nrow(data_species) > 0) {
    # Perform Dunn's Test using the filtered data
    dunn_test <- dunn.test(
      x = log2(data_species$value + 0.01),  # Log-transformed values for the current species
      g = data_species$Type_2,              # Grouping variable for the current species
      method = "bonferroni"
    )
    
    # Store results in a list, using the species name as the key
    dunn_results[[speciess]] <- dunn_test
  } else {
    # Handle cases where there is no data for a species
    cat("No data for species:", speciess, "\n")
  }
}


# Step 2: Extract pairwise p-values (and format them)
# We will loop over the Dunn test results for each species

# Step 1: Initialize an empty data frame for the results
pval_labels <- data.frame()

# Step 2: Loop through each species in the Dunn's test results
for(species in names(dunn_results)) {
  # Get the adjusted p-values and comparisons
  pvals <- dunn_results[[species]]$P.adjusted
  comparisons <- dunn_results[[species]]$comparisons
  
  # Create a data frame for the current species
  temp_df <- data.frame(
    comparison = comparisons,     # Comparison names
    pval = pvals,                 # Adjusted p-values
    Species = species             # Add the species name
  )
  
  # Bind the temporary data frame to the main pval_labels data frame
  pval_labels <- rbind(pval_labels, temp_df)
}


y_positions <- c(16.5, 15, 14)

my_comparisons <- list(
  c("PC genes", "LincRNA"),  
  c("PC genes", "LncRNA AS"),
  c("LincRNA", "LncRNA AS")
)
expression_plot_1 <- ggplot(balanced_data_expression, aes(x=Type_2, y=log2(as.numeric(as.character(value))+0.01),fill=Type_2))+
  geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 18) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 2, 
              vjust = 0.5, 
              size = 0.5) +
  guides(fill = guide_legend(title = "")) +
  labs(title = "Expression", 
       y = "Log2(mean TPM+0.01)", x="") +
  labs(title = "Expression", 
       y = "Log2(mean TPM+0.01)", x="") + 
  scale_fill_manual(values =  c(my_palette[6], my_palette[2], "#1b9e77"))+  #c(my_palette[6], my_palette[2],my_palette[1], my_palette[8], my_palette[4])
  theme(legend.position = "none")+
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() )+  # Remove minor gridlines)  +
  panel_border()
  
ggsave("../../paper/plots/expression_plot.svg", width = 8, height = 4, dpi = 600, device = "svg")
ggsave("../../paper/plots/expression_plot.tiff", width = 8, height = 4, dpi = 600, device = "tiff", compression = "lzw")


# only intergenic conservation
# data for plot with conservation info
expression_filtered_data <- expression %>%
  mutate(Type_2= case_when(
    Type_2 == "protein" ~ "PC genes",
    Type_2 == 'intergenic' & ID %in% conserved_ids ~ "LincRNA conserved",
    Type_2 == "intergenic" & !(ID %in% conserved_ids) ~ "LincRNA non-conserved",
    # # Antisense lncRNA classification
    Type_2 == 'antisense' ~ "LncRNA AS",
    
    # Default case to keep Type_2 unchanged
    TRUE ~ Type_2
  ))
expression_filtered_data$Type_2<-factor(expression_filtered_data$Type_2,
                                        levels = c("PC genes", 
                                                   "LincRNA conserved", 
                                                   "LincRNA non-conserved",
                                                   "LncRNA AS"))
expression_filtered_data$Species <- factor(all_data$Species, levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))


my_comparisons <- list(
  c("PC genes", "LincRNA conserved"),  
  c("LincRNA conserved", "LincRNA non-conserved"),
  c("LncRNA AS", "LincRNA conserved"),
  c("LncRNA AS", "LincRNA non-conserved")
)

expression_plot <- ggplot(expression_filtered_data, aes(x=Type_2, y=log2(as.numeric(as.character(value))+0.01),fill=Type_2))+
  geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 19) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 2, 
              vjust = 0.5, 
              size = 0.5) +
  guides(fill = guide_legend(title = "")) +
  labs(title = "Expression", 
       y = "Log2(mean TPM+0.01)", x="") + 
  scale_fill_manual(values = c(my_palette[6], my_palette[2],my_palette[1], my_palette[8], my_palette[4]))+theme(legend.position = "none")+
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1))  +
  panel_border()

ggsave(filename = "expression_conservation.svg", plot = expression_plot, width = 13)


# boostrap with replacemnet
# Step 1: Find the minimum group size
# Find the minimum group size for each species and Type_2
min_size_per_species <- expression_filtered_data %>%
  group_by(species, Type_2) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(species) %>%
  summarise(min_size = min(n)) # Get the minimum size per species


# Step 2: Perform sampling with replacement for each group
# Merge the minimum sizes back with the original dataset
balanced_data_expression <- expression_filtered_data %>%
  left_join(min_size_per_species, by = "species") %>%
  group_by(species, Type_2) %>%
  sample_n(size = first(min_size), replace = TRUE) %>%  # Use the species-specific min size
  ungroup()


# Step 3: Ensure factor levels remain correct
balanced_data_expression$Type_2 <- factor(balanced_data_expression$Type_2, 
                               levels = c("PC genes", "LincRNA conserved", 
                                          "LincRNA non-conserved", "LncRNA AS"))

balanced_data_expression$Species <- factor(balanced_data_expression$Species, 
                                levels = c("A. flavus", "A. fumigatus", "A. nidulans", "A. niger"))

y_positions <- c(15.5, 14.5, 13, 12)

expression_plot <- ggplot(balanced_data_expression, aes(x=Type_2, y=log2(as.numeric(as.character(value))+0.01),fill=Type_2))+
  geom_boxplot() +
  facet_grid(. ~ Species) +
  # Global Kruskal-Wallis test
  stat_compare_means(method = "kruskal.test", label.y = 17) +  
  # Add pairwise comparisons with ggsignif and custom y positions
  geom_signif(comparisons = my_comparisons, 
              map_signif_level = TRUE, 
              y_position = y_positions,  # Use the defined y positions
              textsize = 2, 
              vjust = 0.5, 
              size = 0.5) +
  guides(fill = guide_legend(title = "")) +
  labs(title = "Expression", 
       y = "Log2(mean TPM+0.01)", x="") + 
  scale_fill_manual(values =  c(my_palette[6], "#d95f02", "#7570b3", "#1b9e77"))+  #c(my_palette[6], my_palette[2],my_palette[1], my_palette[8], my_palette[4])
  theme(legend.position = "none")+
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        legend.position = "top",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) + # Remove minor gridlines)  +
  panel_border()

ggsave("../../paper/plots/expression_conserved_plot.svg", width = 8, height = 4, dpi = 600, device = "svg")
ggsave("../../paper/plots/expression_conserved_plot.tiff", width = 8, height = 4, dpi = 600, device = "tiff", compression = "lzw")

density_plot <- ggplot(balanced_data_expression, aes(x = log2(as.numeric(as.character(value))+0.01), fill = Type_2)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~Species) +
  scale_fill_manual(values = c(my_palette[6], "#d95f02", "#7570b3", "#1b9e77")) +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) +
  labs(title = "Distribution of normalized counts", x = "Log2(mean TPM+0.01)", fill = "Gene Type") +
  panel_border()


ggsave("../../paper/plots/expression_conserved_density.svg", width = 9, height = 5, dpi = 600, device = "svg")
ggsave("../../paper/plots/expression_conserved_density.tiff", width = 9, height = 5, dpi = 600, device = "tiff", compression = "lzw")


