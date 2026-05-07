library(tidyverse)
library(igraph)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(amelia)


# Define species and base paths
species_list <- c("afla", "afum", "anid", "anig")
base_path <- "~/data_drive/Aspergilli_paper/data/processed/wgcna/network_R_data"

all_spp_data <- list()

for (spp in species_list) {
  message("Parsing connectivity for: ", spp)
  
  # Load your pre-calculated iK file
  ik_path <- sprintf("%s/%s/%s_connectivity.txt", base_path, spp, spp)
  
  if(file.exists(ik_path)){
    ik_data <- read.table(ik_path, header = TRUE, sep = "\t")
    
    # Organize data
    ik_data <- ik_data %>%
      mutate(
        Species = spp,
        # kTotal is the Weighted Degree (Global)
        # kWithin is the Intramodular Connectivity (Local)
        SimpleType = ifelse(Type == "Protein", "Protein", "lncRNA")
      )
    
    all_spp_data[[spp]] <- ik_data
  } else {
    warning(paste("File missing for:", spp))
  }
}

# Combine all species
master_df <- bind_rows(all_spp_data)
master_df$Species <- recode(master_df$Species,
                            "anid" = "A. nidulans",
                            "anig" = "A. niger",
                            "afla" = "A. flavus",
                            "afum" = "A. fumigatus")

master_df$Species <- factor(master_df$Species,
                            levels = c("A. nidulans",
                                       "A. niger",
                                       "A. flavus",
                                       "A. fumigatus"))
# -------------------------------------------------------------------
# PLOT A: Barplot (Chi-square test for composition)
# -------------------------------------------------------------------
# Since Barplots show counts, we usually report a global Chi-square 
# in the text rather than on the bars.
species_list <- c("afla", "afum", "anid", "anig")
genome_counts_list <- list()

for (spp in species_list) {
  # Path to your raw expression table
  raw_file_path <- sprintf("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/%s_expression.txt", spp)
  
  if (file.exists(raw_file_path)) {
    # Read only the first column (gene IDs) to save memory/time
    raw_ids <- read.table(raw_file_path, header = TRUE, sep = "\t", check.names = FALSE)
    gene_names <- as.character(raw_ids[[1]]) # Assuming first column is the ID
    
    # Identify types based on MSTRG prefix
    n_lncrna <- sum(grepl("^MSTRG", gene_names))
    n_protein <- sum(!grepl("^MSTRG", gene_names))
    
    genome_counts_list[[spp]] <- data.frame(
      Species = spp,
      SimpleType = c("lncRNA", "Protein"),
      TotalAnnotated = c(n_lncrna, n_protein)
    )
    message(sprintf("%s: Found %s proteins and %s lncRNAs in raw file.", spp, n_protein, n_lncrna))
  } else {
    warning(sprintf("File not found: %s", raw_file_path))
  }
}

genome_counts <- bind_rows(genome_counts_list)
genome_counts <- genome_counts %>%
  mutate(Species = case_when(
    Species == "afla" ~ "A. flavus",
    Species == "afum" ~ "A. fumigatus",
    Species == "anid" ~ "A. nidulans",
    Species == "anig" ~ "A. niger",
    TRUE ~ Species # Keep as is if already changed
  ))

plot_a <- master_df %>%
  group_by(Species, SimpleType) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Species) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ggplot(aes(x = Species, y = Percentage, fill = SimpleType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("lncRNA" = "#c83737ff", "Protein" = "#62a0caFF"), name = "Gene type") +
  labs(title = "Network Composition (Expressed Genes)", y = "% of genes in network", x = NULL) +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 14),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        # panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) + # Remove minor gridlines)  +
  panel_border()

# 1. Get counts of genes that made it into the network
network_counts <- master_df %>%
  distinct(Species, gene, .keep_all = TRUE) %>%
  group_by(Species, SimpleType) %>%
  summarise(InNetwork = n(), .groups = 'drop')

# 2. Merge and calculate percentage
# 2. Merge, calculate percentage, and reorder Species
plot_data_genome <- network_counts %>%
  mutate(Species = as.character(Species)) %>%
  left_join(genome_counts, by = c("Species", "SimpleType")) %>%
  mutate(PercOfGenome = (InNetwork / TotalAnnotated) * 100) %>%
  # Correct way to reorder the Species factor inside a pipe:
  mutate(Species = factor(Species, 
                          levels = c("A. nidulans", 
                                     "A. niger", 
                                     "A. flavus", 
                                     "A. fumigatus")))


# 3. Create the barplot
plot_a_2 <- ggplot(plot_data_genome, aes(x = Species, y = PercOfGenome, fill = SimpleType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(PercOfGenome, 1), "%")), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("lncRNA" = "#c83737ff", "Protein" = "#62a0caFF"), 
                    name = "Gene type") +
  labs(title = "Genome representation in co-expression network", 
       subtitle = "Proportion of annotated genes retained after network construction filters",
       y = "% of total annotated genes", 
       x = NULL) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 14),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        # panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) + # Remove minor gridlines)  +
  panel_border()

# -------------------------------------------------------------------
# PLOT B: Connectivity with Wilcoxon Stats
# -------------------------------------------------------------------

# Define the comparisons we want to test
# Comparing each lncRNA type against Protein
# Define comparisons
my_comparisons <- list(
  c("LincRNA", "Protein"),
  c("LncRNA AS", "Protein")
)

plot_b <- master_df %>%
  select(Species, Type, kTotal, kWithin) %>%
  pivot_longer(cols = c(kTotal, kWithin),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(MetricLabel = ifelse(Metric == "kTotal",
                              "Global weighted\ndegree (kTotal)",
                              "Intramodular\nconnectivity (kWithin)"))%>%
  ggplot(aes(x = Type, y = Value, fill = Type)) +
  geom_violin(alpha = 1, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               fill = "white", alpha = 0.6) +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     step.increase = 0.12,
                     tip.length = 0.02,
                     hide.ns = FALSE) +
  facet_grid(MetricLabel ~ Species, scales = "free_y") +
  scale_fill_manual(values = c("LincRNA" = "#70bc6bFF",
                               "LncRNA AS" = "#edbc4cFF",
                               "Protein" = "#62a0caFF")) +
  labs(title = "Topological feature distribution",
       y = "Connectivity",
       x = NULL) +
  theme_minimal_grid() +
  theme(strip.text = element_text(face = "italic", size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
       # panel.background = element_rect(fill = "white", color = "black"),  # White panel, black border
        plot.background = element_rect(fill = "white"),  # White background
        panel.spacing = unit(1, "lines"),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank() ) + # Remove minor gridlines)  +
  panel_border()

# -------------------------------------------------------------------
# PLOT C: histogram of nodes vs. degree
# -------------------------------------------------------------------
plot_c <- master_df %>%
  ggplot(aes(x = kTotal, fill = SimpleType)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Species, scales = "free", ncol = 4) +
  scale_fill_manual(
    values = c("lncRNA" = "#c83737", "Protein" = "#62a0ca"),
    name = "Gene type"
  ) +
  labs(
    title = "Connectivity Distribution",
    x = "Soft Connectivity (kTotal)",
    y = "Density"
  ) +
  theme_minimal_grid() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "italic", size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  panel_border()

# # 1. Calculate the distribution data
# dist_data <- master_df %>%
#   # Rounding connectivity to integers to create frequency bins
#   mutate(k_bin = round(kTotal)) %>%
#   group_by(Species, k_bin) %>%
#   summarise(Frequency = n(), .groups = 'drop') %>%
#   filter(k_bin > 0) # Log(0) is undefined
# 
# # 2. Plotting the Power Law
# plot_scale_free <- ggplot(dist_data, aes(x = k_bin, y = Frequency)) +
#   geom_point(alpha = 0.5, color = "midnightblue") +
#   # Use log-log scales to check for a straight line (linear decay)
#   scale_x_log10() + 
#   scale_y_log10() +
#   geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
#   facet_wrap(~Species, ncol = 4) +
#   labs(title = "Scale-Free Topology Check",
#        x = "Log10(Global Weighted Degree)",
#        y = "Log10(Frequency)") +
#   theme_bw()
# 
# print(plot_scale_free)
# -------------------------------------------------------------------
# COMBINE AND SAVE
# -------------------------------------------------------------------
final_fig <- (plot_a_2 / plot_c / plot_b) + 
  plot_layout(heights = c(1.5, 1.2, 3)) + # Adjusting ratios to give violins more space
  plot_annotation(tag_levels = 'A') # Automatically adds A, B, C labels
final_fig
ggsave("results/figure_draft/supplementary/Supp_Fig_Topology_coexpression_with_Stats2.pdf", final_fig, width = 10, height = 12)
ggsave("results/figure_draft/supplementary/Supp_Fig_Topology_coexpression_with_Stats2.svg", final_fig, width = 10, height = 12)
ggsave("results/figure_draft/supplementary/Supp_Fig_Topology_coexpression_with_Stats2.png", final_fig, width = 10, height = 12)


# 1. Calculate the distribution data
dist_data <- master_df %>%
  # Rounding connectivity to integers to create frequency bins
  mutate(k_bin = round(kTotal)) %>%
  group_by(Species, k_bin) %>%
  summarise(Frequency = n(), .groups = 'drop') %>%
  filter(k_bin > 0) # Log(0) is undefined

# 2. Plotting the Power Law
plot_scale_free <- ggplot(dist_data, aes(x = k_bin, y = Frequency)) +
  geom_point(alpha = 0.5, color = "midnightblue") +
  # Use log-log scales to check for a straight line (linear decay)
  scale_x_log10() + 
  scale_y_log10() +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
  facet_wrap(~Species, ncol = 4) +
  labs(title = "Scale-Free Topology Check",
       x = "Log10(Global Weighted Degree)",
       y = "Log10(Frequency)") +
  theme_bw()

print(plot_scale_free)

library(broom)

# Calculate R2 for each species
scale_free_stats <- dist_data %>%
  group_by(Species) %>%
  do(glance(lm(log10(Frequency) ~ log10(k_bin), data = .))) %>%
  select(Species, r.squared) %>%
  mutate(label = paste0("R² = ", round(r.squared, 2)))

print(scale_free_stats)

plot_scale_free_with_stats <- plot_scale_free +
  geom_text(data = scale_free_stats, 
            aes(x = Inf, y = Inf, label = label), 
            hjust = 1.1, vjust = 1.5, 
            inherit.aes = FALSE, 
            fontface = "bold", size = 4)

print(plot_scale_free_with_stats)

# Define your species and file paths
library(data.table)
library(dplyr)

species_list <- c("anid", "anig", "afla", "afum")
sample_summary <- list()

for (spp in species_list) {
  file_path <- sprintf("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/%s_expression.txt", spp)
  
  if (file.exists(file_path)) {
    # fread automatically detects if it is tab, comma, or space separated
    # fill = TRUE handles rows with inconsistent lengths
    temp_data <- fread(file_path, header = TRUE, nrows = 1)
    
    n_samples <- ncol(temp_data) - 1
    sample_names <- colnames(temp_data)[-1]
    
    sample_summary[[spp]] <- data.frame(
      SpeciesShort = spp,
      TotalSamples = n_samples,
      SampleList = paste(sample_names, collapse = ", ")
    )
  } else {
    message("File not found for: ", spp)
  }
}

final_sample_counts <- bind_rows(sample_summary)
print(final_sample_counts)

