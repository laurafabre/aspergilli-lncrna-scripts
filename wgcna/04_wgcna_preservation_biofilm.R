# ===========================================================================
# MODULE PRESERVATION ANALYSIS
# Original data: log2(TPM+1)  (what you used for network)
# Biofilm data: vst from DESeq2 (normalized counts)
# ===========================================================================

library(WGCNA)

# 1. PREPARE ORIGINAL DATASET (log2 TPM - same as network construction)
# ----------------------------------------------------------------------

# Your original TPM data
orig_TPM <- get(paste0(spp, "_TPM"))  # e.g., afum_TPM

# Apply same filtering as network construction
orig_expr <- orig_TPM[rowSums(orig_TPM >= 0.1) >= (0.8 * ncol(orig_TPM)), ]

# Log2 transform (EXACTLY as you did for network)
orig_expr_log <- log2(orig_expr + 1)

# Transpose to samples x genes
orig_data <- t(orig_expr_log)

# Get module colors (must match gene order)
genes_orig <- colnames(orig_data)
module_colors_orig <- moduleColors[genes_orig]

# Remove any genes without module assignment
if(any(is.na(module_colors_orig))) {
  keep <- !is.na(module_colors_orig)
  orig_data <- orig_data[, keep]
  module_colors_orig <- module_colors_orig[keep]
}

# 2. PREPARE BIOFILM DATASET (vst from DESeq2)
# ---------------------------------------------

# Use your existing vst-normalized data
# norm_expr is already vst-transformed (from your earlier code)
biofilm_data <- t(norm_expr)  # samples x genes

# 3. MATCH GENES BETWEEN DATASETS
# --------------------------------

common_genes <- intersect(colnames(orig_data), colnames(biofilm_data))
message(paste("Genes in common:", length(common_genes)))

# Subset both datasets
orig_data_common <- orig_data[, common_genes]
biofilm_data_common <- biofilm_data[, common_genes]
module_colors_common <- module_colors_orig[common_genes]

# 4. CREATE MULTI-SET FORMAT
# ---------------------------

multiExpr <- list(
  Original = list(data = orig_data_common),
  Biofilm = list(data = biofilm_data_common)
)

multiColor <- list(
  Original = module_colors_common,
  Biofilm = NULL
)

# 5. RUN PRESERVATION ANALYSIS
# -----------------------------
library(parallel)  # for detectCores

set.seed(12345)
nCores <- detectCores() - 1
enableWGCNAThreads(nCores)

message("Starting module preservation analysis...")
message("Original data: log2(TPM+1)")
message("Biofilm data: vst from DESeq2")

preservation <- modulePreservation(
  multiExpr,
  multiColor,
  referenceNetworks = 1,
  nPermutations = 200,
  verbose = 3,
  parallelCalculation = TRUE,
  quickCor = 0
)

# 6. EXTRACT RESULTS
# -------------------

# Extract the Zsummary preservation data
Z_preservation <- preservation$preservation$Z$ref.Original$inColumnsAlsoPresentIn.Biofilm

# Convert to data frame
Z_df <- as.data.frame(Z_preservation)

# Add module names from row names
Z_df$module <- rownames(Z_df)

# Remove gold module if present
Z_df <- Z_df[Z_df$module != "gold", ]

# View all modules with Zsummary.pres
message("\n=== ALL MODULES - PRESERVATION STATISTICS ===")
print(Z_df[, c("module", "moduleSize", "Zsummary.pres")] %>% 
        arrange(desc(Zsummary.pres)))

# Your specialist modules (from earlier)
specialist_modules <- c("purple", "bisque4", "darkseagreen3", "paleturquoise", 
                        "darkgreen", "lightsteelblue1", "darkred", "darkseagreen4",
                        "ivory", "grey60", "lightcoral", "turquoise", "orange",
                        "brown2", "saddlebrown", "darkviolet", "darkslateblue",
                        "magenta4", "coral2")

# Filter to specialist modules
specialist_preservation <- Z_df %>%
  filter(module %in% specialist_modules) %>%
  dplyr::select(module, moduleSize, Zsummary.pres) %>%
  arrange(desc(Zsummary.pres))

message("\n=== SPECIALIST MODULES - PRESERVATION STATISTICS ===")
print(specialist_preservation)

# 7. FOCUS ON SPECIALIST MODULES
# --------------------------------
# Categorize preservation strength
specialist_preservation <- specialist_preservation %>%
  mutate(
    preservation_category = case_when(
      Zsummary.pres > 10 ~ "High preservation",
      Zsummary.pres > 2 ~ "Moderate preservation",
      TRUE ~ "Low preservation"
    )
  )

# View by category
message("\n=== PRESERVATION BY CATEGORY ===")
print(specialist_preservation %>%
        group_by(preservation_category) %>%
        summarise(
          count = n(),
          modules = paste(module, collapse = ", ")
        ))

# Print actual Z values for specialist modules
specialist_preservation %>% 
  dplyr::select(module, Zsummary.pres, preservation_category) %>%
  arrange(desc(Zsummary.pres))


# ===========================================================================
# CREATE PRESERVATION SCATTER PLOT (Standard in WGCNA papers)
# ===========================================================================

library(ggplot2)
library(ggrepel)

# Get preservation data for ALL modules (not just specialists)
Z_all <- as.data.frame(preservation$preservation$Z$ref.Original$inColumnsAlsoPresentIn.Biofilm)
Z_all$module <- rownames(Z_all)
Z_all <- Z_all[Z_all$module != "gold", ]

# Add module sizes if not already present
if(!"moduleSize" %in% colnames(Z_all)) {
  module_sizes <- table(module_colors_common)
  Z_all$moduleSize <- as.numeric(module_sizes[Z_all$module])
}

# Add category for coloring
Z_all$category <- case_when(
  Z_all$Zsummary.pres > 10 ~ "High preservation",
  Z_all$Zsummary.pres > 2 ~ "Moderate preservation",
  TRUE ~ "Low preservation"
)

# Create the standard preservation plot
p_preservation <- ggplot(Z_all, aes(x = moduleSize, y = Zsummary.pres)) +
  # Add threshold lines (standard in the field)
  geom_hline(yintercept = 10, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue", size = 0.8) +
  
  # Add points for all modules
  geom_point(aes(color = category), size = 3, alpha = 0.7) +
  
  # Highlight your specialist modules with labels
  geom_point(data = Z_all %>% filter(module %in% specialist_modules),
             aes(x = moduleSize, y = Zsummary.pres),
             color = "black", size = 3.5, shape = 21, stroke = 1.5) +
  
  # Add labels for specialist modules
  geom_text_repel(data = Z_all %>% filter(module %in% specialist_modules),
                  aes(label = module, color = category),
                  size = 3.5, fontface = "bold", 
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0.1, max.overlaps = Inf) +
  
  # Annotate threshold meanings (standard labels)
  annotate("text", x = max(Z_all$moduleSize, na.rm = TRUE) * 0.7, 
           y = 11, label = "Highly preserved (Z > 10)", 
           color = "red", hjust = 0, size = 4) +
  annotate("text", x = max(Z_all$moduleSize, na.rm = TRUE) * 0.7, 
           y = 3, label = "Moderately preserved (Z > 2)", 
           color = "blue", hjust = 0, size = 4) +
  
  # Labels and theme
  labs(x = "Module size (number of genes)",
       y = "Preservation Zsummary",
       title = "Module preservation in independent biofilm dataset",
       color = "Preservation level") +
  scale_color_manual(values = c("High preservation" = "darkgreen",
                                "Moderate preservation" = "orange",
                                "Low preservation" = "grey50")) +
  theme_minimal_grid() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14)) +
  panel_border() +
  scale_y_continuous(limits = c(0, max(Z_all$Zsummary.pres, na.rm = TRUE) * 1.1))

print(p_preservation)

#figure caption
#Scatter plot of Zsummary preservation statistics vs. module size. Each point represents a module, colored by preservation category. Dashed red line indicates Z=10 (high preservation threshold), dashed blue line indicates Z=2 (moderate preservation threshold). Modules with Z>10 (turquoise, darkgreen, paleturquoise, lightcoral, saddlebrown, ivory) show strong preservation across datasets. Modules with 2<Z<10 (bisque4, orange, darkred, darkseagreen3, purple, magenta4, darkslateblue, coral2, darkviolet) show moderate preservation. Modules with Z<2 (lightsteelblue1, brown2, grey60, darkseagreen4) show weak preservation, suggesting condition-specific regulation.

# Save
ggsave("results/figure_draft/Fig5D_Module_Preservation.pdf", 
       p_preservation, width = 10, height = 6)
ggsave("results/figure_draft/Fig5D_Module_Preservation.png", 
       p_preservation, width = 10, height = 6, dpi = 300)

#===============================================================================
# MedianRank plot (alternative visualization)
# Get median rank data from the observed statistics
medianRank_data <- as.data.frame(preservation$preservation$observed$ref.Original$inColumnsAlsoPresentIn.Biofilm)
medianRank_data$module <- rownames(medianRank_data)
medianRank_data <- medianRank_data[medianRank_data$module != "gold", ]
# Add module sizes (you may need to add these if not present)
if(!"moduleSize" %in% colnames(medianRank_data)) {
  module_sizes <- table(module_colors_common)
  medianRank_data$moduleSize <- as.numeric(module_sizes[medianRank_data$module])
}

# Add preservation category based on Zsummary (from your earlier analysis)
# You'll need to merge with Z_all to get categories
medianRank_data <- medianRank_data %>%
  left_join(Z_all %>% dplyr::select(module, category), by = "module")

# Create median rank plot
p_medianRank <- ggplot(medianRank_data, aes(x = moduleSize, y = medianRank.pres)) +
  geom_point(aes(color = category), size = 3, alpha = 0.7) +
  geom_text_repel(data = medianRank_data %>% filter(module %in% specialist_modules),
                  aes(label = module, color = category),
                  size = 3.5, fontface = "bold", 
                  box.padding = 0.5, point.padding = 0.3) +
  # Lower median rank = better preserved (so we want points at the bottom)
  labs(x = "Module size (number of genes)",
       y = "Median rank (lower = better preserved)",
       title = "Module preservation - Median rank",
       color = "Preservation level") +
  scale_color_manual(values = c("High preservation" = "darkgreen",
                                "Moderate preservation" = "orange",
                                "Low preservation" = "grey50")) +
  theme_minimal_grid() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  panel_border()

print(p_medianRank)
