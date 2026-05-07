# 4_filter_coexpression.R
library(dplyr)

# Load precomputed raw correlation results
coexp_raw_df <- readRDS("../Aspergilli_paper/data/processed/correlation_function/coexpression_raw.rds")

# Set thresholds
cor_cutoff <- 0.80
pval_cutoff <- 0.01

# Apply filtering
coexp_filtered <- coexp_raw_df %>%
  filter(!is.na(cor), abs(cor) > cor_cutoff, padj < pval_cutoff)

# Save filtered results
output_file <- paste0("../Aspergilli_paper/data/processed/correlation_function/coexpressed_pairs_corr", cor_cutoff * 100, "_padj", pval_cutoff, ".csv")
write.csv(coexp_filtered, output_file, row.names = FALSE)
