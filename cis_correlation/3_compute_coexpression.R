# 3_compute_coexpression.R
library(Hmisc)     # for rcorr()
library(dplyr)     # for data wrangling
library(tibble)    # for tidy output
library(purrr)     # for functional programming

load("/home/laura/Documents/Aspergillus_lncRNA/Aspergilli_paper/data/processed/correlation_function/normalized_tpm_all_species.RData")

#-------------------
# All vs All
#-------------------
# Function to compute all-vs-all Spearman correlation with adjusted p-values
compute_spearman_all_vs_all <- function(expr_matrix, species_name) {
  # Transpose: rcorr wants samples as rows, genes as columns
  transposed <- t(expr_matrix)
  
  # Calculate correlation
  cor_result <- rcorr(transposed, type = "spearman")
  
  # Extract correlation and p-values
  corr_mat <- cor_result$r
  pval_mat <- cor_result$P
  
  # Get upper triangle indexes (avoid redundancy)
  upper_idx <- upper.tri(corr_mat, diag = FALSE)
  
  # Create a tidy data frame
  result_df <- tibble(
    gene1 = rownames(corr_mat)[row(corr_mat)[upper_idx]],
    gene2 = colnames(corr_mat)[col(corr_mat)[upper_idx]],
    spearman_rho = corr_mat[upper_idx],
    p_value = pval_mat[upper_idx]
  ) %>%
    mutate(
      padj = p.adjust(p_value, method = "BH"),
      species = species_name
    )
  
  return(result_df)
}

# Apply to all species in normalized_list
all_species_cor <- imap_dfr(normalized_list, compute_spearman_all_vs_all)

save(all_species_cor, file = '/home/laura/Documents/Aspergillus_lncRNA/Aspergilli_paper/data/processed/correlation_function/all_vs_all_spearman.RData')
# #-------------------
# # CIS
# #-------------------
cis_pairs <- read.csv("Aspergilli_paper/data/processed/correlation_function/cis_targets_all_species_filtered_100bp.csv")

# Define safe correlation function
safe_cor_test <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) {
    return(list(cor = NA, pval = NA))
  } else {
    test <- cor.test(x, y, method = "spearman") #PEARSON
    return(list(cor = test$estimate, pval = test$p.value))
  }
}

# Compute all correlations and store unfiltered results
coexp_all <- list()

for (sp in unique(cis_pairs$species)) {
  expr <- normalized_list[[sp]]
  cis_sp <- filter(cis_pairs, species == sp)

  results <- cis_sp %>%
    rowwise() %>%
    mutate(
      tmp = list(safe_cor_test(expr[lncRNA_ID, ], expr[gene_ID, ])),
      cor = tmp$cor,
      pval = tmp$pval
    ) %>%
    ungroup() %>%
    dplyr::select(-tmp) %>%
    mutate(species = sp)
  # Add adjusted p-values (within each species)
  results <- results %>%
    mutate(padj = p.adjust(pval, method = "BH"))  # Benjamini-Hochberg method
  coexp_all[[sp]] <- results
}

# Combine and save raw results
coexp_raw_df <- bind_rows(coexp_all)
saveRDS(coexp_raw_df, "Aspergilli_paper/data/processed/correlation_function/coexpression_raw_spearman\`.rds")
# 
