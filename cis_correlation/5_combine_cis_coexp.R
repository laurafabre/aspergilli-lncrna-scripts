# 4_combine_cis_coexp.R
library(dplyr)

cis_df <- read.csv("../Aspergilli_paper/data/processed/correlation_function/cis_targets_all_species_filtered_100bp.csv")
coexp_df <- read.csv("../Aspergilli_paper/data/processed/correlation_function/coexpressed_pairs_corr80_padj0.01.csv")

merged <- inner_join(cis_df, coexp_df, by = c("lncRNA_ID", "gene_ID", "species"))
write.csv(merged, "final_cis_coexpressed_pairscorr80_padj0.01.csv", row.names = FALSE)
