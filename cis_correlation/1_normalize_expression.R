library(tidyverse)

setwd("~/Documents/Aspergillus_lncRNA/Coefficient_variation")
# analysis of the coefficient variation , conservation and evolution
# input list:
# expression counts table:

afum_exp <- as.matrix(read.table("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/afum_expression.txt",check.names = FALSE, header = TRUE, ))[,-1]
afla_exp <- as.matrix(read.table("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/afla_expression.txt",check.names = FALSE, header = TRUE, ))[,-1]
anid_exp <- as.matrix(read.table("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/anid_expression.txt",check.names = FALSE, header = TRUE, ))[,-1]
anig_exp <- as.matrix(read.table("~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/anig_expression.txt",check.names = FALSE, header = TRUE, ))[,-1]

# Function to normalize the expression matrix and log-transform it to TPM
normalize_expression <- function(expression) {
  # Extract gene lengths (assuming 'length' data frame is provided)
  length <- data.frame(do.call('rbind', strsplit(as.character(row.names(expression)), "|", fixed = TRUE)))
  
  # Normalize the expression values by gene length
  norm_exp_len <- expression / as.numeric(length[, 2])  # Dividing by gene length
  
  # Convert to TPM (Transcripts Per Million)
  TPM_initial <- t(t(norm_exp_len) * 1e6 / colSums(norm_exp_len))  # Normalize across samples
  
  # Log2-transform the TPM values
  log2_TPM <- log2(TPM_initial + 1)  # Adding 1 to avoid log(0)
  
  return(log2_TPM)
}

# Normalize the expression data
afum_normalized <- normalize_expression(afum_exp)
afla_normalized <- normalize_expression(afla_exp)
anid_normalized <- normalize_expression(anid_exp)
anig_normalized <- normalize_expression(anig_exp)

normalized_list <- list(afum =afum_normalized,
                        afla= afla_normalized,
                        anid =anid_normalized,
                        anig =anig_normalized
)
save(normalized_list, file = "../Aspergilli_paper/data/processed/correlation_function/normalized_tpm_all_species.RData")
