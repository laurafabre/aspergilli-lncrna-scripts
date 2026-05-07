# =============================================================================
# WGCNA Co-expression Network Analysis
# Aspergillus lncRNA Project
# Disclaimer: This script was modified from https://github.com/Gabaldonlab/lncRNAs

# Description: Builds co-expression networks and modules from RNA-seq TPM data,
#              performs GO, KEGG, and PFAM enrichment per module, and exports
#              results and plots.
#
# Usage: Rscript wgcna_network_analysis.R [spp]
#        e.g. Rscript wgcna_network_analysis.R anid
#        Or source() interactively, setting SPP_TO_RUN in the CONFIG block below.
# =============================================================================

# =============================================================================
# SECTION 1: CONFIG
# Edit all paths, parameters, and species settings here.
# =============================================================================

# --- Paths ---
PATHS <- list(
  working_dir       = "~/Aspergilli_paper/data/processed/wgcna/network_R_data/",
  expression_dir    = "~/Documents/Aspergillus_lncRNA/WGCNA/expression_tables/",
  prediction_dir    = "~/Documents/Aspergillus_lncRNA/WGCNA/../../prediction/",  # for GO/PFAM input files
  output_base       = "~/Aspergilli_paper/data/processed/wgcna/network_R_data/"
)

# --- Species to run (can be overridden by command-line arg) ---
SPP_TO_RUN <- c("anid")   # change to e.g. c("afla","afum","anid","anig") to run multiple

# --- Per-species parameters ---
# Each entry: list(soft_threshold, TPM_filt, percent_sample_filt, MEDissThres, minClusterSize, kegg_code, pfam_file)
SPECIES_CONFIG <- list(
  afla = list(
    soft_threshold      = 10,
    TPM_filt            = 0.1,
    percent_sample_filt = 0.8,
    MEDissThres         = 0.25,
    minClusterSize      = 20,
    kegg_code           = "afv",
    pfam_file           = "afla_iprscan.out"
  ),
  afum = list(
    soft_threshold      = 6,
    TPM_filt            = 0.1,
    percent_sample_filt = 0.8,
    MEDissThres         = 0.25,
    minClusterSize      = 20,
    kegg_code           = "afm",
    pfam_file           = "afum_iprscan.out"
  ),
  anid = list(
    soft_threshold      = 12,
    TPM_filt            = 0.1,
    percent_sample_filt = 0.8,
    MEDissThres         = 0.25,
    minClusterSize      = 20,
    kegg_code           = "ani",
    pfam_file           = "anid_iprscan.out"
  ),
  anig = list(
    soft_threshold      = 7,
    TPM_filt            = 0.1,
    percent_sample_filt = 0.8,
    MEDissThres         = 0.25,
    minClusterSize      = 20,
    kegg_code           = "ang",
    pfam_file           = "anig_iprscan.out"
  )
)

# =============================================================================
# SECTION 2: LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(stringr)
  library(WGCNA)
  library(ggrepel)
  library(GO.db)
  library(clusterProfiler)
  library(patchwork)
  library(ggpubr)
  library(ggh4x)
  library(tidyr)
  library(openxlsx)
  library(colorspace)
  library(gridExtra)
  library(svglite)
  library(ggVennDiagram)
  library(RColorBrewer)
  library(ggridges)
})

disableWGCNAThreads()

# =============================================================================
# SECTION 3: REPRODUCIBILITY
# =============================================================================

set.seed(42)

log_session <- function(spp, out_dir) {
  log_path <- file.path(out_dir, sprintf("%s_session_info.txt", spp))
  sink(log_path)
  cat("=== Session Info ===\n")
  cat(sprintf("Date/Time: %s\n", Sys.time()))
  cat(sprintf("Species: %s\n", spp))
  cat(sprintf("R version: %s\n", R.version.string))
  print(sessionInfo())
  sink()
  message(sprintf("[%s] Session info saved to: %s", spp, log_path))
}

# =============================================================================
# SECTION 4: HELPER FUNCTIONS
# =============================================================================

#' Compute TPM from a raw count/length-combined expression matrix.
#' Assumes rownames encode gene|length (pipe-separated, length in position 2).
compute_TPM <- function(expression_matrix) {
  len_info  <- data.frame(do.call("rbind", strsplit(rownames(expression_matrix), "|", fixed = TRUE)))
  norm_len  <- expression_matrix / as.numeric(len_info[, 2])
  TPM       <- t(t(norm_len) * 1e6 / colSums(norm_len))
  return(TPM)
}

#' Load and compute TPM for a species.
load_expression <- function(spp, expression_dir) {
  path <- file.path(expression_dir, sprintf("%s_expression.txt", spp))
  message(sprintf("[%s] Loading expression from: %s", spp, path))
  raw <- as.matrix(read.table(path, check.names = FALSE, header = TRUE))[, -1]
  TPM <- compute_TPM(raw)
  return(TPM)
}

#' Filter TPM table: keep genes expressed >= TPM_filt in >= percent_sample_filt of samples.
filter_TPM <- function(TPM_table, TPM_filt, percent_sample_filt) {
  keep <- rowSums(TPM_table >= TPM_filt) >= percent_sample_filt * ncol(TPM_table)
  message(sprintf("  Filtering: %d / %d genes pass TPM >= %.2f in >= %.0f%% of samples",
                  sum(keep), nrow(TPM_table), TPM_filt, percent_sample_filt * 100))
  TPM_table[keep, , drop = FALSE]
}

#' Log-transform helper (log2(x + 1)), transposed for WGCNA (samples x genes).
log2_t <- function(mat) t(log2(mat + 1))

#' Classify gene type from gene ID patterns.
classify_gene_type <- function(gene_ids) {
  dplyr::case_when(
    grepl("\\|u", gene_ids) ~ "LincRNA",
    grepl("\\|x", gene_ids) ~ "LncRNA AS",
    TRUE                     ~ "Protein"
  )
}

#' Build output directory for a species, creating subdirs as needed.
make_out_dirs <- function(spp, output_base) {
  dirs <- list(
    base    = file.path(output_base, spp),
    go      = file.path(output_base, spp, "go_terms"),
    network = file.path(output_base, spp)
  )
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
  dirs
}

# =============================================================================
# SECTION 5: ENRICHMENT FUNCTIONS
# =============================================================================

#' Prepare GO universe data frames split by ontology.
prepare_GO_universe <- function(spp, prediction_dir) {
  all_go  <- as.data.frame(GOTERM)
  all_go  <- all_go[!(all_go$go_id %in% c("GO:0003674", "GO:0008150", "GO:0005575")), ]
  goterms <- Term(GOTERM)
  go_names <- cbind(rownames(as.data.frame(goterms)), as.data.frame(goterms))

  yeast_go <- read.table(
    file.path(prediction_dir, spp, sprintf("%s_go.txt", spp))
  )

  list(
    MF       = yeast_go[yeast_go$V1 %in% all_go$go_id[all_go$Ontology == "MF"], ],
    BP       = yeast_go[yeast_go$V1 %in% all_go$go_id[all_go$Ontology == "BP"], ],
    CC       = yeast_go[yeast_go$V1 %in% all_go$go_id[all_go$Ontology == "CC"], ],
    go_names = go_names
  )
}

#' Run GO enrichment for one module and ontology; return enrich result or NULL.
run_GO_enrichment <- function(gene_ids_no_lnc, universe_df, go_names, ontology,
                               module, spp, out_dir) {
  tryCatch({
    ego <- enricher(
      gene_ids_no_lnc,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      universe      = as.character(universe_df$V2),
      minGSSize     = 2,
      maxGSSize     = NA,
      TERM2GENE     = universe_df,
      TERM2NAME     = go_names
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      p <- dotplot(ego, showCategory = 5)
      ggplot2::ggsave(
        file.path(out_dir, sprintf("%s_module_%s_%s.png", spp, module, ontology)),
        p, units = "in", width = 10, height = 7, dpi = 600
      )
      write.table(ego,
        file.path(out_dir, sprintf("%s_module_%s_%s.txt", spp, module, ontology)),
        sep = "\t")
    }
    ego
  }, error = function(e) {
    message(sprintf("  [GO %s] Module %s: %s", ontology, module, e$message))
    NULL
  })
}

#' Run KEGG enrichment for one module.
run_KEGG_enrichment <- function(gene_ids_no_lnc, kegg_code) {
  tryCatch({
    enrichKEGG(gene = gene_ids_no_lnc, organism = kegg_code, pvalueCutoff = 0.05)
  }, error = function(e) {
    message(sprintf("  [KEGG] %s", e$message))
    NULL
  })
}

#' Run PFAM enrichment for one module.
run_PFAM_enrichment <- function(gene_ids_no_lnc, pfam_df) {
  universe_pfam <- unique(pfam_df[pfam_df$V12 != "NULL", c(12, 1)])
  pfam_name_df  <- unique(pfam_df[pfam_df$V12 != "NULL", c(12, 13)])
  tryCatch({
    enricher(
      gene_ids_no_lnc,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      universe      = as.character(universe_pfam$V1),
      minGSSize     = 2,
      maxGSSize     = 10000,
      TERM2GENE     = universe_pfam,
      TERM2NAME     = pfam_name_df
    )
  }, error = function(e) {
    message(sprintf("  [PFAM] %s", e$message))
    NULL
  })
}

#' Safely extract top N descriptions from an enrichResult.
top_descriptions <- function(ego, n = 5) {
  df <- tryCatch(as.data.frame(ego), error = function(e) data.frame())
  if (nrow(df) == 0) return("No enrichment")
  toString(as.character(df[seq_len(min(n, nrow(df))), 2]))
}

# =============================================================================
# SECTION 6: MAIN ANALYSIS FUNCTION
# =============================================================================

#' Full WGCNA pipeline for one species.
#'
#' @param spp          Character. Species code (must match key in SPECIES_CONFIG).
#' @param TPM_table    Matrix. Genes x samples TPM matrix.
#' @param cfg          List.   Species config from SPECIES_CONFIG.
#' @param paths        List.   Global PATHS.
build_network_and_modules <- function(spp, TPM_table, cfg, paths) {

  message(sprintf("\n========== [%s] Starting analysis ==========", spp))

  # --- Output dirs ---
  out_dirs <- make_out_dirs(spp, paths$output_base)

  # --- Log session ---
  log_session(spp, out_dirs$base)

  # --- Filter ---
  TPM_filt <- filter_TPM(TPM_table, cfg$TPM_filt, cfg$percent_sample_filt)

  # --- Build network ---
  message(sprintf("[%s] Computing adjacency (power=%d)...", spp, cfg$soft_threshold))
  adjacency <- adjacency(log2_t(TPM_filt), type = "unsigned", power = cfg$soft_threshold)

  message(sprintf("[%s] Computing TOM...", spp))
  TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
  rownames(TOM) <- colnames(TOM) <- rownames(TPM_filt)

  saveRDS(TOM, file.path(out_dirs$base, sprintf("%s_TOM_2026.rds", spp)))

  dissTOM  <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  # --- Dynamic tree cut ---
  message(sprintf("[%s] Detecting modules (minClusterSize=%d)...", spp, cfg$minClusterSize))
  dynamicMods   <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3,
                                  pamRespectsDendro = FALSE, minClusterSize = cfg$minClusterSize)
  dynamicColors <- labels2colors(dynamicMods)

  # --- Merge similar modules ---
  MEList <- moduleEigengenes(log2_t(TPM_filt), colors = dynamicColors,
                              subHubs = TRUE, excludeGrey = TRUE)
  MEs <- MEList$eigengenes

  merge        <- mergeCloseModules(log2_t(TPM_filt), dynamicColors,
                                     cutHeight = cfg$MEDissThres, verbose = 3)
  mergedColors <- merge$colors
  mergedMEs    <- merge$newMEs

  # --- Plot dendrogram ---
  pdf_path <- file.path(out_dirs$base,
    sprintf("%s_modules_TPM>=%.1f_in_%.0fpct.pdf", spp, cfg$TPM_filt, cfg$percent_sample_filt * 100))
  pdf(pdf_path, width = 12, height = 5)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
    c("Gene\nclustering", "Eigengene\nclustering"),
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    main = sprintf("Cluster Dendrogram (%s)", spp))
  dev.off()

  # --- Finalize module assignments ---
  moduleColors <- mergedColors
  names(moduleColors) <- rownames(TPM_filt)
  MEs <- mergedMEs

  # --- Gene type classification ---
  lncRNA_genes  <- grep("\\|u|\\|x", rownames(TPM_filt), value = TRUE)
  protein_genes <- setdiff(rownames(TPM_filt), lncRNA_genes)

  module_df <- data.frame(
    gene   = rownames(TPM_filt),
    module = moduleColors,
    type   = ifelse(rownames(TPM_filt) %in% lncRNA_genes, "lncRNA", "protein")
  )

  # --- Save network objects ---
  colorOrder   <- c("grey", standardColors(50))
  moduleLabels <- match(moduleColors, colorOrder) - 1

  save(MEs, moduleLabels, moduleColors, module_df, TPM_filt,
       file = file.path(out_dirs$base, sprintf("%s_network_2026.RData", spp)))
  write.csv(module_df, file.path(out_dirs$base, "module_labels_df.csv"))

  message(sprintf("[%s] Network saved. Modules detected: %d (excl. grey)",
                  spp, length(setdiff(unique(moduleColors), "grey"))))

  # --- Hub genes ---
  top_hubs <- chooseTopHubInEachModule(log2_t(TPM_filt), colorh = moduleColors,
                                        power = cfg$soft_threshold, type = "signed hybrid")
  top_hubs_df <- data.frame(module = names(top_hubs), gene = top_hubs)
  write.table(top_hubs_df,
    file.path(out_dirs$base, sprintf("%s_top_hubs.txt", spp)),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  # --- Intramodular connectivity ---
  iK <- intramodularConnectivity.fromExpr(log2_t(TPM_filt),
                                           power = cfg$soft_threshold, colors = moduleColors)
  rownames(iK) <- rownames(TPM_filt)
  iK$Type <- factor(classify_gene_type(rownames(iK)),
                    levels = c("LincRNA", "LncRNA AS", "Protein"))
  iK$gene <- rownames(iK)

  write.table(iK,
    file.path(out_dirs$base, sprintf("%s_connectivity.txt", spp)),
    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  p_conn <- ggplot(iK, aes(x = Type, y = log2(kWithin + 0.01), fill = Type)) +
    geom_boxplot() +
    scale_fill_manual(values = c("darkorange", "darkgreen", "darkblue")) +
    labs(y = "Log2(intramodular connectivity)", title = spp,
         fill = "Class code") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 22), axis.title = element_text(size = 22),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 22), legend.title = element_text(size = 22))

  ggplot2::ggsave(
    file.path(out_dirs$base, sprintf("%s_connectivity_TPM>=%.1f_%.0fpct.png",
                                      spp, cfg$TPM_filt, cfg$percent_sample_filt * 100)),
    p_conn, width = 10, height = 7, dpi = 300)

  # --- Prepare enrichment data ---
  message(sprintf("[%s] Loading GO/PFAM data...", spp))
  go_universe <- prepare_GO_universe(spp, paths$prediction_dir)

  pfam_raw <- read.table(
    file.path(paths$prediction_dir, spp, cfg$pfam_file),
    header = FALSE, fill = TRUE, sep = "\t", quote = "\""
  )

  # --- lncRNA module summary ---
  all_lncRNAs_module <- NULL

  # --- Per-module enrichment loop ---
  go_module_table <- NULL

  for (module in setdiff(unique(moduleColors), "grey")) {

    gene_ids      <- names(moduleColors)[moduleColors == module]
    gene_ids_prot <- sapply(strsplit(gene_ids[!grepl("^MSTRG", gene_ids)], "\\|"), "[", 1)
    gene_ids_lnc  <- gene_ids[grepl("^MSTRG", gene_ids)]

    lncRNA_str <- if (length(gene_ids_lnc) > 0) toString(gene_ids_lnc) else "NA"

    if (length(gene_ids_lnc) > 0) {
      all_lncRNAs_module <- rbind(all_lncRNAs_module,
                                   data.frame(gene = gene_ids_lnc, module = module))
    }

    # GO
    bp_desc <- "No enrichment"
    for (ont in c("MF", "CC", "BP")) {
      ego <- run_GO_enrichment(gene_ids_prot, go_universe[[ont]], go_universe$go_names,
                                ont, module, spp, out_dirs$go)
      if (ont == "BP") bp_desc <- top_descriptions(ego)
    }

    # KEGG
    kegg_desc <- top_descriptions(run_KEGG_enrichment(gene_ids_prot, cfg$kegg_code))

    # PFAM
    pfam_desc <- top_descriptions(run_PFAM_enrichment(gene_ids_prot, pfam_raw))

    go_module_table <- rbind(go_module_table, data.frame(
      Module   = module,
      GO_BP    = bp_desc,
      PFAM     = pfam_desc,
      KEGG     = kegg_desc,
      Gene_ids = toString(gene_ids_prot),
      lncRNA_ids = lncRNA_str,
      stringsAsFactors = FALSE
    ))
  }

  # --- lncRNA-module file ---
  if (!is.null(all_lncRNAs_module)) {
    write.table(all_lncRNAs_module,
      file.path(out_dirs$base, sprintf("%s_all_lncRNAs_in_modules_id_module.txt", spp)),
      quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }

  # --- Module stats ---
  module_stats <- do.call(rbind, lapply(setdiff(unique(moduleColors), "grey"), function(module) {
    gene_ids <- names(moduleColors)[moduleColors == module]
    n_lnc    <- length(grep("MSTRG", gene_ids))
    data.frame(
      Module       = module,
      N_genes      = length(gene_ids),
      N_lncRNAs    = n_lnc,
      prcnt_lncRNAs = round(n_lnc / length(gene_ids) * 100, 2),
      Spp          = spp,
      stringsAsFactors = FALSE
    )
  }))

  module_stats_final <- merge(module_stats, go_module_table, by = "Module")
  module_stats_final <- module_stats_final[order(module_stats_final$N_genes), ]
  module_stats_final$Module <- factor(module_stats_final$Module,
                                       levels = module_stats_final$Module)

  write.csv(module_stats_final,
    file.path(out_dirs$base, sprintf("%s_module_stats_final.csv", spp)),
    row.names = FALSE)

  # --- Module size barplot ---
  p_bar <- ggplot(module_stats_final,
    aes(x = Module, y = N_genes, fill = Module,
        label = paste0(N_lncRNAs, " (", prcnt_lncRNAs, "%)"))) +
    geom_bar(stat = "identity", width = 0.9) +
    scale_fill_manual(values = as.character(module_stats_final$Module)) +
    geom_text(hjust = 0, size = 5) +
    coord_flip() +
    scale_y_continuous(
      breaks = seq(0, max(module_stats_final$N_genes) + 300, 200),
      limits = c(0, max(module_stats_final$N_genes) + 300)
    ) +
    labs(y = "Number of genes in module", x = sprintf("Modules (%s)", spp)) +
    theme_bw() +
    theme(legend.position = "none")

  ggplot2::ggsave(
    file.path(out_dirs$base, sprintf("%s_N_genes_in_modules.pdf", spp)),
    p_bar, width = 12, height = 6, dpi = 300)

  message(sprintf("[%s] Done. Results in: %s", spp, out_dirs$base))

  # Return key objects invisibly for interactive inspection
  invisible(list(
    moduleColors      = moduleColors,
    MEs               = MEs,
    module_stats      = module_stats_final,
    iK                = iK,
    TPM_filt          = TPM_filt
  ))
}

# =============================================================================
# SECTION 7: RUN
# =============================================================================

# Allow overriding species via command-line: Rscript script.R anid
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) SPP_TO_RUN <- args

setwd(PATHS$working_dir)

results <- list()

for (spp in SPP_TO_RUN) {

  if (!spp %in% names(SPECIES_CONFIG)) {
    warning(sprintf("Species '%s' not found in SPECIES_CONFIG. Skipping.", spp))
    next
  }

  cfg        <- SPECIES_CONFIG[[spp]]
  TPM_table  <- load_expression(spp, PATHS$expression_dir)

  results[[spp]] <- build_network_and_modules(
    spp       = spp,
    TPM_table = TPM_table,
    cfg       = cfg,
    paths     = PATHS
  )
}

message("\n===== All species complete =====")
