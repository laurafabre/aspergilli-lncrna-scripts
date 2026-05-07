# =============================================================================
# WGCNA lncRNA Module Summary, Enrichment Tables & Heatmap
# Aspergillus lncRNA Project
#
# Description: Loads per-species WGCNA module stats, computes lncRNA functional
#              summaries, exports multi-tab Excel tables, and generates a
#              lncRNA regulatory investment heatmap.
#
# Usage: Rscript wgcna_lncrna_summary_tables.R
#        Or source() interactively.
# =============================================================================

# =============================================================================
# SECTION 1: CONFIG
# =============================================================================

PATHS <- list(
  working_dir  = "~/data_drive/Aspergilli_paper",
  base_dir     = "~/data_drive/Aspergilli_paper/data/processed/wgcna/network_R_data",
  output_dir   = "results/tables"
)

SPECIES_CONFIG <- list(
  afla = list(name = "A. flavus",    total_lncRNAs = 2741, role = "Crop Pest"),
  afum = list(name = "A. fumigatus", total_lncRNAs = 1165, role = "Pathogen"),
  anid = list(name = "A. nidulans",  total_lncRNAs = 2291, role = "Model"),
  anig = list(name = "A. niger",     total_lncRNAs = 2356, role = "Industry")
)

SPECIES_LIST <- names(SPECIES_CONFIG)

# Module filters for "high-confidence" lncRNA-dense modules
HC_MODULE_MAX_SIZE    <- 200   # max module size OR ...
HC_MODULE_MIN_PCT_LNC <- 10    # ... min % lncRNA to be included

HEATMAP_LNC_PCT_THRESHOLD <- 15  # for the simple heatmap (prcnt_lncRNAs > X)

# =============================================================================
# SECTION 2: LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tools)
  library(openxlsx)
  library(pheatmap)
  library(RColorBrewer)
  library(gt)
})

# =============================================================================
# SECTION 3: REPRODUCIBILITY
# =============================================================================

set.seed(42)

log_session <- function(out_dir) {
  log_path <- file.path(out_dir, "session_info_summary_tables.txt")
  sink(log_path)
  cat("=== Session Info ===\n")
  cat(sprintf("Date/Time : %s\n", Sys.time()))
  cat(sprintf("R version : %s\n", R.version.string))
  print(sessionInfo())
  sink()
  message(sprintf("Session info saved to: %s", log_path))
}

# =============================================================================
# SECTION 4: HELPER — species lookup tables
# =============================================================================

total_lncRNAs_per_species <- data.frame(
  Species      = SPECIES_LIST,
  Total_lncRNAs = sapply(SPECIES_CONFIG, `[[`, "total_lncRNAs"),
  Species_name  = sapply(SPECIES_CONFIG, `[[`, "name"),
  Role          = sapply(SPECIES_CONFIG, `[[`, "role"),
  stringsAsFactors = FALSE
)

species_name_map <- setNames(
  sapply(SPECIES_CONFIG, `[[`, "name"),
  SPECIES_LIST
)

# =============================================================================
# SECTION 5: DATA LOADING FUNCTIONS
# =============================================================================

#' Read module stats CSV for one species and standardise column names.
read_module_stats <- function(spp, base_dir) {
  path <- file.path(base_dir, spp, sprintf("%s_module_stats_sort_final.csv", spp))
  df   <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- make.names(colnames(df))
  
  df %>%
    select(
      Module, N_genes, N_lncRNAs, prcnt_lncRNAs,
      Putative_Function = any_of(c("GO.terms",   "GO_terms")),
      lncRNA_full_id    = any_of(c("lncRNA.ids", "lncRNA_ids"))
    ) %>%
    mutate(
      Species           = spp,
      Putative_Function = if_else(
        is.na(Putative_Function) | Putative_Function %in% c("", "NA"),
        "No enrichment", Putative_Function
      )
    )
}

#' Load all species module stats into one data frame.
load_all_module_stats <- function(species_list, base_dir) {
  bind_rows(lapply(species_list, read_module_stats, base_dir = base_dir))
}

# =============================================================================
# SECTION 6: lncRNA-LEVEL PROCESSING
# =============================================================================

#' Expand module stats to one row per lncRNA, with enrichment flags.
expand_module_lncRNAs <- function(all_modules) {
  all_modules %>%
    mutate(
      Has_GO          = Putative_Function != "No enrichment",
      Has_Enrichment  = Has_GO
    ) %>%
    filter(!is.na(lncRNA_full_id) & lncRNA_full_id != "" & lncRNA_full_id != "NA") %>%
    separate_rows(lncRNA_full_id, sep = ",\\s*") %>%
    mutate(lncRNA_id = sub("\\|.*", "", lncRNA_full_id)) %>%
    select(Species, Module, lncRNA_id, Putative_Function,
           N_genes, N_lncRNAs, prcnt_lncRNAs, Has_GO, Has_Enrichment)
}

# =============================================================================
# SECTION 7: FUNCTIONAL CATEGORISATION
# =============================================================================

#' Assign a broad functional category from a GO term string.
#' Priority order matters: first match wins.
categorise_function <- function(function_string) {
  if (is.na(function_string) || function_string %in% c("", "No enrichment")) {
    return("Unannotated")
  }
  terms <- tolower(function_string)
  
  case_when_vec <- list(
    "Pathogenesis & Nutrient Scavenging" =
      "iron|siderophore|fusarinine|triacetylfusarinine|copper|zinc|homeostasis|virulence|nutritional immunity|host",
    "Specialized Toxins & SM" =
      "aflatoxin|sterigmatocystin|mycotoxin|toxin|melanin|pks|nrps|bgc|secondary metabolite",
    "Adaptive Metabolic Versatility" =
      "hemicellulose|oxaloacetate|singlet oxygen|sequestering|metabolic versatility",
    "Industrial & Carbon Metabolism" =
      "carbohydrate|glycolytic|glucan|xylan|galactose|pectin|cellulose|citric|acid|protease|amylase|secretion",
    "Growth & Development" =
      "conidium|hyphal|spore|cleistothecia|mating|sexual|asexual|development|spindle|morphogenesis|cell cycle|division|actin|microtubule",
    "Stress & Resilience" =
      "stress|heat|oxid|peroxide|thermal|unfolded protein|eisosome|cell wall|chitin",
    "Core Housekeeping" =
      "ribosom|translation|mitochondrial|atp|dna repair|replication"
  )
  
  for (category in names(case_when_vec)) {
    if (grepl(case_when_vec[[category]], terms)) return(category)
  }
  "Other Biological Processes"
}

#' Add Broad_Category column to a data frame that has a Putative_Function column.
add_broad_category <- function(df) {
  df %>% mutate(Broad_Category = sapply(Putative_Function, categorise_function))
}

# =============================================================================
# SECTION 8: SUMMARY STATISTICS
# =============================================================================

compute_analysis_summary <- function(module_lncRNAs, total_lncRNAs_per_species) {
  module_lncRNAs %>%
    group_by(Species) %>%
    summarise(
      Module_lncRNAs  = n_distinct(lncRNA_id),
      With_GO_info    = n_distinct(lncRNA_id[Has_GO]),
      With_Enrichment = n_distinct(lncRNA_id[Has_Enrichment]),
      .groups = "drop"
    ) %>%
    left_join(total_lncRNAs_per_species, by = "Species") %>%
    mutate(
      Pct_in_modules               = round(Module_lncRNAs  / Total_lncRNAs * 100, 1),
      Pct_with_GO                  = round(With_GO_info    / Total_lncRNAs * 100, 1),
      Pct_with_enrichment          = round(With_Enrichment / Total_lncRNAs * 100, 1),
      Pct_of_module_with_enrichment = round(With_Enrichment / Module_lncRNAs * 100, 1)
    )
}

# =============================================================================
# SECTION 9: TABLE BUILDERS
# =============================================================================

build_table1_summary <- function(analysis_summary) {
  analysis_summary %>%
    select(Species_name, Total_lncRNAs, Module_lncRNAs, Pct_in_modules,
           With_Enrichment, Pct_with_enrichment) %>%
    rename(
      Species                  = Species_name,
      `Total lncRNAs`          = Total_lncRNAs,
      `In modules (N)`         = Module_lncRNAs,
      `In modules (%)`         = Pct_in_modules,
      `With GO enrichment (N)` = With_Enrichment,
      `With GO enrichment (%)` = Pct_with_enrichment
    ) %>%
    bind_rows(tibble(
      Species                  = "TOTAL / AVERAGE",
      `Total lncRNAs`          = sum(analysis_summary$Total_lncRNAs),
      `In modules (N)`         = sum(analysis_summary$Module_lncRNAs),
      `In modules (%)`         = round(sum(analysis_summary$Module_lncRNAs) /
                                         sum(analysis_summary$Total_lncRNAs) * 100, 1),
      `With GO enrichment (N)` = sum(analysis_summary$With_Enrichment),
      `With GO enrichment (%)` = round(sum(analysis_summary$With_Enrichment) /
                                         sum(analysis_summary$Module_lncRNAs) * 100, 1)
    ))
}

build_table2_high_conf <- function(module_lncRNAs, species_name_map,
                                   max_size, min_pct_lnc) {
  module_lncRNAs %>%
    filter(Has_Enrichment) %>%
    group_by(Species, Module) %>%
    summarise(
      Module_Size     = dplyr::first(N_genes),
      LncRNAs         = dplyr::first(N_lncRNAs),
      Pct_lncRNA      = dplyr::first(prcnt_lncRNAs),
      Enriched_Functions = dplyr::first(Putative_Function),
      .groups = "drop"
    ) %>%
    filter(Module_Size < max_size | Pct_lncRNA > min_pct_lnc) %>%
    mutate(Species = species_name_map[Species]) %>%
    arrange(Species, desc(Pct_lncRNA)) %>%
    rename(
      `Module size`       = Module_Size,
      `% lncRNA`          = Pct_lncRNA,
      `Enriched functions` = Enriched_Functions
    )
}

build_table3_secondary_metabolism <- function(module_lncRNAs, species_name_map,
                                              max_size, min_pct_lnc) {
  sm_pattern <- paste(c(
    "secondary metabolite", "aflatoxin", "sterigmatocystin", "melanin",
    "asperthecin", "terpenoid", "siderophore", "endocrocin", "fumigaclavine", "chitin"
  ), collapse = "|")
  
  high_conf <- module_lncRNAs %>%
    filter(Has_Enrichment) %>%
    group_by(Species, Module) %>%
    summarise(
      Module_Size = dplyr::first(N_genes),
      pct_lnc     = dplyr::first(prcnt_lncRNAs),
      Functions   = dplyr::first(Putative_Function),
      .groups = "drop"
    ) %>%
    filter((Module_Size < max_size | pct_lnc > min_pct_lnc) &
             grepl(sm_pattern, Functions, ignore.case = TRUE))
  
  high_conf %>%
    mutate(
      Category = case_when(
        str_detect(Functions, "(?i)aflatoxin|sterigmatocystin|asperthecin|secondary|biosynthetic|microperfuranone|fumigaclavine") ~
          "I. Specialized Metabolism (Toxins)",
        str_detect(Functions, "(?i)siderophore|iron|triacetylfusarinine") ~
          "II. Iron Homeostasis & Scavenging",
        str_detect(Functions, "(?i)melanin|chitin|pigment|spore wall") ~
          "III. Physical Defense & Pigmentation",
        TRUE ~ "IV. Other Specialized Processes"
      ),
      Species = species_name_map[Species]
    ) %>%
    select(Category, Species, Module,
           `Size (Genes)` = Module_Size,
           `% lncRNA`     = pct_lnc,
           `Functional Enrichment` = Functions) %>%
    arrange(Category, desc(`% lncRNA`))
}

build_table4_aflatoxin <- function(module_lncRNAs) {
  module_lncRNAs %>%
    filter(Species == "afla", Has_Enrichment,
           grepl("aflatoxin", Putative_Function, ignore.case = TRUE)) %>%
    group_by(Module) %>%
    summarise(
      `Module size`       = dplyr::first(N_genes),
      `LncRNAs (N)`       = dplyr::first(N_lncRNAs),
      `% lncRNA`          = dplyr::first(prcnt_lncRNAs),
      `LncRNA IDs`        = paste(unique(lncRNA_id), collapse = "; "),
      `Enriched functions` = dplyr::first(Putative_Function),
      .groups = "drop"
    )
}

build_table5_all_modules <- function(all_modules) {
  all_modules %>%
    select(Species, Module, N_genes, N_lncRNAs, prcnt_lncRNAs, Putative_Function)
}

build_table6_all_lncrnas <- function(module_lncRNAs, all_lnc_categorized,
                                     species_name_map) {
  module_lncRNAs %>%
    left_join(
      all_lnc_categorized %>% select(lncRNA_id, Broad_Category) %>% distinct(),
      by = "lncRNA_id"
    ) %>%
    mutate(
      Species           = species_name_map[Species],
      Annotation_status = case_when(
        !Has_Enrichment & is.na(Putative_Function)              ~ "No GO terms",
        !Has_Enrichment & Putative_Function == "No enrichment"  ~ "No enrichment",
        TRUE                                                     ~ "Annotated"
      )
    ) %>%
    select(
      Species,
      lncRNA_id,
      Module,
      `Annotation status`    = Annotation_status,
      `Enriched GO terms`    = Putative_Function,
      `Module genes`         = N_genes,
      `Module lncRNAs`       = N_lncRNAs,
      `% lncRNA in module`   = prcnt_lncRNAs
    ) %>%
    distinct() %>%
    arrange(Species, lncRNA_id)
}

# =============================================================================
# SECTION 10: EXCEL EXPORT
# =============================================================================

#' Write a data frame to a worksheet with a bold header row.
add_bold_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data)
  addStyle(wb, sheet_name,
           createStyle(textDecoration = "bold"),
           rows = 1, cols = seq_len(ncol(data)))
}

export_excel <- function(tables, out_path) {
  wb <- createWorkbook()
  
  for (nm in names(tables)) {
    add_bold_sheet(wb, nm, tables[[nm]])
  }
  
  readme <- data.frame(Description = c(
    "Table1_Summary          : Overall statistics by species",
    "Table2_HighConfModules  : Modules with Size < 200 OR % lncRNA > 10",
    "Table3_SecondaryMetab   : Modules with secondary metabolism GO terms",
    "Table4_Flavus_Aflatoxin : Aflatoxin-associated modules in A. flavus",
    "Table5_AllModules       : Complete list of all modules, all species",
    "Table6_AllLncRNAs       : All lncRNAs with module and annotation info"
  ))
  addWorksheet(wb, "README")
  writeData(wb, "README", readme)
  addStyle(wb, "README", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("Excel saved: %s", out_path))
}

# =============================================================================
# SECTION 11: HEATMAP
# =============================================================================

#' Compute per-species lncRNA regulatory investment (weighted by module size)
#' for each functional category.
build_heatmap_matrix <- function(all_modules, species_name_map) {
  all_modules %>%
    add_broad_category() %>%
    filter(!Broad_Category %in% c("Unannotated", "Other Biological Processes")) %>%
    group_by(Species, Broad_Category) %>%
    summarise(
      Investment = (sum(N_lncRNAs) / sum(N_genes)) * 100,
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = Species, values_from = Investment, values_fill = 0) %>%
    column_to_rownames("Broad_Category") %>%
    rename_with(~ species_name_map[.x]) %>%
    as.matrix()
}

plot_investment_heatmap <- function(heatmap_mat, col_anno, ann_colors, out_path) {
  pheatmap(
    heatmap_mat,
    main            = "lncRNA Regulatory Investment Across Species Niches",
    annotation_col  = col_anno,
    annotation_colors = ann_colors,
    color           = colorRampPalette(
      c("#f7fbff", "#deebf7", "#c6dbef", "#4292c6", "#084594"))(100),
    display_numbers = TRUE,
    number_format   = "%.1f",
    fontsize_number = 9,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    angle_col       = 45,
    border_color    = "white",
    cellwidth       = 50,
    cellheight      = 30,
    filename        = out_path
  )
  message(sprintf("Heatmap saved: %s", out_path))
}

# =============================================================================
# SECTION 12: RUN
# =============================================================================

setwd(PATHS$working_dir)
dir.create(PATHS$output_dir, showWarnings = FALSE, recursive = TRUE)
log_session(PATHS$output_dir)

# --- Load data ---
message("Loading module stats...")
all_modules    <- load_all_module_stats(SPECIES_LIST, PATHS$base_dir)
module_lncRNAs <- expand_module_lncRNAs(all_modules)

# --- Summary stats ---
analysis_summary <- compute_analysis_summary(module_lncRNAs, total_lncRNAs_per_species)
message("=== Summary ===")
print(analysis_summary)

# --- Categorised lncRNAs (for Table 6 join) ---
all_lnc_categorized <- module_lncRNAs %>% add_broad_category()

# --- Build tables ---
message("Building tables...")
tables <- list(
  Table1_Summary         = build_table1_summary(analysis_summary),
  Table2_HighConfModules = build_table2_high_conf(
    module_lncRNAs, species_name_map, HC_MODULE_MAX_SIZE, HC_MODULE_MIN_PCT_LNC),
  Table3_SecondaryMetab  = build_table3_secondary_metabolism(
    module_lncRNAs, species_name_map, HC_MODULE_MAX_SIZE, HC_MODULE_MIN_PCT_LNC),
  Table4_Flavus_Aflatoxin = build_table4_aflatoxin(module_lncRNAs),
  Table5_AllModules       = build_table5_all_modules(all_modules),
  Table6_AllLncRNAs       = build_table6_all_lncrnas(
    module_lncRNAs, all_lnc_categorized, species_name_map)
)

# --- Export Excel ---
export_excel(
  tables,
  file.path(PATHS$output_dir, "Suppl_Tables_GOTerm_enrichment_wgcna.xlsx")
)

# --- Heatmap ---
message("Building heatmap...")
heatmap_mat <- build_heatmap_matrix(all_modules, species_name_map)

col_anno <- data.frame(
  Role = sapply(SPECIES_CONFIG, `[[`, "role"),
  row.names = sapply(SPECIES_CONFIG, `[[`, "name")
)[colnames(heatmap_mat), , drop = FALSE]

ann_colors <- list(
  Role = c(
    "Crop Pest" = "#D95F02",
    "Pathogen"  = "#E7298A",
    "Model"     = "#7570B3",
    "Industry"  = "#1B9E77"
  )
)

plot_investment_heatmap(
  heatmap_mat, col_anno, ann_colors, file.path(PATHS$output_dir, "teste_Fig_lncRNA_investment_heatmap.pdf")
)

# Save heatmap backbone table for supplementary
heatmap_backbone <- all_modules %>%
  add_broad_category() %>%
  filter(!Broad_Category %in% c("Unannotated", "Other Biological Processes")) %>%
  select(Species, Broad_Category, Module, N_genes, N_lncRNAs, prcnt_lncRNAs, Putative_Function) %>%
  arrange(Broad_Category, Species, desc(prcnt_lncRNAs))

write.csv(
  heatmap_backbone,
  file.path(PATHS$output_dir, "Suppl_table_WGCNA_Heatmap_Module_Categorization.csv"),
  row.names = FALSE
)

message("\n===== All done =====")
