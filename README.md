# Aspergilli lncRNA Study — Scripts Compendium

> ⚠️ **Note:** This repository is a compendium of scripts used in our study and is **not designed as a plug-and-play pipeline**. Scripts are provided for transparency and reproducibility purposes, but running them end-to-end will require adapting paths, environments, and inputs to your own setup.

---

## Overview

This repository contains the scripts used for the analysis of long non-coding RNAs (lncRNAs) in *Aspergillus* species. The workflow is split across two GitHub repositories:
---

## Repository Structure

```
.
├── prediction/          # lncRNA prediction pipeline
├── snakemake/           # RNA-Seq mapping & transcriptome assembly
├── cis_correlation/     # Cis-regulatory co-expression analysis
├── synteny_conservation/# lncRNA family classification by synteny
├── TE/                  # Repeat/transposable element analysis
└── wgcna/               # Weighted gene co-expression network analysis
```

---

## Module Descriptions

### `prediction/`
Scripts for lncRNA prediction, largely adapted from [Gabaldonlab/lncRNAs](https://github.com/Gabaldonlab/lncRNAs). Includes subtelomeric enrichment analysis.

Key scripts:
- `lncrna_prediction.py` — Main prediction script
- `filter_lncrna.sbatch` — SLURM batch job for filtering
- `non_cod.py`, `select_longer_200.py`, `select_only_prot_cod.py` — Filtering utilities
- `subtelomeres/` — Subtelomeric enrichment analysis

Environment: `lncrna_prediction.yml`

---

### `snakemake/`
Snakemake workflow for QC, trimming and mapping RNA-Seq reads to reference genomes and performing genome-guided transcriptome assembly (among other analyses).

Key files:
- `snakefile_paired.py` — Main Snakemake workflow (paired-end reads)
- `config.yaml` — Configuration file (paths, samples, parameters)
- `create_folder.sh` — Helper to set up directory structure

Environment: `snakemake.yml`

---

### `cis_correlation/`
R scripts for normalizing expression data and predicting cis-regulatory targets based on co-expression.

Scripts (run in order):
1. `1_normalize_expression.R`
2. `2_predict_cis_targets.R`
3. `3_compute_coexpression.R`
4. `4_filter_coexpression.R`
5. `5_combine_cis_coexp.R`

See `session_info.txt` for R package versions.

---

### `synteny_conservation/`
Scripts and data for classifying lncRNAs into families based on synteny conservation across species.

Scripts (run in order):
1. `00_genespace.R` — GENESPACE synteny analysis
2. `01_parse_orthologs.R` — Parse orthologs
3. `02_run_synteny_analysis.py` — Synteny analysis
4. `03_conservation_analysis.R` — Conservation scoring
5. `04_get_family_anchors.py` — Identify family anchors
6. `05_get_families_proteins_id.py` — Map to protein IDs

Includes `Pegueroles_etal_2019/` — modified scripts from [Pegueroles et al. 2019](https://doi.org/10.1093/molbev/msz108) adapted for fungal species.

Example input/output files are provided in `examples/`.

See `session_info.txt` for R package versions.

---

### `TE/`
Scripts for repeat calling with [EarlGrey](https://github.com/TobyBaril/EarlGrey) and quantifying the overlap between lncRNAs and transposable elements.

Key scripts:
- `earlgrey/00_earlgrey_mycomobilome.sbatch` — EarlGrey SLURM job
- `00_TE_density.sh` — Compute TE density
- `01_parallel_innnovation_analysis.sh` — Parallel innovation analysis
- `02_blastx.sh` — BLASTx for TE classification

Example input/output files are provided in `examples/`.

Environment: `earlgrey/earlgrey.yml`

---

### `wgcna/`
R scripts for weighted gene co-expression network analysis (WGCNA), including lncRNA summary tables and network topology and preservation analyses.

Scripts (run in order):
1. `01_wgcna_network_analysis.R`
2. `02_wgcna_lncrna_summary_tables.R`
3. `03_wgcna_topology_2026.R`
4. `04_wgcna_preservation_biofilm.R`

Example expression, GO annotation, and InterProScan files are provided in `examples/`.

See `session_info.txt` and `session_info_summary_tables.txt` for R package versions.

---

## Dependencies

Each module has its own environment or `session_info.txt`. There is no single unified environment for the full workflow.

| Module | Environment file |
|--------|-----------------|
| `prediction/` | `prediction/lncrna_prediction.yml` |
| `snakemake/` | `snakemake/snakemake.yml` |
| `TE/earlgrey/` | `TE/earlgrey/earlgrey.yml` |
| `cis_correlation/` | `cis_correlation/session_info.txt` |
| `synteny_conservation/` | `synteny_conservation/session_info.txt` |
| `wgcna/` | `wgcna/session_info.txt` |

---

## Citation

If you use these scripts, please cite our paper:

> *[Citation to be added upon publication]*

For scripts adapted from Gabaldonlab:
> Hovhannisyan H, Gabaldón T. (2021) *Nat Commun. *https://doi.org/10.1038/s41467-021-27635-4

For synteny scripts adapted from Pegueroles et al.:
> Pegueroles et al. (2019) *Mol. Biol. Evol.* https://doi.org/10.1093/molbev/msz108

---

## Contact

For questions about specific scripts, please open an issue in this repository.
