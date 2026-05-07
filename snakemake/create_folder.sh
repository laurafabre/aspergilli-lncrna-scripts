#!/bin/bash

cd /vast/sa58ged/snakemake_afla/snakemake_afla
# Configuration
BASE_DIR="/vast/sa58ged/snakemake_afla/snakemake_afla"

# Define directories
RAW_READS_DIR="$BASE_DIR/rawReads"
TRIMMED_READS_DIR="$BASE_DIR/trimmedReads"
TRIMMED_QC_DIR="$BASE_DIR/trimmedQC"
MULTIQC="$BASE_DIR/trimmedQC/multiqc"
STAR_INDEX_DIR="$BASE_DIR/starindex"
STAR_ALIGNED_DIR="$BASE_DIR/starAligned"
STRINGTIE_DIR="$BASE_DIR/stringtieAssembly"

# Create directories if they don't exist
mkdir -p $TRIMMED_READS_DIR $STAR_INDEX_DIR $STAR_ALIGNED_DIR $STRINGTIE_DIR $MULTIQC

#