#!/bin/bash
# Phase 2.5: Parallel Innovation Analysis (Strict TE Filter)

cd /vast/sa58ged/TE_analysis

SPECIES=("afla" "afum" "anid" "anig")

# Define the noise filter: removes Simple Repeats, Microsatellites, Low Complexity, and RNAs
# We do NOT include "Unclassified" here so they remain in the analysis.
NOISE_FILTER="Simple_repeat|Low_complexity|Microsatellite|RNA|Satellite"

for SP in "${SPECIES[@]}"; do
    echo "Categorizing TRUE TE innovation modes for $SP..."

    # # 1. Capture "Promoter-Donation" (TSS overlap)
    # # We intersect, then immediately filter out the "Other" categories
    # bedtools intersect -a results/${SP}_TSS_50bp.bed \
    # -b /vast/sa58ged/earlgrey_all/curate_mycoMobilome/${SP}_toby/${SP}_EarlGrey/${SP}_summaryFiles/${SP}.filteredRepeats.bed \
    # -wa -wb | grep -vE "$NOISE_FILTER" > results/${SP}_mode_Promoter_Donation.txt

    # 2. Capture "Structural Exonization" (Internal Exons)
    # First, get the raw overlaps excluding the TSS hits
    # Then, filter out the "Other" categories from column 10
    bedtools intersect -a results/${SP}_exon_te_raw.txt -b results/${SP}_TSS_te_intersect.txt -v | \
    grep -vE "$NOISE_FILTER" > results/${SP}_mode_Exonization_Internal.txt

    # 3. Stats Collection (Unique lncRNA IDs)
    # prom_count=$(cut -f4 results/${SP}_mode_Promoter_Donation.txt | sort -u | wc -l)
    exon_count=$(cut -f4 results/${SP}_mode_Exonization_Internal.txt | sort -u | wc -l)

    echo "--- $SP Statistics (True TEs + Unclassified) ---"
    # echo "  - lncRNAs with TE Promoters: $prom_count"
    echo "  - lncRNAs with Internal TE Exons: $exon_count"
done