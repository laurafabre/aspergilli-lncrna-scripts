#!/bin/bash

cd /vast/sa58ged/TE_analysis

SPECIES=("afla" "afum" "anid" "anig")

echo "--- Starting Correct TE Density Calculation ---"

for SP in "${SPECIES[@]}"; do
    echo "Processing $SP..."

    RAW="results/${SP}_exon_te_raw.txt"
    EXONS="results/${SP}_exons.bed"

    # 1. Extract REAL TE coordinates per transcript
    # Columns from bedtools -wo:
    # 4 = transcript ID
    # 8 = TE start
    # 9 = TE end
    # Filter out simple repeats etc.
    awk -F'\t' '
        BEGIN{OFS="\t"}
        $10 !~ /Simple_repeat|Low_complexity|Microsatellite|RNA|Satellite/ {
            print $4, $8, $9
        }
    ' "$RAW" | sort -k1,1 -k2,2n > results/${SP}_te_intervals.bed

    # 2. Merge TE intervals per transcript (non-redundant TE coverage)
    # Output: transcriptID start end
    if [ -s results/${SP}_te_intervals.bed ]; then
        bedtools merge -i results/${SP}_te_intervals.bed -c 1 -o distinct \
            > results/${SP}_te_merged.bed
    else
        > results/${SP}_te_merged.bed
    fi

    # 3. Sum merged TE lengths per transcript
    awk '
        BEGIN{FS=OFS="\t"}
        {len[$4] += ($3-$2)}
        END {for (id in len) print id, len[id]}
    ' results/${SP}_te_merged.bed \
        > results/${SP}_te_bp.txt

    # 4. Compute total exon length per transcript
    awk -F'\t' '
        {len[$4] += ($3-$2)}
        END {for (id in len) print id, len[id]}
    ' "$EXONS" | sort -k1,1 > results/${SP}_exon_len.txt

    # 5. Join TE bp + exon length and compute TE density
    sort -k1,1 results/${SP}_te_bp.txt > results/${SP}_te_bp_sorted.txt

    join -a2 -1 1 -2 1 results/${SP}_te_bp_sorted.txt results/${SP}_exon_len.txt \
    | awk '
        BEGIN{OFS="\t"; print "Transcript_ID","TE_BP","Total_Len","TE_Percentage"}
        {
            id=$1
            if (NF==3) {te=$2; total=$3}
            else {te=0; total=$2}
            if (total==0) total=1
            perc=(te/total)*100
            if (perc>100) perc=100
            print id, te, total, perc
        }
    ' > results/${SP}_final_te_density_corrected.tsv

    echo "Finished $SP"
done

echo "--- All species processed ---"
