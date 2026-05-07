#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <name> <directory>"
    exit 1
fi


# Assign the provided name to a variable
name="$1"
directory="$2" 

# Activate conda environment
source /home/sa58ged/.bashrc
source /home/sa58ged/mambaforge/etc/profile.d/conda.sh
mamba activate lncrna_prediction

# Merge annotations
agat_sp_merge_annotations.pl --gff "$directory/$name.gff3" --gff "$name/lncRNA_prediction/lncrna.gtf" --out "$name/lncRNA_prediction/${name}_gene_lncrna"

# Find and remove non-stranded transcripts
grep -P "\\.\t\\." "$name/lncRNA_prediction/${name}_gene_lncrna.gff" > "$name/lncRNA_prediction/to_delete.gff"
echo "Removing non-stranded transcripts..."
grep -v -x -f "$name/lncRNA_prediction/to_delete.gff" "$name/lncRNA_prediction/${name}_gene_lncrna.gff" > "$name/lncRNA_prediction/${name}_validate.gff3"

# Re-checking
echo "Re-checking..."
result=$(grep -P "\\.\t\\." "$name/lncRNA_prediction/${name}_validate.gff3")
echo "$result"

# Convert GFF to GTF
#gffread -E "$name/lncRNA_prediction/${name}_validate.gff3" -T -o "$name/lncRNA_prediction/${name}_validate_final.gtf"
