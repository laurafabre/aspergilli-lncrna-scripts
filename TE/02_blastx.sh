#!/bin/bash
#SBATCH --job-name=lnc_blastx
#SBATCH --output=lnc_blastx_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3   # This creates 4 tasks (0, 1, 2, 3)
#SBATCH --partition=standard
# Load environment
module purge
mamba activate lncrna_prediction

# Define the species array
SPECIES=("afum" "afla" "anid" "anig")
SP=${SPECIES[$SLURM_ARRAY_TASK_ID]}  # Picks the species based on the array index

DB_PATH="/vast/sa58ged/TE_analysis/databases/MycoMobilome_prot_db"
RESULT_DIR="/vast/sa58ged/TE_analysis/results/blastx"

echo "Processing species: $SP"

INPUT_FASTA="/vast/sa58ged/prediction/${SP}/lncRNA_prediction/${SP}_lncRNAs.fasta"
CLEAN_FASTA="${RESULT_DIR}/${SP}_lncRNAs_clean.fasta"
OUTPUT_BLAST="${RESULT_DIR}/${SP}_blastx_results.tsv"

# 1. Clean header
sed 's/|.*//' "$INPUT_FASTA" > "$CLEAN_FASTA"

# 2. Run BLASTX
blastx -query "$CLEAN_FASTA" \
       -db "$DB_PATH" \
       -out "$OUTPUT_BLAST" \
       -evalue 1e-5 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
       -num_threads 8 \
       -max_target_seqs 1

echo "Finished $SP"