#!/bin/bash

#SBATCH --job-name=filter_lncRNA
#SBATCH --output=out/hmmscan_tasks_%j.out
#SBATCH --error=out/hmmscan_tasks_%j.err
#SBATCH --time=1-00:00:00          # Adjust time limit as needed
#SBATCH --partition=standard       # Adjust partition as needed
#SBATCH --ntasks=1                 # Number of tasks (one for parallel)
#SBATCH --cpus-per-task=12         # Number of CPU cores per task
#SBATCH --mem=50G                  # Memory per node (adjust as needed)

# Load necessary modules (adjust based on your cluster's environment)
module purge

source /home/sa58ged/.bashrc
source /home/sa58ged/mambaforge/etc/profile.d/conda.sh

mamba activate lncrna_prediction

cd /vast/sa58ged/prediction

#rfam
cmscan_cmds=(
  "cmscan --rfam --cut_ga --nohmmonly --tblout afla/filter_lncRNA/mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm afla/lncRNA_prediction/afla_lncRNAs.fasta"
  "cmscan --rfam --cut_ga --nohmmonly --tblout afum/filter_lncRNA/mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm afum/lncRNA_prediction/afum_lncRNAs.fasta"
  "cmscan --rfam --cut_ga --nohmmonly --tblout anid/filter_lncRNA/mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm anid/lncRNA_prediction/anid_lncRNAs.fasta"
  "cmscan --rfam --cut_ga --nohmmonly --tblout anig/filter_lncRNA/mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm anig/lncRNA_prediction/anig_lncRNAs.fasta"
)

# # Define the commands to run
# hmmscan_cmds=(
#   "hmmscan --domtblout afla/filter_lncRNA/pfam_hits.out Pfam-A.hmm afla/lncRNA_prediction/afla_lncRNAs.fasta"
#   "hmmscan --domtblout afum/filter_lncRNA/pfam_hits.out Pfam-A.hmm afum/lncRNA_prediction/afum_lncRNAs.fasta"
#   "hmmscan --domtblout anid/filter_lncRNA/pfam_hits.out Pfam-A.hmm anid/lncRNA_prediction/anid_lncRNAs.fasta"
#   "hmmscan --domtblout anig/filter_lncRNA/pfam_hits.out Pfam-A.hmm anig/lncRNA_prediction/anig_lncRNAs.fasta"
# )

blastx_cmds=(
  "blastx -query afla/lncRNA_prediction/afla_lncRNAs.fasta -db /home/sa58ged/Aspergillus/transdecoder/swissprot -out afla/filter_lncRNA/blastx.outfmt6Pfam -evalue 1e-5 -num_threads 10 -strand plus -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\""
  "blastx -query afum/lncRNA_prediction/afum_lncRNAs.fasta -db /home/sa58ged/Aspergillus/transdecoder/swissprot -out afum/filter_lncRNA/blastx.outfmt6Pfam  -evalue 1e-5 -num_threads 10 -strand plus -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\""
  "blastx -query anid/lncRNA_prediction/anid_lncRNAs.fasta -db /home/sa58ged/Aspergillus/transdecoder/swissprot -out anid/filter_lncRNA/blastx.outfmt6Pfam  -evalue 1e-5 -num_threads 10 -strand plus -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\""
  "blastx -query anig/lncRNA_prediction/anig_lncRNAs.fasta -db /home/sa58ged/Aspergillus/transdecoder/swissprot -out anig/filter_lncRNA/blastx.outfmt6Pfam  -evalue 1e-5 -num_threads 10 -strand plus -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\""
)

# Export commands to a file
#printf "%s\n" "${hmmscan_cmds[@]}" > hmmscan_cmds.txt
printf "%s\n" "${cmscan_cmds[@]}" > cmscan_cmds.txt
printf "%s\n" "${blastx_cmds[@]}" > blastx_cmds.txt

# Run the commands in parallel
cat cmscan_cmds.txt | parallel --eta -j 8 --load 80% --noswap '{}'
cat blastx_cmds.txt | parallel --eta -j 8 --load 80% --noswap '{}'

echo 'processing files...'
# Extract gene IDs and create discard files
species=("afla" "afum" "anid" "anig")

extract_gene_ids() {
    local species=$1
    local tblout_file="${species}/filter_lncRNA/mrum-genome.tblout"
    local blastx_file="${species}/filter_lncRNA/blastx.outfmt6Pfam"
    #local output_file="${species}/filter_lncRNA/${species}_to_discard.txt"
    local output_file="${species}/filter_lncRNA/${species}_to_discard_2.txt"
    #local output_file="${species}/filter_lncRNA/${species}_to_discard_3.txt" # without blastx

    awk 'NR>2 && $0 !~ /^#/ {print $4}' "$tblout_file" > "$output_file"
    # Extract gene IDs from blastx file (1st column) if query coverage (13th column) is >= 50%
    #awk '$13 >= 50 && $3 > 30 {print $1}' "$blastx_file" >> "$output_file"
    # Extract gene IDs from blastx file (1st column) if query coverage (13th column) is >= 50%, 
    # percent identity (3rd column) is > 30, and the ID contains '|u|' but not '|x|'
    awk '$13 >= 50 && $3 > 30 && $1 ~ /\|u\|/ && $1 !~ /\|x\|/ {print $1}' "$blastx_file" >> "$output_file"
    sort -u -o "$output_file" "$output_file"

    echo "Created $output_file with $(wc -l < "$output_file") unique gene IDs"
}

for spp in "${species[@]}"; do
    extract_gene_ids "$spp"
done