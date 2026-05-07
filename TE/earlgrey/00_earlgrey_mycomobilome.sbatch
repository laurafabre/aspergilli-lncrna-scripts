#!/bin/bash

#SBATCH --job-name=TE_earl
#SBATCH --output=out/earlgrey_task_%A_%a.out
#SBATCH --error=out/earlgrey_task_%A_%a.err
#SBATCH --time=2-00:00:00          # Adjust time limit as needed
#SBATCH --partition=standard           # Adjust partition as needed
#SBATCH --cpus-per-task=8         # Number of CPU cores per task
#SBATCH --mem=25G                  # Memory per node (adjust as needed)
#SBATCH --array=0-3               # Array range for 4 parallel tasks (0 to 3)

# Load necessary modules (adjust based on your cluster's environment)
module purge

source /home/sa58ged/.bashrc
source /home/sa58ged/mambaforge/etc/profile.d/conda.sh

mamba activate earlgrey

# Define arrays for input parameters
GENOMES=(/home/sa58ged/Aspergillus/reference_genome/afla/afla.fa \
         /home/sa58ged/Aspergillus/reference_genome/afum/afum.fa \
         /home/sa58ged/Aspergillus/reference_genome/anid/anid.fa \
         /home/sa58ged/Aspergillus/reference_genome/anig/anig.fa )

SPECIES=(afla afum anid anig)

OUTPUTS=(/vast/sa58ged/earlgrey_all/curate_mycoMobilome/afla_toby \
         /vast/sa58ged/earlgrey_all/curate_mycoMobilome/afum_toby \
         /vast/sa58ged/earlgrey_all/curate_mycoMobilome/anid_toby \
         /vast/sa58ged/earlgrey_all/curate_mycoMobilome/anig_toby )

# Run the corresponding task for the array index
genome=${GENOMES[$SLURM_ARRAY_TASK_ID]}
species=${SPECIES[$SLURM_ARRAY_TASK_ID]}
output=${OUTPUTS[$SLURM_ARRAY_TASK_ID]}

earlGreyAnnotationOnly -g "$genome" -s "$species" -o "$output" -l /vast/sa58ged/earlgrey_all/curate_mycoMobilome/MycoMobilome_db/MycoMobilome_v1.0/MycoMobilome_v1.0-allConsensus_TE_library.fasta -t 16
