#!/bin/bash

# Base directory containing the count files
BASE_DIR="/vast/sa58ged/prediction"

# Output directory for discard files
DISCARD_DIR="/vast/sa58ged/data_analysis/samples_discard"

# Run the Python script to process all CSV files
cd /vast/sa58ged/data_analysis
python3 discard_samples_by_mapping.py "$DISCARD_DIR" mapping_percentages_afla.csv mapping_percentages_afum.csv mapping_percentages_anid.csv mapping_percentages_anig.csv

# Check if the Python script ran successfully
if [ $? -ne 0 ]; then
    echo "Error: Python script failed to run."
    exit 1
fi

# Process each discard file
for discard_file in "$DISCARD_DIR"/*_discard.txt; do
    species=$(basename "$discard_file" _discard.txt)
    echo "Processing $species..."

    # Construct the counts directory path for this species
    COUNTS_DIR="${BASE_DIR}/${species}/lncRNA_prediction/counts"

    # Check if the counts directory exists
    if [ ! -d "$COUNTS_DIR" ]; then
        echo "Warning: Counts directory for $species does not exist: $COUNTS_DIR"
        continue
    fi

    # Read the samples to discard
    while IFS= read -r sample
    do
        # Remove the count files for each discarded sample
        if rm -f "${COUNTS_DIR}/${sample}_counts.txt" "${COUNTS_DIR}/${sample}_counts.txt.summary"; then
            echo "Removed files for sample: $sample"
        else
            echo "Warning: Failed to remove files for sample: $sample"
        fi
    done < "$discard_file"

    echo "Finished processing $species"
    echo
done

echo "All species processed. Discarded samples have been removed from their respective count directories."