#!/bin/bash

# Define variables
name="afum"  # Replace with your project name
ids_file="$name/filter_lncRNA/afum_to_discard_2.txt"  # File with IDs to remove
input_file="$name/lncRNA_prediction/${name}_validate.gff3"  # Input GFF/GTF file
output_file="$name/lncRNA_prediction/${name}_validate_filtered.gff3"  # Output filtered file

# Create a file for transformed IDs
transformed_ids_file="$name/filter_lncRNA/transformed_ids.txt"

# Transform IDs: Extract the part after 'MSTRG.' for the GFF3 ID format
awk -F '|' '{print $1}' "$ids_file" | sed 's/^MSTRG\.[0-9]*\.[0-9]*//' > "$transformed_ids_file"

# Ensure there are no extra spaces or empty lines
sed -i '/^$/d' "$transformed_ids_file"

# Debug: Print the transformed IDs
echo "Transformed IDs:"
cat "$transformed_ids_file"

# Create a temporary file for deletions
temp_delete_file="$name/lncRNA_prediction/to_delete.gff"

# Extract lines with IDs to delete
grep -F -f "$transformed_ids_file" "$input_file" > "$temp_delete_file"

# Debug: Print the temporary delete file content
echo "Temporary delete file content:"
cat "$temp_delete_file"

# Check if temporary delete file has content
if [[ ! -s "$temp_delete_file" ]]; then
    echo "No lines matched for deletion."
    rm "$temp_delete_file"
    cp "$input_file" "$output_file"
    echo "Filtered file created without any deletions."
    exit 0
fi

# Remove specified entries
grep -v -F -f "$temp_delete_file" "$input_file" > "$output_file"

# Cleanup
rm "$temp_delete_file" "$transformed_ids_file"

echo "Filtering completed. Check the output file."
