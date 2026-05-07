import re
import sys

input_file_path = sys.argv[1]
ids_file_path = sys.argv[2]
output_file_path = sys.argv[3]


# Step 1: Extract gene IDs from your list (bed lncrna file)
gene_ids = set()
with open(ids_file_path, 'r') as ids_file:
    for line in ids_file:
        match = re.search(r'^([^|]+)', line)
        if match:
            gene_id = match.group(1)
            gene_ids.add(gene_id)

# Step 2: Iterate through the original GTF file and write matching lines to a new GTF file
with open(input_file_path, 'r') as original_file, open(output_file_path, 'w') as novel_genes_file:
    for line in original_file:
        if line.startswith('#'):
            # Skip comments
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue  # Skip lines that do not have sufficient fields
        feature_type = fields[2]
        attributes = fields[8]
        transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
        if transcript_id_match:
            transcript_id = transcript_id_match.group(1)
            
            if transcript_id in gene_ids:
                # If the gene ID matches, write the line to the new GTF file
                novel_genes_file.write(line)

print("New GTF file 'lncrna.gtf' containing matching transcripts created.")


