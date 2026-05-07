import pandas as pd

# Load your files (replace 'family_file.txt' and 'protein_file.txt' with your actual filenames)
# Assuming they are tab-separated as in your example
df_families = pd.read_csv('family_summary_table.txt', sep='\t')
df_proteins = pd.read_csv('orthologs_genespace.txt', sep='\t')

# 1. Filter for CORE families
core_df = df_families[df_families['Core_Status'] == 'CORE']

# 2. Extract all myIDs from the Protein_Anchor_Orthogroups column
# We split the comma-separated strings into a flat list
core_myids = []
for entry in core_df['Protein_Anchor_Orthogroups']:
    core_myids.extend([item.strip() for item in entry.split(',')])

# 3. Filter the protein file for these IDs
result = df_proteins[df_proteins['myID'].isin(core_myids)]

# Save the results
result.to_csv('lnc_conserved_core_proteins_results.csv', index=False)
print(f"Extraction complete. {len(result)} matching proteins found.")
