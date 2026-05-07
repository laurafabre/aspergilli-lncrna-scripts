import sys
import re
from collections import defaultdict, Counter

# Parameters from your original run
GENES_NEARBY = 3 

def load_orthologs(ortho_file):
    """Maps gene IDs to their Orthogroup ID."""
    gene_to_ortho = {}
    try:
        with open(ortho_file, 'r') as f:
            header = next(f)
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 5: continue
                og_id = cols[4] 
                for gene_id in cols[0:4]:
                    if gene_id and gene_id != "NA":
                        gene_to_ortho[gene_id] = og_id
    except FileNotFoundError:
        print(f"Error: {ortho_file} not found.")
        sys.exit(1)
    return gene_to_ortho

def load_gene_orders(species_list):
    """Loads the physical order of genes for each species."""
    orders = {}
    for spp in species_list:
        path = f"gene_lists/{spp}_protcod_and_u.txt"
        try:
            with open(path, 'r') as f:
                orders[spp] = [line.strip() for line in f]
        except FileNotFoundError:
            print(f"Warning: {path} not found. Skipping {spp}.")
    return orders

def get_anchors(linc_id, gene_list, gene_to_ortho):
    """Finds the orthogroup IDs of protein-coding neighbors."""
    try:
        idx = gene_list.index(linc_id)
        start = max(0, idx - GENES_NEARBY)
        end = min(len(gene_list), idx + GENES_NEARBY + 1)
        neighbors = gene_list[start:idx] + gene_list[idx+1:end]
        return set([gene_to_ortho[g] for g in neighbors if g in gene_to_ortho])
    except ValueError:
        return set()

# 1. Setup
species = ['afla', 'afum', 'anid', 'anig']
gene_to_ortho = load_orthologs("orthologs_genespace.txt")
gene_orders = load_gene_orders(species)

# 2. Map Genes to Families
family_map = defaultdict(list)
fam_file = "/vast/sa58ged/data_analysis/synteny/intergenic_last_ortho/intergenic.fam"

with open(fam_file, "r") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            fam_id, gene_id = parts
            family_map[fam_id].append(gene_id)

# 3. Process Families and Prepare Tables
all_anchor_counts = Counter()

with open("family_summary_table.txt", "w") as out_summary:
    # Write Header
    out_summary.write("Family_ID\tGene_Count\tCore_Status\tProtein_Anchor_Orthogroups\n")
    
    # Sort families numerically (fam1, fam2...)
    sorted_fams = sorted(family_map.keys(), key=lambda x: int(re.search(r'\d+', x).group()))

    for fam in sorted_fams:
        genes = family_map[fam]
        family_anchors = set()
        found_spp = set()
        
        for gene in genes:
            spp = next((s for s in species if s in gene), None)
            if spp and spp in gene_orders:
                anchors = get_anchors(gene, gene_orders[spp], gene_to_ortho)
                family_anchors.update(anchors)
                found_spp.add(spp)
        
        # Update global counts for frequency table
        all_anchor_counts.update(family_anchors)
        
        unique_anchors_str = ",".join(sorted(family_anchors))
        is_core = "CORE" if len(found_spp) == 4 else "PARTIAL"
        
        # Save row to Summary Table
        out_summary.write(f"{fam}\t{len(genes)}\t{is_core}\t{unique_anchors_str}\n")

# 4. Save Anchor Frequency Table
with open("anchor_frequency_table.txt", "w") as out_freq:
    out_freq.write("Protein_Orthogroup_ID\tLincRNA_Family_Count\n")
    # Sort by frequency (most common anchors first)
    for og, count in all_anchor_counts.most_common():
        out_freq.write(f"{og}\t{count}\n")

print("Success! Created two files:")
print("- family_summary_table.txt  (Detailed family info)")
print("- anchor_frequency_table.txt (Protein orthogroup popularity)")