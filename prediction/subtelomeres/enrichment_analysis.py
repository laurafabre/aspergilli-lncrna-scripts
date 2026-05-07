# to run mamba activate snakemake

# to run mamba activate snakemake

import os
import subprocess
from scipy.stats import hypergeom
import csv

# Function to calculate the total genome size (N) from the chrNameLength.txt file
def calculate_total_genome_size(chr_length_file):
    total_size = 0
    with open(chr_length_file, 'r') as infile:
        for line in infile:
            # Only process lines with the correct number of columns
            columns = line.strip().split()
            if len(columns) < 3:
                continue  # Skip malformed lines
            _, start, end, _ = columns
            total_size += int(end) - int(start)
    return total_size


# Root directory containing the subfolders
root_dir = "/vast/sa58ged/data_analysis/subtelomers"  # Adjust path relative to the script folder

# Open the output CSV file to write the results
output_file = "/vast/sa58ged/data_analysis/subtelomers/scripts/results.csv"  # Adjust the path as needed
with open(output_file, 'w', newline='') as csvfile:
    # Create a CSV writer
    fieldnames = ['Folder Name', 'Total Genome Size (N)', 'Subtelomeric Region Size (M)', 
                  'Overlapping lncRNA Regions (K)', 'Observed Overlaps (k)', 'total lncrna (j)','P-value']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    # Write the header to the file
    writer.writeheader()

    # Iterate through subfolders
    for folder in os.listdir(root_dir):
        folder_path = os.path.join(root_dir, folder)
        if os.path.isdir(folder_path):  # Ensure it's a directory
            # Step 1: Set up the paths for the lncRNA file and chrNameLength file
            lncrna_bed = os.path.join(folder_path, f"{folder}_lncRNAs.bed")  # Adjust the naming convention if needed
            chr_length_file = os.path.join(folder_path, f"{folder}_chrNameLength.txt")
            
            # Step 2: Check if subtelomeric_regions.bed exists in the folder
            subtelomer_bed_file = os.path.join(folder_path, "subtelomeric_regions.bed")
        

            # Step 3: Calculate total genome size (N) from the chrNameLength.txt file
            N = calculate_total_genome_size(chr_length_file)
            print(f"Total genome size (N) for {folder}: {N}")

            # Step 4: Bedtools intersect - lncRNA with subtelomeric regions
            intersect_bed = os.path.join(folder_path, "lncrna_in_subtelomeres.bed")
            intersect_command = f"bedtools intersect -a {subtelomer_bed_file} -b {lncrna_bed} > {intersect_bed}"
            subprocess.run(intersect_command, shell=True)

            # Step 5: Calculate M (subtelomeric region size) and K (lncRNA region size that overlap with subtelomeric regions)
            M = sum([int(end) - int(start) for chrom, start, end in (line.strip().split() for line in open(subtelomer_bed_file))])
            print(f"Total subtelomeric regions (M) for {folder}: {M}")

            # Calculate K as the total length of overlapping regions
            K = 0
            with open(intersect_bed, 'r') as f:
                for line in f:
                    columns = line.strip().split()
                    start = int(columns[1])
                    end = int(columns[2])
                    K += (end - start)

            print(f"Total overlapping lncrna regions (K) for {folder}: {K}")

            # Count the observed overlap (k) which is the number of lines in the intersect file
            with open(intersect_bed, 'r') as f:
                k = sum(1 for _ in f)
            # Count the lncrna regions which is the number of lines in the lncrna bed file
            with open(lncrna_bed, 'r') as f:
                j = sum(1 for _ in f)

            # Hypergeometric test
            p_value = hypergeom.sf(k - 1, N, M, K)  # sf(k-1) gives the upper tail
            print(f"P-value for enrichment in {folder}: {p_value}")

            # Write the results for the current folder to the CSV file
            writer.writerow({
                'Folder Name': folder,
                'Total Genome Size (N)': N,
                'Subtelomeric Region Size (M)': M,
                'Overlapping lncRNA Regions (K)': K,
                'Observed Overlaps (k)': k,
                'total lncrna (j)': j,
                'P-value': p_value
            })
