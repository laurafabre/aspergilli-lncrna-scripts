import os

# Root directory containing the subfolders
root_dir = ".."  # Adjust path relative to the script folder

# Iterate through subfolders
for folder in os.listdir(root_dir):
    folder_path = os.path.join(root_dir, folder)
    if os.path.isdir(folder_path):  # Ensure it's a directory
        # Find chromosome length file in the folder
        for file in os.listdir(folder_path):
            if file.endswith("_chrNameLength.txt"):  # Identify chromosome length file
                chrom_length_file = os.path.join(folder_path, file)
                output_bed_file = os.path.join(folder_path, "subtelomeric_regions.bed")
                
                print(f"Processing {chrom_length_file}...")

                # Generate subtelomeric regions
                with open(chrom_length_file, "r") as infile, open(output_bed_file, "w") as outfile:
                    for line in infile:
                        # Skip empty lines or lines that don't have enough columns
                        if not line.strip():  # Skip empty lines
                            continue
                        
                        columns = line.strip().split()
                        
                        # Skip lines with less than 4 columns
                        if len(columns) < 4:
                            print(f"Skipping malformed line: {line.strip()}")
                            continue
                        
                        # Unpack columns (chrom, start, end, _)
                        chrom, start, end, _ = columns
                        start, end = int(start), int(end)
                        
                        # Calculate subtelomeric regions (first and last 10% of chromosome)
                        ten_percent = int((end - start) * 0.1)
                        
                        # Write the first 10% region
                        outfile.write(f"{chrom}\t{start}\t{start + ten_percent}\n")
                        
                        # Write the last 10% region
                        outfile.write(f"{chrom}\t{end - ten_percent}\t{end}\n")

                print(f"Subtelomeric regions saved to {output_bed_file}")
