import os
import sys
import subprocess

species_u = ["afla_u", "afum_u", "anid_u", "anig_u"]

with open("orthologs_genespace.txt", "r") as orthologs:
    orth = []
    header = next(orthologs).strip().split("\t")  # Read header and split
    for line in orthologs:
        line = line.strip()  # Remove extra spaces and newlines
        if not line:
            continue  # Skip empty lines

        cols = line.split("\t")  # Split the line by tabs

        if len(cols) >= 4:  # Ensure the line has enough columns
            orth.append([cols[0], cols[1], cols[2], cols[3]])  # Store the first four columns
        else:
            print("Skipping malformed line:", line)  # Debugging message


### Run synteny script

def make_family_command(spp_list, directory):
    N_seqs_list = []

    for spp in spp_list:
        if "_u" in spp:
            gene_list_file = "gene_lists/%s_protcod_and_u.txt" % (spp.split("_")[0])
        else:
            gene_list_file = "gene_lists/%s_protcod_and_u_x.txt" % (spp)

        N_seqs = 0
        with open(gene_list_file, "r+") as in_file:
            for line in in_file:
                if "MSTRG" in line:
                    N_seqs += 1
        N_seqs_list.append(str(N_seqs))

    if "_u" in spp:
        run_synteny = "python ../synteny_nematodesv4GH_mod_for_candidas_5spp.py gene_lists/afla_protcod_and_u.txt gene_lists/afum_protcod_and_u.txt gene_lists/anid_protcod_and_u.txt gene_lists/anig_protcod_and_u.txt orthologs_genespace.txt %s/synteny_info_output_intergenic.txt 3 3 1 no > %s/output_synteny_intergenic" % (directory, directory)
        classify = "python ../classifyFamiliesv5_VennGH_mod_for_candidas_5spp.py " + str(" ".join(N_seqs_list)) + " %s/synteny_info_output_intergenic.txt %s/intergenic.fam %s/intergenic.txt %s/intergenic_FAM_VENN.R %s/intergenic_GENES_VENN.R > %s/output_intergenic_families" % (directory, directory, directory, directory, directory, directory)

        print(run_synteny)
        print(classify)
        subprocess.call(run_synteny, shell=True)
        subprocess.call(classify, shell=True)

    ### run r scripts to make venn diagrams
    run_r_scripts = "Rscript %s/intergenic_FAM_VENN.R && Rscript %s/intergenic_GENES_VENN.R" % (directory, directory)
    subprocess.call(run_r_scripts, shell=True)   

make_family_command(species_u, "intergenic")

