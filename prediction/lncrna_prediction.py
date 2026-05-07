### See annotations above each command to see its purpose.
### Each command can be executed by uncommenting it at the end of the script.

import sys
import os
import subprocess
from Bio import SeqIO

spp=sys.argv[1]
assemblies = sys.argv[2]

###########################################################
######                 Modify                       #######

# Define the common part of the path
reference_genome_path = "/home/sa58ged/Aspergillus/reference_genome/%s"%(spp)
# Define Bam path
BAM_path = "/home/sa58ged/Aspergillus/bam_file/%s"%(spp)
# Define the path to the strand info file
strand_file_path = "/home/sa58ged/Aspergillus/data_fetching_and_QC_and_strand_detection/%s/strand_detection/pseudomapping/strand_info.txt"%(spp)
###########################################################

# Function to create directory if it doesn't exist
def create_dir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

# Create directories if they don't exist
create_dir_if_not_exists("%s/lncRNA_prediction/feelnc" % (spp))
create_dir_if_not_exists("%s/lncRNA_prediction/cpc" % (spp))
create_dir_if_not_exists("%s/lncRNA_prediction/counts" % (spp))
create_dir_if_not_exists("%s/lncRNA_prediction/DE_analysis" % (spp))

### extract genes from fasta and gff
extract_genes = "gffread -w %s/lncRNA_prediction/%s_genes.fasta -W -F -g %s/%s.fa \
%s/%s.gff3"%(spp,spp,reference_genome_path,spp,reference_genome_path,spp)

### renamed entries in fasta file
rename_genes = "sed 's/ /|/g' %s/lncRNA_prediction/%s_genes.fasta > %s/lncRNA_prediction/%s_genes_renamed.fasta"%(spp,spp,spp,spp)

### selects only protein coding gene ids
select_prot_cod_ids = "grep -A 1 $'\tprotein_coding_gene\t' %s/%s.gff3 | \
grep $'\tmRNA\t'| cut -f 9| cut -f 2 -d ';'| sed 's/Parent=//g' > %s/lncRNA_prediction/%s_prot_cod_ids.txt"%(reference_genome_path,spp,spp,spp)

### selects only protein coding sequences
select_prot_cod_fasta = "python scripts/select_only_prot_cod.py %s/lncRNA_prediction/%s_genes_renamed.fasta %s/lncRNA_prediction/%s_prot_cod_ids.txt %s/lncRNA_prediction/%s_only_prot_cod.fasta"%(spp,spp,spp,spp,spp,spp)

### remove prot. coding seqeunces with ambiguous nucleotides
remove_amb_nucl_prot_cod = "python scripts/remove_seqs_with_amb_nucl.py  \
%s/lncRNA_prediction/%s_only_prot_cod.fasta %s/lncRNA_prediction/feelnc/%s_prot_cod_genes_no_amb_nucl.fasta 1> %s/lncRNA_prediction/feelnc/amb_nuclt_output" %(spp,spp,spp,spp,spp)

### run stringtie merge
#### added -F option following Cristina Zamora-Ballesteros et. al > 0.1 FPKM
merge = "stringtie --merge -o %s/lncRNA_prediction/%s_merged.gtf -g 50 \
                -v -G %s/%s.gff3 -c 1.5 -f 0.05\
                %s"%(spp,spp,reference_genome_path,spp,assemblies)


### run gffcompare        
compare = "gffcompare -V %s/lncRNA_prediction/%s_merged.gtf -o %s/lncRNA_prediction/%s_merged_compared \
                -r %s/%s.gff3"%(spp, spp, spp, spp,reference_genome_path,spp)

### select u and x class codes                
non_code = "python scripts/non_cod.py %s/lncRNA_prediction/%s_merged_compared.annotated.gtf %s/lncRNA_prediction/unknown_and_antisense_ids.gtf"%(spp,spp,spp)
        
### select only longest transcirpts with cgat
cgat = "source /home/sa58ged/mambaforge/bin/activate cgat && cgat gtf2gtf \
        --method=filter --filter-method=longest-transcript -I %s/lncRNA_prediction/unknown_and_antisense_ids.gtf \
        > %s/lncRNA_prediction/unknown_and_antisense_ids_longest_transcripts.gtf && conda deactivate"%(spp,spp)

### make fasta file of u and x class codes        
gffread = "gffread -w %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts.fasta -W -F -g %s/%s.fa %s/lncRNA_prediction/unknown_and_antisense_ids_longest_transcripts.gtf"%(spp,spp,reference_genome_path,spp,spp)
    
### rename 
sed = "sed 's/ /_/g' %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts.fasta > %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed.fasta"%(spp,spp,spp,spp)
     
### select transcirpts longer than 200 bp        
select_200 = "python scripts/select_longer_200.py %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed.fasta \
        %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta"%(spp,spp,spp,spp)
        
### remove sequences with ambigouos nucleotides        
remove_amb_nucl_lncRNA = "python scripts/remove_seqs_with_amb_nucl.py %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta \
        %s/lncRNA_prediction/feelnc/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta > %s/lncRNA_prediction/removed_trascripts_with_amb_nuclt.txt"%(spp,spp,spp,spp,spp)
        
### make a gtf file including coding genes. 
mRNA_gtf = "python scripts/helper_mrna_gtf_afum.py %s/%s.gtf %s/lncRNA_prediction/%s_prot_cod_ids.txt %s/lncRNA_prediction/feelnc/%s_cod_genes.gtf"%(reference_genome_path,spp,spp,spp,spp,spp)
#mRNA_gtf = "python scripts/helper_mrna_gtf_anid.py %s/%s.gtf %s/lncRNA_prediction/%s_prot_cod_ids.txt %s/lncRNA_prediction/feelnc/%s_cod_genes.gtf"%(reference_genome_path,spp,spp,spp,spp,spp)
#mRNA_gtf = "python scripts/helper_mrna_gtf_anig.py %s/%s.gtf %s/lncRNA_prediction/%s_prot_cod_ids.txt %s/lncRNA_prediction/feelnc/%s_cod_genes.gtf"%(reference_genome_path,spp,spp,spp,spp,spp)
# for afla I use the modified helper_mrna_gtf_anig.py

#####################################
###    If problem with fasta    #####
### run FEELnc
#feelnc = "source /home/sa58ged/mambaforge/bin/activate lncrna_prediction &&\
# FEELnc_codpot.pl --outdir='%s/lncRNA_prediction/feelnc/feelnc_codpot_out/' -i %s/lncRNA_prediction/unknown_and_antisense_ids_longest_transcripts.gtf -a %s/%s.gtf -g %s/%s.fa --mode=shuffle \
# && conda deactivate"%(spp,spp,reference_genome_path,spp,reference_genome_path,spp)

### rename
#rename_rf="grep -Ff %s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc.txt %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta | awk '/^>/ {header=\"\"; getline header; if (header != \"\") { gsub(/>/, \"\", $0); gsub(/>/, \"\", header); print $0 ORS header }}' > %s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc_rename.txt"%(spp, spp, spp, spp)
#####################################

#feelnc = "source /home/sa58ged/mambaforge/bin/activate lncrna_prediction &&\
# FEELnc_codpot.pl --outdir='%s/lncRNA_prediction/feelnc/feelnc_codpot_out/' -i %s/lncRNA_prediction/feelnc/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta  -a %s/lncRNA_prediction/feelnc/%s_prot_cod_genes_no_amb_nucl.fasta  --mode=shuffle \
# && conda deactivate"%(spp,spp,spp, spp,spp)

feelnc = "source /home/sa58ged/mambaforge/bin/activate lncrna_prediction &&\
 FEELnc_codpot.pl --outdir='%s/lncRNA_prediction/feelnc/feelnc_codpot_out/' -i %s/lncRNA_prediction/feelnc/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta -a /vast/sa58ged/prediction/%s/lncRNA_prediction/feelnc/%s_cod_genes.gtf -g %s/%s.fa --mode=shuffle \
 && conda deactivate"%(spp,spp,spp, spp, spp, reference_genome_path,spp)
 
### select non-coding sequences from FEELnc result
select_non_cod_feelnc= "awk '$11==\"0\"' %s/lncRNA_prediction/feelnc/feelnc_codpot_out/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta_RF.txt| cut -f 1 > %s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc.txt"%(spp,spp,spp)

### run cpc 
cpc = "~/programs/CPC2-beta/bin/CPC2.py -i %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta -o %s/lncRNA_prediction/cpc/%s_results.tab"%(spp,spp,spp,spp)

### select non-coding sequences from cpc result
select_non_cod_cpc = "awk '$8==\"noncoding\"' %s/lncRNA_prediction/cpc/%s_results.tab.txt | cut -f 1 > %s/lncRNA_prediction/cpc/noncoding_ids_cpc.txt"%(spp,spp,spp)

def find_ovelap():
	with open("%s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc.txt"%(spp), "r") as feelnc_ids, open("%s/lncRNA_prediction/cpc/noncoding_ids_cpc.txt"%(spp),"r") as cpc_ids, open("%s/lncRNA_prediction/removed_trascripts_with_amb_nuclt.txt"%(spp), "r") as removed_transcipts , open("%s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt"%(spp),"w") as output:
		cpc_noncod=[]
		feelnc_noncod=[]
		cpat_noncod=[]
		
		for line in feelnc_ids:
			line=line.rstrip()
			feelnc_noncod.append(line.upper())
		
		for line in cpc_ids:
			line=line.rstrip()
			cpc_noncod.append(line.upper())
			

		overlap = list(set(feelnc_noncod) & set(cpc_noncod)) 
		for line in overlap:
			output.write("%s\n"%(line))
			
			
		for transcript in removed_transcipts:
			#print transcript
			transcript=transcript.rstrip().upper()
			if transcript in cpc_noncod:
				#print transcript
				output.write("%s\n"%(transcript))

### make a gtf file including coding genes, x and u transcripts. 
helper_gtf = "python scripts/helper_generate_saf.py %s/lncRNA_prediction/%s_merged_compared.annotated.gtf %s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt %s/lncRNA_prediction/genes_and_noncod_u_and_x.gtf"%(spp,spp,spp,spp)

### again select only the longest transcripts
# dont use this gtf for generating any final GTF becuase it is removing transcript that are in the orginal gtf
cgat2 = "source /home/sa58ged/mambaforge/bin/activate cgat && cgat gtf2gtf\
        --method=filter --filter-method=longest-transcript -I %s/lncRNA_prediction/genes_and_noncod_u_and_x.gtf  \
        > %s/lncRNA_prediction/genes_and_noncod_u_and_x_longest.gtf  && conda deactivate"%(spp,spp)

### code below generates the saf and bed files using the gtf file obtained above. This whole procedure ensures that all non-coding x and u ids are included in final saf file, and that it does not contain k or p class codes.
generate_saf = "python scripts/generate_saf_mstrg.py \
%s/lncRNA_prediction/genes_and_noncod_u_and_x_longest.gtf %s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt %s %s/lncRNA_prediction/counts/%s_unsorted.saf %s/lncRNA_prediction/%s_unsorted.bed \
%s/%s.gff3"%(spp,spp,spp,spp,spp,spp,spp,reference_genome_path,spp)
sort_saf = "sort -k2,2 -k3,3n %s/lncRNA_prediction/counts/%s_unsorted.saf > %s/lncRNA_prediction/counts/%s.saf"%(spp,spp,spp,spp)
sort_bed = "sort -k1,1 -k2,2n %s/lncRNA_prediction/%s_unsorted.bed > %s/lncRNA_prediction/%s.bed"%(spp,spp,spp,spp)

###------------------
# replace the 2 rules upstream fr the new ones to keep exon information
helper_gtf_exon = "python scripts/helper_to_get_exons_lncrna.py %s/lncRNA_prediction/%s_merged_compared.annotated.gtf %s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt %s/lncRNA_prediction/lncRNAs_with_exons.gtf"%(spp,spp,spp,spp)
cgat2_exon = "source /home/sa58ged/mambaforge/bin/activate cgat && cgat gtf2gtf\
        --method=filter --filter-method=longest-transcript -I %s/lncRNA_prediction/lncRNAs_with_exons.gtf  \
        > %s/lncRNA_prediction/%s_lncRNAs_with_exons_longest.gtf  && conda deactivate"%(spp,spp,spp)
#--------------------

### generates fasta and bed files for lncRNAs
generate_lncRNA_fasta = "grep 'MSTRG' %s/lncRNA_prediction/%s.bed > %s/lncRNA_prediction/%s_lncRNAs.bed && gffread -w %s/lncRNA_prediction/%s_lncRNAs.fasta -W -F -g %s/%s.fa %s/lncRNA_prediction/%s_lncRNAs.bed"%(spp,spp,spp,spp,spp,spp,reference_genome_path,spp,spp,spp)

######################### ADD TO THE PIPELINE
get_lncrna_ids = "awk '{print $4}' %s/lncRNA_prediction/%s_lncRNAs.bed > %s/lncRNA_prediction/lncrna_ids.txt"%(spp,spp, spp)

### generates the gtf file for lncrna novel genes
gtf_make = "python scripts/helper_gtf.py %s/lncRNA_prediction/%s_lncRNAs_with_exons_longest.gtf %s/lncRNA_prediction/lncrna_ids.txt %s/lncRNA_prediction/lncrna.gtf"%(spp, spp, spp, spp)
##############################################

### Run Feelnc classifier
feelnc_class = "source /home/sa58ged/mambaforge/bin/activate lncrna_prediction &&\
 FEELnc_classifier.pl -i %s/lncRNA_prediction/lncrna.gtf -a %s/%s.gtf > %s/lncRNA_prediction/feelnc/feelnc_codpot_out/lncRNA_classes.txt \
 && conda deactivate"%(spp,reference_genome_path,spp, spp)

### Merged novel lncrna genes with the reference annotation and delete all non stranded lncRNA 
#merge_delete = "bash scripts/helper_step2.sh %s %s"%(spp, reference_genome_path)


### generate counts for S dataset
#bams_for_mapping_bigRNAseq = "ls %s/*subset*bam > %s/lncRNA_prediction/counts/files_for_counting.txt"%(spp,spp)
#generate_counts_bigRNAseq = 'while read filename; do prefix_filename=\"$(echo ${filename} | rev | cut -f 2 -d"/" | rev)" && featureCounts -F SAF -p -T 2 -s 2 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/${prefix_filename}_counts.txt ${filename}; done < %s/lncRNA_prediction/counts/files_for_counting.txt'%(spp,spp,spp,spp)


def generate_counts_public():
    with open(strand_file_path, "r") as strand_file:
        next(strand_file)
        for line in strand_file:
            line = line.rstrip().split("\t")
            if line[3] == "used for transcript reconstruction":
                if "SR" in line[2]  or "ISR" in line[2]:
                    featurecounts = "featureCounts -F SAF -p -T 2 -s 2 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/%s_counts.txt %s/%s.out.bam"%(spp, spp, spp, line[0], BAM_path, line[0])
                elif "SF" in line[2]  or "ISF" in line[2]:
                    featurecounts = "featureCounts -F SAF -p -T 2 -s 1 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/%s_counts.txt %s/%s.out.bam"%(spp, spp, spp, line[0], BAM_path, line[0])
                subprocess.call(featurecounts, shell=True)

def generate_counts_without_stradness():
    with open(strand_file_path, "r") as strand_file:
        next(strand_file)
        for line in strand_file:
            line = line.rstrip().split("\t")
            if line[3] == "used for transcript reconstruction":
                featurecounts_2 = "featureCounts -F SAF -p -T 2 -s 0 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts_wo_stradness/%s_counts_without_stradness.txt %s/%s.out.bam"%(spp, spp, spp, line[0], BAM_path, line[0])
            subprocess.call(featurecounts_2, shell=True)

def generate_gc_content_data(output_file, input_file):
	concat_fasta="cat %s/lncRNA_prediction/%s_genes.fasta %s/lncRNA_prediction/%s_lncRNAs.fasta > %s/lncRNA_prediction/%s_genes_and_lncRNAs.fasta"%(spp,spp,spp,spp,spp,spp)
	subprocess.call(concat_fasta, shell=True)
	with open("%s"%(output_file),"w") as gc_output:
		for seq_record in SeqIO.parse("%s"%(input_file), "fasta"):
			name = str(seq_record.id)
			seq = str(seq_record.seq)
			seq=seq.upper()
			gc=float(seq.count("G") + seq.count("C")) / float(len(seq))
			if "|u|" in name:
				class_code="u"
			elif "|x|" in name:
				class_code="x"
			elif "intergenic" in name:
				class_code="inter"
			else:
				class_code="pc"
			
			gc_output.write("%s\t%s\t%s\t%s\n"%(name,class_code,gc,spp))



#### Executing

#subprocess.call(extract_genes, shell=True)
#subprocess.call(rename_genes, shell=True)
#subprocess.call(select_prot_cod_ids, shell=True)
#subprocess.call(select_prot_cod_fasta, shell=True)
#subprocess.call(remove_amb_nucl_prot_cod, shell=True)
#subprocess.call(merge, shell=True)
#subprocess.call(compare, shell=True)
#subprocess.call(non_code, shell=True)
#subprocess.call(cgat, shell=True)
#subprocess.call(gffread, shell=True)
#subprocess.call(sed, shell=True)
#subprocess.call(select_200, shell=True)
#subprocess.call(remove_amb_nucl_lncRNA, shell=True)
#subprocess.call(mRNA_gtf, shell=True)
#subprocess.call(feelnc, shell=True)
#subprocess.call(rename_rf, shell=True) # If problem with fasta
#subprocess.call(select_non_cod_feelnc,shell=True)
#subprocess.call(cpc,shell=True)
#subprocess.call(select_non_cod_cpc,shell=True)
#find_ovelap()
#subprocess.call(helper_gtf, shell=True)
#subprocess.call(cgat2, shell=True)
#subprocess.call(generate_saf,shell=True)
#subprocess.call(sort_saf,shell=True)
#subprocess.call(sort_bed,shell=True)
subprocess.call(helper_gtf_exon,shell=True)
subprocess.call(cgat2_exon,shell=True)
#subprocess.call(generate_lncRNA_fasta,shell=True)
#subprocess.call(get_lncrna_ids, shell=True)
subprocess.call(gtf_make,shell=True)
#subprocess.call(feelnc_class, shell=True) 
#subprocess.call(merge_delete,shell=True)
#subprocess.call(bams_for_mapping_bigRNAseq,shell=True)
#subprocess.call(generate_counts_bigRNAseq,shell=True)
#generate_counts_public()
#generate_counts_without_stradness()

# this runs gc content for all species
#for spp in ["afla","afum", "anid", "anig"]: #add more species if needed
#    generate_gc_content_data("%s/lncRNA_prediction/%s_intergenic_gc.tsv" % (spp, spp), \
#                             "%s/lncRNA_prediction/%s_genes_and_lncRNAs.fasta" % (spp, spp))
