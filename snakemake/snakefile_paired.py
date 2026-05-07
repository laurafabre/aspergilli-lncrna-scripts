import os
from snakemake.io import glob_wildcards
import glob

configfile: "/vast/sa58ged/snakemake_afla/snakemake_afla/config.yaml"
#The rule copy_bam_files needs to be specified the OUTPUT PATH to where COPY the BAM files needed for lncRNA_prediction in the step 2

#change the extension accordingly (.fq, .fastq) in the entire script, especially for the first rules
SAMPLE, READ, = glob_wildcards("rawReads/{sample}_{read}.fastq.gz")


rule all:
	input:
		#expand("rawQC/{sample}_{read}{extention}", sample=SAMPLE, read=READ, extention=["_fastqc.zip", "_fastqc.html"]),
		#"rawQC/multiqc/multiqc_report.html",
		#expand("trimmedQC/{sample}_{ext}{extention}", sample=SAMPLE, ext=["1P","2P"], extention=["_fastqc.zip", "_fastqc.html"]),
		#"trimmedQC/multiqc/multiqc_report.html",
		expand("stringtieAssembly/{sample}.gtf", sample=SAMPLE),
		expand("/home/sa58ged/Aspergillus/bam_file/afla/{sample}.out.bam", sample=SAMPLE),
		"stringtieAssembly/stringtie_assemblies.txt"

	
rule fastp:
	input: 
		read1 = "rawReads/{sample}_1.fastq.gz",
		read2 = "rawReads/{sample}_2.fastq.gz"
	output:
		forwardPaired="trimmedReads/{sample}_1P.fq",
		reversePaired="trimmedReads/{sample}_2P.fq",
		fastp_json="trimmedQC/multiqc/{sample}_fastp.json"
	threads:
		8
	params:
		log="trimmedReads/{sample}.log"
	shell:
		"""
		OMP_NUM_THREADS={threads}
		fastp -i {input.read1} -I {input.read2} --length_required 49 -5 3 -3 3 -W 4 -M 5 --trim_poly_g --trim_poly_x --thread {threads} -o {output.forwardPaired} -O {output.reversePaired} -j {output.fastp_json} 2>{params.log} 
		"""

rule starindex:
	input:
		genome_fasta= config["configs"][config["myspecies"]]["reference_genome"],
		annotation_gff3= config["configs"][config["myspecies"]]["annotation_gff3"]
	output:
		sa="starindex/SA",
		genome_parameter="starindex/genomeParameters.txt"
	threads:
		12
	params:
		path="starindex/"
	shell:
		"""
		STAR --runThreadN 12 --runMode genomeGenerate --genomeDir {params.path} --genomeFastaFiles {input.genome_fasta} --sjdbGTFfile {input.annotation_gff3}  --genomeSAindexNbases 11
		"""

def get_mem_mb(wildcards, attempt):
    return attempt * 50000

rule star:
	input:
		read1=rules.fastp.output.forwardPaired,
		read2=rules.fastp.output.reversePaired,
		genome_parameter=rules.starindex.output.genome_parameter
	output:
		bam="starAligned/{sample}Aligned.sortedByCoord.out.bam",
		log="starAligned/{sample}Log.final.out"
	threads:
		16
	resources:
		mem_mb = get_mem_mb,
		time="02:00:00"
	params:
		prefix="starAligned/{sample}",
		annotation_gff3=config["configs"][config["myspecies"]]["annotation_gff3"]
	shell:
		"""
		OMP_NUM_THREADS={threads}
		STAR --runThreadN {threads} --genomeDir starindex --readFilesIn {input.read1} {input.read2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {params.annotation_gff3} --alignIntronMax 2000
		"""

rule stringtie:
	input:
		bam= rules.star.output.bam,
		annotation_gff3= config["configs"][config["myspecies"]]["annotation_gff3"]
	output:
		assembly_gtf="stringtieAssembly/{sample}.gtf",
		assembly_tab="stringtieAssembly/{sample}_expression_fus.tab"
	params:
		log="stringtieAssembly/{sample}.log"
	shell:
		"""
		stringtie -v -G {input.annotation_gff3} -o {output.assembly_gtf} --rf -A {output.assembly_tab} {input.bam} 2>{params.log}
		"""

rule copy_bam_files:
    input:
        bam=rules.star.output.bam,
    output:
        bam="/home/sa58ged/Aspergillus/bam_file/afla/{sample}.out.bam"
    shell:
        """
        cp {input.bam} {output.bam}
        """
rule create_file_list:
    input:
        expand("stringtieAssembly/{sample}.gtf", sample=SAMPLE)
    output:
        "stringtieAssembly/stringtie_assemblies.txt"
    shell:
        """
        ls -1 stringtieAssembly/*.gtf > {output}
        """