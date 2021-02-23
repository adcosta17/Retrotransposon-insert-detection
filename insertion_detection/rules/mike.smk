import os

# To run on SGE:
#   snakemake --jobs 500 --keep-going --latency-wait 120 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_stack=32M -l h_rt={params.run_time} -P simpsonlab -b y"

configfile: "config.yaml"

# Get path information for each program
conda_dir = "/.mounts/labs/simpsonlab/users/mmolnar/code/miniconda3/bin"
scripts_dir = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/scripts"
root_dir = os.getcwd()
flye_dir = "/.mounts/labs/simpsonlab/users/mmolnar/code/Flye/bin"

# Declare variables
genome_length = 3100000000
reference = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_no_alt_analysis_set.fna"
chromosome_sizes = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_chromosome_sizes.tsv"
genes = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_genes.bed"
promoters = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_promoters.bed"
ctcf = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_CTCF_binding_sites.bed"
repeats = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_repeats.bed"
calls_header = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/headers/nanopolish_call_header.tsv"
frequency_header = "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/headers/nanopolish_frequency_header.tsv"


chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
               'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']

def get_regions():
    f = open("/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/reference/GRCh38_regions_list.txt")
    regions = list()
    for line in f:
        regions.append(line.rstrip())
    return regions

regions = get_regions()

# Default rule to make everything
rule all:
    input:
        expand("{sample}/mapped/{sample}.bam.bai", sample=config["samples"]),
        expand("{sample}/mapped/{sample}.methylation.bam.bai", sample=config["samples"]),
        expand("{sample}/mapped/{sample}.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/mapped/{sample}.phased.methylation.bam.bai", sample=config["samples"]),
        
        expand("{sample}/structural_variants/{sample}.cuteSV.bedpe", sample=config["samples"]),
                
        expand("{sample}/methylation/{sample}.nanopolish_frequency.tsv", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.nanopolish_calls.tsv.gz", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.pycoMeth_Interval_Aggregate.bed", sample=config["samples"]),
        expand("{sample}/methylation/{sample}.pycoMeth_Interval_Aggregate.tsv", sample=config["samples"]),
        
        expand("{tumor}/dmrs/{tumor}.{normal}.dmrs.annotated.bed", normal=config["normal"], tumor=config["tumor"]),
        expand("{tumor}/dmrs/{tumor}.{normal}.dmrs.promoters.bed", normal=config["normal"], tumor=config["tumor"]),
        expand("{tumor}/dmrs/{tumor}.{normal}.dmrs.ctcf.bed", normal=config["normal"], tumor=config["tumor"]),
        
        expand("{sample}/analysis/assembly/quast_flye/report.tsv", sample=config["samples"]),
        expand("{sample}/analysis/assembly/quast_racon/report.tsv", sample=config["samples"]),
#        expand("{sample}/analysis/assembly/quast_medaka/report.tsv", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{sample}.mosdpeth.html", sample=config["samples"]),
        expand("{tumor}/analysis/coverage/plots/{normal}/{tumor}.depth.10000_window.pdf", normal=config["normal"], tumor=config["tumor"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{tumor}/analysis/dmrs/plots/{normal}", normal=config["normal"], tumor=config["tumor"]),
        expand("{sample}/analysis/nanoplot/{sample}LengthvsQualityScatterPlot_kde.png", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.annotated_structural_variants.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.structural_variant_primers.txt", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.igv.bat", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/plots/{sample}.sniffles_DEL_small.pdf", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/plots/{sample}.sniffles_INS_small.pdf", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/plots/{sample}.sniffles_DEL_large.pdf", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/plots/{sample}.sniffles_INS_large.pdf", sample=config["samples"])
		

include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/analysis.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/assembly.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/dmrs.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/mapping.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/methylation.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/phasing.smk"
include: "/.mounts/labs/simpsonlab/users/mmolnar/code/nanopore_workflow/files/rules/structural_variants.smk"
