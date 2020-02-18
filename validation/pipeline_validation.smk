import glob
import os

#
# Helper functions
#

def get_base_reference_for_test(wildcards):
    return config["base_ref"]

def get_bed_for_test(wildcards):
    return config["bed_"+wildcards.test]

def get_hc_bed_for_test(wildcards):
    return config["hc_bed_"+wildcards.test]

# Search fastq folder from config
def get_sample_fastq(wildcards):
    path = config['base_dir'] + wildcards.sample + "/*.fastq.gz"
    return glob.glob(path)

def get_sample_fastq_folder(wildcards):
    return config['base_dir'] + wildcards.sample + "/fastq/"

def get_chrom_for_test(wildcards):
    return wildcards.chrom

def get_sample(wildcards):
    return wildcards.sample

def get_test(wildcards):
    return wildcards.test

#
# Top level targets
#



rule all_spanning:
    input:
        expand("{s}/results/{s}.{t}.read_insertions.truth_set.{c}.tsv", t=config["tests"], s=config["samples"], c=config["chromosomes"]),
        expand("{s}/results/{s}.{t}.read_soft_clipped.truth_set.{c}.tsv", t=config["tests"], s=config["samples"], c=config["chromosomes"])



# Rules

rule minimap2_align:
    input:
        get_sample_fastq
    output:
        "{sample}/ref_mapped/{sample}.{test}.sorted.bam"
    params:
        memory_per_thread="10G",
        ref_to_use= get_base_reference_for_test,
        input_fastq_folder=get_sample_fastq_folder,
        tmp_loc="{sample}/ref_mapped/{sample}.{test}.tmp",
        tmp_output="{sample}/ref_mapped/{sample}.{test}.sorted.tmp.bam"
    threads: 20
    shell:
        """
        i=0
        for filename in {params.input_fastq_folder}*.fastq.gz; do
            echo $filename
            {config[minimap_dir]} -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} $filename | {config[samtools_dir]} sort -T {params.tmp_loc} -o {params.tmp_output}_$i
            ((i=i+1))
        done
        {config[samtools_dir]} merge {output} {params.tmp_output}*
        rm {params.tmp_output}*
        """

rule split_by_chrom:
    input:
        "{sample}/ref_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/ref_mapped/{sample}.{test}.sorted.bam.bai"
    output:
        "{sample}/ref_mapped/{sample}.{test}.sorted.{chrom}.bam"
    params:
        memory_per_thread="10G",
        chrom=get_chrom_for_test
    threads: 1
    shell:
        """
        {config[samtools_dir]} view -b {input} "{params.chrom}" > {output}
        """


# Index a bam
rule make_bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    params:
        memory_per_thread="16G"
    threads: 1
    shell:
        "{config[samtools_dir]} index {input}"

rule get_true_insertions_and_soft_clip_tsv:
    input:
        insert="{sample}/results/{sample}.{test}.spanning_only_inserts.{chrom}.txt",
        sc="{sample}/results/{sample}.{test}.spanning_only_soft_clips.{chrom}.txt"
    output:
        insert_out="{sample}/results/{sample}.{test}.read_insertions.truth_set.{chrom}.tsv",
        sc_out="{sample}/results/{sample}.{test}.read_soft_clipped.truth_set.{chrom}.tsv"
    threads: 1
    params:
        true_insert_script = srcdir("get_true_insert_tsv.py"),
        memory_per_thread="24G",
        input_fastq_folder=get_sample_fastq_folder
    shell:
        "{config[python_dir]} {params.true_insert_script} --insert-out {output.insert_out} --soft-clip-out {output.sc_out} --spanning-insert-reads {input.insert} --spanning-sc-reads {input.sc} --fastq-folder {params.input_fastq_folder}"


rule get_inserts_from_spanning_only:
    input:
        bam="{sample}/ref_mapped/{sample}.{test}.sorted.{chrom}.bam",
        bam_index="{sample}/ref_mapped/{sample}.{test}.sorted.{chrom}.bam.bai"
    output:
        insert="{sample}/results/{sample}.{test}.spanning_only_inserts.{chrom}.txt",
        sc="{sample}/results/{sample}.{test}.spanning_only_soft_clips.{chrom}.txt"
    params:
        spanning_script = srcdir("get_all_spanning_reads.py"),
        wrapper = srcdir("bash_wrapper.sh"),
        memory_per_thread="10G",
        bed_to_use = get_bed_for_test
    threads: 20
    shell:
        "{params.wrapper} '{config[python_dir]} {params.spanning_script} --bed {params.bed_to_use} --threads {threads} --read-to-reference-bam {input.bam} --output-insert-list {output.insert} --output-sc-list {output.sc}'"

rule combine_spanning_inserts:
    output:
        insert="{sample}/results/{sample}.{test}.spanning_only_inserts.txt",
        sc="{sample}/results/{sample}.{test}.spanning_only_soft_clips.txt"
    threads: 1
    params:
        memory_per_thread="10G",
        sample=get_sample,
        test=get_test
    shell:
        """
        cat {params.sample}/results/{params.sample}.{params.test}.spanning_only_inserts.c* >> {output.insert}
        cat {params.sample}/results/{params.sample}.{params.test}.spanning_only_soft_clips.c* >> {output.sc}
        """