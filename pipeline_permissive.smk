import glob
import os

#
# Helper functions
#

# return a list of runs for the sample


def get_reference_for_test(wildcards):
    return config["reference_"+wildcards.test]

def get_hc_bed_for_test(wildcards):
    return config["hc_bed_"+wildcards.test]

# Search fastq folder from config
def get_sample_fastq(wildcards):
    path = config['base_dir'] + wildcards.sample + "/fastq/*.fastq.gz"
    return glob.glob(path)

def get_sample_fastq_folder(wildcards):
    return config['base_dir'] + wildcards.sample + "/fastq/"

def get_run_fastq_folder(wildcards):
    return config['base_dir'] + wildcards.sample + "/run_fastq/"

def get_plot_dir(wildcards):
    return config['base_dir'] + wildcards.sample + "/plots/"

#
# Top level targets
#

rule all_plots:
    input:
        expand("{s}/plots/{s}.{t}.LINE.png", t=config["tests"], s=config["samples"])

rule all_normalized_counts:
    input:
        expand("{s}/read_analysis_permissive/{s}.{t}.normalized_counts.txt", t=config["tests"], s=config["samples"])

rule all_sample_bams:
    input:
        expand("{s}/filtered_mapped/{s}.{t}.sorted.bam.bai", t=config["tests"], s=config["samples"]),
        expand("{s}/filtered_mapped/{s}.{t}.sorted.bam", t=config["tests"], s=config["samples"])

rule all_repbase_annotated_tsv:
    input:
        expand("{s}/read_analysis_permissive/{s}.{t}.read_insertions.repbase_annotated.tsv", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis_permissive/{s}.{t}.read_soft_clipped.repbase_annotated.tsv", t=config["tests"], s=config["samples"])

rule all_high_confidence_tsv:
    input:
        expand("{s}/read_analysis_permissive/{s}.{t}.read_insertions.repbase_annotated.high_confidence.tsv", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis_permissive/{s}.{t}.read_soft_clipped.repbase_annotated.high_confidence.tsv", t=config["tests"], s=config["samples"])

rule all_filtered_tsv:
    input:
        expand("{s}/read_analysis_permissive/{s}.{t}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.reference_filtered.tsv", t=config["tests"], s=config["samples"])

#
# Mapping rules
#

# Map a run to the reference genome

#rule minimap2_align:
#    input:
#        get_sample_fastq
#    output:
#        "{sample}/mapped/{sample}.sorted.bam"
#    params:
#        memory_per_thread="10G",
#        ref_to_use= get_reference_for_test,
#        input_fastq_folder=get_sample_fastq_folder,
#        tmp_loc="{sample}/ref_mapped/{sample}.{test}.tmp",
#        tmp_output="{sample}/ref_mapped/{sample}.{test}.sorted.tmp.bam"
#    threads: 20
#    shell:
#        """
#        i=0
#        for filename in {params.input_fastq_folder}*.fastq.gz; do
#            echo $filename
            #{config[minimap_dir]} -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} $filename | {config[samtools_dir]} sort -T {params.tmp_loc} -o {params.tmp_output}_$i
#            ((i=i+1))
#        done
#        {config[samtools_dir]} merge {output} {params.tmp_output}*
#        rm {params.tmp_output}*
#        """

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

#
# Extract reads containing candidate insertions
#


rule run_get_filtered_bams:
    input:
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
    output:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        merged_reads="{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt"
    threads: 1
    params:
        merge_alignments_script = srcdir("merge_split_alignments.py"),
        input_fastq_folder=get_sample_fastq_folder,
        memory_per_thread="200G",
        tmp = "{sample}/filtered_mapped/{sample}.{test}.tmp_output.bam"
    shell:
        """
        {config[python_dir]} {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.bam} --indel-size {config[min_insertion_length]} --fastq-folder {params.input_fastq_folder} --bam-output {params.tmp} --merged {output.merged_reads}
        {config[samtools_dir]} sort {params.tmp} -T {output.bam}.tmp -o {output.bam}
        rm {params.tmp}
        """

rule run_get_filtered_candidate_insertions:
    input:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/filtered_mapped/{sample}.{test}.sorted.bam.bai",
        full_bam="{sample}/mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
        merged_reads="{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("get_candidate_insertions.py"),
        memory_per_thread="128G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_script} --bam {input.bam} --original-bam {input.full_bam} --merged {input.merged_reads} --min-insertion-length {config[min_insertion_length]} > {output}"


rule run_get_filtered_candidate_soft_clipped:
    input:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/filtered_mapped/{sample}.{test}.sorted.bam.bai"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.tsv"
    threads: 1
    params:
        candidate_insertion_script = srcdir("get_candidate_insertions_soft_clipped.py"),
        memory_per_thread="128G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_script} --bam {input.bam} --min-insertion-length {config[min_insertion_length]} --sc-length-filter {config[soft_clip_length]} > {output}"


rule run_convert_candidate_insertions_to_fasta:
    input:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.tsv"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.fa"
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_convert_candidate_soft_clipped_to_fasta:
    input:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.tsv"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.fa"
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_conversion_script} --input {input} > {output}"


rule run_candidate_soft_clipped_annotation:
    input:
        tsv="{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.tsv",
        tab="{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.mapped_to_repbase.last.tab"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[soft_clipped_mapping_fraction]} > {output}"

rule run_candidate_insertion_annotation:
    input:
        tsv="{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.tsv",
        tab="{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"

rule get_high_confidence_soft_clipped_tsv:
    input:
        tsv="{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.repbase_annotated.tsv",
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.tsv"
    threads: 20
    params:
        memory_per_thread="3G",
        hc_soft_clip_script = srcdir("get_high_confidence_inserts_and_soft_clips.py"),
        hc_bed=get_hc_bed_for_test
    shell:
        "{config[python_dir]} {params.hc_soft_clip_script} --threads {threads} --sc 1 --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} > {output}"

rule get_high_confidence_inserts_tsv:
    input:
        tsv="{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.repbase_annotated.tsv",
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv"
    threads: 20
    params:
        memory_per_thread="3G",
        hc_insert_script = srcdir("get_high_confidence_inserts_and_soft_clips.py"),
        hc_bed=get_hc_bed_for_test
    shell:
        "{config[python_dir]} {params.hc_insert_script} --threads {threads} --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} > {output}"

rule filter_by_reference_sample:
    input:
        file="{sample}/read_analysis_permissive/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv",
        sc="{sample}/read_analysis_permissive/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.tsv"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.reference_filtered.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("remove_insertions_in_reference_samples.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.reference_filter_script} --reference-sample {config[reference_tsv_to_filter]} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --sc {input.sc} > {output}"

rule normalize_counts:
    input:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.reference_filtered.tsv"
    output:
        "{sample}/read_analysis_permissive/{sample}.{test}.normalized_counts.txt"
    threads: 1
    params:
        reference_filter_script = srcdir("normalize_calls.py"),
        memory_per_thread="150G",
        fastq_folder=get_run_fastq_folder
    shell:
        "{config[python_dir]} {params.reference_filter_script} --input {input} --fastq-folder {params.fastq_folder} > {output}"


rule make_plots:
    input:
        "{sample}/read_analysis_permissive/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.reference_filtered.tsv"
    output:
        "{sample}/plots/{sample}.{test}.LINE.png"
    threads: 1
    params:
        reference_filter_script = srcdir("plot_family_coverage.py"),
        memory_per_thread="150G",
        plot_dir=get_plot_dir
    shell:
        "{config[python_dir]} {params.reference_filter_script} --input {input} --output-dir {params.plot_dir} --repbase-fasta {config[repbase]}"


#
# Gather insertion sequences from sniffles output
#

rule make_lastdb:
    input:
        config["repbase"]
    output:
        config["repbase"] + ".lastdb.suf"
    threads: 1
    params:
        memory_per_thread="16G"
    shell:
        "lastdb {config[repbase]}.lastdb {input}"

rule map_insertion_sequences_last:
    input:
        fa="{base}.fa",
        db=config["repbase"] + ".lastdb.suf"
    output:
        "{base}.mapped_to_repbase.last.tab"
    threads: 8
    params:
        memory_per_thread="12G"
    shell:
        "{config[lastal_dir]} -P {threads} -f tab -r1 -a1 -b1 {config[repbase]}.lastdb {input.fa} > {output}"

