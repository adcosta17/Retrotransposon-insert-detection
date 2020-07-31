import glob
import os

#
# Helper functions
#

# return a list of runs for the sample


def get_reference_for_test(wildcards):
    return config["reference_"+wildcards.test]

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

def get_sample_name(wildcards):
    return wildcards.sample

def get_suffix(wildcards):
    return "all.read_insertions.repbase_annotated.high_confidence.read_length.tsv"

def get_old_tsv(wildcards):
    return "Oct12_tsvs/"+wildcards.sample+".read_insertions.repbase_annotated.reference_filtered.tsv"

def get_regions():
    f = open("/.mounts/labs/simpsonlab/users/adcosta/code/Retrotransposon-insert-detection/insertion_detection/split_regions.bed")
    regions = list()
    for line in f:
        regions.append(line.rstrip())
    return regions

regions = get_regions()

def get_reference_default(wildcards):
    return config["reference_all"]

#
# Top level targets
#

rule all_plots:
    input:
        expand("{s}/plots/{s}.{t}.LINE.png", t=config["tests"], s=config["samples"])

rule all_normalized_counts:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_counts.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts_updated_annotation:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts_combined:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_sample_bams:
    input:
        expand("{s}/filtered_mapped/{s}.{t}.sorted.bam.bai", t=config["tests"], s=config["samples"]),
        expand("{s}/filtered_mapped/{s}.{t}.sorted.bam", t=config["tests"], s=config["samples"])

rule all_repbase_annotated_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.tsv", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.read_soft_clipped.repbase_annotated.tsv", t=config["tests"], s=config["samples"])

rule all_high_confidence_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.high_confidence.tsv", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.read_soft_clipped.repbase_annotated.high_confidence.tsv", t=config["tests"], s=config["samples"])

rule all_read_length_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.high_confidence.read_length.tsv", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.read_soft_clipped.repbase_annotated.high_confidence.read_length.tsv", t=config["tests"], s=config["samples"])

rule all_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.read_length.reference_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_chimeric_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_chimeric_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_insert_high_confidence_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.high_confidence.tsv", t=config["tests"], s=config["samples"])

rule all_repbase_annotated_insert_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.tsv", t=config["tests"], s=config["samples"])

rule all_old_tsv:
    input:
        expand("{s}/old/{s}.{t}.old_mapped_new.tsv", t=config["tests"], s=config["samples"])

rule all_sample_filter_tsv:
    input:
        expand("multi_sample_tsv/combined_sample.tsv")

rule all_vs_all_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_insert_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.tsv", t=config["tests"], s=config["samples"])

rule all_insert_bam:
    input:
        expand("{s}/mapped_inserts/{s}.{t}.sorted.bam", t=config["tests"], s=config["samples"])

rule all_chimeric_filtered:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_phased_bams:
    input:
        expand("{s}/phased/{s}.sorted.phased.bam", s=config["samples"])

rule all_combined_multi_sample:
    input:
        "combined/combined_multi_sample_x_large.tsv",
        "combined/combined_multi_sample_large.tsv",
        "combined/combined_multi_sample_medium.tsv",
        "combined/combined_multi_sample_small.tsv"

#
# Mapping rules
#

# Map a run to the reference genome

rule minimap2_align:
    input:
        get_sample_fastq
    output:
        protected("{sample}/mapped/{sample}.{test}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_for_test,
        input_fastq_folder=get_sample_fastq_folder,
        tmp_loc="{sample}/mapped/{sample}.{test}.tmp",
        tmp_output="{sample}/mapped/{sample}.{test}.sorted.tmp.bam"
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


## Index a bam
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
        bam=protected("{sample}/filtered_mapped/{sample}.{test}.sorted.bam"),
        merged_reads=protected("{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt")
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
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.tsv"
    threads: 1
    params:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/filtered_mapped/{sample}.{test}.sorted.bam.bai",
        full_bam="{sample}/mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
        merged_reads="{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt",
        candidate_insertion_script = srcdir("get_candidate_insertions.py"),
        memory_per_thread="12G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_script} --bam {params.bam} --merged {params.merged_reads} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"


rule run_get_filtered_candidate_soft_clipped:
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_soft_clipped.tsv")
    threads: 1
    params:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/filtered_mapped/{sample}.{test}.sorted.bam.bai",
        candidate_insertion_script = srcdir("get_candidate_insertions_soft_clipped.py"),
        memory_per_thread="128G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_script} --bam {params.bam} --min-insertion-length {config[min_insertion_length]} --sc-length-filter {config[soft_clip_length]} > {output}"


rule run_convert_candidate_insertions_to_fasta:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.tsv"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_convert_candidate_soft_clipped_to_fasta:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_soft_clipped.tsv"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_soft_clipped.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_conversion_script} --input {input} > {output}"


rule run_candidate_soft_clipped_annotation:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_soft_clipped.tsv",
        tab="{sample}/read_analysis/{sample}.{test}.read_soft_clipped.mapped_to_repbase.last.tab"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.tsv")
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[soft_clipped_mapping_fraction]} > {output}"

rule run_candidate_insertion_annotation:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions.tsv",
        tab="{sample}/read_analysis/{sample}.{test}.read_insertions.mapped_to_repbase.last.tab"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.tsv")
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"

rule get_high_confidence_soft_clipped_tsv:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.tsv"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_soft_clip_script = srcdir("get_high_confidence_inserts_and_soft_clips.py"),
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai"
    shell:
        "{config[python_dir]} {params.hc_soft_clip_script} --threads {threads} --sc 1 --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {params.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]}> {output}"

rule get_high_confidence_inserts_tsv:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_insert_script = srcdir("get_high_confidence_inserts_and_soft_clips.py"),
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai"
    shell:
        "{config[python_dir]} {params.hc_insert_script} --threads {threads} --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {params.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"

rule run_convert_hc_insertions_to_fasta:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("convert_hc_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule map_hc_inserts:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.fa"
    output:
        bam=protected("{sample}/mapped_inserts/{sample}.{test}.sorted.bam")
    params:
        memory_per_thread="3G",
        tmp_loc="{sample}/mapped_inserts/{sample}.{test}.tmp",
        ref_to_use= get_reference_for_test
    threads: 20
    shell:
        """
        {config[minimap_dir]} -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | {config[samtools_dir]} sort -T {params.tmp_loc} -o {output}
        """

rule remove_chimeric_inserts_tsv:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv",
        bam="{sample}/mapped_inserts/{sample}.{test}.sorted.bam",
        bam_index="{sample}/mapped_inserts/{sample}.{test}.sorted.bam.bai",
        phased_bam="{sample}/phased/{sample}.sorted.phased.bam",
        phased_bam_index="{sample}/phased/{sample}.sorted.phased.bam.bai"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv"
    threads: 1
    params:
        memory_per_thread="48G",
        chimeric_insert_script = srcdir("remove_chimeric_and_add_haplotype.py")
    shell:
        "{config[python_dir]} {params.chimeric_insert_script} --tsv {input.tsv} --insert-bam {input.bam} --full-bam {input.phased_bam} > {output}"


rule get_read_length_filtered_soft_clipped_tsv:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.read_length.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        filter_script = srcdir("read_length_filter.py"),
        input_fastq_folder=get_sample_fastq_folder
    shell:
        "{config[python_dir]} {params.filter_script} --tsv {input.tsv} --fastq-folder {params.input_fastq_folder} > {output}"

rule get_read_length_filtered_inserts_tsv:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.read_length.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        filter_script = srcdir("read_length_filter.py"),
        input_fastq_folder=get_sample_fastq_folder 
    shell:
        "{config[python_dir]} {params.filter_script} --tsv {input.tsv} --fastq-folder {params.input_fastq_folder} > {output}"

rule filter_by_reference_sample:
    input:
        file="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.read_length.tsv",
        sc="{sample}/read_analysis/{sample}.{test}.read_soft_clipped.repbase_annotated.high_confidence.read_length.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.read_length.reference_filtered.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("remove_insertions_in_reference_samples.py"),
        memory_per_thread="96G"
    shell:
        "{config[python_dir]} {params.reference_filter_script} --reference-sample {config[reference_tsv_to_filter]} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --sc {input.sc} --reference-repeat-bed {config[repeat_bed]} > {output}"

rule filter_by_reference_sample_chimeric:
    input:
        file="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_filtered.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("remove_insertions_in_reference_samples.py"),
        memory_per_thread="36G"
    shell:
        "{config[python_dir]} {params.reference_filter_script} --reference-sample {config[reference_tsv_to_filter]} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --reference-repeat-bed {config[repeat_bed]} > {output}"

rule get_multi_sample_filter_xl:
    input:
        expand("{sample}/phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", sample=config["samples"])
    output:
        x_large="combined/combined_multi_sample_x_large.tsv"
    threads: 1
    params:
        gen_script = srcdir("generate_multi_sample_tsv.py"),
        memory_per_thread="64G"
    shell:
        """
        {config[python_dir]} {params.gen_script} --sample {config[samples_csv]} --max-distance {config[xlarge_window]} > {output.x_large}
        """

rule get_multi_sample_filter_large:
    input:
        expand("{sample}/phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", sample=config["samples"])
    output:
        large="combined/combined_multi_sample_large.tsv"
    threads: 1
    params:
        gen_script = srcdir("generate_multi_sample_tsv.py"),
        memory_per_thread="64G"
    shell:
        """
        {config[python_dir]} {params.gen_script} --sample {config[samples_csv]} --max-distance {config[large_window]} > {output.large}
        """

rule get_multi_sample_filter_medium:
    input:
        expand("{sample}/phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", sample=config["samples"])
    output:
        medium="combined/combined_multi_sample_medium.tsv"
    threads: 1
    params:
        gen_script = srcdir("generate_multi_sample_tsv.py"),
        memory_per_thread="64G"
    shell:
        """
        {config[python_dir]} {params.gen_script} --sample {config[samples_csv]} --max-distance {config[medium_window]} > {output.medium}
        """

rule get_multi_sample_filter_small:
    input:
        expand("{sample}/phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", sample=config["samples"])
    output:
        small="combined/combined_multi_sample_small.tsv"
    threads: 1
    params:
        gen_script = srcdir("generate_multi_sample_tsv.py"),
        memory_per_thread="64G"
    shell:
        """
        {config[python_dir]} {params.gen_script} --sample {config[samples_csv]} --max-distance {config[small_window]} > {output.small}
        """

rule ava_filter_by_reference_sample:
    input:
        file="{sample}/read_analysis/{sample}.{test}.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv",
        x_large="combined/combined_multi_sample_x_large.tsv",
        large="combined/combined_multi_sample_large.tsv",
        medium="combined/combined_multi_sample_medium.tsv",
        small="combined/combined_multi_sample_small.tsv",
        insert_filter="combined/combined_multi_sample_insert_location.tsv",
        low_mapq_filter="combined/combined_multi_sample_low_mapq.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.tsv"
    threads: 2
    params:
        reference_filter_script = srcdir("remove_insertions_in_reference_samples_multi_filter.py"),
        memory_per_thread="8G"
    shell:
        "{config[python_dir]} {params.reference_filter_script} --multi-reference-filter-xlarge-tsv {input.x_large} --multi-reference-filter-large-tsv {input.large} --multi-reference-filter-medium-tsv {input.medium} --multi-reference-filter-small-tsv {input.small} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --sample {wildcards.sample} --reference-repeat-bed {config[repeat_bed]} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} --low-mapq-filter {input.low_mapq_filter} --insert-filter {input.insert_filter} > {output}"


rule filter_inserts_within_other_inserts:
    input:
        expand("{sample}/phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv", sample=config["samples"]),
        expand("{sample}/filtered_mapped/{sample}.all.merged_reads.txt", sample=config["samples"])
    output:
        inserts="combined/combined_multi_sample_insert_location.tsv",
        low_mapq="combined/combined_multi_sample_low_mapq.tsv"
    threads: 20
    params:
        script = srcdir("filter_inserts_within_insert_sequence.py"),
        region_list = srcdir("split_regions.bed"),
        memory_per_thread="10G"
    shell:
        "{config[python_dir]} {params.script} --sample {config[samples_csv]} --insert-location-output {output.inserts} --low-mapq-output {output.low_mapq} --regions-list {params.region_list} --threads {threads}"

rule filter_old_tsvs:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.read_length.reference_filtered.tsv"
    output:
        tsv = "{sample}/old/{sample}.{test}.old_mapped_new.tsv",
        txt = "{sample}/old/{sample}.{test}.old_mapped_new_counts.txt"
    threads: 1
    params:
        filter_script = srcdir("get_diff_versions.py"),
        memory_per_thread="96G",
        old_tsv=get_old_tsv
    shell:
        "{config[python_dir]} {params.filter_script} --new-tsv {input.tsv} --output-tsv {output.tsv} --sample {wildcards.sample} --old-tsv {params.old_tsv} > {output.txt}"

rule get_multi_sample_filter:
    output:
        "multi_sample_tsv/combined_sample.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("generate_multi_sample_tsv.py"),
        suffix=get_suffix,
        memory_per_thread="64G"
    shell:
        "{config[python_dir]} {params.reference_filter_script} --sample {config[sample_csv]} --suffix {params.suffix} > {output}"

rule normalize_counts:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.read_length.reference_filtered.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.normalized_counts.txt"
    threads: 1
    params:
        reference_filter_script = srcdir("normalize_calls.py"),
        memory_per_thread="64G",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "{config[python_dir]} {params.reference_filter_script} --input {input} --sample {params.sample_name} --fastq-folder {params.fastq_folder} > {output}"

rule normalize_counts_ava:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.tsv"
    output:
        all="{sample}/read_analysis/{sample}.{test}.normalized_ava_counts.txt",
        mapped="{sample}/read_analysis/{sample}.{test}.normalized_ava_mapped_counts.txt",
        distributions="{sample}/read_analysis/{sample}.{test}.read_len_distributions.txt"
    threads: 2
    params:
        reference_filter_script = srcdir("normalize_ava_calls.py"),
        memory_per_thread="8G",
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "{config[python_dir]} {params.reference_filter_script} --input {input.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions}"


rule run_candidate_insertion_updated_annotation:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.tsv",
        tab="{sample}/read_analysis/{sample}.{test}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("update_ambiguous.py"),
        memory_per_thread="24G"
    shell:
        "{config[python_dir]} {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"

rule normalize_counts_ava_updated_annotation:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv"
    output:
        all="{sample}/read_analysis/{sample}.{test}.normalized_ava_counts.updated_annoation.txt",
        mapped="{sample}/read_analysis/{sample}.{test}.normalized_ava_mapped_counts.updated_annoation.txt",
        distributions="{sample}/read_analysis/{sample}.{test}.read_len_distributions.updated_annoation.txt"
    threads: 2
    params:
        reference_filter_script = srcdir("normalize_ava_calls.py"),
        memory_per_thread="8G",
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "{config[python_dir]} {params.reference_filter_script} --input {input.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions}"




rule make_plots:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.read_length.reference_filtered.tsv"
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



#
# Phasing and haplotyping rules
#

# Split bam by region
rule split_bam:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/phased/{sample}.{regions}.sorted.bam")
    params:
        memory_per_thread="2G",
        wrapper = srcdir("bash_wrapper.sh")
    threads: 1
    shell:
        "{params.wrapper} 'samtools view -b {input.bam} \"{wildcards.regions}\" > {output}'"

# Combined region bams for each sample
rule combined_split_bams:
    input:
        bam = expand("{sample}/phased/{sample}.{{regions}}.sorted.bam", sample=config["samples"]),
        bai = expand("{sample}/phased/{sample}.{{regions}}.sorted.bam.bai", sample=config["samples"])
    output:
        temp("combined/combined_{regions}.sorted.bam")
    params:
        memory_per_thread="2G",
        wrapper = srcdir("bash_wrapper.sh")
    threads: 1
    shell:
        "{params.wrapper} 'samtools merge {output} {input.bam}'"

# Genotype the sample with longshot
rule phase_bam:
    input:
        bam = "combined/combined_{regions}.sorted.bam",
        bai = "combined/combined_{regions}.sorted.bam.bai"
    output:
        temp("combined/combined_{regions}.vcf")
    params:
        memory_per_thread="30G",
        ref_to_use= get_reference_default,
        wrapper = srcdir("bash_wrapper.sh")
    threads: 2
    shell:
        "{params.wrapper} 'longshot -c 3 -C 5000 -r {wildcards.regions} --bam {input.bam} --ref {params.ref_to_use} --out {output}'"

# Zip and index the phased bams
rule phase_bam_index:
    input:
        "combined/combined_{regions}.vcf"
    output:
        zip = temp("combined/combined_{regions}.vcf.gz"),
        index = temp("combined/combined_{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="6G",
        wrapper = srcdir("bash_wrapper.sh")
    threads: 2
    shell:
        """
        {params.wrapper} 'bgzip -c {input} > {output.zip}'
        {params.wrapper} 'tabix -p vcf {output.zip}'
        """

# Merge the vcf files
rule merge_vcf:
    input:
        zip = expand("combined/combined_{regions}.vcf.gz", regions=regions),
        index = expand("combined/combined_{regions}.vcf.gz.tbi", regions=regions)
    output:
        zip = protected("combined/combined.vcf.gz"),
        index = protected("combined/combined.vcf.gz.tbi")
    params:
        memory_per_thread="4G",
        wrapper = srcdir("bash_wrapper.sh")
    threads: 2
    shell:
        """
        {params.wrapper} 'bcftools concat -a -O v {input.zip} | bcftools sort -O z -o {output.zip}'
        {params.wrapper} 'tabix -p vcf {output.zip}'
        """

# Tag the bams with the haplotype information
rule haplotag_bam:
    input:
        vcf = "combined/combined.vcf.gz",
        tbi = "combined/combined.vcf.gz.tbi",
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        bam = protected("{sample}/phased/{sample}.sorted.phased.bam")
    params:
        memory_per_thread="16G",
        wrapper = srcdir("bash_wrapper.sh"),
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "{params.wrapper} 'whatshap haplotag --ignore-read-groups --reference {params.ref_to_use} {input.vcf} {input.bam} | samtools sort -o {output.bam}'"
