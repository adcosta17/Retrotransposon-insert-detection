import os
import glob

rule all:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_minimap:
    input:
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts_combined_lra:
    input:
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts_combined_winnow:
    input:
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_normalized_ava_counts_combined_lra_winnow:
    input:
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/lra_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.txt", t=config["tests"], s=config["samples"]),
        expand("{s}/winnow_read_analysis/{s}.{t}.normalized_ava_counts.updated_annoation.txt", t=config["tests"], s=config["samples"])

rule all_lra_phased:
    input:
        expand("{s}/lra_phased/{s}.sorted.phased.bam",s=config["samples"])

rule all_tra_filtered:
    input:
        expand("{s}/structural_variants/bp_assembly/combined_corrected_small_window.fa",s=config["samples"])

#
# Top level targets
#

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
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.tsv", t=config["tests"], s=config["samples"])

rule all_mapq_ct_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_ma_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.reference_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_chimeric_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", t=config["tests"], s=config["samples"])

rule all_insert_mapq_ct_filtered_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv", t=config["tests"], s=config["samples"])

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
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv", t=config["tests"], s=config["samples"])

rule all_insert_tsv:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.tsv", t=config["tests"], s=config["samples"])

rule all_insert_bam:
    input:
        expand("{s}/mapped_inserts/{s}.{t}.sorted.bam", t=config["tests"], s=config["samples"])

rule all_ma_filtered:
    input:
        expand("{s}/read_analysis/{s}.{t}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", t=config["tests"], s=config["samples"])


rule all_combined_multi_sample:
    input:
        "combined/combined_multi_sample_x_large.tsv",
        "combined/combined_multi_sample_large.tsv",
        "combined/combined_multi_sample_medium.tsv",
        "combined/combined_multi_sample_small.tsv"

rule all_insert_mapped_bams:
    input:
        "combined/combined_inserts_mapped.sorted.bam"

rule all_sample_bams_lra:
    input:
        expand("{s}/lra_mapped/{s}.sorted.bam.bai", s=config["samples"]),
        expand("{s}/lra_mapped/{s}.sorted.bam", s=config["samples"])

rule all_sample_bams_winnow:
    input:
        expand("{s}/winnow_mapped/{s}.sorted.bam.bai", s=config["samples"]),
        expand("{s}/winnow_mapped/{s}.sorted.bam", s=config["samples"])

rule all_sample_bams_lra_winnow:
    input:
        expand("{s}/lra_mapped/{s}.sorted.bam.bai", s=config["samples"]),
        expand("{s}/lra_mapped/{s}.sorted.bam", s=config["samples"]),
        expand("{s}/winnow_mapped/{s}.sorted.bam.bai", s=config["samples"]),
        expand("{s}/winnow_mapped/{s}.sorted.bam", s=config["samples"])


include: "rules/helpers.smk"
include: "rules/mapping.smk"
include: "rules/phasing.smk"
include: "rules/rt_analysis.smk"
include: "rules/structural_variaions.smk"
include: "rules/rt_lra_analysis.smk"
include: "rules/rt_winnowmap_analysis.smk"

