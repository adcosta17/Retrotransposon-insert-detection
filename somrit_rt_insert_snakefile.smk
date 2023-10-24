import os
import glob

rule all:
    input:
        expand("combined_realign_all/Realigned_classified_filtered.tsv")

rule all_extract:
    input:
        expand("{s}/{s}.tsv",s=config["samples"])

rule all_mapped:
    input:
        expand("{s}/mapped/{s}.bam",s=config["samples"]),
        expand("{s}/mapped/{s}.bam.bai",s=config["samples"])

rule all_polymorphic_mapped:
    input:
        expand("{s}/polymorphic/{s}_inserts.hap1.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap1.bam.bai",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam.bai",s=config["samples"]),

rule all_mapped_both:
    input:
        expand("{s}/polymorphic/{s}_inserts.hap1.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap1.bam.bai",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam.bai",s=config["samples"]),
        expand("{s}/mapped/{s}.bam",s=config["samples"]),
        expand("{s}/mapped/{s}.bam.bai",s=config["samples"])

rule all_assembly_filtered:
    input:
        expand("{s}/polymorphic/{s}_inserts.tsv",s=config["samples"]),
        expand("{s}/polymorphic/{s}_normalized_counts.txt",s=config["samples"])

rule all_assembly_filtered_chrom:
    input:
        expand("{s}/polymorphic/{s}_inserts_{c}.tsv",s=config["samples"], c=config["chroms"]),
        expand("{s}/polymorphic/{s}_normalized_counts_{c}.txt",s=config["samples"], c=config["chroms"])

rule all_assembly_filtered_both:
    input:
        expand("{s}/polymorphic/{s}_inserts_{c}.tsv",s=config["samples"], c=config["chroms"]),
        expand("{s}/polymorphic/{s}_normalized_counts_{c}.txt",s=config["samples"], c=config["chroms"]),
        expand("{s}/polymorphic/{s}_inserts.tsv",s=config["samples"]),
        expand("{s}/polymorphic/{s}_normalized_counts.txt",s=config["samples"])

rule all_assembly_filtered_fpr:
    input:
        expand("{s}/polymorphic_fpr/{s}_inserts.tsv",s=config["samples"]),
        expand("{s}/polymorphic_fpr/{s}_normalized_counts.txt",s=config["samples"])

rule all_sniffles:
    input:
        "sniffles_multi_sample.vcf",
        expand("{s}/polymorphic/{s}_inserts.hap1.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap1.bam.bai",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam.bai",s=config["samples"]),

rule all_sniffles_norm:
    input:
         expand("{s}/sniffles/{s}_sniffles_norm.txt",s=config["samples"])

rule all_inserts_sniffles_hc:
    input:
         expand("{s}/results/{s}.hc_norm.txt",s=config["samples"])

rule all_sniffles_tra:
    input:
        expand("{s}/polymorphic/{s}.tra_both_hap.txt",s=config["samples"])

rule all_assembly_filtered_sniffles:
    input:
        expand("{s}/polymorphic/{s}_inserts.tsv",s=config["samples"]),
        expand("{s}/polymorphic/{s}_normalized_counts.txt",s=config["samples"]),
        "combined_sniffles.vcf"

include: "rules/somrit_ddts.smk"
