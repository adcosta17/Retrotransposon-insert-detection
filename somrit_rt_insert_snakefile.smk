import os
import glob

rule all:
    input:
        expand("combined_realign_all/Realigned_classified_filtered.tsv"),
        expand("{s}/polymorphic/{s}_inserts.hap1.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap1.bam.bai",s=config["samples"]),
        expand("{s}/polymorphic/{s}_inserts.hap2.bam.bai",s=config["samples"]),

rule all_assembly_filtered:
    input:
        expand("{s}/results/{s}.norm.txt",s=config["samples"])

include: "rules/somrit_ddts.smk"
