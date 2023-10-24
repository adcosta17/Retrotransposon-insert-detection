#
# Helper functions used throughout 
#

import os
import glob

def get_reference_for_test(wildcards):
    return config["reference_"+wildcards.test]

def get_reference_base(wildcards):
    return config["reference_all"]

def get_reference_base_lra_index(wildcards):
    return config["reference_all"]+".gli"

def get_winnow_kmers(wildcards):
    if config["minimap2_preset"] == "asm20":
        return config["reference_all"]+"_repetitive_k19.txt"
    else:
        return config["reference_all"]+"_repetitive_k15.txt"
# Search fastq folder from config
def get_sample_fastq(wildcards):
    return config['base_dir'] + wildcards.sample + "/fastq/" + wildcards.sample + ".fastq.gz"

def get_sample_fastq_folder(wildcards):
    return config['base_dir'] + "../" +wildcards.sample + "/fastq/"

def get_run_fastq_folder(wildcards):
    return config['base_dir'] + wildcards.sample + "/run_fastq/"

def get_plot_dir(wildcards):
    return config['base_dir'] + wildcards.sample + "/plots/"

def get_sample_name(wildcards):
    return wildcards.sample

def get_old_tsv(wildcards):
    return "Oct12_tsvs/"+wildcards.sample+".read_insertions.repbase_annotated.reference_filtered.tsv"

def get_regions():
    f = open(config['regions'])
    regions = list()
    for line in f:
        regions.append(line.rstrip())
    return regions

regions = get_regions()

def get_chr22_regions():
    f = open(config['chr22_regions'])
    chr22_regions = list()
    for line in f:
        chr22_regions.append(line.rstrip())
    return chr22_regions

chr22_regions = get_chr22_regions()

def get_reference_default(wildcards):
    return config["reference_all"]


