#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["base_ref"]

def get_tldr_bams(wildcards):
    bam_list = ""                                                                                                                                                     
    for sample in config["samples_all"]:                                                                                                                              
        bam_list += config['base_dir']+sample+"/winnow_mapped/"+sample+".sorted.bam,"                                                                                 
    bam_list = bam_list[:-1]                                                                                                                                          
    return bam_list 

rule xtea_before:
    input:
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bai="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        tsv="{sample}/xtea_before/{sample}/classified_results.txt.Combined.txt"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        echo {wildcards.sample} > {wildcards.sample}_before_sample_id.txt
        echo {wildcards.sample} {input.bam} > {wildcards.sample}_before_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i {wildcards.sample}_before_sample_id.txt -b {wildcards.sample}_before_long_read_bam_list.txt -p {wildcards.sample}/xtea_before -o {wildcards.sample}/xtea_before/{wildcards.sample}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190
        rm {wildcards.sample}/xtea_before/{wildcards.sample}_submit_jobs.sh
        cd {wildcards.sample}/xtea_before/{wildcards.sample}
        chmod +x run_xTEA_pipeline.sh
        cd ../../../
        ./{wildcards.sample}/xtea_before/{wildcards.sample}/run_xTEA_pipeline.sh
        """

rule xtea_after:
    input:
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bai="{sample}/winnow_realign/{sample}.sorted.bam.bai"
    output:
        tsv="{sample}/xtea_after/{sample}/classified_results.txt.Combined.txt"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref
    shell:
        """
        echo {wildcards.sample} > {wildcards.sample}_after_sample_id.txt
        echo {wildcards.sample} {input.bam} > {wildcards.sample}_after_long_read_bam_list.txt
        {config[xtea_dir]}/bin/xtea_long -i {wildcards.sample}_after_sample_id.txt -b {wildcards.sample}_after_long_read_bam_list.txt -p {wildcards.sample}/xtea_after -o {wildcards.sample}/xtea_after/{wildcards.sample}_submit_jobs.sh --clean --rmsk {config[xtea_dir]}/LINE/hg38/hg38_L1_larger_500_with_all_L1HS.out -r {params.ref} --cns {config[xtea_dir]}/consensus/LINE1.fa --rep {config[xtea_dir]} --xtea {config[xtea_dir]}/xtea_long/ -f 31 -y 15 -n 10 -m 190
        rm {wildcards.sample}/xtea_after/{wildcards.sample}_submit_jobs.sh
        cd {wildcards.sample}/xtea_after/{wildcards.sample}
        chmod +x run_xTEA_pipeline.sh
        cd ../../
        ./{wildcards.sample}/xtea_after/{wildcards.sample}/run_xTEA_pipeline.sh
        """

rule filter_tldr_bams:
    input:
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bai="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/tldr_bams/{sample}.sorted.bam")
    threads: 1
    params:
        memory_per_thread="12G",
        script = srcdir("../scripts/filter_bam_coverage.py")
    shell:
        """
        python {params.script} --input {input.bam} --output {output}
        """

rule tldr_before_chrom:
    input:
        bam=expand("{sample}/tldr_bams/{sample}.sorted.bam", sample=config["samples_all"]),
        bai=expand("{sample}/tldr_bams/{sample}.sorted.bam.bai", sample=config["samples_all"])
    output:
        tsv="tldr_combined/DDTS_combined.{chrom}.table.txt"
    threads: 20
    params:
        memory_per_thread="15G",
        ref=get_ref,
        bams=get_tldr_bams
    shell:
        """
        cd tldr_combined
        echo {wildcards.chrom} > DDTS_combined.{wildcards.chrom}.txt
        {config[tldr_dir]}/tldr/tldr -b {params.bams} -r {params.ref} -c DDTS_combined.{wildcards.chrom}.txt -e {config[tldr_dir]}/ref/teref.ont.human.fa -p {threads} -m 1 -o DDTS_combined.{wildcards.chrom}
        """

rule combined_split_bams_winnow_realign_tldr:
    input:
        bam = expand("{sample}/winnow_phased/{sample}.{{chr22_regions}}.sorted.bam", sample=config["samples"]),
        bai = expand("{sample}/winnow_phased/{sample}.{{chr22_regions}}.sorted.bam.bai", sample=config["samples"])
    output:
        temp("combined_winnow_tldr/combined_{chr22_regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools merge {output} {input.bam}"

rule tldr_before_chrom_regions:
    input:
        bam="combined_winnow_tldr/combined_{chr22_regions}.sorted.bam",
        bai="combined_winnow_tldr/combined_{chr22_regions}.sorted.bam.bai"
    output:
        tsv="tldr_combined_regions/DDTS_combined.{chr22_regions}.table.txt"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_ref
    shell:
        """
        cd tldr_combined_regions
        {config[tldr_dir]}/tldr/tldr -b ../{input.bam} -r {params.ref} -c DDTS_combined.chr22.txt -e {config[tldr_dir]}/ref/teref.ont.human.fa -p {threads} -m 1 -o DDTS_combined.{wildcards.chr22_regions}
        cd ..
        rm combined_winnow_tldr/combined_{wildcards.chr22_regions}.sorted.bam.bai
        """

