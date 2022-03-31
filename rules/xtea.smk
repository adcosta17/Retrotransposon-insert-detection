#
# Data & Fastq related rules
#

def get_ref(wildcards):
    return config["base_ref"]


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