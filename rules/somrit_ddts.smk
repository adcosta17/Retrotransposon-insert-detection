
# Delete sequence to simulate somatic insertions

def get_sample(wildcards):
    return config["full_sample"]

def get_fastq(wildcards):
    return config["fastq"]

def get_ref(wildcards):
    return config["reference"]

def get_all_samples(wildcards):
    return config["samples_all_csv"]

def max_depth(wildcards):
    return 35

def get_base_dir(wildcards):
    return config["base_dir"]

def get_repbase(wildcards):
    return config["repbase"]

def get_centromeres(wildcards):
    return config["centromeres"]

def get_telomeres(wildcards):
    return config["telomeres"]

def get_chrom_list(wildcards):
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

def get_tsv_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/"+sample+".tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_sample_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += sample+","
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_fastq_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/fastq/"+sample+".fastq.gz,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_bam_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/mapped/"+sample+".bam,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_output_bam_list(wildcards):
    bam_list = config['base_dir']+wildcards.sample+"/realign/"+wildcards.sample+".tmp.bam"
    return bam_list

def get_output_tsv_list(wildcards):
    bam_list = config['base_dir']+wildcards.sample+"/realign/"+wildcards.sample+".tsv"
    return bam_list

def get_asm_hap1(wildcards):
    return config['asm_hap1']

def get_asm_hap2(wildcards):
    return config['asm_hap2']

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
        "samtools index {input}"



rule map_fastqs:
    output:
        bam="{sample}/mapped/{sample}.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        ref=get_ref,
        fastq="{sample}/fastq/{sample}.fastq.gz"
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.ref} {params.fastq} | samtools sort > {output.bam}
        """

rule extract_inserts:
    input:
        bam="{sample}/mapped/{sample}.bam",
        bai="{sample}/mapped/{sample}.bam.bai",
        fastq="{sample}/fastq/{sample}.fastq.gz"
    output:
        tsv="{sample}/{sample}.tsv",
        merged=temp("{sample}/{sample}.merged.txt")
    threads: 10
    params:
        memory_per_thread="8G",
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-merged {params.base_dir}/{output.merged} --output-tsv {params.base_dir}/{output.tsv} --fastq-file {params.base_dir}/{input.fastq} --threads {threads} --min-read-len 1000 --min-flank-size 500
        cd {params.base_dir}
        """

rule realign_inserts:
    input:
        bam=expand("{sample}/mapped/{sample}.bam",sample=config["samples"]),
        tsv=expand("{sample}/{sample}.tsv",sample=config["samples"]),
        fastq=expand("{sample}/fastq/{sample}.fastq.gz",sample=config["samples"])
    output:
        bam="combined_realign_all/{chrom}.Realigned.bam",
        tsv="combined_realign_all/{chrom}.Realigned.tsv"
    threads: 10
    params:
        memory_per_thread="25G",
        base_dir=get_base_dir,
        depth=max_depth,
        bam_list=get_bam_list,
        tsv_list=get_tsv_list,
        fastq_list=get_fastq_list,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.bam_list} --tsv-list {params.tsv_list} --fastq-list {params.fastq_list} --output-dir {params.base_dir}/combined_realign_all --tsv-prefix Realigned --bam-prefix Realigned --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 7000 --max-depth {params.depth} --chromosome-list {wildcards.chrom} --only-realign --filter-pass
        cd {params.base_dir}
        """


rule merge_realign_inserts:
    input:
        bam= expand("{sample}/mapped/{sample}.bam",sample=config["samples"]),
        cbams = expand("combined_realign_all/{chrom}.Realigned.bam",chrom=config["chroms"]),
        cbams_bai = expand("combined_realign_all/{chrom}.Realigned.bam.bai",chrom=config["chroms"]),
        tsv = expand("combined_realign_all/{chrom}.Realigned.tsv",chrom=config["chroms"])
    output:
        bam="combined_realign_all/Realigned.tmp.bam",
        tsv="combined_realign_all/Realigned.tsv"
    threads: 1
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        depth=max_depth,
        ref=get_ref,
        bams=get_bam_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/combined_realign_all --tsv-prefix Realigned --bam-prefix Realigned --bam-list {params.bams}
        cd {params.base_dir}
        mv combined_realign_all/Realigned.bam combined_realign_all/Realigned.tmp.bam
        """

rule sort_realign_bams:
    input:
        bam="combined_realign_all/Realigned.tmp.bam"
    output:
        bam="combined_realign_all/Realigned.bam",
    threads: 1
    params:
        memory_per_thread="24G"
    shell:
        """
        samtools sort {input.bam} > {output.bam}
        rm {input.bam}
        """

rule classify_inserts:
    input:
        bam="combined_realign_all/Realigned.bam",
        bai="combined_realign_all/Realigned.bam.bai",
        realign_tsv="combined_realign_all/Realigned.tsv"
    output:
        tsv=temp("combined_realign_all/Realigned_classified.tsv")
    threads: 1
    params:
        memory_per_thread="64G",
        base_dir=get_base_dir,
        repbase=get_repbase,
        tsv_list=get_tsv_list,
        fastq_list=get_fastq_list,
        sample_list=get_sample_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py classify --sample-list {params.sample_list} --tsv-list {params.tsv_list} --realign-tsv {params.base_dir}/{input.realign_tsv} --annotation-file {params.repbase} --fastq-list {params.fastq_list} --output-tsv {params.base_dir}/{output.tsv}
        cd {params.base_dir}
        """

rule filter_inserts_chrom:
    input:
        tsv="combined_realign_all/Realigned_classified.tsv",
        bam="combined_realign_all/Realigned.bam",
        bai="combined_realign_all/Realigned.bam.bai"
    output:
        tsv="combined_realign_all/Realigned_classified_filtered.{chrom}.tsv"
    threads: 10
    params:
        memory_per_thread="20G",
        base_dir=get_base_dir,
        ref=get_ref,
        repbase=get_repbase,
        centromeres=get_centromeres,
        telomeres=get_telomeres,
        fastq_list=get_fastq_list
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py filter --threads {threads} --input-tsv {params.base_dir}/{input.tsv} --bam {params.base_dir}/{input.bam} --fastq-list {params.fastq_list} --reference-genome {params.ref} --centromeres {params.centromeres} --telomeres {params.telomeres} --output-tsv {params.base_dir}/{output.tsv} --chromosome-list {wildcards.chrom}
        cd {params.base_dir}
        """

rule merge_filter_inserts:
    input:
        expand("combined_realign_all/Realigned_classified_filtered.{chrom}.tsv",chrom=config["chroms"])
    output:
        "combined_realign_all/Realigned_classified_filtered.tsv"
    threads: 1
    params:
        memory_per_thread="10G"
    shell:
        """
        head -1 combined_realign_all/Realigned_classified_filtered.chr1.tsv > {output}
        for file in {input}
        do
          tail -n+2 "$file" >> {output}
        done
        """

rule inserts_to_fasta:
    input:
        tsv="combined_realign_all/Realigned_classified_filtered.tsv",
        fastq="{sample}/fastq/{sample}.fastq.gz"
    output:
        "{sample}/polymorphic/{sample}_inserts.fa"
    threads: 1
    params:
        memory_per_thread="12G",
        script=srcdir("../scripts/get_insert_fasta_somrit.py")
    shell:
        """
        python {params.script} --tsv {input.tsv} --fastq {input.fastq} --sample {wildcards.sample} > {output}
        """


rule map_hap1_inserts:
    output:
        hap1="{sample}/polymorphic/{sample}_inserts.hap1.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        hap1=get_asm_hap1,
        fastq="{sample}/fastq/{sample}.fastq.gz"
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.hap1} {params.fastq} | samtools sort > {output.hap1}
        """

rule map_hap2_inserts:
    output:
        hap2="{sample}/polymorphic/{sample}_inserts.hap2.bam"
    threads: 10
    params:
        memory_per_thread="20G",
        hap2=get_asm_hap2,
        fastq="{sample}/fastq/{sample}.fastq.gz"
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {params.hap2} {params.fastq} | samtools sort > {output.hap2}
        """

rule filter_polymorphic_inserts:
    output:
        tsv="{sample}/polymorphic/{sample}_inserts.tsv",
        normalized="{sample}/polymorphic/{sample}_normalized_counts.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_polymorphic_somrit.py"),
        bam="{sample}/mapped/{sample}.bam",
        tsv="combined_realign_all/Realigned_classified_filtered.tsv",
        hap1="{sample}/polymorphic/{sample}_inserts.hap1.bam",
        hap2="{sample}/polymorphic/{sample}_inserts.hap2.bam",
        hap1_bai="{sample}/polymorphic/{sample}_inserts.hap1.bam.bai",
        hap2_bai="{sample}/polymorphic/{sample}_inserts.hap2.bam.bai",
    shell:
        """
        python {params.script} --tsv {params.tsv} --hap1 {params.hap1} --hap2 {params.hap2} --bam {params.bam} --bed {params.bed} --sample {wildcards.sample} --normalized {output.normalized} > {output.tsv}
        """

rule filter_polymorphic_inserts_chrom:
    output:
        reads="{sample}/polymorphic/{sample}_reads_{chrom}.txt"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/get_polymorphic_somrit_split.py"),
        bam="{sample}/mapped/{sample}.bam",
        tsv="combined_realign_all/Realigned_classified_filtered.tsv",
        hap1="{sample}/polymorphic/{sample}_inserts.hap1.bam",
        hap2="{sample}/polymorphic/{sample}_inserts.hap2.bam",
        hap1_bai="{sample}/polymorphic/{sample}_inserts.hap1.bam.bai",
        hap2_bai="{sample}/polymorphic/{sample}_inserts.hap2.bam.bai",
    shell:
        """
        python {params.script} --tsv {params.tsv} --bam {params.bam} --sample {wildcards.sample} --chrom {wildcards.chrom} --read-lengths {output.reads}
        """

rule filter_polymorphic_inserts_combined:
    input:
        expand("{{sample}}/polymorphic/{{sample}}_reads_{chrom}.txt",chrom=config["chroms"])
    output:
        norm="{sample}/results/{sample}.norm.txt",
        inserts="{sample}/results/{sample}.inserts.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        script=srcdir("../scripts/normalize_effective_all.py"),
        bam="{sample}/mapped/{sample}.bam",
        tsv="combined_realign_all/Realigned_classified_filtered.tsv",
        hap1="{sample}/polymorphic/{sample}_inserts.hap1.bam",
        hap2="{sample}/polymorphic/{sample}_inserts.hap2.bam",
        hap1_bai="{sample}/polymorphic/{sample}_inserts.hap1.bam.bai",
        hap2_bai="{sample}/polymorphic/{sample}_inserts.hap2.bam.bai",
        samples=get_all_samples,
        chroms=get_chrom_list
    shell:
        """
        python {params.script} --tsv {params.tsv} --samples {params.samples} --sample {wildcards.sample} --hap1 {params.hap1} --hap2 {params.hap2} --chroms {params.chroms} --folder polymorphic --bam {params.bam}  --normalized {output.norm}  > {output.inserts}
        """




rule sniffles_vcf:
    input:
        bam="combined_realign_all/Realigned.bam",
        bai="combined_realign_all/Realigned.bam.bai"
    output:
        "combined_sniffles.vcf"
    params:
        memory_per_thread="20G",
        min_reads = "1",
    threads: 10
    shell:
        """
        sniffles --input {input.bam} --vcf {output}
        """


rule sniffles_sample_snf:
    input:
        bam="{sample}/mapped/{sample}.bam",
        bai="{sample}/mapped/{sample}.bam.bai"
    output:
        "{sample}/sniffles/{sample}_sniffles.snf"
    params:
        memory_per_thread="12G",
        min_reads = "1",
        bam="{sample}/mapped/{sample}.bam",
        bai="{sample}/mapped/{sample}.bam.bai"
    threads: 5
    shell:
        """
        sniffles --input {params.bam} --snf {output} --non-germline --output-rnames
        """


rule sniffles_sample_combined:
    input:
        expand("{sample}/sniffles/{sample}_sniffles.snf", sample=config["samples"])
    output:
        "sniffles_multi_sample.vcf"
    params:
        memory_per_thread="12G",
        min_reads = "1",
    threads: 5
    shell:
        """
        sniffles --input {input} --vcf {output} --output-rnames
        """

rule sniffles_summarize_per_sample:
    output:
        "{sample}/sniffles/{sample}_sniffles_norm.txt"
    params:
        memory_per_thread="24G",
        min_reads = "1",
        script=srcdir("../scripts/sniffles_comp.py"),
        vcf="sniffles_multi_sample.vcf",
        bam="{sample}/mapped/{sample}.bam",
        samples=get_all_samples,
        hap1="{sample}/polymorphic/{sample}_inserts.hap1.bam",
        hap2="{sample}/polymorphic/{sample}_inserts.hap2.bam",
    threads: 2
    shell:
        """
        python {params.script} --sample {wildcards.sample} --samples {params.samples} --vcf {params.vcf} --bam {params.bam} --bed {params.bed} --hap1 {params.hap1} --hap2 {params.hap2} > {output}
        """

rule sniffles_hap1:
    input:
        bam="{sample}/polymorphic/{sample}_inserts.hap1.bam",
        bam_index="{sample}/polymorphic/{sample}_inserts.hap1.bam.bai"
    output:
        vcf="{sample}/hap1/{sample}.sniffles.vcf",
        snf="{sample}/hap1/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --snf {output.snf} --output-rnames --minsupport 1
        """

rule sniffles_hap2:
    input:
        bam="{sample}/polymorphic/{sample}_inserts.hap2.bam",
        bam_index="{sample}/polymorphic/{sample}_inserts.hap2.bam.bai"
    output:
        vcf="{sample}/hap2/{sample}.sniffles.vcf",
        snf="{sample}/hap2/{sample}.sniffles.snf"
    threads: 10
    params:
        memory_per_thread="10G",
        ref=get_ref
    shell:
        """
        sniffles -t {threads} -i {input.bam} -v {output.vcf} --snf {output.snf} --output-rnames --minsupport 1
        """

rule sniffles_tra_both_haps:
    input:
        vcf1="{sample}/hap1/{sample}.sniffles.vcf",
        vcf2="{sample}/hap2/{sample}.sniffles.vcf"
    output:
        "{sample}/polymorphic/{sample}.tra_both_hap.txt"
    threads: 1
    params:
        memory_per_thread="12G",
        script=srcdir("../scripts/identify_translocations_haps.py"),
        fastq="{sample}/fastq/{sample}.fastq.gz"
    shell:
        """
        python {params.script} --fastq {params.fastq} --hap1 {input.vcf1} --hap2 {input.vcf2} > {output}
        """
