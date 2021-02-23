#
# Rules for SV calling
#
import os
import re

include: "helpers.smk"

chromosomes = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13',
               'chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
# Return minimum threshold for translocation calls
def get_TRA_coverage_threshold(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    default_threshold = 4
#    if os.path.exists(file_name):
#        f = open(file_name)
#        pattern = "Total bases"    
#        for line in f:
#            if re.search(pattern, line):
#                total_bases = line.split(":")[1].strip().replace(',', '')
#                new_threshold = int(0.4 * (float(total_bases)/3100000000))
#                if new_threshold > default_threshold:
#                    return new_threshold
#                else:
#                    return default_threshold
#    else:
    return default_threshold


def get_normalization_factor(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    default_threshold = 40
#    if os.path.exists(file_name):
#        f = open(file_name)
#        pattern = "Total bases"
#        for line in f:
#            if re.search(pattern, line):
#                total_bases = line.split(":")[1].strip().replace(',', '')
#                return int(2 * (int(float(total_bases))/3100000000))
#    else:
    return default_threshold


rule all_sniffles:
    input:
        expand("{s}/structural_variants/{s}.sniffles.bedpe", s=config["samples"])

# Sniffles

rule sniffles:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.sniffles.bedpe")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        min_reads = "3",
        sniffles_dir = config["sniffles_dir"]
    threads: 10
    shell:
        "{params.sniffles_dir} -s {params.min_reads} --max_num_splits 1 -t {threads} -m {input.bam} -b {output} -n -1"

# Get FASTQ stats from NanoPlot
rule nanoplot:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        protected("{sample}/analysis/nanoplot/{sample}LengthvsQualityScatterPlot_kde.png"),
        protected("{sample}/analysis/nanoplot/{sample}NanoStats.txt")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        scale="--loglength"
    threads: 8
    shell:
        "NanoPlot {params.scale} -t {threads} -o {wildcards.sample}/analysis/nanoplot \
        -p {wildcards.sample} --title {wildcards.sample} --fastq {input}"

# Filter sniffles calls for translocations 
rule filter_TRA:
    input:
        variants = "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.filter_TRA.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        grep_TRA="TRA",
        grep_chr="-v -E 'random|chrUn|alt|fix|chrM|chrEBV'",
        threshold = lambda wildcards: get_TRA_coverage_threshold(wildcards.sample),
        awk="'$2>=0 && $3>=0 && $5>=0 && $6>=0 && $6>=$5 && $3>=$2'"
    threads: 1
    shell:
        "grep {params.grep_TRA} {input.variants} | grep {params.grep_chr} | awk {params.awk} | \
        awk '$12>={params.threshold}' > {output}"


# Find large insertions in the dataset
rule find_large_INS:
    input:
        "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.large_INS.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        grep_INS="INS",
        grep_chr="-v -E 'random|chrUn|alt|fix|chrM|chrEBV'",
        awk="'$17>200 && $17<8000'",
        cut="-f1-6",
        slop="5000",
        chromosome_sizes = srcdir("../utils/GRCh38_chromosome_sizes.tsv")
    threads: 1
    shell:
        "grep {params.grep_INS} {input} | grep {params.grep_chr} | awk {params.awk} | cut {params.cut} | \
        bedtools slop -i stdin -g {params.chromosome_sizes} -b {params.slop} > {output}"


# Filter regions with large insertions
rule filter_large_INS:
    input:
        bedpe = "{sample}/structural_variants/{sample}.filter_TRA.bedpe",
        large_INS = "{sample}/structural_variants/{sample}.large_INS.bed"
    output:
        protected("{sample}/structural_variants/{sample}.filter_large_INS.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk="-v OFS='\t' -F'\t' '{{print $4, $5, $6, $1, $2, $3}}'"
    threads: 1
    shell:
        """
        bedtools intersect -v -wa -a {input.bedpe} -b {input.large_INS} | awk {params.awk} | \
        bedtools intersect -v -wa -a stdin -b {input.large_INS} | awk {params.awk} > {output}
        """


# Find reads with low mapping quality
rule find_low_mapped_reads:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.sam")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        quality="20",
        max="61"
    threads: 1
    shell:
        "samtools view -h {input.bam} | awk '$5<{params.quality} || $5>={params.max} {{print $0}}' > {output}"

# Find the coverage of regions with low mapping quality reads
rule find_low_mapping_coverage:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.sam"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.cov"),
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        grep="-v -E 'random|chrUn|alt|fix|chrM|chrEBV'"
    threads: 1
    shell:
        "samtools view -S -b -h {input} | samtools depth - | grep {params.grep} > {output}"
      
# Make a bed file of the regions with low mapping quality      
rule make_low_mapped_bed:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.cov"
    output:
        protected("{sample}/structural_variants/{sample}.low_mapped_reads.bed")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        slop="5000",
        distance="1000",
        min_coverage="2",
        chromosome_sizes = srcdir("../utils/GRCh38_chromosome_sizes.tsv")
    threads: 1
    shell:
        "SURVIVOR bincov {input} {params.distance} {params.min_coverage} | bedtools slop -i stdin -g {params.chromosome_sizes} -b {params.slop} > {output}"

# Filter regions with low mapping quality reads
rule filter_low_mapped_regions:
    input:
        bedpe = "{sample}/structural_variants/{sample}.filter_large_INS.bedpe",
        bed = "{sample}/structural_variants/{sample}.low_mapped_reads.bed"       
    output:
        protected("{sample}/structural_variants/{sample}.filter_low_mapped_reads.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk="-v OFS='\t' -F'\t' '{{print $4, $5, $6, $1, $2, $3}}'",
    threads: 1
    shell:
        """
        bedtools intersect -v -wa -a {input.bedpe} -b {input.bed} | awk {params.awk} | bedtools intersect -v -wa -a stdin -b {input.bed} | awk {params.awk} > {output}
        """

rule samtools_depth:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    threads: 1
    shell:
        "samtools depth -r {wildcards.chromosomes} {input.bam} > {output}"


rule coverage_by_window_100:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv"
    output:
        temp("{sample}/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="100",
        threshold = lambda wildcards: get_normalization_factor(wildcards.sample),
        coverage_by_window_script = srcdir("../scripts/coverage_by_window.py")
    threads: 1
    shell:
        "python {params.coverage_by_window_script} -i {input.depth} \
        -o {output} -w {params.window} -f {params.threshold}"

rule cat_coverage_by_window_100:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_100/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    threads: 1
    shell:
        "cat {input} > {output}"


# Make a bed file of regions with high or low coverage
rule find_abnormal_coverage:
    input:
        "{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed"
    output:
        protected("{sample}/structural_variants/{sample}.abnormal_coverage.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk="'$5==0 || $5==1'",
        merge_length="1",
        slop="5000",
        chromosome_sizes = srcdir("../utils/GRCh38_chromosome_sizes.tsv")
    threads: 1
    shell:
        "awk {params.awk} {input} | bedtools merge -i stdin -d {params.merge_length} | \
        bedtools slop -i stdin -g {params.chromosome_sizes} -b {params.slop} > {output}"


# Filter translocations in regions with high or low coverage
rule filter_abnormal_coverage:
    input:
        translocations = "{sample}/structural_variants/{sample}.filter_low_mapped_reads.bedpe",
        cov_bed = "{sample}/structural_variants/{sample}.abnormal_coverage.bed"
    output:
        protected("{sample}/structural_variants/{sample}.filter_abnormal_coverage.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk="-v OFS='\t' -F'\t' '{{print $4, $5, $6, $1, $2, $3}}'"
    threads: 1
    shell:
        """
        bedtools intersect -v -wa -a {input.translocations} -b {input.cov_bed} | awk {params.awk} | \
        bedtools intersect -v -wa -a stdin -b {input.cov_bed} | awk {params.awk} > {output}
        """


# Find structural variants with cuteSV
rule cuteSV:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.cuteSV.vcf")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        working_dir="{sample}/structural_variants/"
    threads: 8
    shell:
        "cuteSV --min_support 3 --max_split_parts 2 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
        --threads {threads} {input.bam} {config[reference_all]} {output} {params.working_dir}"


rule vcf2bedpe:
    input:
        "{sample}/structural_variants/{sample}.cuteSV.vcf"
    output:
        protected("{sample}/structural_variants/{sample}.cuteSV.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        script = srcdir("../scripts/vcf2bedpe.py")
    threads: 1
    shell:
        "python {params.script} {input} {output}"


# Break point assembly
rule breakpoint_assembly:
    input:
        fastq = "{sample}/fastq/{sample}.fastq.gz",
        gzi = "{sample}/fastq/{sample}.fastq.gz.gzi",
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai",
        SVs = "{sample}/structural_variants/{sample}.filter_abnormal_coverage.bedpe"
    output:
        protected("{sample}/structural_variants/bp_assembly/combined_corrected_small_window.fa")
    params:
        memory_per_thread="64G",
        run_time="1:0:0:0",
        bp_assemble_script = srcdir("../scripts/bp_assemble.py")
    threads: 1
    shell:
        "python {params.bp_assemble_script} \
        --sniffles-input {input.SVs} --input-bam {input.bam} --input-fastq {input.fastq} \
        --output-folder {wildcards.sample}/structural_variants/bp_assembly \
        --reference-genome {config[reference_all]} --racon {config[racon_dir]} --small-window"

