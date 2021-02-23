import os
import re

# Return minimum threshold for translocation calls
def get_TRA_coverage_threshold(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    default_threshold = 4
    if os.path.exists(file_name):
        f = open(file_name)
        pattern = "Total bases"    
        for line in f:
            if re.search(pattern, line):
                total_bases = line.split(":")[1].strip().replace(',', '')
                new_threshold = int(0.4 * (float(total_bases)/3100000000))
                if new_threshold > default_threshold:
                    return new_threshold
                else:
                    return default_threshold
    else:
        return default_threshold

# Find structural variants with Sniffles
rule sniffles:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.sniffles.bedpe")
    params:
        memory_per_thread="32G",
        run_time="1:0:0:0",
        min_reads = "3",
        max_splits = "1"
    log:
        "{sample}/logs/structural_variants/{sample}.sniffles.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.sniffles.txt"
    threads: 1
    shell:
        "{conda_dir}/sniffles -s {params.min_reads} --max_num_splits {params.max_splits} \
        -t {threads} -m {input.bam} -b {output} &> {log}"
                
# Filter sniffles calls for translocations 
rule filter_TRA:
    input:
        variants = "{sample}/structural_variants/{sample}.sniffles.bedpe",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        protected("{sample}/structural_variants/{sample}.filter_TRA.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        grep_TRA="TRA",
        grep_chr="-v -E 'random|chrUn|alt|fix|chrM|chrEBV'",
        threshold = lambda wildcards: get_TRA_coverage_threshold(wildcards.sample),
        awk="'$2>=0 && $3>=0 && $5>=0 && $6>=0 && $6>=$5 && $3>=$2'"
    log:
        "{sample}/logs/structural_variants/{sample}.filter_TRA.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.filter_TRA.txt"
    threads: 1
    shell:
        "grep {params.grep_TRA} {input.variants} 2> {log} | \
        grep {params.grep_chr} 2>> {log} | awk {params.awk} 2>> {log} | \
        awk '$12>={params.threshold}' > {output} 2>> {log}"
      
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
        slop="5000"
    log:
        "{sample}/logs/structural_variants/{sample}.find_large_INS_regions.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.find_large_INS_regions.txt"
    threads: 1
    shell:
        "grep {params.grep_INS} {input} 2> {log} | grep {params.grep_chr} 2>> {log} | \
        awk {params.awk} 2>> {log} | cut {params.cut} 2>> {log} | \
        {conda_dir}/bedtools slop -i stdin -g {chromosome_sizes} -b {params.slop} \
        > {output} 2>> {log}"

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
    log:
        "{sample}/logs/structural_variants/{sample}.filter_large_INS.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.filter_large_INS.txt"
    threads: 1
    shell:
        """
        {conda_dir}/bedtools intersect -v -wa -a {input.bedpe} -b {input.large_INS} \
        2> {log} | awk {params.awk} 2>> {log} | \
        {conda_dir}/bedtools intersect -v -wa -a stdin -b {input.large_INS} 2>> {log} | \
        awk {params.awk} > {output} 2>> {log}
        """

# Find reads with low mapping quality
rule find_low_mapped_reads:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.sam")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        quality="20",
        max="255"
    log:
        "{sample}/logs/structural_variants/{sample}.find_low_mapped_reads.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.find_low_mapped_reads.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -h {input.bam} 2> {log} | \
        awk '$5<{params.quality} || $5>={params.max} {{print $0}}' > {output} 2>> {log}"

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
    log:
        "{sample}/logs/structural_variants/{sample}.find_low_map_regions.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.find_low_map_regions.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -S -b -h {input} 2> {log}| \
        {conda_dir}/samtools depth - 2>> {log} | grep {params.grep} > {output} 2>> {log}"
      
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
        min_coverage="2"
    log:
        "{sample}/logs/structural_variants/{sample}.make_low_mapped_bed.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.make_low_mapped_bed.txt"
    threads: 1
    shell:
        "{conda_dir}/SURVIVOR bincov {input} {params.distance} {params.min_coverage} \
        2> {log} | {conda_dir}/bedtools slop -i stdin -g {chromosome_sizes} \
        -b {params.slop} > {output} 2>> {log}"

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
    log:
        "{sample}/logs/structural_variants/{sample}.filter_low_mapped_regions.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.filter_low_mapped_regions.txt"
    threads: 1
    shell:
        """
        {conda_dir}/bedtools intersect -v -wa -a {input.bedpe} -b {input.bed} 2> {log} | \
        awk {params.awk} 2>> {log} | \
        {conda_dir}/bedtools intersect -v -wa -a stdin -b {input.bed} 2>> {log} | \
        awk {params.awk} 2>> {log} > {output}
        """


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
        slop="5000"
    log:
        "{sample}/logs/structural_variants/{sample}.find_abnormal_coverage.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.find_abnormal_coverage.txt"
    threads: 1
    shell:
        "awk {params.awk} {input} 2> {log} | \
        {conda_dir}/bedtools merge -i stdin -d {params.merge_length} 2>> {log} | \
        {conda_dir}/bedtools slop -i stdin -g {chromosome_sizes} -b {params.slop} \
        > {output} 2>> {log}"

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
    log:
        "{sample}/logs/structural_variants/{sample}.filter_abnormal_coverage.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.filter_abnormal_coverage.txt"
    threads: 1
    shell:
        """
        {conda_dir}/bedtools intersect -v -wa -a {input.translocations} -b \
        {input.cov_bed} 2> {log} | awk {params.awk} 2>> {log} | \
        {conda_dir}/bedtools intersect -v -wa -a stdin -b {input.cov_bed} 2>> {log} | \
        awk {params.awk} > {output} 2>> {log}
        """
          
# Find structural variants with cuteSV
rule cuteSV:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.cuteSV.vcf")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        working_dir="{sample}/structural_variants/"
    log:
        "{sample}/logs/structural_variants/{sample}.cuteSV.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.cuteSV.txt"
    threads: 8
    shell:
        "{conda_dir}/cuteSV --min_support 3 --max_split_parts 2 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
        --threads {threads} {input.bam} {reference} {output} {params.working_dir} &> {log}"

rule vcf2bedpe:
    input:
        "{sample}/structural_variants/{sample}.cuteSV.vcf"
    output:
        protected("{sample}/structural_variants/{sample}.cuteSV.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/structural_variants/{sample}.vcf2bedpe.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.vcf2bedpe.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/vcf2bedpe.py {input} {output} &> {log}"

# Break point assembly
rule breakpoint_assembly:
    input:
        fastq = "{sample}/fastq/{sample}.fastq.gz",
        gzi = "{sample}/fastq/{sample}.fastq.gz.gzi",
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai",
        SVs = "{sample}/structural_variants/{sample}.filter_abnormal_coverage.bedpe"
    output:
        protected("{sample}/structural_variants/bp_assembly/combined_corrected_small_window.fa")
    params:
        memory_per_thread="64G",
        run_time="1:0:0:0"
    log:
        "{sample}/logs/structural_variants/{sample}.breakpoint_assembly.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.breakpoint_assembly.txt"
    threads: 1
    shell:
        "{scripts_dir}/bash_wrapper.sh '{conda_dir}/python {scripts_dir}/bp_assemble.py \
        --sniffles-input {input.SVs} --input-bam {input.bam} --input-fastq {input.fastq} \
        --output-folder {wildcards.sample}/structural_variants/bp_assembly \
        --reference-genome {reference} --racon {conda_dir}/racon --small-window &> {log}'"

