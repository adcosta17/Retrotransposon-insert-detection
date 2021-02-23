import os
# Coverage used to set thresholds
def get_normalization_factor(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    default_threshold = 40
    if os.path.exists(file_name):
        f = open(file_name)
        pattern = "Total bases"
        for line in f:
            if re.search(pattern, line):
                total_bases = line.split(":")[1].strip().replace(',', '')
                return int(2 * (int(float(total_bases))/3100000000))
    else:
        return default_threshold

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
    log:
        "{sample}/logs/analysis/{sample}.nanoplot.log"
    benchmark:
        "{sample}/benchmarks/analysis/{sample}.nanoplot.txt"
    threads: 8
    shell:
        "{conda_dir}/NanoPlot {params.scale} -t {threads} -o {wildcards.sample}/analysis/nanoplot \
        -p {wildcards.sample} --title {wildcards.sample} --fastq {input} &> {log}"
        
# Calculate assembly statistics with QUAST
rule quast_flye:
    input:
        "{sample}/assembly/assembly.fasta"
    output:
        protected("{sample}/analysis/assembly/quast_flye/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_flye"
    log:
        "{sample}/logs/analysis/assembly/quast_flye.log"
    benchmark:
        "{sample}/benchmarks/analysis/assembly/quast_flye.txt"
    threads: 8
    shell:
        "{conda_dir}/python {conda_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"

rule quast_racon:
    input:
        "{sample}/assembly/{sample}.racon.fasta"
    output:
        protected("{sample}/analysis/assembly/quast_racon/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_racon"
    log:
        "{sample}/logs/analysis/assembly/quast_racon.log"
    benchmark:
        "{sample}/benchmarks/analysis/assembly/quast_racon.txt"
    threads: 8
    shell:
        "{conda_dir}/python {conda_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"

rule quast_medaka:
    input:
        "{sample}/assembly/medaka_consensus/consensus.fasta"
    output:
        protected("{sample}/analysis/assembly/quast_medaka/report.tsv")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        genome_size="--large",
        del_temp="--space-efficient",
        output_folder="{sample}/analysis/assembly/quast_medaka"
    log:
        "{sample}/logs/analysis/assembly/quast_medaka.log"
    benchmark:
        "{sample}/benchmarks/analysis/assembly/quast_medaka.txt"
    threads: 8
    shell:
        "{conda_dir}/python {conda_dir}/quast {params.genome_size} -r {reference} \
		-o {params.output_folder} -t {threads} {params.del_temp} {input} &> {log}"
        
rule SV_circos_plots:
    input:
        "{sample}/structural_variants/{sample}.filter_low_mapped_reads.bedpe"
    output:
        protected("{sample}/analysis/structural_variants/plots/{sample}.SVs.pdf")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/analysis/structural_variants/plots/SV_circos_plot.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/plots/SV_circos_plot.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/RCircosSV.r {input} {output} &> {log}"

rule samtools_depth:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp("{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/temp_files/samtools_depth/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/benchmarks/temp_files/samtools_depth/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools depth -r {wildcards.chromosomes} {input.bam} > {output} 2> {log}"

rule coverage_by_window_10000:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="10000",
        threshold = lambda wildcards: get_normalization_factor(wildcards.sample)
    log:
        "{sample}/logs/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/benchmarks/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -f {params.threshold} &> {log}"

rule coverage_by_window_100:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="100",
        threshold = lambda wildcards: get_normalization_factor(wildcards.sample)
    log:
        "{sample}/logs/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/benchmarks/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -f {params.threshold} &> {log}"

rule cat_coverage_by_window_100:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_100/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/analysis/coverage/samtools_depth/{sample}.depth.100_window.log"
    benchmark:
        "{sample}/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.100_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

rule cat_coverage_by_window_10000:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_10000/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/analysis/coverage/samtools_depth/{sample}.depth.10000_window.log"
    benchmark:
        "{sample}/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.10000_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

rule plot_genome_coverage:
    input:
        normal = "{normal}/analysis/coverage/samtools_depth/{normal}.depth.10000_window.bed",
        tumor = "{tumor}/analysis/coverage/samtools_depth/{tumor}.depth.10000_window.bed"
    output:
        protected("{tumor}/analysis/coverage/plots/{normal}/{tumor}.depth.pdf")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        genome="hg38"
    log:
        "{tumor}/logs/analysis/coverage/plots/{normal}/{tumor}.depth.log"
    benchmark:
        "{tumor}/benchmarks/analysis/coverage/plots/{normal}/{tumor}.depth.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploteR_coverage.R {input.normal} {input.tumor} {output} {params.genome} &> {log}"
        
rule sniffles_INS_plots_large:
    input:
        "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/analysis/structural_variants/plots/{sample}.sniffles_INS_large.pdf")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        type="INS",
        min="1000",
        max="10000",
        bin_size="100"
    log:
        "{sample}/logs/analysis/structural_variants/plots/{sample}.sniffles_INS_large.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/plots/{sample}.sniffles_INS_large.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/frequency_histogram.R {input} {output} {params.type} \
        {params.min} {params.max} {params.bin_size} '{wildcards.sample} - Sniffles Large {params.type}' &> {log}"

rule sniffles_INS_plots_small:
    input:
        "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/analysis/structural_variants/plots/{sample}.sniffles_INS_small.pdf")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        type="INS",
        min="50",
        max="1000",
        bin_size="10"
    log:
        "{sample}/logs/analysis/structural_variants/plots/{sample}.sniffles_INS_small.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/plots/{sample}.sniffles_INS_small.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/frequency_histogram.R {input} {output} {params.type} \
        {params.min} {params.max} {params.bin_size} '{wildcards.sample} - Sniffles Small {params.type}' &> {log}"

rule sniffles_DEL_plots_large:
    input:
        "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/analysis/structural_variants/plots/{sample}.sniffles_DEL_large.pdf")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        type="DEL",
        min="1000",
        max="10000",
        bin_size="100"
    log:
        "{sample}/logs/analysis/structural_variants/plots/{sample}.sniffles_DEL_large.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/plots/{sample}.sniffles_DEL_large.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/frequency_histogram.R {input} {output} {params.type} \
        {params.min} {params.max} {params.bin_size} '{wildcards.sample} - Sniffles Large {params.type}' &> {log}"

rule sniffles_DEL_plots_small:
    input:
        "{sample}/structural_variants/{sample}.sniffles.bedpe"
    output:
        protected("{sample}/analysis/structural_variants/plots/{sample}.sniffles_DEL_small.pdf")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        type="DEL",
        min="1",
        max="1000",
        bin_size="10"
    log:
        "{sample}/logs/analysis/structural_variants/plots/{sample}.sniffles_DEL_small.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/plots/{sample}.sniffles_DEL_small.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/frequency_histogram.R {input} {output} {params.type} \
        {params.min} {params.max} {params.bin_size} '{wildcards.sample} - Sniffles Small {params.type}' &> {log}"

rule mosdepth:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        bed = temp("{sample}/analysis/coverage/mosdepth/{sample}.{chromosomes}.regions.bed.gz"),
        dist = protected("{sample}/analysis/coverage/mosdepth/{sample}.{chromosomes}.mosdepth.global.dist.txt"),
        summary = protected("{sample}/analysis/coverage/mosdepth/{sample}.{chromosomes}.mosdepth.summary.txt")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        prefix="{sample}/analysis/coverage/mosdepth/{sample}.{chromosomes}"
    log:
        "{sample}/logs/analysis/coverage/mosdepth/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/benchmarks/analysis/coverage/mosdepth/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/mosdepth -n --fast-mode --by 100 --chrom {wildcards.chromosomes} {params.prefix} {input.bam} &> {log}"

rule mosdepth_plots:
    input:
        expand("{{sample}}/analysis/coverage/mosdepth/{{sample}}.{chromosomes}.mosdepth.global.dist.txt", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/plots/{sample}.mosdpeth.html")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/logs/analysis/coverage/mosdepth/{sample}.mosdepth_plots.log"
    benchmark:
        "{sample}/benchmarks/analysis/coverage/mosdepth/{sample}.mosdepth_plots.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/plot-dist.py -o {output} {input} &> {log}"

rule merge_dmrs:
    input:
        "{tumor}/dmrs/{tumor}.{normal}.dmrs.bed"
    output:
        protected("{tumor}/analysis/dmrs/{tumor}.{normal}.merged_dmrs.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk="'$4>1000 && $5>20'",
        merge_length="10000"
    log:
        "{tumor}/logs/analysis/dmrs/{tumor}.{normal}.merged_dmrs.log"
    benchmark:
        "{tumor}/benchmarks/analysis/dmrs/{tumor}.{normal}.merged_dmrs.txt"
    threads: 1
    shell:
        "awk {params.awk} {input} | {conda_dir}/bedtools merge -i stdin \
        -d {params.merge_length} &> {log} > {output}"
                
rule plot_dmrs:
    input:
        normal = "{normal}/methylation/{normal}.nanopolish_frequency.tsv",
        tumor = "{tumor}/methylation/{tumor}.nanopolish_frequency.tsv",
        dmrs = "{tumor}/dmrs/{tumor}.{normal}.dmrs.bed",
        merged_dmrs = "{tumor}/analysis/dmrs/{tumor}.{normal}.merged_dmrs.bed"
    output:
        directory("{tumor}/analysis/dmrs/plots/{normal}")
    params:
        memory_per_thread="32G",
        run_time="1:0:0:0",
        extend="12000",
        output_dir="{tumor}/analysis/dmrs/plots/{normal}"
    log:
        "{tumor}/logs/analysis/dmrs/plots/{normal}/{tumor}.plot_dmrs.log"
    benchmark:
        "{tumor}/benchmarks/analysis/dmrs/plots/{normal}/{tumor}.plot_dmrs.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploteR_methylation.R {input.normal} {input.tumor} \
        {input.dmrs} {ctcf} {input.merged_dmrs} {wildcards.normal} {wildcards.tumor} \
        {params.extend} {root_dir}/{params.output_dir} &> {log}"

rule annotate_SVs:
    input:
        "{sample}/structural_variants/{sample}.filter_abnormal_coverage.bedpe" 
    output:
        protected("{sample}/analysis/structural_variants/{sample}.annotated_structural_variants.bedpe")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        awk1="-v 'OFS=\t' '-F\t' '{{print $4, $5, $6, $1, $2, $3, $11, $12}}'",
        awk2="-v 'OFS=\t' '-F\t' '{{print $4, $5, $6, $7, $8, $1, $2, $3, $13, $14'}}"
    log:
        "{sample}/logs/analysis/structural_variants/{sample}.annotated_SVs.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/{sample}.annotated_SVs.txt"
    threads: 1
    shell:
        """
        {conda_dir}/bedtools intersect -loj -a {input} -b {genes} 2>> {log} | \
        awk {params.awk1} 2>> {log} | \
        {conda_dir}/bedtools intersect -loj -a stdin -b {genes} 2>> {log} | \
        awk {params.awk2} > {output} 2>> {log}
        """

rule make_igv_bed:
    input:
        "{sample}/structural_variants/{sample}.filter_abnormal_coverage.bedpe"
    output:
        temp1 = temp("{sample}/analysis/structural_variants/{sample}.temp1.bed"),
        temp2 = temp("{sample}/analysis/structural_variants/{sample}.temp2.bed"),
        bed = protected("{sample}/analysis/structural_variants/{sample}.structural_variants.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        cut_1="-f1-3",
        cut_2="-f4-6",
        sort="-k 1,1 -k2,2n -k3,3n"
    log:
        "{sample}/logs/analysis/structural_variants/{sample}.make_igv_bed.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/{sample}.make_igv_bed.txt"
    threads: 1
    shell:
        """
        cut {params.cut_1} {input} > {output.temp1} 2>{log}
        cut {params.cut_2} {input} > {output.temp2} 2>>{log}
        cat {output.temp1} {output.temp2} 2>>{log} | sort {params.sort} > {output.bed} 2>> {log}
        """

rule make_igv_batch:
    input:
        bed = "{sample}/analysis/structural_variants/{sample}.structural_variants.bed",
        bam = "{sample}/mapped/{sample}.phased.bam"
    output:
        protected("{sample}/analysis/structural_variants/{sample}.igv.bat")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        output_folder="{sample}/analysis/structural_variants/plots/igv",
        extend="5000"
    log:
        "{sample}/logs/analysis/structural_variants/{sample}.make_igv_batch.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/{sample}.make_igv_batch.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/make_igv_batch.py -i {input.bed} -b {root_dir}/{input.bam} \
        -o {root_dir}/{params.output_folder} -e {params.extend} > {output} 2>{log}"

# Convert .fa output to BoulderIO format for Primer3
rule convert_to_boulder_io:
    input:
        "{sample}/structural_variants/bp_assembly/combined_corrected_small_window.fa"
    output:
        protected("{sample}/structural_variants/{sample}_SV_boulderIO.txt")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0"
    log:
        "{sample}/logs/structural_variants/{sample}.boulderIO.log"
    benchmark:
        "{sample}/benchmarks/structural_variants/{sample}.boulderIO.txt"
    threads: 1
    shell:
        "{scripts_dir}/convert_to_boulder_io.sh {input} > {output}"
		
# Design primers for SVs
rule primer3:
    input:
        "{sample}/structural_variants/{sample}_SV_boulderIO.txt"
    output:
        protected("{sample}/analysis/structural_variants/{sample}.structural_variant_primers.txt")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0"
    log:
        "{sample}/logs/analysis/structural_variants/{sample}.primer3.log"
    benchmark:
        "{sample}/benchmarks/analysis/structural_variants/{sample}.primer3.txt"
    threads: 1
    shell:
        "{conda_dir}/primer3_core --output={output} --error={log} {input}"
		
