#
# Rules for SV calling
#
import os
import re

include: "helpers.smk"

# Sniffles

rule sniffles:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.md.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.md.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.sniffles.bedpe")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        min_reads = "1",
        sniffles_dir = config["sniffles_dir"]
    threads: 10
    shell:
        "{params.sniffles_dir} -s {params.min_reads} --max_num_splits 1 -t {threads} -m {input.bam} -b {output} -n -1"


rule sniffles_vcf:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.md.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.md.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.sniffles.vcf")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        min_reads = "1",
        sniffles_dir = config["sniffles_dir"]
    threads: 10
    shell:
        "{params.sniffles_dir} -s {params.min_reads} --max_num_splits 1 -t {threads} -m {input.bam} -v {output} -n -1"


rule cute_sv:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.md.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.md.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.cuteSV.vcf")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        min_reads="1",
        ref=get_reference_base
    threads: 10
    shell:
        "cuteSV {input.bam} {params.ref} {output} {wildcards.sample}/structural_variants/ --min_support {params.min_reads} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads {threads} --report_readid --retain_work_dir"


rule nanovar:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.md.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.md.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.nanovar.pass.vcf")
    params:
        memory_per_thread="10G",
        run_time="1:0:0:0",
        ref=get_reference_base
    threads: 10
    shell:
        "nanovar -t {threads} -f hg38 {input.bam} {params.ref} {wildcards.sample}/structural_variants/"


rule sniffles_vcf_to_fasta:
    input:
        "{sample}/structural_variants/{sample}.sniffles.vcf"
    output:
        temp("{sample}/structural_variants/{sample}.sniffles.fa")
    threads: 1
    params:
        script = srcdir("../scripts/convert_vcf_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.script} --input {input} > {output}"

rule cutesv_vcf_to_fasta:
    input:
        "{sample}/structural_variants/{sample}.cuteSV.vcf"
    output:
        temp("{sample}/structural_variants/{sample}.cuteSV.fa")
    threads: 1
    params:
        script = srcdir("../scripts/convert_vcf_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.script} --input {input} > {output}"

rule annotate_cuteSV:
    input:
        vcf="{sample}/structural_variants/{sample}.cuteSV.vcf",
        tab="{sample}/structural_variants/{sample}.cuteSV.mapped_to_repbase.last.tab"
    output:
        "{sample}/structural_variants/{sample}.cuteSV.annotated.vcf"
    threads: 1
    params:
        annotation_script = srcdir("../scripts/annotate_vcf_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.annotation_script} --input {input.vcf} --last {input.tab} > {output}"

rule annotate_sniffles:
    input:
        vcf="{sample}/structural_variants/{sample}.sniffles.vcf",
        tab="{sample}/structural_variants/{sample}.sniffles.mapped_to_repbase.last.tab"
    output:
        "{sample}/structural_variants/{sample}.sniffles.annotated.vcf"
    threads: 1
    params:
        annotation_script = srcdir("../scripts/annotate_vcf_from_repbase.py"),
        memory_per_thread="48G"
    shell:
        "python {params.annotation_script} --input {input.vcf} --last {input.tab} > {output}"


rule get_mapq_ct_filtered_inserts_vcf_sniffles:
    input:
        vcf="{sample}/structural_variants/{sample}.sniffles.annotated.vcf",
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.vcf"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_insert_script = srcdir("../scripts/filter_low_mapq_centromeric_telomeric_inserts_vcf.py"),
    shell:
        "python {params.hc_insert_script} --threads {threads} --vcf {input.vcf} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"



rule get_mapq_ct_filtered_inserts_vcf_cutesv:
    input:
        vcf="{sample}/structural_variants/{sample}.cuteSV.annotated.vcf",
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.vcf"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_insert_script = srcdir("../scripts/filter_low_mapq_centromeric_telomeric_inserts_vcf.py"),
    shell:
        "python {params.hc_insert_script} --threads {threads} --vcf {input.vcf} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"


rule run_convert_mapq_ct_insertions_to_fasta_sniffles:
    input:
        "{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.vcf"
    output:
        temp("{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_vcf_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule map_mapq_ct_inserts_sniffles:
    input:
        "{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.fa"
    output:
        bam=protected("{sample}/mapped_inserts_sniffles/{sample}.sniffles.sorted.bam")
    params:
        memory_per_thread="3G",
        tmp_loc="{sample}/mapped_inserts_sniffles/{sample}.sniffles.tmp",
        ref_to_use= get_reference_base
    threads: 20
    shell:
        """
        minimap2 -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | samtools sort -T {params.tmp_loc} -o {output}
        """

rule remove_ma_inserts_sniffles:
    input:
        vcf="{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.vcf",
        bam="{sample}/mapped_inserts_sniffles/{sample}.sniffles.sorted.bam",
        bam_index="{sample}/mapped_inserts_sniffles/{sample}.sniffles.sorted.bam.bai",
        full_bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/structural_variants/{sample}.sniffles.annotated.mapq_ct_filtered.ma_filtered.vcf"
    threads: 1
    params:
        memory_per_thread="48G",
        ma_filter_script = srcdir("../scripts/remove_mapping_artifacts_add_haplotype_tags_vcf.py")
    shell:
        "python {params.ma_filter_script} --vcf {input.vcf} --insert-bam {input.bam} --full-bam {input.full_bam} > {output}"


rule run_convert_mapq_ct_insertions_to_fasta_cuteSV:
    input:
        "{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.vcf"
    output:
        temp("{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_vcf_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule map_mapq_ct_inserts_cuteSV:
    input:
        "{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.fa"
    output:
        bam=protected("{sample}/mapped_inserts_cuteSV/{sample}.cuteSV.sorted.bam")
    params:
        memory_per_thread="3G",
        tmp_loc="{sample}/mapped_inserts_sniffles/{sample}.cuteSV.tmp",
        ref_to_use= get_reference_base
    threads: 20
    shell:
        """
        minimap2 -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | samtools sort -T {params.tmp_loc} -o {output}
        """

rule remove_ma_inserts_cuteSV:
    input:
        vcf="{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.vcf",
        bam="{sample}/mapped_inserts_cuteSV/{sample}.cuteSV.sorted.bam",
        bam_index="{sample}/mapped_inserts_cuteSV/{sample}.cuteSV.sorted.bam.bai",
        full_bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/structural_variants/{sample}.cuteSV.annotated.mapq_ct_filtered.ma_filtered.vcf"
    threads: 1
    params:
        memory_per_thread="48G",
        ma_filter_script = srcdir("../scripts/remove_mapping_artifacts_add_haplotype_tags_vcf.py")
    shell:
        "python {params.ma_filter_script} --vcf {input.vcf} --insert-bam {input.bam} --full-bam {input.full_bam} > {output}"


rule graph_align:
    input:
        get_sample_fastq
    output:
        "{sample}/structural_variants/{sample}.gaf"
    threads: 20
    params:
        memory_per_thread="10G",
        input_fastq_folder=get_sample_fastq_folder,
        output_prefix="{sample}/structural_variants/{sample}.tmp"
    shell:
        """
        i=0
        for filename in {params.input_fastq_folder}*.fastq.gz; do
            GraphAligner -g {config[illumina_assembly]} -f $filename -a {params.output_prefix}_$i.gaf -x dbg -t {threads}
            cat {params.output_prefix}_$i.gaf >> {output}
            rm {params.output_prefix}_$i.gaf
            ((i=i+1))
        done
        """


rule path_sv:
    input:
        gaf="{sample}/structural_variants/{sample}.gaf",
        tsv="{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv"
    output:
        "{sample}/structural_variants/{sample}.path_sv.tsv"
    threads: 1
    params:
        memory_per_thread="64G",
        tmp="{sample}/structural_variants/{sample}.tmp"
    shell:
        """
        grep PASS {input.tsv} > {params.tmp}
        {config[path_sv]}/PathSV -i {input.gaf} -g {config[illumina_assembly]} -n {input.gaf} -t dbg -m filter -p {params.tmp} -w 250 > {output}
        """

rule path_sv_filter:
    input:
        "{sample}/structural_variants/{sample}.path_sv.tsv"
    output:
        "{sample}/structural_variants/{sample}.path_sv.filtered.tsv"
    threads: 1
    params:
        memory_per_thread="12G"
    shell:
        """
        python {config[path_sv]}/GroupPathCalls.py --input-tsv {input} > {output}
        """

