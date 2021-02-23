#
# Rules for insert detection analysis
#


rule run_get_filtered_candidate_insertions:
    input:
        bam="{sample}/filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/filtered_mapped/{sample}.{test}.sorted.bam.bai"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.tsv"
    threads: 1
    params:
        full_bam="{sample}/mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
        merged_reads="{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt",
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="12G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --merged {params.merged_reads} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"


rule run_convert_candidate_insertions_to_fasta:
    input:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.tsv"
    output:
        temp("{sample}/read_analysis/{sample}.{test}.read_insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_candidate_insertion_annotation:
    input:
        tsv="{sample}/read_analysis/{sample}.{test}.read_insertions.tsv",
        tab="{sample}/read_analysis/{sample}.{test}.read_insertions.mapped_with_parasail.tsv"
    output:
        "{sample}/read_analysis/{sample}.{test}.read_insertions.parasail_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_insertions_from_parasail.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --alignment {input.tab} > {output}"


