#
# Rules for insert detection analysis
#


rule run_get_filtered_candidate_insertions_winnow:
    input:
        bam="{sample}/winnow_filtered_mapped/{sample}.{test}.sorted.bam",
        bam_index="{sample}/winnow_filtered_mapped/{sample}.{test}.sorted.bam.bai"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.tsv"
    threads: 1
    params:
        full_bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        full_bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai",
        merged_reads="{sample}/winnow_filtered_mapped/{sample}.{test}.merged_reads.txt",
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="12G"
    shell:
        "python {params.candidate_insertion_script} --bam {input.bam} --merged {params.merged_reads} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]} > {output}"


rule run_convert_candidate_insertions_to_fasta_winnow:
    input:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.tsv"
    output:
        temp("{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_candidate_insertion_annotation_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.tsv",
        tab="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"

rule get_mapq_ct_filtered_inserts_tsv_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.tsv",
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_insert_script = srcdir("../scripts/filter_low_mapq_centromeric_telomeric_inserts.py"),
    shell:
        "python {params.hc_insert_script} --threads {threads} --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"

rule run_convert_mapq_ct_insertions_to_fasta_winnow:
    input:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv"
    output:
        temp("{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_mapq_ct_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule map_mapq_ct_inserts_winnow:
    input:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.fa"
    output:
        bam=protected("{sample}/mapped_inserts_winnow/{sample}.{test}.sorted.bam")
    params:
        memory_per_thread="3G",
        tmp_loc="{sample}/mapped_inserts_winnow/{sample}.{test}.tmp",
        ref_to_use= get_reference_for_test
    threads: 20
    shell:
        """
        minimap2 -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | samtools sort -T {params.tmp_loc} -o {output}
        """

rule remove_ma_inserts_tsv_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv",
        bam="{sample}/mapped_inserts_winnow/{sample}.{test}.sorted.bam",
        bam_index="{sample}/mapped_inserts_winnow/{sample}.{test}.sorted.bam.bai",
        phased_bam="{sample}/winnow_phased/{sample}.sorted.phased.bam",
        phased_bam_index="{sample}/winnow_phased/{sample}.sorted.phased.bam.bai",
        insert_filter="combined_winnow/combined_multi_sample_insert_location.tsv"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv"
    threads: 1
    params:
        memory_per_thread="48G",
        ma_filter_script = srcdir("../scripts/remove_mapping_artifacts_add_haplotype_tags.py")
    shell:
        "python {params.ma_filter_script} --tsv {input.tsv} --insert-bam {input.bam} --full-bam {input.phased_bam} --insert-filter {input.insert_filter} > {output}"


rule filter_by_reference_sample_single_winnow:
    input:
        file="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.reference_filtered.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("../scripts/remove_insertions_in_reference_samples.py"),
        memory_per_thread="36G"
    shell:
        "python {params.reference_filter_script} --reference-sample {config[reference_tsv_to_filter]} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --reference-repeat-bed {config[repeat_bed]} > {output}"

rule get_multi_sample_filter_xl_winnow:
    input:
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        x_large="combined_winnow/combined_multi_sample_x_large.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_read_analysis --bam-folder winnow_phased --merged-folder winnow_filtered_mapped --max-distance {config[xlarge_window]} --reads-to-exclude {config[exclude]} > {output.x_large}
        """

rule get_multi_sample_filter_large_winnow:
    input:
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        large="combined_winnow/combined_multi_sample_large.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_read_analysis --bam-folder winnow_phased --merged-folder winnow_filtered_mapped --max-distance {config[large_window]} --reads-to-exclude {config[exclude]} > {output.large}
        """

rule get_multi_sample_filter_medium_winnow:
    input:
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        medium="combined_winnow/combined_multi_sample_medium.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_read_analysis --bam-folder winnow_phased --merged-folder winnow_filtered_mapped --max-distance {config[medium_window]} --reads-to-exclude {config[exclude]} > {output.medium}
        """

rule get_multi_sample_filter_small_winnow:
    input:
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        small="combined_winnow/combined_multi_sample_small.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_read_analysis --bam-folder winnow_phased --merged-folder winnow_filtered_mapped --max-distance {config[small_window]} --reads-to-exclude {config[exclude]} > {output.small}
        """

rule reference_filter_and_haplotype_analysis_winnow:
    input:
        file="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv",
        x_large="combined_winnow/combined_multi_sample_x_large.tsv",
        large="combined_winnow/combined_multi_sample_large.tsv",
        medium="combined_winnow/combined_multi_sample_medium.tsv",
        small="combined_winnow/combined_multi_sample_small.tsv",
        low_mapq="combined_winnow/combined_multi_sample_low_mapq.tsv"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv"
    threads: 2
    params:
        reference_filter_script = srcdir("../scripts/ref_filter_and_haplotype_check.py"),
        memory_per_thread="8G"
    shell:
        "python {params.reference_filter_script} --mapq0-filter {input.low_mapq} --multi-reference-filter-xlarge-tsv {input.x_large} --multi-reference-filter-large-tsv {input.large} --multi-reference-filter-medium-tsv {input.medium} --multi-reference-filter-small-tsv {input.small} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --sample {wildcards.sample} --ref-sample {config[reference_sample]} --reference-repeat-bed {config[repeat_bed]} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"


rule filter_inserts_within_other_inserts_winnow:
    input:
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.tsv", sample=config["samples"]),
        expand("{sample}/winnow_filtered_mapped/{sample}.all.merged_reads.txt", sample=config["samples"])
    output:
        inserts="combined_winnow/combined_multi_sample_insert_location.tsv",
        low_mapq="combined_winnow/combined_multi_sample_low_mapq.tsv"
    threads: 20
    params:
        script = srcdir("../scripts/filter_inserts_within_insert_sequence.py"),
        region_list = srcdir("../utils/split_regions.bed"),
        memory_per_thread="10G"
    shell:
        "python {params.script} --sample {config[samples_csv]} --folder winnow_read_analysis --merged-folder winnow_filtered_mapped --bam-folder winnow_mapped --bam-suffix .sorted.bam --insert-location-output {output.inserts} --low-mapq-output {output.low_mapq} --regions-list {params.region_list} --threads {threads}"


rule normalize_counts_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv"
    output:
        all="{sample}/winnow_read_analysis/{sample}.{test}.normalized_ava_counts.txt",
        mapped="{sample}/winnow_read_analysis/{sample}.{test}.normalized_ava_mapped_counts.txt",
        distributions="{sample}/winnow_read_analysis/{sample}.{test}.read_len_distributions.txt"
    threads: 1
    params:
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {input.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions}"


rule run_candidate_insertion_updated_annotation_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv",
        tab="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    threads: 5
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/update_ambiguous.py"),
        memory_per_thread="10G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last-tab {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} --sub-family-fasta {config[l1_filter]} > {output}"

rule normalize_counts_ava_updated_annotation_winnow:
    input:
        tsv="{sample}/winnow_read_analysis/{sample}.{test}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    output:
        all="{sample}/winnow_read_analysis/{sample}.{test}.normalized_ava_counts.updated_annoation.txt",
        mapped="{sample}/winnow_read_analysis/{sample}.{test}.normalized_ava_mapped_counts.updated_annoation.txt",
        distributions="{sample}/winnow_read_analysis/{sample}.{test}.read_len_distributions.updated_annoation.txt"
    threads: 1
    params:
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {input.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions} --subfamily"


rule filter_tsvs:
    input:
        old="{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.tsv",
        tsv="{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    output:
        "{sample}/winnow_read_analysis/{sample}.filter_tsvs.tsv"
    threads: 1
    params:
        script= srcdir("../scripts/filter_tsvs.py"),
        memory_per_thread="8G",
        ref_file="Sample1/winnow_read_analysis/Sample1.all.read_insertions.tsv"
    shell:
        "python {params.script} --input {input.old} --reference-sample {params.ref_file} > {output}"



rule find_differences:
    input:
        filtered_tsv="{sample}/winnow_read_analysis/{sample}.filter_tsvs.tsv",
        old_tsv="{sample}/winnow_read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    output:
        "{sample}/winnow_read_analysis/{sample}.differences.tsv"
    threads: 1
    params:
        memory_per_thread="8G",
        script= srcdir("../scripts/find_differences.py")
    shell:
        "python {params.script} --without-filters {input.filtered_tsv} --with-filters {input.old_tsv} > {output}"