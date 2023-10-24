#
# Rules for insert detection analysis
#

def get_ref(wildcards):
    return config["reference_all"]

def get_base_dir(wildcards):
    return config['base_dir']

def get_fastq(wildcards):
   return config['base_dir'] + "/"+ wildcards.sample + "/fastq/" + wildcards.sample + ".fastq.gz"

def get_all_fastq(wildcards):
    fastq_list = ""
    for sample in config["samples_all"]:
        fastq_list += config['base_dir']+sample+"/fastq/"+sample+".fastq.gz,"
    return fastq_list[:-1]

def get_bam_list(wildcards):
    bam_list = config['base_dir']+wildcards.sample+"/winnow_mapped/"+wildcards.sample+".sorted.bam"
    return bam_list
    
def get_all_bam_list(wildcards):
    bam_list = ""
    for sample in config["samples_all"]:
        bam_list += config['base_dir']+sample+"/winnow_mapped/"+sample+".sorted.bam,"
    bam_list = bam_list[:-1]
    return bam_list

def get_output_bam_list(wildcards):
    bam_list = config['base_dir']+wildcards.sample+"/winnow_realign/"+wildcards.sample+".tmp.bam"
    return bam_list

def get_tsv_list(wildcards):
    tsv_list = ""
    for sample in config["samples_all"]:
        tsv_list += config['base_dir']+sample+"/winnow_realign/"+sample+".tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_final_tsv_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/winnow_realign_read_analysis/"+sample+".read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list

def get_final_tsv_list_novel(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+".read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list


def get_no_polymorphic_tsv_list(wildcards):
    tsv_list = ""
    for sample in config["samples"]:
        tsv_list += config['base_dir']+sample+"/winnow_realign_read_analysis/"+sample+".read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv,"
    tsv_list = tsv_list[:-1]
    return tsv_list


# Add rules for ReAlignment

rule extract_inserts:
    input:
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        tsv="{sample}/winnow_realign/{sample}.tsv",
        merged_reads="{sample}/winnow_realign/{sample}.merged_reads.txt"
    threads: 20
    params:
        memory_per_thread="15G",
        fastq=get_fastq,
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-tsv {params.base_dir}/{output.tsv} --output-merged {params.base_dir}/{output.merged_reads} --fastq-file {params.fastq} --threads {threads}
        cd {params.base_dir}
        """
 

rule realign_inserts:
    input:
        bam=expand("{sample}/winnow_mapped/{sample}.sorted.bam", sample=config["samples_all"]),
        tsv=expand("{sample}/winnow_realign/{sample}.tsv", sample=config["samples_all"])
    output:
        bam="combined_realign_winnow_all/{chrom}.sorted.bam",
        tsv="combined_realign_winnow_all/{chrom}.realigned.tsv"
    threads: 10
    params:
        memory_per_thread="25G",
        fastq=get_all_fastq,
        bams=get_all_bam_list,
        tsvs=get_tsv_list,
        base_dir=get_base_dir,
        ref=get_ref
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py realign --bam-list {params.bams} --tsv-list {params.tsvs} --fastq-list {params.fastq} --output-dir {params.base_dir}/combined_realign_winnow_all --tsv-prefix realigned --bam-prefix sorted --reference-genome {params.ref} --threads {threads} --filter-depth --max-insert-size 8000 --max-depth 35 --chromosome-list {wildcards.chrom} --only-realign
        cd {params.base_dir}
        """

rule merge_realign_inserts:
    input:
        bam= "{sample}/winnow_mapped/{sample}.sorted.bam",
        bai= "{sample}/winnow_mapped/{sample}.sorted.bam.bai",
        cbams = expand("combined_realign_winnow_all/{chrom}.sorted.bam",chrom=config["chroms"]),
        cbams_index = expand("combined_realign_winnow_all/{chrom}.sorted.bam.bai",chrom=config["chroms"]),
        tsv = expand("combined_realign_winnow_all/{chrom}.realigned.tsv",chrom=config["chroms"])
    output:
        bam="{sample}/winnow_realign/{sample}.tmp.bam"
    threads: 1
    params:
        memory_per_thread="24G",
        base_dir=get_base_dir,
        out_bam_list=get_output_bam_list,
        bams=get_bam_list,
        ref=get_ref,
        tsv="{sample}/winnow_realign/{sample}.all.realigned.tsv"
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py merge --output-dir {params.base_dir}/combined_realign_winnow_all --bam-list {params.bams} --tsv-prefix realigned --bam-prefix sorted --output-bam-list {params.out_bam_list} --output-tsv-list {params.base_dir}/{params.tsv}
        cd {params.base_dir}
        """

rule sort_realign_bams:
    input:
        bam="{sample}/winnow_realign/{sample}.tmp.bam"
    output:
        bam="{sample}/winnow_realign/{sample}.sorted.bam"
    threads: 1
    params:
        memory_per_thread="24G"
    shell:
        """
        samtools sort {input.bam} > {output.bam}
        rm {input.bam}
        """

rule run_get_filtered_candidate_insertions_winnow_realign:
    input:
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.tsv"
    threads: 20
    params:
        merged_reads="{sample}/winnow_realign/{sample}.merged_reads.after.txt",
        candidate_insertion_script = srcdir("../scripts/get_candidate_insertions.py"),
        memory_per_thread="15G",
        fastq=get_fastq,
        base_dir=get_base_dir
    shell:
        """
        cd {config[somrit_dir]}
        python somrit.py extract --bam {params.base_dir}/{input.bam} --output-tsv {params.base_dir}/{output} --output-merged {params.base_dir}/{params.merged_reads} --fastq-file {params.fastq} --threads {threads} --min-insertion-length {config[min_insertion_length]} --min-mapq {config[min_mapq]} --min-detected-inclusion-length {config[min_detected_inclusion_length]}
        cd {params.base_dir}
        """


rule run_convert_candidate_insertions_to_fasta_winnow_realign:
    input:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.tsv"
    output:
        temp("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_candidate_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule run_candidate_insertion_annotation_winnow_realign:
    input:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.tsv",
        tab="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/annotate_insertions_from_repbase.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} > {output}"

rule get_mapq_ct_filtered_inserts_tsv_winnow_realign:
    input:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.tsv",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv"
    threads: 20
    params:
        memory_per_thread="6G",
        hc_insert_script = srcdir("../scripts/filter_low_mapq_centromeric_telomeric_inserts.py"),
    shell:
        "python {params.hc_insert_script} --threads {threads} --tsv {input.tsv} --window-size {config[hc_window_size]} --min-mapq-fraction {config[min_mapq_fraction]} --read-to-reference-bam {input.bam} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"

rule run_convert_mapq_ct_insertions_to_fasta_winnow_realign:
    input:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv"
    output:
        temp("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_mapq_ct_insertions_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"

rule map_mapq_ct_inserts_winnow_realign:
    input:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.fa"
    output:
        bam=protected("{sample}/mapped_inserts_realign_winnow/{sample}.sorted.bam")
    params:
        memory_per_thread="3G",
        tmp_loc="{sample}/mapped_inserts_realign_winnow/{sample}.tmp",
        ref_to_use= get_ref
    threads: 20
    shell:
        """
        minimap2 -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | samtools sort -T {params.tmp_loc} -o {output}
        """

rule remove_ma_inserts_tsv_winnow_realign:
    input:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv",
        bam="{sample}/mapped_inserts_realign_winnow/{sample}.sorted.bam",
        bam_index="{sample}/mapped_inserts_realign_winnow/{sample}.sorted.bam.bai",
        phased_bam="{sample}/winnow_realign_phased/{sample}.sorted.phased.bam",
        phased_bam_index="{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai",
        insert_filter="combined_realign_winnow/combined_multi_sample_insert_location.tsv"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv"
    threads: 1
    params:
        memory_per_thread="48G",
        ma_filter_script = srcdir("../scripts/remove_mapping_artifacts_add_haplotype_tags.py")
    shell:
        "python {params.ma_filter_script} --tsv {input.tsv} --insert-bam {input.bam} --full-bam {input.phased_bam} --insert-filter {input.insert_filter} > {output}"


rule filter_by_reference_sample_single_winnow_realign:
    input:
        file="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.reference_filtered.tsv"
    threads: 1
    params:
        reference_filter_script = srcdir("../scripts/remove_insertions_in_reference_samples.py"),
        memory_per_thread="36G"
    shell:
        "python {params.reference_filter_script} --reference-sample {config[reference_tsv_to_filter]} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --reference-repeat-bed {config[repeat_bed]} > {output}"

rule get_multi_sample_filter_xl_winnow_realign:
    input:
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        x_large="combined_realign_winnow/combined_multi_sample_x_large.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_realign_read_analysis --suffix .read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv --bam-folder winnow_realign_phased --merged-folder winnow_realign --merged-suffix .merged_reads.txt --max-distance {config[xlarge_window]} --reads-to-exclude {config[exclude]} > {output.x_large}
        """

rule get_multi_sample_filter_large_winnow_realign:
    input:
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        large="combined_realign_winnow/combined_multi_sample_large.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_realign_read_analysis --suffix .read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv --bam-folder winnow_realign_phased --merged-folder winnow_realign --merged-suffix .merged_reads.txt --max-distance {config[large_window]} --reads-to-exclude {config[exclude]} > {output.large}
        """

rule get_multi_sample_filter_medium_winnow_realign:
    input:
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        medium="combined_realign_winnow/combined_multi_sample_medium.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_realign_read_analysis --suffix .read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv --bam-folder winnow_realign_phased --merged-folder winnow_realign --merged-suffix .merged_reads.txt --max-distance {config[medium_window]} --reads-to-exclude {config[exclude]} > {output.medium}
        """

rule get_multi_sample_filter_small_winnow_realign:
    input:
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv", sample=config["samples"])
    output:
        small="combined_realign_winnow/combined_multi_sample_small.tsv"
    threads: 12
    params:
        gen_script = srcdir("../scripts/generate_multi_sample_tsv.py"),
        memory_per_thread="6G"
    shell:
        """
        python {params.gen_script} --sample {config[samples_csv]} --threads {threads} --folder winnow_realign_read_analysis --suffix .read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv --bam-folder winnow_realign_phased --merged-folder winnow_realign --merged-suffix .merged_reads.txt --max-distance {config[small_window]} --reads-to-exclude {config[exclude]} > {output.small}
        """

rule reference_filter_and_haplotype_analysis_winnow_realign:
    input:
        file="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv",
        x_large="combined_realign_winnow/combined_multi_sample_x_large.tsv",
        large="combined_realign_winnow/combined_multi_sample_large.tsv",
        medium="combined_realign_winnow/combined_multi_sample_medium.tsv",
        small="combined_realign_winnow/combined_multi_sample_small.tsv",
        low_mapq="combined_realign_winnow/combined_multi_sample_low_mapq.tsv"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv"
    threads: 2
    params:
        reference_filter_script = srcdir("../scripts/ref_filter_and_haplotype_check.py"),
        memory_per_thread="8G"
    shell:
        "python {params.reference_filter_script} --mapq0-filter {input.low_mapq} --multi-reference-filter-xlarge-tsv {input.x_large} --multi-reference-filter-large-tsv {input.large} --multi-reference-filter-medium-tsv {input.medium} --multi-reference-filter-small-tsv {input.small} --max-distance {config[max_distance_to_reference_insertion]} --input {input.file} --sample {wildcards.sample} --ref-sample {config[reference_sample]} --reference-repeat-bed {config[repeat_bed]} --telomere-filter {config[telomere_filter]} --centromere-filter {config[centromere_filter]} > {output}"


rule filter_inserts_within_other_inserts_winnow_realign:
    input:
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam", sample=config["samples"]),
        expand("{sample}/winnow_realign_phased/{sample}.sorted.phased.bam.bai", sample=config["samples"]),
        expand("{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.tsv", sample=config["samples"]),
        expand("{sample}/winnow_realign/{sample}.merged_reads.txt", sample=config["samples"])
    output:
        inserts="combined_realign_winnow/combined_multi_sample_insert_location.tsv",
        low_mapq="combined_realign_winnow/combined_multi_sample_low_mapq.tsv"
    threads: 20
    params:
        script = srcdir("../scripts/filter_inserts_within_insert_sequence.py"),
        region_list = srcdir("../utils/split_regions.bed"),
        memory_per_thread="10G"
    shell:
        "python {params.script} --sample {config[samples_csv]} --folder winnow_realign_read_analysis --suffix .read_insertions.repbase_annotated.mapq_ct_filtered.tsv --merged-reads .merged_reads.txt --merged-folder winnow_realign --bam-folder winnow_realign --bam-suffix .sorted.bam --insert-location-output {output.inserts} --low-mapq-output {output.low_mapq} --regions-list {params.region_list} --threads {threads}"


rule normalize_counts_winnow_realign:
    output:
        all="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_counts.txt",
        mapped="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_mapped_counts.txt",
        distributions="{sample}/winnow_realign_read_analysis/{sample}.read_len_distributions.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv",
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions}"


rule normalize_counts_winnow_realign_effective_bases:
    output:
        counts="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_mapped_counts_effective_bases.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv",
        normalize_script = srcdir("../scripts/normalize_effective_bases.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name,
        tsv_list=get_final_tsv_list
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --tsv-list {params.tsv_list} > {output.counts}"

rule run_candidate_insertion_updated_annotation_winnow_realign:
    input:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv",
        tab="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    threads: 1
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/update_ambiguous.py"),
        memory_per_thread="24G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last-tab {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} --sub-family-fasta {config[l1_filter]} --alu-tab {config[alu_tab]} > {output}"

rule run_candidate_insertion_updated_annotation_winnow_realign_no_polymorphic:
    input:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv",
        tab="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.mapped_to_repbase.last.tab"
    output:
        "{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.updated_annotation.tsv"
    threads: 5
    params:
        candidate_insertion_annotation_script = srcdir("../scripts/update_ambiguous.py"),
        memory_per_thread="10G"
    shell:
        "python {params.candidate_insertion_annotation_script} --input {input.tsv} --last-tab {input.tab} --min-mapped-fraction {config[insert_mapping_fraction]} --sub-family-fasta {config[l1_filter]} --alu-tab {config[alu_tab]} > {output}"

rule normalize_counts_ava_updated_annotation_winnow_realign:
    output:
        all="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_counts.updated_annoation.txt",
        mapped="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_mapped_counts.updated_annoation.txt",
        distributions="{sample}/winnow_realign_read_analysis/{sample}.read_len_distributions.updated_annoation.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv",
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions} --subfamily"


rule normalize_counts_winnow_realign_updated_annotation_effective_bases:
    output:
        counts="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_mapped_counts_updated_annotation_effective_bases.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv",
        normalize_script = srcdir("../scripts/normalize_effective_bases.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name,
        tsv_list=get_final_tsv_list
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --tsv-list {params.tsv_list} --resolved 1 > {output.counts}"


rule normalize_counts_ava_updated_annotation_winnow_realign_no_polymorphic:
    output:
        all="{sample}/winnow_realign_read_analysis/{sample}.normalized_no_polymorphic_ava_counts.updated_annoation.txt",
        mapped="{sample}/winnow_realign_read_analysis/{sample}.normalized_no_polymorphic_ava_mapped_counts.updated_annoation.txt",
        distributions="{sample}/winnow_realign_read_analysis/{sample}.read_len_distributions_no_polymorphic.updated_annoation.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.updated_annotation.tsv",
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions} --subfamily"


rule normalize_counts_winnow_realign_updated_annotation_effective_bases_no_polymorphic:
    output:
        counts="{sample}/winnow_realign_read_analysis/{sample}.normalized_ava_mapped_counts_no_polymorphic_updated_annotation_effective_bases.txt"
    threads: 1
    params:
        tsv="{sample}/winnow_realign_read_analysis/{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.updated_annotation.tsv",
        normalize_script = srcdir("../scripts/normalize_effective_bases.py"),
        memory_per_thread="16G",
        bam="{sample}/winnow_realign/{sample}.sorted.bam",
        bam_index="{sample}/winnow_realign/{sample}.sorted.bam.bai",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name,
        tsv_list=get_no_polymorphic_tsv_list
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --tsv-list {params.tsv_list} --resolved 1 > {output.counts}"


rule normalize_counts_winnow_realign_updated_annotation_effective_bases_no_polymorphic_new:
    output:
        counts="{sample}.normalized_ava_mapped_counts_no_polymorphic_updated_annotation_effective_bases.novel.txt"
    threads: 1
    params:
        tsv="{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.updated_annotation.tsv",
        normalize_script = srcdir("../scripts/normalize_effective_bases.py"),
        memory_per_thread="16G",
        bam="{sample}.sorted.bam",
        sample_name=get_sample_name,
        tsv_list=get_final_tsv_list_novel
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --tsv-list {params.tsv_list} --resolved 1 > {output.counts}"

rule normalize_counts_winnow_realign_updated_annotation_effective_bases_new:
    output:
        counts="{sample}.normalized_ava_mapped_counts_updated_annotation_effective_bases.novel.txt"
    threads: 1
    params:
        tsv="{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv",
        normalize_script = srcdir("../scripts/normalize_effective_bases.py"),
        memory_per_thread="16G",
        bam="{sample}.sorted.bam",
        sample_name=get_sample_name,
        tsv_list=get_final_tsv_list_novel
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --tsv-list {params.tsv_list} --resolved 1 > {output.counts}"

rule normalize_counts_ava_updated_annotation_winnow_realign_new:
    output:
        all="{sample}.normalized_ava_counts.updated_annoation.txt",
        mapped="{sample}.normalized_ava_mapped_counts.updated_annoation.txt",
        distributions="{sample}.read_len_distributions.updated_annoation.txt"
    threads: 1
    params:
        tsv="{sample}.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv",
        normalize_script = srcdir("../scripts/normalize_ava_calls.py"),
        memory_per_thread="16G",
        bam="{sample}.sorted.bam",
        fastq_folder=get_sample_fastq_folder,
        sample_name=get_sample_name
    shell:
        "python {params.normalize_script} --input {params.tsv} --bam {params.bam} --sample {params.sample_name} --fastq-folder {params.fastq_folder} --output-all {output.all} --output-mapped {output.mapped} --output-distributions {output.distributions} --subfamily"

