#
# Mapping related rules
#

# Add MD to winnowmap calls

rule winnow_add_md:
    input:
        "{sample}/winnow_mapped/{sample}.sorted.bam"
    output:
        "{sample}/winnow_mapped/{sample}.sorted.md.bam"
    params:
        memory_per_thread="12G", 
        ref=get_reference_base
    threads: 1
    shell:
        "samtools calmd {input} {params.ref} -b > {output}"

# Winnowmap Alignmnet

rule winnow_index:
    input:
        get_reference_base
    output:
        config["reference_all"]+"_repetitive_k15.txt"
    params:
        memory_per_thread="48G"
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=15 output merylDB {input}
        {config[meryl_dir]} print greater-than distinct=0.9998 merylDB > {output}
        rm -r merylDB
        """

rule winnow_index_19:
    input:
        get_reference_base
    output:
        config["reference_all"]+"_repetitive_k19.txt"
    params:
        memory_per_thread="48G"
    threads: 1
    shell:
        """
        {config[meryl_dir]} count k=19 output merylDB {input}
        {config[meryl_dir]} print greater-than distinct=0.9998 merylDB > {output}
        rm -r merylDB
        """

rule winnow_align:
    input:
        get_winnow_kmers
    output:
        protected("{sample}/winnow_mapped/{sample}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_base,
        input_fastq=get_sample_fastq,
        tmp_loc="{sample}/winnow_mapped/{sample}.tmp",
        tmp_output="{sample}/winnow_mapped/{sample}.sorted.tmp.bam"
    threads: 20
    shell:
        """
        {config[winnow_dir]} -W {input} -ax {config[minimap2_preset]} -t {threads} {params.ref_to_use} {params.input_fastq} | samtools sort -T {params.tmp_loc} -o {output}
        """

# LRA Alignment

rule lra_index:
    input:
        get_reference_base
    output:
        config["reference_all"]+".gli"
    params:
        memory_per_thread="48G"
    threads: 1
    shell:
        "lra index -{config[lra_preset]}  {input}"

rule lra_align:
    input:
        get_reference_base_lra_index
    output:
        protected("{sample}/lra_mapped/{sample}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_base,
        input_fastq_folder=get_sample_fastq_folder,
        tmp_loc="{sample}/lra_mapped/{sample}.tmp",
        tmp_output="{sample}/lra_mapped/{sample}.sorted.tmp.bam"
    threads: 20
    shell:
        """
        i=0
        for filename in {params.input_fastq_folder}*.fastq.gz; do
            echo $filename
            zcat $filename | lra align -ONT {params.ref_to_use} /dev/stdin -t {threads} -p s | samtools sort -T {params.tmp_loc} -o {params.tmp_output}_$i
            ((i=i+1))
        done
        samtools merge {output} {params.tmp_output}*
        rm {params.tmp_output}*
        """


## Map reads to reference
rule minimap2_align:
    input:
        get_sample_fastq
    output:
        protected("{sample}/mapped/{sample}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_base,
        input_fastq_folder=get_sample_fastq_folder,
        tmp_loc="{sample}/mapped/{sample}.tmp",
        tmp_output="{sample}/mapped/{sample}.sorted.tmp.bam"
    threads: 20
    shell:
        """
        i=0
        for filename in {params.input_fastq_folder}*.fastq.gz; do
            echo $filename
            minimap2 -x {config[minimap2_preset]} -a -2 --MD -t {threads} {params.ref_to_use} $filename | samtools sort -T {params.tmp_loc} -o {params.tmp_output}_$i
            ((i=i+1))
        done
        samtools merge {output} {params.tmp_output}*
        rm {params.tmp_output}*
        """

rule run_get_filtered_bams:
    input:
        bam="{sample}/mapped/{sample}.sorted.bam",
        bam_index="{sample}/mapped/{sample}.sorted.bam.bai",
    output:
        bam=protected("{sample}/filtered_mapped/{sample}.{test}.sorted.bam"),
        merged_reads=protected("{sample}/filtered_mapped/{sample}.{test}.merged_reads.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq_folder=get_sample_fastq_folder,
        memory_per_thread="200G",
        tmp = "{sample}/filtered_mapped/{sample}.{test}.tmp_output.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.bam} --indel-size {config[min_insertion_length]} --fastq-folder {params.input_fastq_folder} --bam-output {params.tmp} --merged {output.merged_reads}
        samtools sort {params.tmp} -T {output.bam}.tmp -o {output.bam}
        rm {params.tmp}
        """

rule run_get_filtered_lra_bams:
    input:
        bam="{sample}/lra_mapped/{sample}.sorted.bam",
        bam_index="{sample}/lra_mapped/{sample}.sorted.bam.bai",
    output:
        bam=protected("{sample}/lra_filtered_mapped/{sample}.{test}.sorted.bam"),
        merged_reads=protected("{sample}/lra_filtered_mapped/{sample}.{test}.merged_reads.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq_folder=get_sample_fastq_folder,
        memory_per_thread="200G",
        tmp = "{sample}/lra_filtered_mapped/{sample}.{test}.tmp_output.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.bam} --indel-size {config[min_insertion_length]} --fastq-folder {params.input_fastq_folder} --bam-output {params.tmp} --merged {output.merged_reads}
        samtools sort {params.tmp} -T {output.bam}.tmp -o {output.bam}
        rm {params.tmp}
        """

rule run_get_filtered_winnow_bams:
    input:
        bam="{sample}/winnow_mapped/{sample}.sorted.bam",
        bam_index="{sample}/winnow_mapped/{sample}.sorted.bam.bai",
    output:
        bam=protected("{sample}/winnow_filtered_mapped/{sample}.{test}.sorted.bam"),
        merged_reads=protected("{sample}/winnow_filtered_mapped/{sample}.{test}.merged_reads.txt")
    threads: 1
    params:
        merge_alignments_script = srcdir("../scripts/merge_split_alignments.py"),
        input_fastq_folder=get_sample_fastq_folder,
        memory_per_thread="200G",
        tmp = "{sample}/winnow_filtered_mapped/{sample}.{test}.tmp_output.bam"
    shell:
        """
        python {params.merge_alignments_script} --reference-gap-minimum 100 --read-to-reference-bam {input.bam} --indel-size {config[min_insertion_length]} --fastq-folder {params.input_fastq_folder} --bam-output {params.tmp} --merged {output.merged_reads}
        samtools sort {params.tmp} -T {output.bam}.tmp -o {output.bam}
        rm {params.tmp}
        """

rule parasail_align:
    input:
        fa="{sample}/read_analysis/{sample}.{test}.read_insertions.fa"
    output:
        mapped="{sample}/read_analysis/{sample}.{test}.read_insertions.mapped_with_parasail.tsv"
    threads: 1
    params:
        script=srcdir("../scripts/align_parasail.py"),
        memory_per_thread="25G"
    shell:
        "python {params.script} --inserts {input.fa} --targets {config[repbase]} > {output.mapped}"

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



rule make_lastdb:
    input:
        config["repbase"]
    output:
        config["repbase"] + ".lastdb.suf"
    threads: 1
    params:
        memory_per_thread="16G"
    shell:
        "lastdb {config[repbase]}.lastdb {input}"


rule map_insertion_sequences_last:
    input:
        fa="{base}.fa",
        db=config["repbase"] + ".lastdb.suf"
    output:
        "{base}.mapped_to_repbase.last.tab"
    threads: 8
    params:
        memory_per_thread="12G"
    shell:
        "lastal -P {threads} -f tab -r1 -a1 -b1 {config[repbase]}.lastdb {input.fa} > {output}"


rule run_convert_passing_insertions_to_fasta:
    input:
        "{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv"
    output:
        temp("{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.fa")
    threads: 1
    params:
        candidate_insertion_conversion_script = srcdir("../scripts/convert_passing_inserts_to_fasta.py"),
        memory_per_thread="8G"
    shell:
        "python {params.candidate_insertion_conversion_script} --input {input} > {output}"


rule minimap2_align_passing_inserts:
    input:
        "{sample}/read_analysis/{sample}.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.fa"
    output:
        protected("{sample}/inserts_mapped/{sample}.sorted.bam")
    params:
        memory_per_thread="10G",
        ref_to_use= get_reference_base,
        tmp_loc="{sample}/mapped/{sample}.tmp"
    threads: 4
    shell:
        "minimap2 -x map-ont -a -2 --MD -t {threads} {params.ref_to_use} {input} | {config[samtools_dir]} sort -T {params.tmp_loc} -o {output} "

rule merge_aligned_passing_inserts:
    input:
        expand("{sample}/inserts_mapped/{sample}.sorted.bam", sample=config["samples"])
    output:
        protected("combined/combined_inserts_mapped.sorted.bam")
    params:
        memory_per_thread="10G"
    threads: 1
    shell:
        """
        samtools merge {output} {input}
        samtools index {output}
        """


