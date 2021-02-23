#
# Rules for Phasing and Haplotyping
#


# Split bam by region

rule split_bam:
    input:
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/phased/{sample}.{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools view -b {input.bam} {wildcards.regions} > {output}"

# Combined region bams for each sample
rule combined_split_bams:
    input:
        bam = expand("{sample}/phased/{sample}.{{regions}}.sorted.bam", sample=config["samples"]),
        bai = expand("{sample}/phased/{sample}.{{regions}}.sorted.bam.bai", sample=config["samples"])
    output:
        temp("combined/combined_{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools merge {output} {input.bam}"

# Genotype the sample with longshot
rule phase_bam:
    input:
        bam = "combined/combined_{regions}.sorted.bam",
        bai = "combined/combined_{regions}.sorted.bam.bai"
    output:
        temp("combined/combined_{regions}.vcf")
    params:
        memory_per_thread="48G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "longshot -c 3 -C 10000 -r {wildcards.regions} --bam {input.bam} --ref {params.ref_to_use} --out {output}"

# Zip and index the phased bams
rule phase_bam_index:
    input:
        "combined/combined_{regions}.vcf"
    output:
        zip = temp("combined/combined_{regions}.vcf.gz"),
        index = temp("combined/combined_{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="6G"
    threads: 2
    shell:
        """
        bgzip -c {input} > {output.zip}
        tabix -p vcf {output.zip}
        """

# Merge the vcf files
rule merge_vcf:
    input:
        zip = expand("combined/combined_{regions}.vcf.gz", regions=regions),
        index = expand("combined/combined_{regions}.vcf.gz.tbi", regions=regions)
    output:
        zip = protected("combined/combined.vcf.gz"),
        index = protected("combined/combined.vcf.gz.tbi")
    params:
        memory_per_thread="4G"
    threads: 2
    shell:
        """
        bcftools concat -a -O v {input.zip} | bcftools sort -O z -o {output.zip}
        tabix -p vcf {output.zip}
        """

# Tag the bams with the haplotype information
rule haplotag_bam:
    input:
        vcf = "combined/combined.vcf.gz",
        tbi = "combined/combined.vcf.gz.tbi",
        bam = "{sample}/mapped/{sample}.sorted.bam",
        bai = "{sample}/mapped/{sample}.sorted.bam.bai"
    output:
        bam = protected("{sample}/phased/{sample}.sorted.phased.bam")
    params:
        memory_per_thread="16G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "whatshap haplotag --ignore-read-groups --reference {params.ref_to_use} {input.vcf} {input.bam} | samtools sort -o {output.bam}"



rule split_bam_winnow:
    input:
        bam = "{sample}/winnow_mapped/{sample}.sorted.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/winnow_phased/{sample}.{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools view -b {input.bam} {wildcards.regions} > {output}"

# Combined region bams for each sample
rule combined_split_bams_winnow:
    input:
        bam = expand("{sample}/winnow_phased/{sample}.{{regions}}.sorted.bam", sample=config["samples"]),
        bai = expand("{sample}/winnow_phased/{sample}.{{regions}}.sorted.bam.bai", sample=config["samples"])
    output:
        temp("combined_winnow/combined_{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools merge {output} {input.bam}"

# Genotype the sample with longshot
rule phase_bam_winnow:
    input:
        bam = "combined_winnow/combined_{regions}.sorted.bam",
        bai = "combined_winnow/combined_{regions}.sorted.bam.bai"
    output:
        temp("combined_winnow/combined_{regions}.vcf")
    params:
        memory_per_thread="48G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "longshot -c 3 -C 10000 -r {wildcards.regions} --bam {input.bam} --ref {params.ref_to_use} --out {output}"

# Zip and index the phased bams
rule phase_bam_index_winnow:
    input:
        "combined_winnow/combined_{regions}.vcf"
    output:
        zip = temp("combined_winnow/combined_{regions}.vcf.gz"),
        index = temp("combined_winnow/combined_{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="6G"
    threads: 2
    shell:
        """
        bgzip -c {input} > {output.zip}
        tabix -p vcf {output.zip}
        """

# Merge the vcf files
rule merge_vcf_winnow:
    input:
        zip = expand("combined_winnow/combined_{regions}.vcf.gz", regions=regions),
        index = expand("combined_winnow/combined_{regions}.vcf.gz.tbi", regions=regions)
    output:
        zip = protected("combined_winnow/combined.vcf.gz"),
        index = protected("combined_winnow/combined.vcf.gz.tbi")
    params:
        memory_per_thread="4G"
    threads: 2
    shell:
        """
        bcftools concat -a -O v {input.zip} | bcftools sort -O z -o {output.zip}
        tabix -p vcf {output.zip}
        """

# Tag the bams with the haplotype information
rule haplotag_bam_winnow:
    input:
        vcf = "combined_winnow/combined.vcf.gz",
        tbi = "combined_winnow/combined.vcf.gz.tbi",
        bam = "{sample}/winnow_mapped/{sample}.sorted.bam",
        bai = "{sample}/winnow_mapped/{sample}.sorted.bam.bai"
    output:
        bam = protected("{sample}/winnow_phased/{sample}.sorted.phased.bam")
    params:
        memory_per_thread="16G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "whatshap haplotag --ignore-read-groups --reference {params.ref_to_use} {input.vcf} {input.bam} | samtools sort -o {output.bam}"



rule split_bam_lra:
    input:
        bam = "{sample}/lra_mapped/{sample}.sorted.bam",
        bai = "{sample}/lra_mapped/{sample}.sorted.bam.bai"
    output:
        temp("{sample}/lra_phased/{sample}.{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools view -b {input.bam} {wildcards.regions} > {output}"

# Combined region bams for each sample
rule combined_split_bams_lra:
    input:
        bam = expand("{sample}/lra_phased/{sample}.{{regions}}.sorted.bam", sample=config["samples"]),
        bai = expand("{sample}/lra_phased/{sample}.{{regions}}.sorted.bam.bai", sample=config["samples"])
    output:
        temp("combined_lra/combined_{regions}.sorted.bam")
    params:
        memory_per_thread="2G"
    threads: 1
    shell:
        "samtools merge {output} {input.bam}"

# Genotype the sample with longshot
rule phase_bam_lra:
    input:
        bam = "combined_lra/combined_{regions}.sorted.bam",
        bai = "combined_lra/combined_{regions}.sorted.bam.bai"
    output:
        temp("combined_lra/combined_{regions}.vcf")
    params:
        memory_per_thread="48G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "longshot -c 3 -C 10000 -r {wildcards.regions} --bam {input.bam} --ref {params.ref_to_use} --out {output}"

# Zip and index the phased bams
rule phase_bam_index_lra:
    input:
        "combined_lra/combined_{regions}.vcf"
    output:
        zip = temp("combined_lra/combined_{regions}.vcf.gz"),
        index = temp("combined_lra/combined_{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="6G"
    threads: 2
    shell:
        """
        bgzip -c {input} > {output.zip}
        tabix -p vcf {output.zip}
        """

# Merge the vcf files
rule merge_vcf_lra:
    input:
        zip = expand("combined_lra/combined_{regions}.vcf.gz", regions=regions),
        index = expand("combined_lra/combined_{regions}.vcf.gz.tbi", regions=regions)
    output:
        zip = protected("combined_lra/combined.vcf.gz"),
        index = protected("combined_lra/combined.vcf.gz.tbi")
    params:
        memory_per_thread="4G"
    threads: 2
    shell:
        """
        bcftools concat -a -O v {input.zip} | bcftools sort -O z -o {output.zip}
        tabix -p vcf {output.zip}
        """

# Update the bit flags of the lra mapped bams to work with whatshap haplotag
rule update_bam_lra:
    input:
        bam = "{sample}/lra_mapped/{sample}.sorted.bam",
        bai = "{sample}/lra_mapped/{sample}.sorted.bam.bai"
    output:
        bam = "{sample}/lra_mapped/{sample}.sorted.updated.bam"
    params:
        memory_per_thread="48G",
        script = srcdir("../scripts/update_lra_bams.py")
    threads: 1
    shell:
        #"python {params.script} --bam {input.bam} --output-bam {output.bam}"
        "cp {input.bam} {output.bam}"

# Tag the bams with the haplotype information
rule haplotag_bam_lra:
    input:
        vcf = "combined_lra/combined.vcf.gz",
        tbi = "combined_lra/combined.vcf.gz.tbi",
        bam = "{sample}/lra_mapped/{sample}.sorted.updated.bam",
        bai = "{sample}/lra_mapped/{sample}.sorted.updated.bam.bai"
    output:
        bam = protected("{sample}/lra_phased/{sample}.sorted.phased.bam")
    params:
        memory_per_thread="16G",
        ref_to_use= get_reference_default
    threads: 2
    shell:
        "whatshap haplotag --ignore-read-groups --reference {params.ref_to_use} {input.vcf} {input.bam} | samtools sort -o {output.bam}"
