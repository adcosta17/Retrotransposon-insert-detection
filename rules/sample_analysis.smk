

rule sample_for_rate:
    output:
        fastq="Num_Samples_Test/Sample_{mean}_{effect}_.txt"
    threads: 1
    params:
        memory_per_thread="12G",
        script = srcdir("../utils/num_samples.py")
    shell:
        """
        python {params.script} --mean {wildcards.mean} --effect {wildcards.effect} > {output}
        """

