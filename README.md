# Retrotransposon-insert-detection

A Snakemake analysis pipeline for detecting novel retrotransposon insertions in long reads. 

Split into two parts.

A. Insertion Detection
  - Detecting novel retrotransposon insertions in long reads
B. Insertion Validation
  - Methods to validate calls made by the Insertion Detection pipeline in A

## Dependencies
- pysam
- snakemake
- pybedtools
- minimap2
- lastal
- samtools
