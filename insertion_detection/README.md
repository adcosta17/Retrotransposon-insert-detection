# Retrotransposon-insert-detection

A Snakemake analysis pipeline for detecting novel retrotransposon insertions in long reads

## Dependencies
- pysam
- snakemake
- pybedtools
- minimap2
- lastal
- samtools
- longshot
- whatshap

## Pipeline Overview

The pipeline takes in one or more fastq's from one or more samples as input. The location of the folder containing the fastq's can be provided in the config. These fastq's are mapped to a provided reference genome, with alignments filtered and altered to recover possible missed insertions. Candidate insertions are computed, with the insert sequence being aligned to a database on known retrotranpsosn concensus sequences. Candidates that map well to a concensus sequence and that are found in an area with high average mapping quality are considered. A control sample can be used to filter hits unique to the sample(s). Reads from the sample can be phased with the assigned haplotypes being used to futher remove mapping artifacts and polymorphic insertions, leaving a set of novel inserts.

## Usage

The pipeline is run using snakemake with a config file specified 

```sh
# Setup:
git clone https://github.com/adcosta17/Retrotransposon-insert-detection.git
cd Retrotransposon-insert-detection/insertion_detection

# Usage: 
snakemake -s pipeline.smk --configfile project_config.yaml <target>

# End-To-End Run:
snakemake -s pipeline.smk --configfile project_config.yaml all_normalized_ava_counts
```
The pipeline requires that the fastq foolder is provided and will generate all required data and files as it proceeds. The final output contains two files for each sample. The first is a list of all insertions that pass our filters plus those that failed annotated with the respective reason. This output is acompanied by a second file that contains the counts of insertions for each repeat family in the sample normalized against the number of reads that mapped with a mapping quality of at least 20. 
