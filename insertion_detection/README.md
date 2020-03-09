# Retrotransposon-insert-detection

A Snakemake analysis pipeline for detecting novel retrotransposon insertions in long reads

## Dependencies
- pysam
- snakemake
- pybedtools
- minimap2
- lastal
- samtools

## Pipeline Overview

The pipeline takes in one or more fastq's from one or more samples as input. The location of the folder containing the fastq's can be provided in the config. These fastq's are mapped to a provided reference genome, with alignments filtered and altered to recover possible missed insertions. Candidate insertions are computed, with the insert sequence being aligned to a database on known retrotranpsosn concensus sequences. Candidates that map well to a concensus sequence and that are found in an area with high average mapping quality are considered. A control sample can be used to filter hits unique to the sample(s) leaving a set of novel inserts.

## Usage

The pipeline is run using snakemake with a config file specified 

```sh
# Setup:
git clone https://github.com/adcosta17/Retrotransposon-insert-detection.git
cd Retrotransposon-insert-detection/insertion_detection

# Usage: 
snakemake -s pipeline.smk --configfile project_config.yaml <target>

# End-To-End Run:
# A. Generate tsv files with candidate inserts and softclips for each sample
snakemake -s pipeline.smk --configfile project_config.yaml all_high_confidence_tsv
# B. Filter against hits in control sample and generate normalized insert and softclip counts
snakemake -s pipeline.smk --configfile project_config.yaml all_normalized_counts
```

Note the end to end run is split into two parts. This is required if the control is one of the samples in the run. The high confidence tsv is from the control sample is used as a filter for the other samples and is generated in the first step. The second step does this filtering and generated normalized counts
