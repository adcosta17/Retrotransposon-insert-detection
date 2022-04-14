# Retrotransposon-insert-detection

A Snakemake analysis pipeline for detecting novel retrotransposon insertions in long reads

## Dependencies
- pysam
- snakemake
- pybedtools
- minimap2
- winnowmap2
- lastal
- samtools
- longshot
- whatshap

## Pipeline Overview

The pipeline takes in one or more fastq's from one or more samples as input. These fastqs are first mapped to a provided reference genome, with insertions extracted from these alignments. A re-alignment process is then run using [somrit](https://github.com/adcosta17/somrit) to reduce the ambiguity of insertion position between reads and recover missed insertions. Insertions are then extracted from the re-aligned BAMs with alignments filtered and altered to recover possible missed insertions. Candidate insertions are computed, with the insert sequence being aligned to a database on known retrotranpsosn concensus sequences. Candidates that map well to a concensus sequence and that are found in an area with high average mapping quality are considered. Candidates are filtered to remove mapping artifacts and chimeric reads. A control sample is to filter hits unique to the sample(s), representing variation between the control sample and the reference genome. Reads from the sample are be phased with the assigned haplotypes being used to futher remove mapping artifacts and polymorphic insertions, leaving a set of novel inserts. 

## Usage

The pipeline is run using snakemake with a config file specified. The recomended aligner is winnowmap2. The pipeline is designed to work with minimap2 if requested, but realignment is only currently supported with winnowmap2.

### Prerequisites

For each sample listed in the config file, a folder should be generated, and a subfolder `fastq` created and populated with the fastqs for that sample. 

```sh
# Setup:
git clone https://github.com/adcosta17/Retrotransposon-insert-detection.git
cd Retrotransposon-insert-detection

# Usage: 
snakemake -s rt_insert_snakefile.smk --configfile project_config.yaml <target>

# End-To-End Run with winnowmap2:
snakemake -s rt_insert_snakefile.smk --configfile project_config.yaml all_normalized_ava_counts_combined_winnow
```
The pipeline requires that the fastq folder is provided and will generate all required data and files as it proceeds. The pipeline outputs a tsv for each sample listing each insertion detected and if it passed the filters and was considered novel, or what filter(s) it failed. It additionally outputs summary stats for each sample, breaking down the passing insertions into their respective repeat families and normalizing them by the number of mapped reads and mapped bases.  
