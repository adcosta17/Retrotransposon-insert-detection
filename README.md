# Retrotransposon Insertion Detection

Two Snakemake analysis pipelines for detecting novel retrotransposon insertions in long reads

The first, RTD-Control, requires two or more samples, one of which is noted as a control sample, while the others are test samples where insertions can be detected from. All samples must be replicates or different conditions of the same underlying genome. 

The second, RTD-Diploid, can be run on one or more samples and uses a phased diploid assembly of the underlying baseline genome.

## Dependencies
- pysam
- snakemake
- bgzip
- minimap2
- winnowmap2
- lastal
- samtools
- longshot
- whatshap
- [somrit](https://github.com/adcosta17/somrit)


## Overview of the Pipelines

### RTD-Control

The pipeline takes in one or more fastq's from two or more samples as input. These fastqs are first mapped to a provided reference genome, with insertions extracted from these alignments. A re-alignment process is then run using [somrit](https://github.com/adcosta17/somrit) to reduce the alignment ambiguity related to read insertion positions and recover missed insertions. Insertions are then extracted from the re-aligned BAMs with alignments filtered and altered to recover possible missed insertions. Candidate insertions are computed, with the insert sequence being aligned to a database on known retrotransposon consensus sequences. Candidates that map well to a consensus sequence and that are found in an area with high average mapping quality are considered. Candidates are filtered to remove mapping artifacts and chimeric reads. A control sample is used to filter hits unique to the sample(s), representing variation between the control sample and the reference genome. Reads from the sample are phased with the assigned haplotypes being used to further remove mapping artifacts and polymorphic insertions, leaving a set of novel inserts.

### RTD-Diploid

The pipeline takes in one or more fastqs from one or more samples as input. As with the RTD-Control pipeline the fastqs are first mapped to a provided reference genome. [Somrit](https://github.com/adcosta17/somrit) is then used to identify insertions, perform local realignment to reduce alignment ambiguity, classify insertions to retrotransposon repeat families and filter mapping artifacts and remove polymorphic insertions based on coverage. Insertions passing these filters have their supporting reads aligned to a phased diploid assembly with these alignments parsed to flag any remaining polymorphic insertions, leaving a set of novel insertion calls. 

Insertions from both pipelines can be normalized using the Effective Bases normalization approach that accounts for sequencing depth and read length while normalizing insertion counts per repeat family. 


## Usage

Both pipelines are run using snakemake with a specified config file. 

### Prerequisites

For each sample listed in the config file, a folder should be generated, and a subfolder `fastq` created and populated with the fastqs for that sample. Each fastq should be compressed with bgzip.
The diploid assembly must be generated before running the RTD-Diploid pipeline with the location of the haplotype asssembly fasta files specified in the config. 

```sh
# Setup:
git clone https://github.com/adcosta17/Retrotransposon-insert-detection.git
cd Retrotransposon-insert-detection

# RTD-Control Usage: 
snakemake -s rt_insert_snakefile.smk --configfile project_config.yaml <target>

# RTD-Control End-To-End Run with winnowmap2, first step to generate the data, second to normalize counts:
snakemake -s rt_insert_snakefile.smk --configfile project_config.yaml all_winnow_realign
snakemake -s rt_insert_snakefile.smk --configfile project_config.yaml all_winnow_realign_normalized

# RTD-Diploid Usage: 
snakemake -s somrit_rt_insert_snakefile.smk --configfile project_config.yaml <target>

# RTD-Diploid End-To-End Run, first step to generate the data, second to normalize counts:
snakemake -s somrit_rt_insert_snakefile.smk --configfile project_config.yaml all
snakemake -s somrit_rt_insert_snakefile.smk --configfile project_config.yaml all_assembly_filtered

```
The pipelines require that the fastq folder is provided and they will generate all required data and files as it proceed unless otherwise noted. Both pipelines output a tsv for each sample listing each insertion detected and if it passed the filters and was considered novel, or what filter(s) it failed. It additionally outputs summary stats for each sample, breaking down the passing insertions into their respective repeat families and normalizing them using the Effective Bases normalization approach.  
