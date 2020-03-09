# Retrotransposon-insert-detection - Validation 

A Snakemake pipeline for assesing the sensitivity and precision of insertion detection

## Dependencies
- pysam
- snakemake
- pybedtools
- minimap2
- lastal
- samtools

## Pipeline Overview

The validation pipeline is designed to be run after high confidence calls have been made using the insertion detection pipeline and a modifided reference genome. Regions of known retrotransposon sequence in the reference genome are deleted to generate a modified genome. The insertion detection pipeline is called using this genome to generate high confidence calls. The positions where deletions have been made in the original genome should be specified in a BED file.
Reads are aligned to the normal reference genome, and the alignments compared ot the BED file of deleted regions to generate a truth set.
These reads are then compared to the results of the Insertion Detection pipeline and inserts classified as true positive, false positive and false negatives based on if an insertion was called or not. A reference tsv can be specified to act as a control sample.

## Usage

The pipeline is run using snakemake with a config file specified 

```sh
# Setup:
git clone https://github.com/adcosta17/Retrotransposon-insert-detection.git
cd Retrotransposon-insert-detection/validation

# Usage: 
snakemake -s pipeline_validation.smk --configfile project_config.yaml <target>

# Get truth set per chromosome
snakemake -s pipeline_validation.smk --configfile project_config.yaml all_spanning
# Combine truth set for sample
snakemake -s pipeline_validation.smk --configfile project_config.yaml all_combined
# Assess Sensitivity and Precision of insert detection
python annotate_tsv_tp_fp.py --insert-tsv <sample>/read_analysis/<sample>.<test>.read_insertions.repbase_annotated.high_confidence.tsv \
							 --soft-clip-tsv <sample>/read_analysis/<sample>.<test>.read_soft_clipped.repbase_annotated.high_confidence.tsv \
							 --insert-truth <sample>/results/<sample>.<test>.read_insertions.truth_set.tsv \
							 --soft-clip-truth <sample>/results/<sample>.<test>.read_soft_clipped.truth_set.tsv \
							 --reference-tsv control_sample.tsv\
							 --output-tsv <sample>.<test>.annotated.output.tsv
```

The reference tsv can be a control sample that has been aligned to the same modified genome using the insertion detection pipeline, or an assembly from the same source genome as the reads. 
