import pysam
import argparse
import sys
import csv
import os
import math
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree

from collections import defaultdict

def get_filters_to_exclude(pf_annotation, filters, repbase_annotation):
    # Base case if pass
    if pf_annotation == "PASS":
        return True
    # pf_annotation may have one or more filters
    # Need to check and elimiate each one found in filters
    # If they are all found in filters, then we can return True.
    # Otherwise there are things in
    if filters == "" or filters == "NA":
        return False
    pf_annotation_split = pf_annotation.split(",")
    count = 0
    for item in pf_annotation_split:
        if item in filters:
            # explicit check to ensure we don't include things that aren't mapped to repbase
            if item == "mapped_fraction":
                if repbase_annotation != "no_repbase_mapping":
                    count += 1
            else:
                count += 1
    if count == len(pf_annotation_split):
        return True
    return False


def seen_before(read, chrom, start, end, seen_positions):
    if read in seen_positions:
        for position in seen_positions[read]:
            if chrom == position.split(':')[0]:
                if abs(int(start) - int(position.split(":")[1].split("-")[0])) <= 100:
                    return True
    return False


parser = argparse.ArgumentParser( description='Compute the number of insertions for each repeat family, raw and normalized by sample size (Gbp and RPM)')
parser.add_argument('--input', required=True)
parser.add_argument('--output-all', required=True)
parser.add_argument('--output-mapped', required=True)
parser.add_argument('--output-distributions', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--fastq-folder', type=str, required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--minimum-mapping-qual', type=int, default=20)
parser.add_argument('--filters-to-exclude', default="NA")
args = parser.parse_args()

reads_to_use = {}
# Open bam and get a list of reads that map with at least mapq >= minimum_mapping_quality
header = pysam.AlignmentFile(args.bam).header
for sq in header['SQ']:
    sam_reader = pysam.AlignmentFile(args.bam)
    tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
    for record in tmp_sam_reader:
        if record.mapping_quality >= args.minimum_mapping_qual:
            reads_to_use[record.query_name] = 1

# get number of reads and bases in fastq
total_read_count = 0.0
total_base_count = 0.0
total_read_lens = defaultdict(int)
mapped_read_count = 0.0
mapped_base_count = 0.0
mapped_read_lens = defaultdict(int)
onlyfiles = [f for f in listdir(args.fastq_folder) if isfile(join(args.fastq_folder, f)) and f.endswith("fastq.gz")]
for file in onlyfiles:
    with pysam.FastaFile(filename = args.fastq_folder+ "/" + file, filepath_index_compressed = args.fastq_folder+ "/" +file + ".gzi") as in_fq:
        for read in in_fq.references:
            seq = in_fq.fetch(read)
            total_read_count += 1
            total_base_count += len(seq)
            total_read_lens[math.floor(len(seq)/1000.0)*1000] += 1
            if read in reads_to_use:
                mapped_read_count += 1
                mapped_base_count += len(seq)
                mapped_read_lens[math.floor(len(seq)/1000.0)*1000] += 1

#print(args.sample + "\t" + str(read_count) + "\t" + str(base_count))
# Print raw, base normalized and read normalized insert counts
# Header:
# LINE  LINE_GBP    LINE_RPM    SINE    SINE_GBP    SINE_RPM    SVA     SVA_GBP     SVA_RPM     ERV     ERV_GBP     ERV_RPM

# Keep track of all the positives we see. If two entries on a read within 100bp of each other combine into one
seen_positions = defaultdict(lambda: defaultdict(int))

counts_pass = defaultdict(int)

#read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if row_args[8] == "PASS":
            if seen_before(row_args[3], row_args[0], row_args[1], row_args[2], seen_positions):
                continue
            seen_positions[row_args[3]][row_args[0]+":"+row_args[1]+"-"+row_args[2]] = 1
            # Care about these entries
            if "LINE" in row_args[9]:
                counts_pass["LINE"] += 1
            elif "SINE" in row_args[9]:
                counts_pass["SINE"] += 1
            elif "SVA" in row_args[9]:
                counts_pass["SVA"] += 1
            elif "ERV" in row_args[9]:
                counts_pass["ERV"] += 1
            elif "ambiguous" in row_args[9]:
                counts_pass["ambiguous"] += 1

mbp_denom = total_base_count/1000000000
rpm_denom = total_read_count/1000000

with open(args.output_all, 'w') as out_all:
    line_mbp = counts_pass["LINE"]/mbp_denom
    line_rpm = counts_pass["LINE"]/rpm_denom
    sine_mbp = counts_pass["SINE"]/mbp_denom
    sine_rpm = counts_pass["SINE"]/rpm_denom
    ambiguous_mbp = counts_pass["ambiguous"]/mbp_denom
    ambiguous_rpm = counts_pass["ambiguous"]/rpm_denom
    sva_mbp = counts_pass["SVA"]/mbp_denom
    sva_rpm = counts_pass["SVA"]/rpm_denom
    erv_mbp = counts_pass["ERV"]/mbp_denom
    erv_rpm = counts_pass["ERV"]/rpm_denom
    out_all.write(args.sample+"\t"+
        str(counts_pass["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+
        str(counts_pass["SINE"])+"\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+
        str(counts_pass["SVA"])+"\t"+str(sva_mbp)+"\t"+str(sva_rpm)+"\t"+
        str(counts_pass["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm)+"\t"+
        str(counts_pass["ambiguous"])+"\t"+str(ambiguous_mbp)+"\t"+str(ambiguous_rpm)+"\n")

mbp_denom = mapped_base_count/1000000000
rpm_denom = mapped_read_count/1000000

with open(args.output_mapped, 'w') as out_mapped:
    line_mbp = counts_pass["LINE"]/mbp_denom
    line_rpm = counts_pass["LINE"]/rpm_denom
    sine_mbp = counts_pass["SINE"]/mbp_denom
    sine_rpm = counts_pass["SINE"]/rpm_denom
    ambiguous_mbp = counts_pass["ambiguous"]/mbp_denom
    ambiguous_rpm = counts_pass["ambiguous"]/rpm_denom
    sva_mbp = counts_pass["SVA"]/mbp_denom
    sva_rpm = counts_pass["SVA"]/rpm_denom
    erv_mbp = counts_pass["ERV"]/mbp_denom
    erv_rpm = counts_pass["ERV"]/rpm_denom
    out_mapped.write(args.sample+"\t"+
        str(counts_pass["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+
        str(counts_pass["SINE"])+"\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+
        str(counts_pass["SVA"])+"\t"+str(sva_mbp)+"\t"+str(sva_rpm)+"\t"+
        str(counts_pass["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm)+"\t"+
        str(counts_pass["ambiguous"])+"\t"+str(ambiguous_mbp)+"\t"+str(ambiguous_rpm)+"\n")

with open(args.output_distributions, 'w') as out_dist:
    #out_dist.write("")
    for i in total_read_lens:
        out_dist.write(str(int(i))+"\t"+str(int(i+999))+"\t"+str(total_read_lens[i])+"\t"+str(mapped_read_lens[i])+"\n")

