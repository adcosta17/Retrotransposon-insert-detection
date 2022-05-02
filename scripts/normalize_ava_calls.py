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

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

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


def filter_poly_AT(insert_seq):
    has_poly_AT = False
    start = insert_seq[:50]
    end = insert_seq[len(insert_seq)-50:]
    pA = "A"*5
    pT = "T"*5
    if pA in start or pT in start or pA in end or pT in end:
        has_poly_AT = True
    return has_poly_AT

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
parser.add_argument('--subfamily', type=str2bool, nargs='?',const=True, default=False)
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
    try:
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
    except OSError:
        pass

#print(args.sample + "\t" + str(read_count) + "\t" + str(base_count))
# Print raw, base normalized and read normalized insert counts
# Header:
# LINE  LINE_GBP    LINE_RPM    SINE    SINE_GBP    SINE_RPM    SVA     SVA_GBP     SVA_RPM     ERV     ERV_GBP     ERV_RPM

# Keep track of all the positives we see. If two entries on a read within 100bp of each other combine into one
seen_positions = defaultdict(lambda: defaultdict(int))

counts_pass = defaultdict(int)
counts_pass_subfamily = defaultdict(int)

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
            seen_annotation_count = 0
            # Care about these entries
            counts_pass["All"] += 1
            if "LINE" in row_args[9]:
                counts_pass["LINE"] += 1
                seen_annotation_count += 1
                if filter_poly_AT(row_args[7]):
                    counts_pass["LINE_PolyA"] += 1
                if "novel" in row:
                    counts_pass["LINE_Novel"] += 1
                else:
                    counts_pass["LINE_HotSpot"] += 1
            if "SINE" in row_args[9]:
                counts_pass["SINE"] += 1
                seen_annotation_count += 1
                if filter_poly_AT(row_args[7]):
                    counts_pass["SINE_PolyA"] += 1
                if "novel" in row:
                    counts_pass["SINE_Novel"] += 1
                else:
                    counts_pass["SINE_HotSpot"] += 1
            if "SVA" in row_args[9]:
                counts_pass["SVA"] += 1
                seen_annotation_count += 1
                if "novel" in row:
                    counts_pass["SVA_Novel"] += 1
                else:
                    counts_pass["SVA_HotSpot"] += 1
            if "ERV" in row_args[9]:
                counts_pass["ERV"] += 1
                seen_annotation_count += 1
                if "novel" in row:
                    counts_pass["ERV_Novel"] += 1
                else:
                    counts_pass["ERV_HotSpot"] += 1
            if "ambiguous" in row_args[9] and seen_annotation_count == 0:
                counts_pass["ambiguous"] += 1
                if "novel" in row:
                    counts_pass["ambiguous_Novel"] += 1
                else:
                    counts_pass["ambiguous_HotSpot"] += 1
            if args.subfamily and "Subfamily" in row_args[9]:
                sub = row_args[9].split("Subfamily:")[1]
                sub_count = 0
                if "L1H" in sub:
                    counts_pass_subfamily["L1H"] += 1
                if "L1M" in sub:
                    counts_pass_subfamily["L1M"] += 1
                if "L1P" in sub:
                    counts_pass_subfamily["L1P"] += 1
                if "L1" in sub:
                    counts_pass_subfamily["L1"] += 1


mbp_denom = total_base_count/1000000000
rpm_denom = total_read_count/1000000

with open(args.output_all, 'w') as out_all:
    all_mbp = counts_pass["All"]/mbp_denom
    all_rpm = counts_pass["All"]/rpm_denom
    line_mbp = counts_pass["LINE"]/mbp_denom
    line_rpm = counts_pass["LINE"]/rpm_denom
    line_novel_mbp = counts_pass["LINE_Novel"]/mbp_denom
    line_novel_rpm = counts_pass["LINE_Novel"]/rpm_denom
    line_hotspot_mbp = counts_pass["LINE_HotSpot"]/mbp_denom
    line_hotspot_rpm = counts_pass["LINE_HotSpot"]/rpm_denom
    line_polya_mbp = counts_pass["LINE_PolyA"]/mbp_denom
    line_polya_rpm = counts_pass["LINE_PolyA"]/rpm_denom
    sine_mbp = counts_pass["SINE"]/mbp_denom
    sine_rpm = counts_pass["SINE"]/rpm_denom
    sine_novel_mbp = counts_pass["SINE_Novel"]/mbp_denom
    sine_novel_rpm = counts_pass["SINE_Novel"]/rpm_denom
    sine_hotspot_mbp = counts_pass["SINE_HotSpot"]/mbp_denom
    sine_hotspot_rpm = counts_pass["SINE_HotSpot"]/rpm_denom
    sine_polya_mbp = counts_pass["SINE_PolyA"]/mbp_denom
    sine_polya_rpm = counts_pass["SINE_PolyA"]/rpm_denom
    ambiguous_mbp = counts_pass["ambiguous"]/mbp_denom
    ambiguous_rpm = counts_pass["ambiguous"]/rpm_denom
    ambiguous_novel_mbp = counts_pass["ambiguous_Novel"]/mbp_denom
    ambiguous_novel_rpm = counts_pass["ambiguous_Novel"]/rpm_denom
    ambiguous_hotspot_mbp = counts_pass["ambiguous_HotSpot"]/mbp_denom
    ambiguous_hotspot_rpm = counts_pass["ambiguous_HotSpot"]/rpm_denom
    sva_mbp = counts_pass["SVA"]/mbp_denom
    sva_rpm = counts_pass["SVA"]/rpm_denom
    sva_novel_mbp = counts_pass["SVA_Novel"]/mbp_denom
    sva_novel_rpm = counts_pass["SVA_Novel"]/rpm_denom
    sva_hotspot_mbp = counts_pass["SVA_HotSpot"]/mbp_denom
    sva_hotspot_rpm = counts_pass["SVA_HotSpot"]/rpm_denom
    erv_mbp = counts_pass["ERV"]/mbp_denom
    erv_rpm = counts_pass["ERV"]/rpm_denom
    erv_novel_mbp = counts_pass["ERV_Novel"]/mbp_denom
    erv_novel_rpm = counts_pass["ERV_Novel"]/rpm_denom
    erv_hotspot_mbp = counts_pass["ERV_HotSpot"]/mbp_denom
    erv_hotspot_rpm = counts_pass["ERV_HotSpot"]/rpm_denom
    l1_mbp = counts_pass_subfamily["L1"]/mbp_denom
    l1_rpm = counts_pass_subfamily["L1"]/rpm_denom
    l1p_mbp = counts_pass_subfamily["L1P"]/mbp_denom
    l1p_rpm = counts_pass_subfamily["L1P"]/rpm_denom
    l1m_mbp = counts_pass_subfamily["L1M"]/mbp_denom
    l1m_rpm = counts_pass_subfamily["L1M"]/rpm_denom
    l1h_mbp = counts_pass_subfamily["L1H"]/mbp_denom
    l1h_rpm = counts_pass_subfamily["L1H"]/rpm_denom
    out_all.write(args.sample+"\t"+
        str(counts_pass["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+
        str(counts_pass["LINE_Novel"])+"\t"+str(line_novel_mbp)+"\t"+str(line_novel_rpm)+"\t"+
        str(counts_pass["LINE_HotSpot"])+"\t"+str(line_hotspot_mbp)+"\t"+str(line_hotspot_rpm)+"\t"+
        str(counts_pass["SINE"])+"\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+
        str(counts_pass["SINE_Novel"])+"\t"+str(sine_novel_mbp)+"\t"+str(sine_novel_rpm)+"\t"+
        str(counts_pass["SINE_HotSpot"])+"\t"+str(sine_hotspot_mbp)+"\t"+str(sine_hotspot_rpm)+"\t"+
        str(counts_pass["SVA"])+"\t"+str(sva_mbp)+"\t"+str(sva_rpm)+"\t"+
        str(counts_pass["SVA_Novel"])+"\t"+str(sva_novel_mbp)+"\t"+str(sva_novel_rpm)+"\t"+
        str(counts_pass["SVA_HotSpot"])+"\t"+str(sva_hotspot_mbp)+"\t"+str(sva_hotspot_rpm)+"\t"+
        str(counts_pass["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm)+"\t"+
        str(counts_pass["ERV_Novel"])+"\t"+str(erv_novel_mbp)+"\t"+str(erv_novel_rpm)+"\t"+
        str(counts_pass["ERV_HotSpot"])+"\t"+str(erv_hotspot_mbp)+"\t"+str(erv_hotspot_rpm)+"\t"+
        str(counts_pass["ambiguous"])+"\t"+str(ambiguous_mbp)+"\t"+str(ambiguous_rpm)+"\t"+
        str(counts_pass["ambiguous_Novel"])+"\t"+str(ambiguous_novel_mbp)+"\t"+str(ambiguous_novel_rpm)+"\t"+
        str(counts_pass["ambiguous_HotSpot"])+"\t"+str(ambiguous_hotspot_mbp)+"\t"+str(ambiguous_hotspot_rpm)+"\t")
    #if args.subfamily:
    #    out_all.write(str(counts_pass_subfamily["L1"])+"\t"+str(l1_mbp)+"\t"+str(l1_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1H"])+"\t"+str(l1h_mbp)+"\t"+str(l1h_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1M"])+"\t"+str(l1m_mbp)+"\t"+str(l1m_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1P"])+"\t"+str(l1p_mbp)+"\t"+str(l1p_rpm)+"\t")
    out_all.write(str(counts_pass["LINE_PolyA"])+"\t"+str(line_polya_mbp)+"\t"+str(line_polya_rpm)+"\t"+
        str(counts_pass["SINE_PolyA"])+"\t"+str(sine_polya_mbp)+"\t"+str(sine_polya_rpm)+"\n")

mbp_denom = mapped_base_count/1000000000
rpm_denom = mapped_read_count/1000000

with open(args.output_mapped, 'w') as out_mapped:
    all_mbp = counts_pass["All"]/mbp_denom
    all_rpm = counts_pass["All"]/rpm_denom
    line_mbp = counts_pass["LINE"]/mbp_denom
    line_rpm = counts_pass["LINE"]/rpm_denom
    line_novel_mbp = counts_pass["LINE_Novel"]/mbp_denom
    line_novel_rpm = counts_pass["LINE_Novel"]/rpm_denom
    line_hotspot_mbp = counts_pass["LINE_HotSpot"]/mbp_denom
    line_hotspot_rpm = counts_pass["LINE_HotSpot"]/rpm_denom
    line_polya_mbp = counts_pass["LINE_PolyA"]/mbp_denom
    line_polya_rpm = counts_pass["LINE_PolyA"]/rpm_denom
    sine_mbp = counts_pass["SINE"]/mbp_denom
    sine_rpm = counts_pass["SINE"]/rpm_denom
    sine_novel_mbp = counts_pass["SINE_Novel"]/mbp_denom
    sine_novel_rpm = counts_pass["SINE_Novel"]/rpm_denom
    sine_hotspot_mbp = counts_pass["SINE_HotSpot"]/mbp_denom
    sine_hotspot_rpm = counts_pass["SINE_HotSpot"]/rpm_denom
    sine_polya_mbp = counts_pass["SINE_PolyA"]/mbp_denom
    sine_polya_rpm = counts_pass["SINE_PolyA"]/rpm_denom
    ambiguous_mbp = counts_pass["ambiguous"]/mbp_denom
    ambiguous_rpm = counts_pass["ambiguous"]/rpm_denom
    ambiguous_novel_mbp = counts_pass["ambiguous_Novel"]/mbp_denom
    ambiguous_novel_rpm = counts_pass["ambiguous_Novel"]/rpm_denom
    ambiguous_hotspot_mbp = counts_pass["ambiguous_HotSpot"]/mbp_denom
    ambiguous_hotspot_rpm = counts_pass["ambiguous_HotSpot"]/rpm_denom
    sva_mbp = counts_pass["SVA"]/mbp_denom
    sva_rpm = counts_pass["SVA"]/rpm_denom
    sva_novel_mbp = counts_pass["SVA_Novel"]/mbp_denom
    sva_novel_rpm = counts_pass["SVA_Novel"]/rpm_denom
    sva_hotspot_mbp = counts_pass["SVA_HotSpot"]/mbp_denom
    sva_hotspot_rpm = counts_pass["SVA_HotSpot"]/rpm_denom
    erv_mbp = counts_pass["ERV"]/mbp_denom
    erv_rpm = counts_pass["ERV"]/rpm_denom
    erv_novel_mbp = counts_pass["ERV_Novel"]/mbp_denom
    erv_novel_rpm = counts_pass["ERV_Novel"]/rpm_denom
    erv_hotspot_mbp = counts_pass["ERV_HotSpot"]/mbp_denom
    erv_hotspot_rpm = counts_pass["ERV_HotSpot"]/rpm_denom
    l1_mbp = counts_pass_subfamily["L1"]/mbp_denom
    l1_rpm = counts_pass_subfamily["L1"]/rpm_denom
    l1p_mbp = counts_pass_subfamily["L1P"]/mbp_denom
    l1p_rpm = counts_pass_subfamily["L1P"]/rpm_denom
    l1m_mbp = counts_pass_subfamily["L1M"]/mbp_denom
    l1m_rpm = counts_pass_subfamily["L1M"]/rpm_denom
    l1h_mbp = counts_pass_subfamily["L1H"]/mbp_denom
    l1h_rpm = counts_pass_subfamily["L1H"]/rpm_denom
    out_mapped.write(args.sample+"\t"+
        str(counts_pass["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+
        str(counts_pass["LINE_Novel"])+"\t"+str(line_novel_mbp)+"\t"+str(line_novel_rpm)+"\t"+
        str(counts_pass["LINE_HotSpot"])+"\t"+str(line_hotspot_mbp)+"\t"+str(line_hotspot_rpm)+"\t"+
        str(counts_pass["SINE"])+"\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+
        str(counts_pass["SINE_Novel"])+"\t"+str(sine_novel_mbp)+"\t"+str(sine_novel_rpm)+"\t"+
        str(counts_pass["SINE_HotSpot"])+"\t"+str(sine_hotspot_mbp)+"\t"+str(sine_hotspot_rpm)+"\t"+
        str(counts_pass["SVA"])+"\t"+str(sva_mbp)+"\t"+str(sva_rpm)+"\t"+
        str(counts_pass["SVA_Novel"])+"\t"+str(sva_novel_mbp)+"\t"+str(sva_novel_rpm)+"\t"+
        str(counts_pass["SVA_HotSpot"])+"\t"+str(sva_hotspot_mbp)+"\t"+str(sva_hotspot_rpm)+"\t"+
        str(counts_pass["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm)+"\t"+
        str(counts_pass["ERV_Novel"])+"\t"+str(erv_novel_mbp)+"\t"+str(erv_novel_rpm)+"\t"+
        str(counts_pass["ERV_HotSpot"])+"\t"+str(erv_hotspot_mbp)+"\t"+str(erv_hotspot_rpm)+"\t"+
        str(counts_pass["ambiguous"])+"\t"+str(ambiguous_mbp)+"\t"+str(ambiguous_rpm)+"\t"+
        str(counts_pass["ambiguous_Novel"])+"\t"+str(ambiguous_novel_mbp)+"\t"+str(ambiguous_novel_rpm)+"\t"+
        str(counts_pass["ambiguous_HotSpot"])+"\t"+str(ambiguous_hotspot_mbp)+"\t"+str(ambiguous_hotspot_rpm)+"\t")
    #if args.subfamily:
    #    out_mapped.write(str(counts_pass_subfamily["L1"])+"\t"+str(l1_mbp)+"\t"+str(l1_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1H"])+"\t"+str(l1h_mbp)+"\t"+str(l1h_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1M"])+"\t"+str(l1m_mbp)+"\t"+str(l1m_rpm)+"\t"+
    #    str(counts_pass_subfamily["L1P"])+"\t"+str(l1p_mbp)+"\t"+str(l1p_rpm)+"\t")
    out_mapped.write(str(counts_pass["LINE_PolyA"])+"\t"+str(line_polya_mbp)+"\t"+str(line_polya_rpm)+"\t"+
        str(counts_pass["SINE_PolyA"])+"\t"+str(sine_polya_mbp)+"\t"+str(sine_polya_rpm)+"\n")

with open(args.output_distributions, 'w') as out_dist:
    #out_dist.write("")
    for i in total_read_lens:
        out_dist.write(str(int(i))+"\t"+str(int(i+999))+"\t"+str(total_read_lens[i])+"\t"+str(mapped_read_lens[i])+"\n")

