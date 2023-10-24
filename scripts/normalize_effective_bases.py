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


parser = argparse.ArgumentParser( description='Take in a summary counts file for a sample and normalize it using effective bases')
parser.add_argument('--input', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--tsv-list', required=True)
parser.add_argument('--min-mapq', default=20, type=int)
parser.add_argument('--flank', default=500, type=int)
parser.add_argument('--fraction', default=0.8, type=float)
parser.add_argument('--filter-size', type=str2bool, nargs='?',const=True, default=False)
parser.add_argument('--resolved', type=str2bool, nargs='?',const=True, default=False)
args = parser.parse_args()

reads_to_use = {}
effective_bases_line = 0
effective_bases_alu = 0
effective_bases_sva = 0
effective_bases_erv = 0
# Open bam and get a list of reads that map with at least mapq >= minimum_mapping_quality
header = pysam.AlignmentFile(args.bam).header
for sq in header['SQ']:
    #print(sq['SN'])
    sam_reader = pysam.AlignmentFile(args.bam)
    tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
    for record in tmp_sam_reader:
        if record.mapping_quality >= args.min_mapq:
            if record.query_name not in reads_to_use and not record.is_secondary and not record.is_supplementary:
                reads_to_use[record.query_name] = record.query_alignment_length

line_sizes = []
alu_sizes = []
erv_sizes = []
sva_sizes = []
# Compute averaage insert sizes
for tsv in args.tsv_list.split(','):
    with open(tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if "PASS" in row_args[8] and "ambiguous" not in row:
                count += 1
                if "LINE" in row:
                    line_sizes.append(len(row_args[7]))
                if "SINE" in row:
                    alu_sizes.append(len(row_args[7]))
                if "SVA" in row:
                    sva_sizes.append(len(row_args[7]))
                if "ERV" in row:
                    erv_sizes.append(len(row_args[7]))

#print("Count: "+str(count))
avg_line_size = 0
if len(line_sizes) > 0:
    avg_line_size = sum(line_sizes)/len(line_sizes)
avg_alu_size = 0
if len(alu_sizes) > 0:
    avg_alu_size = sum(alu_sizes)/len(alu_sizes)
avg_erv_size = 0
if len(erv_sizes) > 0:
    avg_erv_size = sum(erv_sizes)/len(erv_sizes)
avg_sva_size = 0
if len(sva_sizes) > 0:
    avg_sva_size = sum(sva_sizes)/len(sva_sizes)

#print(avg_line_size)
#print(avg_alu_size)
#print(avg_sva_size)
#print(avg_erv_size)

for read in reads_to_use:
    e_line = reads_to_use[read] - (2*args.flank + avg_line_size)
    if e_line > 0:
        effective_bases_line += e_line
    e_alu = reads_to_use[read] - (2*args.flank + avg_alu_size)
    if e_alu > 0:
        effective_bases_alu += e_alu
    e_sva = reads_to_use[read] - (2*args.flank + avg_sva_size)
    if e_sva > 0:
        effective_bases_sva += e_sva
    e_erv = reads_to_use[read] - (2*args.flank + avg_erv_size)
    if e_erv > 0:
        effective_bases_erv += e_erv

effective_bases_line = effective_bases_line/1000000000
effective_bases_alu = effective_bases_alu/1000000000
effective_bases_sva = effective_bases_sva/1000000000
effective_bases_erv = effective_bases_erv/1000000000

#print(effective_bases_line)
#print(effective_bases_alu)
#print(effective_bases_sva)
#print(effective_bases_erv)

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
            if args.resolved:
                if "ambiguous" in row_args[9]:
                    continue
            #if seen_before(row_args[3], row_args[0], row_args[1], row_args[2], seen_positions):
            #    continue
            seen_positions[row_args[3]][row_args[0]+":"+row_args[1]+"-"+row_args[2]] = 1
            seen_annotation_count = 0
            # Care about these entries
            counts_pass["All"] += 1
            if "LINE" in row_args[9]:
                if args.filter_size:
                    if len(row_args[7])/args.line_size < args.fraction:
                        continue
                counts_pass["LINE"] += 1
                seen_annotation_count += 1
                if filter_poly_AT(row_args[7]):
                    counts_pass["LINE_PolyA"] += 1
                if "novel" in row:
                    counts_pass["LINE_Novel"] += 1
                else:
                    counts_pass["LINE_HotSpot"] += 1
            if "SINE" in row_args[9]:
                if args.filter_size:
                    if len(row_args[7])/args.alu_size < args.fraction:
                        continue
                counts_pass["SINE"] += 1
                seen_annotation_count += 1
                if filter_poly_AT(row_args[7]):
                    counts_pass["SINE_PolyA"] += 1
                if "novel" in row:
                    counts_pass["SINE_Novel"] += 1
                else:
                    counts_pass["SINE_HotSpot"] += 1
            if "SVA" in row_args[9]:
                if args.filter_size:
                    if len(row_args[7])/args.sva_size < args.fraction:
                        continue
                counts_pass["SVA"] += 1
                seen_annotation_count += 1
                if "novel" in row:
                    counts_pass["SVA_Novel"] += 1
                else:
                    counts_pass["SVA_HotSpot"] += 1
            if "ERV" in row_args[9]:
                if args.filter_size:
                    if len(row_args[7])/args.erv_size < args.fraction:
                        continue
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



#print("Sample\tLINE\tLINE_Effective\tLINE_Novel\tLINE_Novel_Effective\tLINE_HotSpot\tLINE_HotSpot_Effective\tALU\tALU_Effective\tALU_Novel\tALU_Novel_Effective\tALU_HotSpot\tALU_HotSpot_Effective\tSVA\tSVA_Effective\tSVA_Novel\tSVA_Novel_Effective\tSVA_HotSpot\tSVA_HotSpot_Effective\tERV\tERV_Effective\tERV_Novel\tERV_Novel_Effective\tERV_HotSpot\tERV_HotSpot_Effective\tLINE_POLY_A\tLINE_POLY_A_Effective\tALU_POLY_A\tALU_POLY_A_Effective")
out = args.sample+"\t"
out += str(counts_pass["LINE"])+"\t"+str(counts_pass["LINE"]/effective_bases_line)+"\t"
out += str(counts_pass["LINE_Novel"])+"\t"+str(counts_pass["LINE_Novel"]/effective_bases_line)+"\t"
out += str(counts_pass["LINE_HotSpot"])+"\t"+str(counts_pass["LINE_HotSpot"]/effective_bases_line)+"\t"
out += str(counts_pass["SINE"])+"\t"+str(counts_pass["SINE"]/effective_bases_alu)+"\t"
out += str(counts_pass["SINE_Novel"])+"\t"+str(counts_pass["SINE_Novel"]/effective_bases_alu)+"\t"
out += str(counts_pass["SINE_HotSpot"])+"\t"+str(counts_pass["SINE_HotSpot"]/effective_bases_alu)+"\t"
out += str(counts_pass["SVA"])+"\t"+str(counts_pass["SVA"]/effective_bases_sva)+"\t"
out += str(counts_pass["SVA_Novel"])+"\t"+str(counts_pass["SVA_Novel"]/effective_bases_sva)+"\t"
out += str(counts_pass["SVA_HotSpot"])+"\t"+str(counts_pass["SVA_HotSpot"]/effective_bases_sva)+"\t"
out += str(counts_pass["ERV"])+"\t"+str(counts_pass["ERV"]/effective_bases_erv)+"\t"
out += str(counts_pass["ERV_Novel"])+"\t"+str(counts_pass["ERV_Novel"]/effective_bases_erv)+"\t"
out += str(counts_pass["ERV_HotSpot"])+"\t"+str(counts_pass["ERV_HotSpot"]/effective_bases_erv)+"\t"
out += str(counts_pass["LINE_PolyA"])+"\t"+str(counts_pass["LINE_PolyA"]/effective_bases_line)+"\t"
out += str(counts_pass["SINE_PolyA"])+"\t"+str(counts_pass["SINE_PolyA"]/effective_bases_alu)
print(out)
