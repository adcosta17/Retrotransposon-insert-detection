import pysam
import argparse
import sys
import csv
import os
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree

from collections import defaultdict

parser = argparse.ArgumentParser( description='Compute the number of insertions for each repeat family, raw and normalized by sample size (Gbp and RPM)')
parser.add_argument('--input', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--fastq-folder', type=str, required=True)
args = parser.parse_args()

# get number of reads and bases in fastq
onlyfiles = [f for f in listdir(args.fastq_folder) if isfile(join(args.fastq_folder, f)) and f.endswith("fastq.gz")]
read_count = 0.0
base_count = 0.0
for file in onlyfiles:
    fq_to_use = args.fastq_folder+"/"+file
    ret = os.popen('zcat %s | paste - - - -| cut -f2 | wc -c -l | awk -v OFS="\n" \'{print "reads: "$1, "bases: "$2-$1}\'' % fq_to_use).read()
    ret_arr = ret.split("\n")
    read_count += float(ret_arr[0].split(":")[1])
    base_count += float(ret_arr[1].split(":")[1])

# Print raw, base normalized and read normalized insert counts
# Header:
# LINE  LINE_GBP    LINE_RPM    SINE    SINE_GBP    SINE_RPM    SVA     SVA_GBP     SVA_RPM     ERV     ERV_GBP     ERV_RPM

counts = defaultdict(int)
# read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if row_args[8] == "PASS" and row_args[19] == "insert":
            # Care about these entries
            if "LINE" in row_args[9]:
                counts["LINE"] += 1
            elif "SINE" in row_args[9]:
                counts["SINE"] += 1
            elif "SVA" in row_args[9]:
                counts["SVA"] += 1
            elif "ERV" in row_args[9]:
                counts["ERV"] += 1

line_mbp = counts["LINE"]/(base_count/1000000000)
line_rpm = counts["LINE"]/(read_count/1000000)
sine_mbp = counts["SINE"]/(base_count/1000000000)
sine_rpm = counts["SINE"]/(read_count/1000000)
sva_mbp = counts["SVA"]/(base_count/1000000000)
sva_rpm = counts["SVA"]/(read_count/1000000)
erv_mbp = counts["ERV"]/(base_count/1000000000)
erv_rpm = counts["ERV"]/(read_count/1000000)

print(args.sample + "\tinsert\t" + str(counts["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+str(counts["SINE"])+
    "\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+str(counts["SVA"])+"\t"+str(sva_mbp)+"\t"+
    str(sva_rpm)+"\t"+str(counts["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm))

counts = defaultdict(int)
# read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if row_args[8] == "PASS" and row_args[19] == "softclip":
            # Care about these entries
            if "LINE" in row_args[9]:
                counts["LINE"] += 1
            elif "SINE" in row_args[9]:
                counts["SINE"] += 1
            elif "SVA" in row_args[9]:
                counts["SVA"] += 1
            elif "ERV" in row_args[9]:
                counts["ERV"] += 1

line_mbp = counts["LINE"]/(base_count/1000000000)
line_rpm = counts["LINE"]/(read_count/1000000)
sine_mbp = counts["SINE"]/(base_count/1000000000)
sine_rpm = counts["SINE"]/(read_count/1000000)
sva_mbp = counts["SVA"]/(base_count/1000000000)
sva_rpm = counts["SVA"]/(read_count/1000000)
erv_mbp = counts["ERV"]/(base_count/1000000000)
erv_rpm = counts["ERV"]/(read_count/1000000)

print(args.sample + "\tsoftclip\t" + str(counts["LINE"])+"\t"+str(line_mbp)+"\t"+str(line_rpm)+"\t"+str(counts["SINE"])+
    "\t"+str(sine_mbp)+"\t"+str(sine_rpm)+"\t"+str(counts["SVA"])+"\t"+str(sva_mbp)+"\t"+
    str(sva_rpm)+"\t"+str(counts["ERV"])+"\t"+str(erv_mbp)+"\t"+str(erv_rpm))

