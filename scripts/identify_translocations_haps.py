import pysam
import argparse
from pysam import VariantFile
import sys
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree


parser = argparse.ArgumentParser( description='Parse Minigraph alignments of reads to an augmented graph. Identify Translocations')
parser.add_argument('--fastq', required=True)
parser.add_argument('--hap1', required=True)
parser.add_argument('--hap2', required=True)
parser.add_argument('--min-reads', required=False, default=3)
args = parser.parse_args()

read_lens = {}
with pysam.FastxFile(args.fastq) as fh:
    for entry in fh:
        if len(entry.sequence) > 10000:
            read_lens[entry.name] = len(entry.sequence)

count = 0
hap1_reads = {}
bcf_in = VariantFile(args.hap1)
for rec in bcf_in.fetch():
     if rec.info["SVTYPE"] == "BND":
        hap1_reads[count] = {}
        for read in rec.info["RNAMES"]:
            if read in read_lens:
                hap1_reads[count][read] = 1
        count += 1

hap2_reads = {}
bcf_in = VariantFile(args.hap2)
for rec in bcf_in.fetch():
     if rec.info["SVTYPE"] == "BND":
        hap2_reads[count] = {}
        for read in rec.info["RNAMES"]:
            if read in read_lens:
                hap2_reads[count][read] = 1
        count += 1

count = 0
for item in hap1_reads:
    read_set = set(hap1_reads[item].keys())
    for item2 in hap2_reads:
        read_set2 = set(hap2_reads[item2].keys())
        combined = read_set.intersection(read_set2)
        if len(combined) >= args.min_reads:
            count += 1
            print(combined)
            break

print(count)
