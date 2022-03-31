import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

parser = argparse.ArgumentParser( description='Filter xTea Calls based on control')
parser.add_argument('--test-samples', required=True)
parser.add_argument('--reference-sample', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=500)
args = parser.parse_args()

# read reference insertions and set up interval trees
intervaltrees = defaultdict(IntervalTree)
with open(args.reference_sample+"/xtea_before/"+args.reference_sample+"/classified_results.txt.Combined.txt") as csvfile:
    for line in csvfile:
        row = line.strip().split("\t")
        chrom = row[0]
        start = int(row[1]) - args.max_distance
        end = int(row[1]) + args.max_distance
        key = chrom + ":" + str(start) + "-" + str(end)
        intervaltrees[chrom][start:end] = key


# read input files and only emit records not near a record in the interval tree
for sample in args.test_samples.split(','):
    if sample == args.reference_sample:
        continue
    with open(sample+"/xtea_before/"+sample+"/classified_results.txt.Combined.txt") as csvfile:
        for line in csvfile:
            row = line.strip().split("\t")
            chrom = row[0]
            start = int(row[1])
            end = int(row[1]) + 1
            nearby = intervaltrees[chrom][start:end]
            if len(nearby) == 0:
                row.append("Insert_PASS")
            else:
                row.append("Insert_In_Control")
            print(sample+"\t"+"\t".join(row))

