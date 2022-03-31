import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--input', required=True)
parser.add_argument('--reference-sample', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=500)
args = parser.parse_args()

# read reference insertions and set up interval trees
intervaltrees = defaultdict(IntervalTree)
with open(args.reference_sample) as csvfile:
    count = 0
    for line in csvfile:
        row = line.strip().split("\t")
        if count == 0:
            count = 1
            continue
        chrom = row[0]
        start = int(row[1]) - args.max_distance
        end = int(row[2]) + args.max_distance
        key = chrom + ":" + str(start) + "-" + str(end)
        intervaltrees[chrom][start:end] = key


# read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for line in csvfile:
        row = line.strip().split("\t")
        if count == 0:
            count = 1
            print("\t".join(row))
            continue
        chrom = row[0]
        start = int(row[1])
        end = int(row[2]) + 1
        nearby = intervaltrees[chrom][start:end]
        if len(nearby) == 0:
            print("\t".join(row))
        else:
            if row[8] == "PASS":
                row[8] = "in_reference_sample"
            else:
                row[8] = row[8] +",in_reference_sample"
            row.append("insert")
            print("\t".join(row))

