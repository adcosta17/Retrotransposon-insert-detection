import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

def get_inserts(all_inserts, chrom, start, end):
    mcount = 0
    for i in range(end-start):
        mcount = max(mcount, all_inserts[chrom][start+i])
    return mcount

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--input', required=True)
parser.add_argument('--reference-sample', type=str, required=True)
parser.add_argument('--reference-repeat-bed', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=100)
args = parser.parse_args()

# Read in all inserts and setup position list
all_inserts = defaultdict(lambda: defaultdict(int))
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        chrom = row_args[0]
        start = int(row_args[1]) - args.max_distance
        end = int(row_args[2]) + args.max_distance
        if start < 0:
            start = 0
        for i in range((end - start)):
            all_inserts[chrom][i+start] = 0


# Read inserts again and count how many fall within max_distance bp of each position in genome
with open(args.input) as csvfile:
    #reader = csv.DictReader(csvfile, delimiter="\t")
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        chrom = row_args[0]
        start = int(row_args[1]) - args.max_distance
        end = int(row_args[2]) + args.max_distance
        if start < 0:
            start = 0
        for i in range((end - start)):
            all_inserts[chrom][i+start] += 1


# read reference insertions and set up interval trees
intervaltrees = defaultdict(IntervalTree)
with open(args.reference_sample) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        # get insertion coordinates, insert into intervaltree allowing for the maximum allowable distance
        if row_args[8] == "PASS":
            chrom = row_args[0]
            start = int(row_args[1]) - args.max_distance
            end = int(row_args[2]) + args.max_distance

            key = chrom + ":" + str(start) + "-" + str(end)
            intervaltrees[chrom][start:end] = key


# read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            print(row.strip()+"\tsource")
            continue
        chrom = row_args[0]
        start = int(row_args[1])
        end = int(row_args[2]) + 1
        key = chrom + ":" + str(start) + "-" + str(end)
        nearby = intervaltrees[chrom][start:end]
        nearby_in_sample = get_inserts(all_inserts, chrom, start, end)
        if len(nearby) == 0:
            if nearby_in_sample >= 3:
                if row_args[8] == "PASS":
                    row_args[8] = "in_reference_sample"
                else:
                    row_args[8] = row_args[8] +",in_reference_sample"
                row_args.append("insert")
                print("\t".join(row_args))
            else:
                row_args.append("insert")
                print("\t".join(row_args))
        else:
            if row_args[8] == "PASS":
                row_args[8] = "in_reference_sample"
            else:
                row_args[8] = row_args[8] +",in_reference_sample"
            row_args.append("insert")
            print("\t".join(row_args))

