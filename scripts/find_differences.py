import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

def get_annotation(item):
    if item == "PASS":
        return item
    ret = []
    for a in item.split(','):
        if a == "mapq<20" or a == "mapped_fraction" or a == "in_centromere":
            continue
        else:
            ret.append(a)
    if len(ret) == 0:
        return "PASS"
    else:
        return ",".join(ret)

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--with-filters', required=True)
parser.add_argument('--without-filters', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=500)
args = parser.parse_args()

# read reference insertions and set up interval trees
inserts_per_read = defaultdict(list)
with open(args.without_filters) as csvfile:
    count = 0
    for line in csvfile:
        row = line.strip().split("\t")
        if count == 0:
            count = 1
            continue
        read = row[3]
        inserts_per_read[read].append(row)


# read input file and only emit records not near a record in the interval tree 
with open(args.with_filters) as csvfile:
    count = 0
    for line in csvfile:
        row = line.strip().split("\t")
        if count == 0:
            count = 1
            print("chrom\tstart\tend\tread\tread_start\tread_end\told_annotation\tnew_annotation\trepeat")
            continue
        read = row[3]
        if read in inserts_per_read:
            for item in inserts_per_read[read]:
                if item[0] == row[0] and item[1] == row[1] and item[2] == row[2] and item[3] == row[3] and item[4] == row[4] and item[5] == row[5]:
                    old = get_annotation(row[8])
                    new = get_annotation(item[8])
                    if old != new:
                        print(item[0]+"\t"+item[1]+"\t"+item[2]+"\t"+item[3]+"\t"+item[4]+"\t"+item[5]+"\t"+old+"\t"+new+"\t"+item[9])

