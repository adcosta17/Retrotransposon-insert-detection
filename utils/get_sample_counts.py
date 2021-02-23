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

def get_string(items):
    new_items = []
    for item in items:
        new_items.append("in_sample_"+item)
    return ",".join(new_items)

def get_sample_string(items):
    new_items = []
    for item in items:
        new_items.append(item)
    return ",".join(new_items)

def get_sample_dict(items):
    new_items = defaultdict(int)
    for item in items.split(","):
        new_items[item] += 1
    return new_items

def is_only_sample(sample_dict, sample):
    count = 0
    for item in sample_dict:
        count += sample_dict[item]
    if count == sample_dict[sample]:
        return True
    return False

def in_multiple_samples(sample_dict, sample_count, sample):
    count = 0
    for item in sample_dict:
        # If there is a single insert in another sample we can allow it
        if item == sample:
            continue
        if sample_dict[item] > 1:
            count += 1
    if count >= sample_count:
        return True
    return False

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in any other sample')
parser.add_argument('--multi-reference-filter-tsv', type=str, required=True)
args = parser.parse_args()

# read reference insertions and set up interval trees

with open(args.multi_reference_filter_tsv) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        # Get the count per sample 
        sample_dict = get_sample_dict(row_args[4])
        chrom = row_args[0]
        start = int(row_args[1])
        end = int(row_args[2])
        annotation = row_args[3]
        print(chrom+"\t"+row_args[1]+"\t"+row_args[2]+"\t"+annotation+"\t"+str(len(sample_dict))+"\t"+get_sample_string(sample_dict))


