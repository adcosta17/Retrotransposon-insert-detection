import pysam
import argparse
import sys
import csv
from collections import defaultdict

parser = argparse.ArgumentParser( description='Convert a .read_insertions.tsv file to fasta')
parser.add_argument('--input', required=True)
args = parser.parse_args()

with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count += 1
            continue
        if "mapq<20" in row_args[9] or "min_insertion_length" in row_args[9]:
            continue
        record_name = row_args[3] + ":" + row_args[4] + "-" + row_args[5]
        print(">%s\n%s" % (record_name, row_args[7]))
