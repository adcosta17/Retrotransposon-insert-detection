import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--min-qual', type=int, default=1)
args = parser.parse_args()


infile = pysam.AlignmentFile(args.input, "rb")
outfile = pysam.AlignmentFile(args.output, "wb", template=infile)
for record in infile:
    if record.mapping_quality >= args.min_qual:
        outfile.write(record)
        