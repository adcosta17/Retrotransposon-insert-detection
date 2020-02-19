import pysam
import argparse
import sys
import csv
import os
from os import listdir
from intervaltree import Interval, IntervalTree

from collections import defaultdict

parser = argparse.ArgumentParser( description='Get a table of data for upset plot')
parser.add_argument('--tsv', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--filter', required=False, default="")
parser.add_argument('--header', nargs='?', const="", default="")
parser.add_argument('--sc', nargs='?', const="", default="")
args = parser.parse_args()

in_type = "insert"
if args.sc:
	in_type = "softclip"

# Go through tsv and write out upset matrix
pass_fail = {}
count = 0
if args.header:
    count = 1
with open(args.tsv) as in_tsv:
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        if args.filter:
            if args.filter != line_args[22] or line_args[21] != in_type:
                continue
        reasons = line_args[8].split(",")
        for reason in reasons:
            pass_fail[reason] = 1

count = 0
if args.header:
    count = 1
with open(args.tsv) as in_tsv:
    with open(args.output, 'w') as out_up:
        out_up.write("Name")
        for reason in pass_fail:
            out_up.write("\t"+reason)
        out_up.write("\n")
        for line in in_tsv:
            if count == 0:
                count = 1
                continue
            line_args = line.strip().split('\t')
            if args.filter:
                if args.filter != line_args[22] or line_args[21] != in_type:
                    continue
            reasons = line_args[8].split(",")
            out_up.write(line_args[3]+":"+line_args[4]+"-"+line_args[5])
            for reason in pass_fail:
                if reason in reasons:
                    out_up.write("\t1")
                else:
                    out_up.write("\t0")
            out_up.write("\n")
