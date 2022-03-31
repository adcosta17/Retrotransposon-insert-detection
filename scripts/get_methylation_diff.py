import pysam
import argparse
import sys
import csv
from collections import defaultdict
import mappy as mp



parser = argparse.ArgumentParser( description='Get L1 Subfamilies')
parser.add_argument('--input', required=True)
args = parser.parse_args()

with pysam.FastxFile(args.input) as fh:
    for entry in fh:
        if "L1" in entry.name:
            print(">"+entry.name)
            print(entry.sequence)

