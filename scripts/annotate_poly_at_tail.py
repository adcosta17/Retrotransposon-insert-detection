import pysam
import argparse
import sys
import csv
from collections import defaultdict

def filter_poly_AT(insert_seq):
    has_poly_AT = False
    start = insert_seq[:50]
    end = insert_seq[len(insert_seq)-50:]
    pA = "A"*5
    pT = "T"*5
    if pA in start or pT in start or pA in end or pT in end:
        has_poly_AT = True
    return has_poly_AT

parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequence to repbase')
parser.add_argument('--input', required=True)
args = parser.parse_args()

with open(args.input, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        if filter_poly_AT(row[7]):
            row[8] = row[8]+',PolyAT'
        print("\t".join(row))
