import sys
import csv
import argparse
import pysam

parser = argparse.ArgumentParser( description='Make regions list for splitting jobs.')
parser.add_argument('-t', '--tsv', type=str, required=True)
parser.add_argument('-f', '--fastq', type=str, required=True)
parser.add_argument('-s', '--sample', type=str, required=True)
args = parser.parse_args()

count = 0
fh = pysam.FastaFile(args.fastq)
seen = {}
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        if row[6] == "PASS" and "Polymorphic" not in line:
            # Need to now get the reads that support
            for item in row[5].split(','):
                if args.sample == item.split(':')[0]:
                    read = item.split(':')[1]
                    if read not in seen:
                        try:
                            seq = fh.fetch(read)
                            print(">"+read+"\n"+seq)
                            seen[read] = 1
                        except:
                            continue


