import sys
import csv
import argparse
import pysam

def get_filters(annotation):
    for item in annotation.split(','):
        if "PASS" in item or "in_sample" in item:
            continue
        else:
            return False
    return True

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
        if get_filters(row[8]):
            # Need to now get the reads that support
            read = row[3]
            count += 1
            if count % 1000 == 0:
                print(count, file=sys.stderr)
            if read not in seen:
                try:
                    seq = fh.fetch(read)
                    print(">"+read+"\n"+seq)
                    seen[read] = 1
                except:
                    continue


