#! /usr/env/python

import parasail
import pysam
import argparse
import sys
import csv

parser = argparse.ArgumentParser( description='Align inserts to a series of target sequences')
parser.add_argument('--inserts', required=True)
parser.add_argument('--targets', required=True)
args = parser.parse_args()

# read in the control sequences
control_seqs = list()
with open(args.targets) as control_fh:
    count = 0
    name = ""
    for line in control_fh:
        count += 1
        if count % 2 == 0:
            control_seqs.append([name, line.strip()])
        else:
            name = line.strip()

scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
print("Read\tReadLen\tReadStart\tReadEnd\tRef\tRefLen\tRefStart\tRefEnd\tScore\tMatches")
for read in pysam.FastxFile(args.inserts):
    # align this read against all oligos
    for (control_name, control_sequence) in control_seqs:
        #result = parasail.sw_stats_table_scan_16(read.sequence, control_sequence, 5, 4, scoring_matrix)
        if len(read.sequence) == 0:
            continue
        result = parasail.ssw(read.sequence, control_sequence, 5, 4, scoring_matrix)
        result2 = parasail.sw_stats_table_striped_16(read.sequence, control_sequence, 5, 4, scoring_matrix)
        print(read.name+"\t"+str(len(read.sequence))+"\t"+str(result.read_begin1)+"\t"+str(result.read_end1)+"\t"+control_name+"\t"+str(len(control_sequence))+"\t"+str(result.ref_begin1)+"\t"+str(result.ref_end1)+"\t"+str(result.score1)+"\t"+str(result2.matches))


