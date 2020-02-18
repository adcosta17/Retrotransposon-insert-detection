import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join

def calculate_sdust_score(seq):
    triplet_counts = defaultdict(int)
    for i in range(0, len(seq) - 2):
        triplet_counts[seq[i:i+3]] += 1

    sum_score = 0
    for triplet in triplet_counts:
        count = triplet_counts[triplet]
        s = float(count * (count - 1)) / 2.0
        sum_score += s
    if (len(seq) - 1) > 0:
        sum_score /= (len(seq) - 1)
        return sum_score
    else:
        return 0

parser = argparse.ArgumentParser( description='Get a list of the true inserts in the tsv')
parser.add_argument('--spanning-insert-reads', required=True)
parser.add_argument('--spanning-sc-reads', required=False)
parser.add_argument('--insert-out', required=True)
parser.add_argument('--soft-clip-out', required=False)
parser.add_argument('--fastq-folder', required=True)
args = parser.parse_args()

supporting_inserts = {}
with open(args.spanning_insert_reads) as in_sr:
    for line in in_sr:
        line_args = line.strip().split('\t')
        if line_args[0] not in supporting_inserts:
            supporting_inserts[line_args[0]] = {}
        supporting_inserts[line_args[0]][line_args[2]+"-"+line_args[3]] = line_args[1]

supporting_sc = {}
if args.spanning_sc_reads:
    with open(args.spanning_sc_reads) as in_sr:
        for line in in_sr:
            line_args = line.strip().split('\t')
            if line_args[0] not in supporting_sc:
                supporting_sc[line_args[0]] = {}
            supporting_sc[line_args[0]][line_args[2]+"-"+line_args[3]] = line_args[1]



onlyfiles = [f for f in listdir(args.fastq_folder) if isfile(join(args.fastq_folder, f)) and f.endswith("fastq.gz")]
with open(args.insert_out, 'w') as out_tsv:
    for file in onlyfiles:
        # open each fastq. If reads we want are in the fastq look them up and get the specific sequence
        with pysam.FastaFile(filename = args.fastq_folder+ "/" + file, filepath_index_compressed = args.fastq_folder+ "/" +file + ".gzi") as fq:
                for read in supporting_inserts:
                    try:
                        seq = fq.fetch(read)
                        for region in supporting_inserts[read]:
                            region_start = int(region.split("-")[0])
                            region_end = int(region.split("-")[1])
                            chrom = supporting_inserts[read][region].split(":")[0]
                            ref_start = supporting_inserts[read][region].split(":")[1].split("-")[0]
                            ref_end = supporting_inserts[read][region].split(":")[1].split("-")[1]
                            dust = str(calculate_sdust_score(seq[region_start:region_end]))
                            out_tsv.write(chrom+"\t"+ref_start+"\t"+ref_end+"\t"+read+"\t"+str(region_start)+"\t"+str(region_end)+"\t"+dust+"\t"+seq[region_start:region_end]+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\tinsert\n")
                    except KeyError:
                        pass
with open(args.soft_clip_out, 'w') as out_tsv:
    for file in onlyfiles:
        # open each fastq. If reads we want are in the fastq look them up and get the specific sequence
        with pysam.FastaFile(filename = args.fastq_folder+ "/" + file, filepath_index_compressed = args.fastq_folder+ "/" +file + ".gzi") as fq:
                for read in supporting_sc:
                    try:
                        seq = fq.fetch(read)
                        for region in supporting_sc[read]:
                            region_start = int(region.split("-")[0])
                            region_end = int(region.split("-")[1])
                            chrom = supporting_sc[read][region].split(":")[0]
                            ref_start = supporting_sc[read][region].split(":")[1].split("-")[0]
                            ref_end = supporting_sc[read][region].split(":")[1].split("-")[1]
                            dust = str(calculate_sdust_score(seq[region_start:region_end]))
                            out_tsv.write(chrom+"\t"+ref_start+"\t"+ref_end+"\t"+read+"\t"+str(region_start)+"\t"+str(region_end)+"\t"+dust+"\t"+seq[region_start:region_end]+"\t0\t0\t0\t0\t0\t0\t0\t0\t0\tsoftclip\n")
                    except KeyError:
                        pass
                