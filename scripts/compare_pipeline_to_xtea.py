import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

parser = argparse.ArgumentParser( description='Filter xTea Calls based on control')
parser.add_argument('--samples', required=True)
parser.add_argument('--max-distance', type=int, default=500)
args = parser.parse_args()

# read xtea insertions and set up interval trees
intervaltrees = defaultdict(IntervalTree)
for sample in args.samples.split(','):
    with open(sample+"/xtea_before/"+sample+"/classified_results.txt.Combined.txt") as csvfile:
        for line in csvfile:
            row = line.strip().split("\t")
            chrom = row[0]
            start = int(row[1]) - args.max_distance
            end = int(row[1]) + args.max_distance
            key = chrom + ":" + str(start) + "-" + str(end)
            intervaltrees[chrom][start:end] = key

# Compare the inserts from each sample to that of xtea over all samples
# Count how many inserts are shared, how many are exclusive to xtea, how many exclusive to pipeline for each sample
# Looking at polymorphic inserts
print("Sample\tShared\tUnique\tSharedPolymorphic\tUniquePolymorphic\tSharedPass\tUniquePass\tSharedLINE\tUniqueLINE\tSharedPolymorphicLINE\tUniquePolymorphicLINE\tSharedPassLINE\tUniquePassLINE\tSharedAlu\tUniqueAlu\tSharedPolymorphicAlu\tUniquePolymorphicAlu\tSharedPassAlu\tUniquePassAlu\tSharedSVA\tUniqueSVA\tSharedPolymorphicSVA\tUniquePolymorphicSVA\tSharedPassSVA\tUniquePassSVA")

for sample in args.samples.split(','):
    shared = 0
    unique = 0
    shared_pass = 0
    unique_pass = 0
    shared_polymorphic = 0
    unique_polymorphic = 0
    shared_line = 0
    unique_line = 0
    shared_pass_line = 0
    unique_pass_line = 0
    shared_polymorphic_line = 0
    unique_polymorphic_line = 0
    shared_sine = 0
    unique_sine = 0
    shared_pass_sine = 0
    unique_pass_sine = 0
    shared_polymorphic_sine = 0
    unique_polymorphic_sine = 0
    shared_sva = 0
    unique_sva = 0
    shared_pass_sva = 0
    unique_pass_sva = 0
    shared_polymorphic_sva = 0
    unique_polymorphic_sva = 0
    with open(sample+"/winnow_realign_read_analysis/"+sample+".read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv") as csvfile:
        count = 0
        for line in csvfile:
            if count == 0:
                count = 1
                continue
            row = line.strip().split("\t")
            if "no_repbase" in line or "min_insert" in line or "mapq" in line:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            nearby = intervaltrees[chrom][start:end]
            if len(nearby) == 0:
                # Unique to sample
                unique += 1
                if "LINE" in row[9]:
                    unique_line += 1
                elif "SINE" in row[9]:
                    unique_sine += 1
                elif "SVA" in row[9]:
                    unique_sva += 1
                if "PASS" in line:
                    unique_pass += 1
                    if "LINE" in row[9]:
                        unique_pass_line += 1
                    elif "SINE" in row[9]:
                        unique_pass_sine += 1
                    elif "SVA" in row[9]:
                        unique_pass_sva += 1
                elif "polymorphic" in line or "ref_control" in line:
                    unique_polymorphic += 1
                    if "LINE" in row[9]:
                        unique_polymorphic_line += 1
                    elif "SINE" in row[9]:
                        unique_polymorphic_sine += 1
                    elif "SVA" in row[9]:
                        unique_polymorphic_sva += 1
            else:
                shared += 1
                if "LINE" in row[9]:
                    shared_line += 1
                elif "SINE" in row[9]:
                    shared_sine += 1
                elif "SVA" in row[9]:
                    shared_sva += 1
                if "PASS" in line:
                    shared_pass += 1
                    if "LINE" in row[9]:
                        shared_pass_line += 1
                    elif "SINE" in row[9]:
                        shared_pass_sine += 1
                    elif "SVA" in row[9]:
                        shared_pass_sva += 1
                elif "polymorphic" in line or "ref_control" in line:
                    shared_polymorphic += 1
                    if "LINE" in row[9]:
                        shared_polymorphic_line += 1
                    elif "SINE" in row[9]:
                        shared_polymorphic_sine += 1
                    elif "SVA" in row[9]:
                        shared_polymorphic_sva += 1
    print(sample+"\t"+str(shared)+"\t"+str(unique)+"\t"+str(shared_polymorphic)+"\t"+str(unique_polymorphic)+"\t"+str(shared_pass)+"\t"+str(unique_pass)+"\t"+str(shared_line)+"\t"+str(unique_line)+"\t"+str(shared_polymorphic_line)+"\t"+str(unique_polymorphic_line)+"\t"+str(shared_pass_line)+"\t"+str(unique_pass_line)+"\t"+str(shared_sine)+"\t"+str(unique_sine)+"\t"+str(shared_polymorphic_sine)+"\t"+str(unique_polymorphic_sine)+"\t"+str(shared_pass_sine)+"\t"+str(unique_pass_sine)+"\t"+str(shared_sva)+"\t"+str(unique_sva)+"\t"+str(shared_polymorphic_sva)+"\t"+str(unique_polymorphic_sva)+"\t"+str(shared_pass_sva)+"\t"+str(unique_pass_sva))