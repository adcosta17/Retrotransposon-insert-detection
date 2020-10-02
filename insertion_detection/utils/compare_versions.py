import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree
from multiprocessing.pool import ThreadPool

def is_hc(intervaltrees, record):
    chrom = record[0]
    start = int(record[1])
    end = int(record[2])
    nearby = intervaltrees[chrom][start:end]
    if len(nearby) == 0:
        return False
    return True

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation

parser = argparse.ArgumentParser( description='Remove inserts where the insertion sequence maps to a location where the read also maps to with sufficent mapq')
parser.add_argument('--sample', required=True)
parser.add_argument('--folder-1', required=True)
parser.add_argument('--folder-2', required=True)
args = parser.parse_args()

pass_folder_1 = {}
pass_folder_2 = {}
total_align = 0
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    with open(args.folder_1+"/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if line_arr[8] == "PASS":
                line_arr.append(sample)
                pass_folder_1[line_arr[3]+":"+line_arr[4]+"-"+line_arr[5]] = line_arr
    with open(args.folder_2+"/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if line_arr[8] == "PASS":
                line_arr.append(sample)
                pass_folder_2[line_arr[3]+":"+line_arr[4]+"-"+line_arr[5]] = line_arr

# Now compare:
#print("Only in Folder 1")
only_in_one = {}
for item in pass_folder_1:
    if item not in pass_folder_2:
        only_in_one[item] = pass_folder_1[item]

#print("Only in Folder 2")
only_in_two = {}
for item in pass_folder_2:
    if item not in pass_folder_1:
        only_in_two[item] = pass_folder_2[item]

# Go back through the lists and output the records that differ
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    with open(args.folder_1+"/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if line_arr[3]+":"+line_arr[4]+"-"+line_arr[5] in only_in_two:
                print("Pass_Folder2\t" +"\t".join(only_in_two[line_arr[3]+":"+line_arr[4]+"-"+line_arr[5]])+"\t"+sample)
                print("Fail_Folder1\t" +"\t".join(line_arr)+"\t"+sample)
    with open(args.folder_2+"/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if line_arr[3]+":"+line_arr[4]+"-"+line_arr[5] in only_in_one:
                print("Pass_Folder1\t" +"\t".join(only_in_one[line_arr[3]+":"+line_arr[4]+"-"+line_arr[5]])+"\t"+sample)
                print("Fail_Folder2\t" +"\t".join(line_arr)+"\t"+sample)