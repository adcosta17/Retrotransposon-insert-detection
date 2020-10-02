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
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv")
args = parser.parse_args()

#read_hap = {}
#total_align = 0
#for sample in args.sample.split(','):
#    print(sample,file=sys.stderr)
#    sam_reader = pysam.AlignmentFile(sample+"/phased/"+sample+".sorted.phased.bam")
#    for record in sam_reader.fetch():
#        if record.is_unmapped:
#            continue
#        if record.query_name not in read_hap:
#            read_hap[record.query_name] = defaultdict(int)
#        tags = record.get_tags()
#        total_align += 1
#        hp = 0
#        for t in tags:
#            if t[0] == "HP":
#                hp = int(t[1])
#        read_hap[record.query_name][hp] += 1
#
## Compute overall haplotype stats
#total_reads = len(read_hap)
#count_0 = 0
#count_1 = 0
#count_2 = 0
#count_1_2 = 0
#count_hap = 0
#count_no_hap = 0
#for read in read_hap:
#    count_0 += read_hap[read][0]
#    count_1 += read_hap[read][1]
#    count_2 += read_hap[read][2]
#    if read_hap[read][1] > 0 and read_hap[read][2] > 0:
#        count_1_2 += 1
#    if read_hap[read][1] > 0 or read_hap[read][2] > 0:
#        count_hap += 1
#    else:
#        count_no_hap += 1
#
## Print stats
#print("Alignment Level:")
#print("Total Aligned "+str(total_align))
#print("No HP Assigned Alignments "+str(count_0)+" "+str(count_0/total_align))
#print("HP 1 Assigned Alignments "+str(count_1)+" "+str(count_1/total_align))
#print("HP 2 Assigned Alignments "+str(count_2)+" "+str(count_2/total_align))
#print("Either HP Assigned Alignments "+str(count_1+count_2)+" "+str((count_1+count_2)/total_align))
#print("Read Level:")
#print("Total Reads "+str(total_reads))
#print("No HP Assigned Reads "+str(count_no_hap)+" "+str(count_no_hap/total_reads))
#print("HP Assigned Reads "+str(count_hap)+" "+str(count_hap/total_reads))
#print("Both HP Assigned Reads "+str(count_1_2)+" "+str(count_1_2/total_reads))

all_insert_hap = defaultdict(int)
all_insert_pass_hap = defaultdict(int)
total_align = 0
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    with open("results_sep_17/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if "min_insertion_length" in line_arr[8] or "mapq<20" in line_arr[8]:
                continue
            if line_arr[8] == "PASS":
                try:
                    all_insert_pass_hap[int(line_arr[18])] += 1
                except ValueError:
                    all_insert_pass_hap[int(line_arr[19])] += 1
            try:
                all_insert_hap[int(line_arr[18])] += 1
            except ValueError:
                all_insert_hap[int(line_arr[19])] += 1

print("Insert Level:")
print("Total Inserts: "+str(all_insert_hap[0]+all_insert_hap[1]+all_insert_hap[2]))
print("No HP Inserts "+str(all_insert_hap[0])+" "+str(all_insert_hap[0]/(all_insert_hap[0]+all_insert_hap[1]+all_insert_hap[2])))
print("HP 1 Inserts "+str(all_insert_hap[1])+" "+str(all_insert_hap[1]/(all_insert_hap[0]+all_insert_hap[1]+all_insert_hap[2])))
print("HP 2 Inserts "+str(all_insert_hap[2])+" "+str(all_insert_hap[2]/(all_insert_hap[0]+all_insert_hap[1]+all_insert_hap[2])))
print("Passing Inserts: "+str(all_insert_pass_hap[0]+all_insert_pass_hap[1]+all_insert_pass_hap[2]))
print("No HP Passing Inserts "+str(all_insert_pass_hap[0])+" "+str(all_insert_pass_hap[0]/(all_insert_pass_hap[0]+all_insert_pass_hap[1]+all_insert_pass_hap[2])))
print("HP 1 Passing Inserts "+str(all_insert_pass_hap[1])+" "+str(all_insert_pass_hap[1]/(all_insert_pass_hap[0]+all_insert_pass_hap[1]+all_insert_pass_hap[2])))
print("HP 2 Passing Inserts "+str(all_insert_pass_hap[2])+" "+str(all_insert_pass_hap[2]/(all_insert_pass_hap[0]+all_insert_pass_hap[1]+all_insert_pass_hap[2])))
