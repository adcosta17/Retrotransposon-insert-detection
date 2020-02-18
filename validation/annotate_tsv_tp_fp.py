import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join
import sklearn.metrics as metrics
import matplotlib.pyplot as plt
import numpy as np
from inspect import signature
from intervaltree import Interval, IntervalTree

def calculate_sdust_score(seq):
    triplet_counts = defaultdict(int)
    for i in range(0, len(seq) - 2):
        triplet_counts[seq[i:i+3]] += 1
    sum_score = 0
    for triplet in triplet_counts:
        count = triplet_counts[triplet]
        s = float(count * (count - 1)) / 2.0
        sum_score += s
    if len(seq) - 1 == 0:
        return 0
    sum_score /= (len(seq) - 1)
    return sum_score    

def is_fp(intervaltrees, record):
    chrom = record[0]
    start = int(record[1])
    end = int(record[2])+1
    nearby = intervaltrees[chrom][start:end]
    if len(nearby) == 0:
        return True
    return False

def is_hc(intervaltrees, record):
    chrom = record[0]
    start = int(record[1])
    end = int(record[2])+1
    nearby = intervaltrees[chrom][start:end]
    if len(nearby) == 0:
        return False
    return True

def get_tp_call(value, truth):
    if value == 1 and truth == 1:
        return 1
    return 0

def get_fp_call(value, truth, ref_insert):
    if value == 1 and truth == 0 and ref_insert == 1:
        return 1
    return 0


parser = argparse.ArgumentParser( description='Output a list of reads that are true positive inserts')
parser.add_argument('--insert-tsv', required=False)
parser.add_argument('--soft-clip-tsv', required=False)
parser.add_argument('--insert-truth', required=False)
parser.add_argument('--soft-clip-truth', required=False)
parser.add_argument('--max-distance', type=int, default=200)
parser.add_argument('--cutoff-insert', type=float, default=0.7)
parser.add_argument('--cutoff-sc', type=float, default=0.95)
parser.add_argument('--fraction', type=float, default=0.5)
parser.add_argument('--dust', type=int, default=5)
parser.add_argument('--length', type=int, default=500)
parser.add_argument('--dust-sc', type=float, default=0.0)
parser.add_argument('--hc-window', type=int, default=25)
parser.add_argument('--reference-tsv', required=True)
parser.add_argument('--output-tsv', required=True)
args = parser.parse_args()


# Read in our truth data
truth_inserts_and_sc = {}
if args.insert_truth:
    with open(args.insert_truth) as in_it:
        for line in in_it:
            line_args = line.strip().split('\t')
            if int(line_args[5]) != 0:
                if line_args[3] not in truth_inserts_and_sc:
                    truth_inserts_and_sc[line_args[3]] = {}
                truth_inserts_and_sc[line_args[3]][line_args[4]+"-"+line_args[5]] = line_args
if args.soft_clip_truth:
    pruned = 0
    total = 0
    with open(args.soft_clip_truth) as in_it:
        for line in in_it:
            line_args = line.strip().split('\t')
            if int(line_args[5]) != 0:
                total+=1
                if line_args[3] not in truth_inserts_and_sc:
                    truth_inserts_and_sc[line_args[3]] = {}
                truth_inserts_and_sc[line_args[3]][line_args[4]+"-"+line_args[5]] = line_args

intervaltrees = defaultdict(IntervalTree)
with open(args.reference_tsv) as in_tsv:
    line = in_tsv.readline()
    line = in_tsv.readline()
    while line != '':
        line_args = line.strip().split('\t')
        chrom = line_args[0]
        start = int(line_args[1]) - args.max_distance
        end = int(line_args[2]) + args.max_distance
        key = chrom + ":" + str(start) + "-" + str(end)
        intervaltrees[chrom][start:end] = key
        line = in_tsv.readline()


with open(args.output_tsv, 'w') as out_tsv:
    tp_count_insert = 0
    tp_count_sc = 0
    tp_count_both = 0
    fp_count_insert = 0
    fp_count_sc = 0
    fp_count_both = 0
    records = defaultdict(list)
    count = 0
    if args.insert_tsv:
        with open(args.insert_tsv) as in_tsv:
            for line in in_tsv:
                if count == 0:
                    count = 1
                    out_tsv.write(line)
                    continue
                line_args = line.strip().split('\t')
                if line_args[8] != "PASS":
                    updated_line = line_args
                    if not is_fp(intervaltrees, updated_line):
                        if line_args[8] == "PASS":
                            updated_line[8] = "in_reference_sample"
                        else:
                            updated_line[8] = updated_line[8]+",in_reference_sample"
                    updated_line.append(0)
                    updated_line.append(0)
                    updated_line.append("insert")
                    records[updated_line[3]].append(updated_line)
                    continue
                if line_args[3] in truth_inserts_and_sc:
                    # know we have a true positive read
                    # Check to see if the positions match
                    insert_match = False
                    for region in truth_inserts_and_sc[line_args[3]]:
                        # Check that we intersect region
                        region_start = int(region.split("-")[0])
                        region_end = int(region.split("-")[1])
                        if ((int(line_args[4]) <= region_start and int(line_args[5]) > region_start) or
                            (int(line_args[4]) < region_end and int(line_args[5]) >= region_end) or
                            (int(line_args[4]) >= region_start and int(line_args[5]) <= region_end)):
                            insert_match = True
                            # Have a TP
                            updated_line = line_args
                            updated_line.append(1)
                            updated_line.append(1)
                            updated_line.append("insert")
                            records[updated_line[3]].append(updated_line)
                            for item in updated_line:
                                out_tsv.write(str(item)+"\t")
                            out_tsv.write("TP\n")
                            tp_count_insert += 1
                    if not insert_match:
                        # Have an insert on a read that didn't intersect with our true postive insert
                        # Have a FP
                        updated_line = line_args
                        updated_line.append(0)
                        if is_fp(intervaltrees, updated_line):
                            updated_line.append(1)
                            updated_line.append("insert")
                            records[updated_line[3]].append(updated_line)
                            fp_count_insert += 1
                            for item in updated_line:
                                out_tsv.write(str(item)+"\t")
                            out_tsv.write("FP\n")
                        else:
                            updated_line[8] = "in_reference_sample"
                            updated_line.append(0)
                            updated_line.append("insert")
                            records[updated_line[3]].append(updated_line)
                else:
                    # Insert called is either not annotated to the repeat we want
                    # And/Or is on a read that doesn't span a region we modified
                    # Have a FP
                    updated_line = line_args
                    updated_line.append(0)
                    if is_fp(intervaltrees, updated_line):
                        updated_line.append(1)
                        updated_line.append("insert")
                        records[updated_line[3]].append(updated_line)
                        fp_count_insert += 1
                        for item in updated_line:
                            out_tsv.write(str(item)+"\t")
                        out_tsv.write("FP\n")
                    else:
                        updated_line[8] = "in_reference_sample"
                        updated_line.append(0)
                        updated_line.append("insert")
                        records[updated_line[3]].append(updated_line)
        count = 0


    if args.soft_clip_tsv:
        with open(args.soft_clip_tsv) as in_tsv:
            count = 0
            for line in in_tsv:
                if count == 0:
                    count = 1
                    if not args.insert_tsv:
                        out_tsv.write(line)
                    continue
                line_args = line.strip().split('\t')
                #impose test filter here:
                if line_args[8] != "PASS":
                    updated_line = line_args
                    if not is_fp(intervaltrees, updated_line):
                        if line_args[8] == "PASS":
                            updated_line[8] = "in_reference_sample"
                        else:
                            updated_line[8] = updated_line[8]+",in_reference_sample"
                    updated_line.append(0)
                    updated_line.append(0)
                    updated_line.append("softclip")
                    records[updated_line[3]].append(updated_line)
                    continue
                if line_args[3] in truth_inserts_and_sc:
                    # know we have a true positive read
                    # Check to see if the positions match
                    insert_match = False
                    for region in truth_inserts_and_sc[line_args[3]]:
                        # Check that we intersect region
                        region_start = int(region.split("-")[0])
                        region_end = int(region.split("-")[1])
                        if ((int(line_args[4]) <= region_start and int(line_args[5]) > region_start) or
                            (int(line_args[4]) < region_end and int(line_args[5]) >= region_end) or
                            (int(line_args[4]) >= region_start and int(line_args[5]) <= region_end)):
                            insert_match = True
                            # Have a TP
                            updated_line = line_args
                            updated_line.append(1)
                            updated_line.append(1)
                            updated_line.append("softclip")
                            records[updated_line[3]].append(updated_line)
                            for item in updated_line:
                                out_tsv.write(str(item)+"\t")
                            out_tsv.write("TP\n")
                            tp_count_sc += 1
                    if not insert_match:
                        # Have an insert on a read that didn't intersect with our true postive insert
                        # Have a FP
                        updated_line = line_args
                        updated_line.append(0)
                        if is_fp(intervaltrees, updated_line):
                            updated_line.append(1)
                            updated_line.append("softclip")
                            records[updated_line[3]].append(updated_line)
                            fp_count_sc += 1
                            for item in updated_line:
                                out_tsv.write(str(item)+"\t")
                            out_tsv.write("FP\n")
                        else:
                            updated_line[8] = "in_reference_sample"
                            updated_line.append(0)
                            updated_line.append("softclip")
                            records[updated_line[3]].append(updated_line)
                else:
                    # Insert called is either not annotated to the repeat we want
                    # And/Or is on a read that doesn't span a region we modified
                    # Have a FP
                    updated_line = line_args
                    updated_line.append(0)
                    if is_fp(intervaltrees, updated_line):
                        updated_line.append(1)
                        updated_line.append("softclip")
                        records[updated_line[3]].append(updated_line)
                        fp_count_sc += 1
                        for item in updated_line:
                            out_tsv.write(str(item)+"\t")
                        out_tsv.write("FP\n")
                    else:
                        updated_line[8] = "in_reference_sample"
                        updated_line.append(0)
                        updated_line.append("softclip")
                        records[updated_line[3]].append(updated_line)


    fn = defaultdict(list)
    fn_count_both = 0
    fn_count_sc = 0
    fn_count_insert = 0
    output_inserts = {}
    for read in truth_inserts_and_sc :
        if read in records:
            # Go through each insert and see if it is in tp
            for region in truth_inserts_and_sc[read]:
                match = []
                region_start = int(region.split("-")[0])
                region_end = int(region.split("-")[1])
                for insert in records[read]:
                    insert_start = int(insert[4])
                    insert_end = int(insert[5])
                    if ((insert_start <= region_start and insert_end > region_start) or
                        (insert_start < region_end and insert_end >= region_end) or
                        (insert_start >= region_start and insert_end <= region_end)):
                        match.append(insert)
                if len(match) == 0:
                    # Read is in our list of inserts but we didn't get the region we were looking for
                    updated_line = truth_inserts_and_sc[read][region]
                    if updated_line[17] == "insert":
                        fn_count_insert += 1
                        updated_line[8] = "insert_not_found"
                        updated_line[9] = "no_repbase_mapping"
                        updated_line[17] = "0"
                        updated_line.extend(["0","0","0", "insert"])
                        for item in updated_line:
                            out_tsv.write(str(item)+"\t")
                        out_tsv.write("FN\n")
                    else:
                        fn_count_sc += 1
                        updated_line[8] = "sc_not_found"
                        updated_line[9] = "no_repbase_mapping"
                        updated_line[17] = "0"
                        updated_line.extend(["0","0","0", "softclip"])
                        for item in updated_line:
                            out_tsv.write(str(item)+"\t")
                        out_tsv.write("FN\n")
                else:
                    # We have an insert that intersects our read. If its a Negative output it, else we've got a TP
                    for insert in match:
                        if int(insert[19]) == 0 and int(insert[20] == 0):
                            # Have a FN insert/sc
                            if insert[21] == "insert":
                                fn_count_insert += 1
                                for item in insert:
                                    out_tsv.write(str(item)+"\t")
                                out_tsv.write("FN\n")
                                output_inserts[insert[3]+":"+insert[4]+"-"+insert[5]] = 1
                            else:
                                fn_count_sc += 1
                                for item in insert:
                                    out_tsv.write(str(item)+"\t")
                                out_tsv.write("FN\n")
                                output_inserts[insert[3]+":"+insert[4]+"-"+insert[5]] = 1

        else:
            # Missed inserts on read
            # Output each of them
            for region in truth_inserts_and_sc[read]:
                updated_line = truth_inserts_and_sc[read][region]
                if updated_line[17] == "insert":
                    fn_count_insert += 1
                    updated_line[8] = "read_not_found"
                    updated_line[9] = "no_repbase_mapping"
                    updated_line[17] = "0"
                    updated_line.extend(["0","0","0", "insert"])
                    for item in updated_line:
                        out_tsv.write(str(item)+"\t")
                    out_tsv.write("FN\n")
                else:
                    fn_count_sc += 1
                    updated_line[8] = "read_not_found"
                    updated_line[9] = "no_repbase_mapping"
                    updated_line[17] = "0"
                    updated_line.extend(["0","0","0", "softclip"])
                    for item in updated_line:
                        out_tsv.write(str(item)+"\t")
                    out_tsv.write("FN\n")

    # Go through the rest of records. 
    # If they haven't been written out yet, write them with a TN

    for read in records:
        for insert in records[read]:
            if insert[3]+":"+insert[4]+"-"+insert[5] in output_inserts:
                continue
            elif int(insert[19]) == 0 and int(insert[20] == 0):
                for item in insert:
                    out_tsv.write(str(item)+"\t")
                out_tsv.write("NC\n")


print("TP " + str(tp_count_insert + tp_count_sc) + "\tFP " + str(fp_count_insert + fp_count_sc) + "\tFN " + str(fn_count_insert + fn_count_sc))
print("TP Insert " + str(tp_count_insert)+ "\tFP Insert " + str(fp_count_insert)+ "\tFN Insert " + str(fn_count_insert))
print("TP SC " + str(tp_count_sc)+ "\tFP SC " + str(fp_count_sc) + "\tFN SC " + str(fn_count_sc))

sensitivity = float(tp_count_insert + tp_count_sc)/(tp_count_insert + tp_count_sc+fn_count_insert + fn_count_sc)
precision = float(tp_count_insert + tp_count_sc)/(tp_count_insert + tp_count_sc+fp_count_insert + fp_count_sc)
print("Sensitivity: " + str(sensitivity) + "\tPrescision: " + str(precision))

sensitivity = float(tp_count_insert)/(tp_count_insert +fn_count_insert )
precision = float(tp_count_insert)/(tp_count_insert +fp_count_insert )
print("Sensitivity Insert: " + str(sensitivity) + "\tPrescision Insert: " + str(precision))

sensitivity = float( tp_count_sc)/( tp_count_sc+ fn_count_sc)
precision = float( tp_count_sc)/( tp_count_sc+ fp_count_sc)
print("Sensitivity SC: " + str(sensitivity) + "\tPrescision SC: " + str(precision))

