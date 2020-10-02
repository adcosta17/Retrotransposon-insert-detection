import pysam
import argparse
import sys
import csv
import os
from os import listdir
from intervaltree import Interval, IntervalTree

from collections import defaultdict

parser = argparse.ArgumentParser( description='Get a table of data for upset plot')
parser.add_argument('--sample', required=True)
args = parser.parse_args()

failure_reason_counts = defaultdict(int)
failure_reason_strings = defaultdict(int)
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    with open("results_sep_17/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.updated_annoation.tsv",'r') as in_tsv:
        count = 0
        for line in in_tsv:
            line_arr = line.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if len(line_arr) < 21:
                continue
            reasons = line_arr[8].split(",")
            if "min_insertion_length" in reasons or "mapq<20" in reasons:
                continue
            #print(line_arr[21],file=sys.stderr)
            for reason in reasons:
                if "in_sample" in reason:
                    # Check the secondary reason
                    if "other_failure" in line_arr[21]:
                        continue
                    elif "ref_control_sample" in line_arr[21]:
                        failure_reason_counts["ref_control"] += 1
                    elif "alignment_mistake" in line_arr[21]:
                        failure_reason_counts["other_nearby_hp_inserts"] += 1
                    elif "No_haplotype_likely_polymorphic" in line_arr[21]:
                        failure_reason_counts["no_hp_polymorphic"] += 1
                    elif "liklely_polymorphic_single_hp" in line_arr[21]:
                        failure_reason_counts["single_hp_polymorphic"] += 1
                    elif "liklely_polymorphic_both_hp" in line_arr[21]:
                        failure_reason_counts["both_hp_polymorphic"] += 1
                else:
                    failure_reason_counts[reason] += 1
            # Get and add the string to the string dict
            reason_string = []
            for reason in reasons:
                if "in_sample" in reason:
                    # Check the secondary reason
                    if "other_failure" in line_arr[21]:
                        continue
                    elif "ref_control_sample" in line_arr[21]:
                        if "ref_control" not in reason_string:
                            reason_string.append("ref_control")
                    elif "alignment_mistake" in line_arr[21]:
                        if "other_nearby_hp_inserts" not in reason_string:
                            reason_string.append("other_nearby_hp_inserts")
                    elif "No_haplotype_likely_polymorphic" in line_arr[21]:
                        if "no_hp_polymorphic" not in reason_string:
                            reason_string.append("no_hp_polymorphic")
                    elif "liklely_polymorphic_single_hp" in line_arr[21]:
                        if "single_hp_polymorphic" not in reason_string:
                            reason_string.append("single_hp_polymorphic")
                    elif "liklely_polymorphic_both_hp" in line_arr[21]:
                        if "both_hp_polymorphic" not in reason_string:
                            reason_string.append("both_hp_polymorphic")
                else:
                    if reason not in reason_string:
                        reason_string.append(reason)
            failure_reason_strings[",".join(reason_string)] += 1

# Print out reasons
print("Reasons Singular")
for reason in failure_reason_counts:
    print(reason+"\t"+str(failure_reason_counts[reason]))
print("\nReasons Combined")
for reason in failure_reason_strings:
    print(reason+"\t"+str(failure_reason_strings[reason]))
