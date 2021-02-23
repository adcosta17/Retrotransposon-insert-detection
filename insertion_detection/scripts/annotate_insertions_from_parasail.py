import pysam
import argparse
import sys
import csv
from collections import defaultdict


parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequences with parasail')
parser.add_argument('--input', required=True)
parser.add_argument('--alignment', required=True)
args = parser.parse_args()

read_to_best_annotation = defaultdict(list)

with open(args.alignment) as align_in:
    count = 0
    for row in align_in:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if row_args[0] not in read_to_best_annotation:
            read_to_best_annotation[row_args[0]].append(row_args)
            continue
        # See if the new alignment is better than the old
        # Sort by score
        if int(row_args[8]) > int(read_to_best_annotation[row_args[0]][0][8]):
            # New row is higher scoring thatn old one
            old = read_to_best_annotation[row_args[0]][0]
            if len(read_to_best_annotation[row_args[0]]) == 1:
                read_to_best_annotation[row_args[0]][0] = row_args;
                read_to_best_annotation[row_args[0]].append(old)
            else:
                # There are 2 entries in the list
                read_to_best_annotation[row_args[0]][0] = row_args;
                read_to_best_annotation[row_args[0]][1] = old;
        elif len(read_to_best_annotation[row_args[0]]) == 1:
            # Score isn't greater but only one entry, add second
            read_to_best_annotation[row_args[0]].append(row_args)
        else:
            # Have at 2 entries here, check to see if row is higher than second
            if int(row_args[8]) >= int(read_to_best_annotation[row_args[0]][1][8]):
                read_to_best_annotation[row_args[0]][1] = row_args;


with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            print(row.strip()+"\t"+"\t".join(["annotation","score","matches","percent_identity_mapped","percent_identity_insert","percent_identity","read_start","read_end","mapping_start","mapping_end"]))
            continue
        annotation = "no_mapping\t0\t0\t0\t0\t0\t0\t0"
        # Get the best annotation for the read an the insert postion
        if row_args[3]+":"+row_args[4]+"-"+row_args[5] in read_to_best_annotation:
            # Check to see if the best mapping is significanly better than the second 
            if (int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][8]) - int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][1][8]) > 10 and
               int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][8])/int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][1][8]) > 1.05):
               # Consider the best score as annotation
               insert_size = int(row_args[5]) - int(row_args[4])
               mapped_size = int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][5])
               pii = int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][9])/insert_size
               pim = int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][9])/mapped_size
               min_size = min(insert_size, mapped_size)
               pi = int(read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][9])/min_size
               annotation = read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][4]+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][8]+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][9]+"\t"\
               +str(pim)+"\t"+str(pii)+"\t"+str(pi)+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][2]+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][3]+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][6]+"\t"\
               +read_to_best_annotation[row_args[3]+":"+row_args[4]+"-"+row_args[5]][0][7]
        pass_fail = row_args[8]
        if "no_mapping" in annotation:
            if "PASS" in pass_fail:
                pass_fail = "mapped_fraction"
            else:
                pass_fail = pass_fail + ",mapped_fraction"
        row_args[8] = pass_fail
        print("\t".join(row_args)+"\t"+annotation)
