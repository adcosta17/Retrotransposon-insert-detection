import pysam
import argparse
import sys
import csv
from collections import defaultdict

def filter_poly_AT(insert_seq):
    has_poly_AT = False
    start = insert_seq[:50]
    end = insert_seq[len(insert_seq)-50:]
    pA = "A"*5
    pT = "T"*5
    if pA in start or pT in start or pA in end or pT in end:
        has_poly_AT = True
    return has_poly_AT

parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequence to repbase')
parser.add_argument('--sample-list', required=True)
args = parser.parse_args()

before_LINE = defaultdict(int)
after_LINE = defaultdict(int)
before_Alu = defaultdict(int)
after_Alu = defaultdict(int)
print("ORIGINAL:")
print("SAMPLE\tLINE_PASS_ALL\tLINE_PASS_POLY_A\tFRACTION\tAlu_PASS_ALL\tAlu_PASS_POLY_A\tFRACTION")
for sample in args.sample_list.split(','):
    with open(sample+"/winnow_read_analysis/"+sample+".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.polyAT.tsv", 'r') as in_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            if "PASS" in row[8] and "LINE" in row[9]:
                before_LINE[sample] += 1
                if 'PolyAT' in row[8]:
                    after_LINE[sample] += 1
            if "PASS" in row[8] and "Alu" in row[9]:
                before_Alu[sample] += 1
                if 'PolyAT' in row[8]:
                    after_Alu[sample] += 1
    print(sample+"\t"+str(before_LINE[sample])+"\t"+str(after_LINE[sample])+"\t"+str(after_LINE[sample]/before_LINE[sample])+"\t"+str(before_Alu[sample])+"\t"+str(after_Alu[sample])+"\t"+str(after_Alu[sample]/before_Alu[sample]))

print("\n")

before_LINE = defaultdict(int)
after_LINE = defaultdict(int)
before_Alu = defaultdict(int)
after_Alu = defaultdict(int)
print("UPDATED:")
print("SAMPLE\tLINE_PASS_ALL\tLINE_PASS_POLY_A\tFRACTION\tAlu_PASS_ALL\tAlu_PASS_POLY_A\tFRACTION")
for sample in args.sample_list.split(','):
    with open(sample+"/winnow_read_analysis/"+sample+".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.polyAT.updated_annotation.tsv", 'r') as in_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            if "PASS" in row[8] and "LINE" in row[9]:
                before_LINE[sample] += 1
                if 'PolyAT' in row[8]:
                    after_LINE[sample] += 1
            if "PASS" in row[8] and "Alu" in row[9]:
                before_Alu[sample] += 1
                if 'PolyAT' in row[8]:
                    after_Alu[sample] += 1
    print(sample+"\t"+str(before_LINE[sample])+"\t"+str(after_LINE[sample])+"\t"+str(after_LINE[sample]/before_LINE[sample])+"\t"+str(before_Alu[sample])+"\t"+str(after_Alu[sample])+"\t"+str(after_Alu[sample]/before_Alu[sample]))