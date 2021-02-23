import pysam
import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def get_annotation(annotation):
    if "LINE" in annotation:
        return "LINE"
    if "SINE" in annotation:
        return "SINE"
    if "SVA" in annotation:
        return "SVA"
    if "ERV" in annotation:
        return "ERV"
    return "NA" 

parser = argparse.ArgumentParser( description='Check polymorphic insertions against external SV calls')
parser.add_argument('--samples', required=True)
parser.add_argument('--suffix', required=False, default=".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.updated_annoation.tsv")
parser.add_argument('--folder', required=False, default="read_analysis")
parser.add_argument('--vcf', required=True)
parser.add_argument('--max-distance', required=False, default=1000)
parser.add_argument('--min-sv-size', required=False, default=50)
parser.add_argument('--genome', required=False, default="")
args = parser.parse_args()

# Read in the SR VCF and note the genomic positions of any Insert calls made > 50bp
# Next compare these to each of the samples to identify if they are within 500 bp and if they're called polymorphic
# Repeat with LR VCF. Print stats for each and the combined version

slr_insert_dict = defaultdict(IntervalTree)
for vcf in args.vcf.split(','):
    slr_inserts = pysam.VariantFile(vcf)
    count = 0
    for variant in slr_inserts.fetch():
        chrom = str(variant.chrom)
        if "chr" not in chrom:
            chrom = "chr"+chrom
        start = variant.start - int(args.max_distance)
        end = variant.stop + int(args.max_distance)
        key = variant.id
        insert = False
        ref = variant.ref
        for alt in variant.alts:
            if len(alt) - len(ref) >= 50:
                insert = True
        if not insert:
            continue
        if args.genome != "":
            if args.genome in key:
                slr_insert_dict[chrom][start:end] = key
        else:
            count += 1
            slr_insert_dict[chrom][start:end] = key

print(count)

tp_pass_no_polymorphic = defaultdict(int)
fp_pass_missed_polymorphic = defaultdict(int)
fn_polymorphic = defaultdict(int)
tn_polymorphic = defaultdict(int)
total = defaultdict(int)

for sample in args.samples.split(','):
    print(sample)
    with open(sample+"/"+args.folder+"/"+sample+args.suffix, 'r') as in_tsv:
        for line in in_tsv:
            # Check to see if the insert is a PASS or polymorphic
            # If PASS check to see that we don't have an SV call here
            # If polymorphic check to see that we do
            if "PASS" not in line and "polymorphic" not in line:
                continue 
            line_arr = line.strip().split('\t')
            annotation = get_annotation(line_arr[9])
            if "PASS" in line:
                # Check to see if we have an insert in the VCFs
                # If so this is a FP that should have been flagged by our 
                nearby = slr_insert_dict[line_arr[0]][int(line_arr[1]):int(line_arr[2])]
                if len(nearby) > 0:
                    fp_pass_missed_polymorphic[annotation] += 1
                else:
                    tp_pass_no_polymorphic[annotation] += 1
            if "polymorphic" in line:
                nearby = slr_insert_dict[line_arr[0]][int(line_arr[1]):int(line_arr[2])]
                if len(nearby) > 0:
                    tn_polymorphic[annotation] += 1
                else:
                    #print(line.strip()+"\t"+sample)
                    fn_polymorphic[annotation] += 1
                total[annotation] += 1

print(tp_pass_no_polymorphic)
print(fp_pass_missed_polymorphic)
print(fn_polymorphic)
print(tn_polymorphic)
print(total)
