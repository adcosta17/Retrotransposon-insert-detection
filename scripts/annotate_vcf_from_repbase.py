import pysam
import argparse
import sys
import csv
from collections import defaultdict

class RepbaseMapping:

    def __init__(self, name, ml, fm, es, st, en, ms, me):
        self.name = name
        self.mapped_length = ml
        self.frac_mapped = fm
        self.escore = es
        self.start = st
        self.end = en
        self.sup = []
        self.mapping_start = ms
        self.mapping_end = me

# Assume the repeat family is the same for the hits
# Takes a repbase mapping and extends it or adds a secondary postion 
def merge_repbase_hits(current_mapping, new_start, new_end):
    mapping_to_ret = current_mapping
    # First see if new positions intersect current ones
    if ((new_start < mapping_to_ret.start and new_end > mapping_to_ret.start) or
        (new_end > mapping_to_ret.end and new_start < mapping_to_ret.end)):
        if new_start < mapping_to_ret.start and new_end > mapping_to_ret.start:
            # Update the start position
            mapping_to_ret.start = new_start
        if new_end > mapping_to_ret.end and new_start < mapping_to_ret.end:
            # Update the end postion
            mapping_to_ret.end = new_end
    # If it doesn't intersec see if the new mapping is fully contained in the old one
    # If it isnt contained then consider it a secondary hit
    elif not (new_start >= mapping_to_ret.start and new_end <= mapping_to_ret.end):
        mapping_to_ret.sup.append((new_start, new_end))

    return mapping_to_ret


def parse_mappings_from_tab(input_tab):
    read_to_best_annotation = {}
    fn = [ "score", "name1", "start1", "alnSize1", "strand1", "seqSize1", "name2", "start2", "alnSize2", "strand2", "seqSize2", "blocks", "EG", "E" ]
    family_count = {}
    with open(input_tab) as tsvfile:
        reader = csv.DictReader(tsvfile, fieldnames=fn, delimiter='\t')
        for row in reader:
            # hack to skip comments
            if row['score'].find("#") != -1:
                continue
            escore = float(row['E'].split('=')[1])
            if escore > args.min_escore:
                continue
            repeat_name = row['name1']
            insert_id = row['name2']
            name_to_use = ""
            if insert_id not in family_count:
                family_count[insert_id] = defaultdict(int)
            if "LINE" in repeat_name:
                family_count[insert_id]["LINE"] += 1
                name_to_use= "LINE"
            elif "SINE" in repeat_name:
                family_count[insert_id]["SINE"] += 1
                name_to_use= "SINE"
            elif "SVA" in repeat_name:
                family_count[insert_id]["SVA"] += 1
                name_to_use= "SVA"
            elif "ERV" in repeat_name:
                family_count[insert_id]["ERV"] += 1
                name_to_use= "ERV"

            if insert_id not in read_to_best_annotation:
                read_to_best_annotation[insert_id] = name_to_use
    return (read_to_best_annotation,family_count)


parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequence to repbase')
parser.add_argument('--input', required=True)
parser.add_argument('--last-tab', required=True)
parser.add_argument('--min-escore', type=float, default=1e-14)
args = parser.parse_args()


read_to_best_annotation,family_count = parse_mappings_from_tab(args.last_tab)

vcf_in = pysam.VariantFile(args.input)
header_lines = str(vcf_in.header).strip().split('\n')
for line in header_lines:
    if "##" in line:
        print(line)
    else:
        print('##INFO=<ID=RTTYPE,Number=1,Type=String,Description="RT Repeat Family">')
        print(line)
for rec in vcf_in.fetch():
    if str(rec.id) in read_to_best_annotation:
        rec_row = str(rec).strip().split('\t')
        rec_row[7] = rec_row[7] + ";RTTYPE="+read_to_best_annotation[str(rec.id)]
        print('\t'.join(rec_row))
        #vcf_out.write(rec)
