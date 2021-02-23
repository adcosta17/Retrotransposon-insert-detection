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

def get_mapped_total(annotation):
    regions = annotation.sup
    regions.append((annotation.start, annotation.end))
    combined_regions = []
    for begin,end in sorted(regions):
        if combined_regions and combined_regions[-1][1] >= begin - 1:
            combined_regions[-1][1] = max(combined_regions[-1][1], end)
        else:
            combined_regions.append([begin, end])
    mapped_total = 0
    for interval in combined_regions:
        mapped_total += (interval[1] - interval[0])
    return mapped_total

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
            read_name = row['name2'].split(":")[0]
            insert_position = row['name2'].split(":")[1]
            # Get start and end postion on 
            read_insert_start = int(row['name2'].split(":")[1].split("-")[0])
            read_alignment_start_pos = read_insert_start+int(row['start2'])
            read_alignment_end_pos = read_alignment_start_pos+int(row['alnSize2'])
            
            if read_name not in family_count:
                family_count[read_name] = {}
            if insert_position not in family_count[read_name]:
                family_count[read_name][insert_position] = defaultdict(int)
            if "LINE" in repeat_name:
                family_count[read_name][insert_position]["LINE"] += 1
            elif "SINE" in repeat_name:
                family_count[read_name][insert_position]["SINE"] += 1
            elif "SVA" in repeat_name:
                family_count[read_name][insert_position]["SVA"] += 1
            elif "ERV" in repeat_name:
                family_count[read_name][insert_position]["ERV"] += 1

            if read_name not in read_to_best_annotation:
                # If read not seen create a dict for it
                read_to_best_annotation[read_name] = defaultdict(RepbaseMapping)
            if insert_position not in read_to_best_annotation[read_name]:
                # If this insert on the read hasn't been seen before, add the hit. 
                # Repbase tab file is sorted by best hit first so first hit should be highest scoring
                read_to_best_annotation[read_name][insert_position] = RepbaseMapping(repeat_name, row['alnSize2'], float(row['alnSize2']) / float(row['seqSize2']), 
                                                                    escore, read_alignment_start_pos, read_alignment_end_pos, int(row['start1']), int(row['start1'])+int(row['alnSize1']))
            elif (("LINE" in repeat_name and "LINE" in read_to_best_annotation[read_name][insert_position].name) or 
                    ("SVA" in repeat_name and "SVA" in read_to_best_annotation[read_name][insert_position].name) or 
                    ("ERV" in repeat_name and "ERV" in read_to_best_annotation[read_name][insert_position].name) or 
                    ("SINE" in repeat_name and "SINE" in read_to_best_annotation[read_name][insert_position].name)):
                # Another hit for an insert we've already seen. Check to see if it is from the same repeat familiy as the best scoring hit, and if it intersects
                read_to_best_annotation[read_name][insert_position] = merge_repbase_hits(read_to_best_annotation[read_name][insert_position], read_alignment_start_pos, read_alignment_end_pos)
    return (read_to_best_annotation,family_count)

def multiple_families(family_count, read_name, insert_position):
    count = 0
    for item in family_count[read_name][insert_position]:
        if family_count[read_name][insert_position][item] > 0:
            count += 1
    if count > 1:
        return True
    return False

parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequence to repbase')
parser.add_argument('--input', required=True)
parser.add_argument('--minimap2-paf', required=False)
parser.add_argument('--last-tab', required=False)
parser.add_argument('--min-mapped-fraction', type=float, default=0.9)
parser.add_argument('--min-mapped-length', type=int, default=100)
parser.add_argument('--min-escore', type=float, default=1e-14)
args = parser.parse_args()

read_to_best_annotation = defaultdict(RepbaseMapping)

if args.minimap2_paf:
    read_to_best_annotation = parse_mappings_from_paf(args.minimap2_paf)
elif args.last_tab:
    read_to_best_annotation,family_count = parse_mappings_from_tab(args.last_tab)

with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            print(row.strip()+"\t"+"annotation"+"\t"+"mapped_length"+"\t"+"fraction_mapped"+"\t"+"alignment_escore"+"\t"+"read_alignment_insert_start"+"\t"+"read_alignment_insert_end"+"\t"+"read_mapped_bases_total"+"\t"+"mapping_start"+"\t"+"mapping_end")
            continue
        annotation = RepbaseMapping("no_repbase_mapping", 0, 0, 0, 0, 0, 0, 0)
        # Get the best annotation for the read an the insert postion
        if row_args[3] in read_to_best_annotation:
            if str(row_args[4]+"-"+row_args[5]) in read_to_best_annotation[row_args[3]]:
                annotation = read_to_best_annotation[row_args[3]][str(row_args[4]+"-"+row_args[5])]
                if multiple_families(family_count, row_args[3], str(row_args[4]+"-"+row_args[5])):
                    annotation.name = "ambiguous_mapping"
        # Get total amount of bases covered by this annotation
        mapped_total = get_mapped_total(annotation)
        pass_fail = row_args[8]
        if not annotation.frac_mapped > args.min_mapped_fraction:
            if "PASS" in pass_fail:
                pass_fail = "mapped_fraction"
            else:
                pass_fail = pass_fail + ",mapped_fraction"
        row_args[8] = pass_fail
        print("\t".join(row_args)+"\t"+str(annotation.name)+"\t"+str(annotation.mapped_length)+"\t"+str(annotation.frac_mapped)+"\t"+str(annotation.escore)+"\t"+str(annotation.start)+"\t"+str(annotation.end)+"\t"+str(mapped_total)+"\t"+str(annotation.mapping_start)+"\t"+str(annotation.mapping_end))
