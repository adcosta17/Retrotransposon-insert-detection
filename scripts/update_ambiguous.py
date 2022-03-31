import pysam
import argparse
import sys
import csv
from collections import defaultdict
import mappy as mp

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

def parse_mappings_from_tab(input_tab, min_escore):
    fn = [ "score", "name1", "start1", "alnSize1", "strand1", "seqSize1", "name2", "start2", "alnSize2", "strand2", "seqSize2", "blocks", "EG", "E" ]
    family_count = {}
    with open(input_tab) as tsvfile:
        reader = csv.DictReader(tsvfile, fieldnames=fn, delimiter='\t')
        for row in reader:

            # hack to skip comments
            if row['score'].find("#") != -1:
                continue

            escore = float(row['E'].split('=')[1])
            if escore > min_escore:
                continue

            repeat_name = row['name1']
            read_name = row['name2'].split(":")[0]
            insert_position = row['name2'].split(":")[1]
            
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
    return family_count

def multiple_families(family_count, read_name, insert_position):
    count = 0
    if read_name not in family_count:
        return False
    if insert_position not in family_count[read_name]:
        return False
    for item in family_count[read_name][insert_position]:
        if family_count[read_name][insert_position][item] > 0:
            count += 1
    if count > 1:
        return True
    return False

def get_family_str(family_count, read_name, insert_position):
    ret = "ambiguous_mapping"
    for item in family_count[read_name][insert_position]:
        if family_count[read_name][insert_position][item] > 0:
            ret = ret +"_"+item
    return ret


parser = argparse.ArgumentParser( description='Annotate .read_insertions.tsv with mappings of the insertion sequence to repbase')
parser.add_argument('--input', required=True)
parser.add_argument('--minimap2-paf', required=False)
parser.add_argument('--last-tab', required=False)
parser.add_argument('--min-mapped-fraction', type=float, default=0.9)
parser.add_argument('--min-mapped-length', type=int, default=100)
parser.add_argument('--min-escore', type=float, default=1e-14)
parser.add_argument('--sub-family-fasta', required=False)
args = parser.parse_args()

read_to_best_annotation = defaultdict(RepbaseMapping)

if args.minimap2_paf:
    read_to_best_annotation = parse_mappings_from_paf(args.minimap2_paf)
elif args.last_tab:
    family_count = parse_mappings_from_tab(args.last_tab, args.min_escore)

if args.sub_family_fasta:
    a = mp.Aligner(args.sub_family_fasta)  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")

with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            print(row.strip())
            continue
        # Get the best annotation for the read an the insert postion
        if multiple_families(family_count, row_args[3], str(row_args[4]+"-"+row_args[5])):
            row_args[9] = get_family_str(family_count, row_args[3], str(row_args[4]+"-"+row_args[5]))
        if "LINE" in row_args[9] and args.sub_family_fasta:
            # Now check for sub_family mapping
            best_score = 0
            best_hit = ""
            for hit in a.map(row_args[7]): # traverse alignments
                a_score = (hit.mlen-hit.NM)/hit.blen
                if a_score > best_score:
                    best_hit = hit.ctg
                    best_score = a_score
            if best_hit != "":
                row_args[9] = row_args[9]+"___Subfamily:_"+best_hit
        print("\t".join(row_args))
