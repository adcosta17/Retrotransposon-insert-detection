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

def parse_mappings_from_tab(input_tab, min_escore, sva_locations):
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
            read_name = row['name2'].split(":")[0]#+":"+row['name2'].split(":")[1]
            insert_position = row['name2'].split(":")[1]
            if read_name not in family_count:
                family_count[read_name] = {}
            if insert_position not in family_count[read_name]:
                family_count[read_name][insert_position] = defaultdict(list)
            if "LINE" in repeat_name:
                family_count[read_name][insert_position]["LINE"].append([repeat_name, int(row["start1"]), int(row["start1"])+int(row["alnSize1"])])
            elif "SINE" in repeat_name:
                family_count[read_name][insert_position]["SINE"].append([repeat_name, int(row["start1"]), int(row["start1"])+int(row["alnSize1"])])
            elif "SVA" in repeat_name:
                family_count[read_name][insert_position]["SVA"].append([repeat_name, int(row["start1"]), int(row["start1"])+int(row["alnSize1"])])
            elif "ERV" in repeat_name:
                family_count[read_name][insert_position]["ERV"].append([repeat_name, int(row["start1"]), int(row["start1"])+int(row["alnSize1"])])
    # Iterate over SVAs and check 
    for read in family_count:
        for pos in family_count[read]:
            new_sva = []
            if len(family_count[read][pos]["SVA"]) > 0:
                # Have an SVA, Check if the position for the SVA is only where the Alu aligns
                #print(read+" "+str(pos),file=sys.stderr)
                for item in family_count[read][pos]["SVA"]:
                    if item[1] <= sva_locations[item[0]][0] and item[2] <= sva_locations[item[0]][0]:
                        # Lines up with the Alu, skip
                        continue
                    elif item[1] >= sva_locations[item[0]][1] and item[2] >= sva_locations[item[0]][1]:
                        continue
                    else:
                        #print(sva_locations[item[0]],file=sys.stderr)
                        #print(item,file=sys.stderr)
                        new_sva.append(item)
                family_count[read][pos]["SVA"] = new_sva
    return family_count

def parse_alu_from_tab(input_tab, min_escore):
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
            read_name = row['name2']
            if repeat_name not in family_count:
                family_count[repeat_name] = [-10, int(row["seqSize1"])]
            middle = int(row["seqSize1"])/2
            start = int(row["start1"]) - 20
            end = int(row["start1"]) + int(row["alnSize1"]) + 20
            if start < middle and end < middle:
                # At start of SVA 
                if end > family_count[repeat_name][0]:
                    family_count[repeat_name][0] = end
            if start > middle and end > middle:
                # At end of SVA
                if start < family_count[repeat_name][1]:
                    family_count[repeat_name][1] = start
    return family_count

def multiple_families(family_count, read_name, insert_position):
    count = 0
    if read_name not in family_count:
        return False
    if insert_position not in family_count[read_name]:
        return False
    for item in family_count[read_name][insert_position]:
        if len(family_count[read_name][insert_position][item]) > 0:
            count += 1
    if count > 1:
        return True
    return False

def get_family_str(family_count, read_name, insert_position):
    if multiple_families(family_count, read_name, insert_position):
        ret = "ambiguous_mapping"
        for item in family_count[read_name][insert_position]:
            if len(family_count[read_name][insert_position][item]) > 0:
                ret = ret +"_"+item
        return ret
    else:
        ret = "resolved"
        for item in family_count[read_name][insert_position]:
            if len(family_count[read_name][insert_position][item]) > 0:
                ret = ret +"_"+item
        return ret


parser = argparse.ArgumentParser( description='Update Insert TSV to resolve ambiguous_mappings')
parser.add_argument('--input', required=True)
parser.add_argument('--alu-tab', required=True)
parser.add_argument('--minimap2-paf', required=False)
parser.add_argument('--last-tab', required=False)
parser.add_argument('--min-mapped-fraction', type=float, default=0.9)
parser.add_argument('--min-mapped-length', type=int, default=100)
parser.add_argument('--min-escore', type=float, default=1e-14)
parser.add_argument('--sub-family-fasta', required=False)
args = parser.parse_args()

read_to_best_annotation = defaultdict(RepbaseMapping)

# Parse Repbase Fasta, Align Alu sequences to SVA, identify all high quality mappings

sva_locations = parse_alu_from_tab(args.alu_tab, args.min_escore)

#exit(0)
if args.minimap2_paf:
    read_to_best_annotation = parse_mappings_from_paf(args.minimap2_paf)
elif args.last_tab:
    family_count = parse_mappings_from_tab(args.last_tab, args.min_escore, sva_locations)

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
        #if multiple_families(family_count, row_args[3], str(row_args[4]+"-"+row_args[5])):
        if "ambiguous_mapping" in row:
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
