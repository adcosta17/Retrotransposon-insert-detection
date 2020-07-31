import pysam
import argparse
import sys
import re
from collections import defaultdict

def calculate_sdust_score(seq):
    if seq == "":
        return 0
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


def get_deletion_pos(cigarstring):
    # Count up the position on the read until we get to a deletion
    deletion_positions = []
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) >= 100:
                deletion_positions.append((count,int(cg[:cg.find("D")])))
    return deletion_positions

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation

parser = argparse.ArgumentParser( description='Extract reads with long insertions')
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--min-insertion-length', type=int, default=100)
parser.add_argument('--min-detected-inclusion-length', type=int, default=50)
parser.add_argument('--min-flank-size', required=False, default=100)
parser.add_argument('--min-mapq', required=False, type=int, default=20)
parser.add_argument('--reference-gap-minimum', type=int, default=100)
parser.add_argument('--merged', type=str, required=True)
args = parser.parse_args()

merged_reads = {}
with open(args.merged) as in_merge:
    for line in in_merge:
        line = line.strip()
        merged_reads[line] = 1

max_ref_gap_at_candidate = args.reference_gap_minimum

mapped_count = {}
records_to_output = defaultdict(list)

sam_reader = pysam.AlignmentFile(args.bam)
print("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "dust_score", "insertion_sequence","pass_fail"]))
for record in sam_reader.fetch():
    if record.is_unmapped:
        continue

    read_annotation = "PASS"
    if record.mapq < args.min_mapq:
        read_annotation = "mapq<20" 

    # check for any long insertions
    aligned_pairs = record.get_aligned_pairs(matches_only=True)
    for idx in range(0, len(aligned_pairs) - 1):
        read_gap = aligned_pairs[idx + 1][0] - aligned_pairs[idx][0]
        ref_gap = aligned_pairs[idx + 1][1] - aligned_pairs[idx][1]
        if read_gap >= args.min_detected_inclusion_length and ref_gap <= max_ref_gap_at_candidate:
            ref_start = aligned_pairs[idx][1]
            ref_end = aligned_pairs[idx+1][1]
            read_start = aligned_pairs[idx][0]
            read_end = aligned_pairs[idx+1][0]
            insertion_sequence = ""
            if record.query_sequence is not None:
                insertion_sequence = record.query_sequence[read_start:read_end]
            sdust_score = calculate_sdust_score(insertion_sequence)
            annotation = read_annotation
            if not ((abs(ref_start - record.reference_start) > int(args.min_flank_size)) and 
                (abs(record.reference_end - ref_end) > int(args.min_flank_size))):
                annotation = update_annotation(annotation, "flank_size")
            if not read_gap >= args.min_insertion_length:
                annotation = update_annotation(annotation, "min_insertion_length")
            # convert read start and end positions if the alignment is rc
            if record.is_reverse and insertion_sequence != "":
                tmp = len(record.query_sequence) - read_start
                read_start = len(record.query_sequence) - read_end
                read_end = tmp
            records_to_output[record.query_name].append(["%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, sdust_score, insertion_sequence), annotation])


for read in records_to_output:
    for item in records_to_output[read]:
        annotation = item[1]
        print(item[0]+"\t"+annotation)
