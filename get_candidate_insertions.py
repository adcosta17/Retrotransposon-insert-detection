import pysam
import argparse
import sys
from collections import defaultdict

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

parser = argparse.ArgumentParser( description='Extract reads with long insertions')
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--original-bam', type=str, required=True)
parser.add_argument('--min-insertion-length', type=int, default=100)
parser.add_argument('--min-flank-size', required=False, default=100)
parser.add_argument('--reference-gap-minimum', type=int, default=100)
parser.add_argument('--merged', type=str, required=True)
args = parser.parse_args()

merged_reads = {}
with open(args.merged) as in_merge:
    for line in in_merge:
        line = line.strip()
        merged_reads[line] = 1

max_ref_gap_at_candidate = args.reference_gap_minimum
sam_reader = pysam.AlignmentFile(args.original_bam)

mapped_count = {}
records_to_output = defaultdict(list)

for record in sam_reader.fetch():
    if record.is_unmapped:
        continue
    if record.query_name not in mapped_count:
        mapped_count[record.query_name] = 0
    mapped_count[record.query_name] += 1

sam_reader = pysam.AlignmentFile(args.bam)
print("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "dust_score", "insertion_sequence","pass_fail"]))
for record in sam_reader.fetch():
    if record.is_unmapped:
        continue

    read_annotation = "PASS"
    if record.mapq < 60:
        read_annotation = "mapq<60" 

    # check for any long insertions
    aligned_pairs = record.get_aligned_pairs(matches_only=True)
    for idx in range(0, len(aligned_pairs) - 1):
        read_gap = aligned_pairs[idx + 1][0] - aligned_pairs[idx][0]
        ref_gap = aligned_pairs[idx + 1][1] - aligned_pairs[idx][1]
        if read_gap >= args.min_insertion_length and ref_gap <= max_ref_gap_at_candidate:
            ref_start = aligned_pairs[idx][1]
            ref_end = aligned_pairs[idx+1][1]
            read_start = aligned_pairs[idx][0]
            read_end = aligned_pairs[idx+1][0]
            insertion_sequence = record.query_sequence[read_start:read_end]
            sdust_score = calculate_sdust_score(insertion_sequence)
            annotation = read_annotation
            if not ((abs(ref_start - record.reference_start) > int(args.min_flank_size)) and 
                (abs(record.reference_end - ref_end) > int(args.min_flank_size))):
                if "PASS" in annotation:
                    annotation = "flank_size"
                else:
                    annotation = annotation + ",flank_size"
            records_to_output[record.query_name].append(("%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, sdust_score, insertion_sequence), annotation))

for read in records_to_output:
    for item in records_to_output[read]:
        annotation = item[1]
        if mapped_count[read] > 1 and read not in merged_reads:
            if "PASS" not in annotation:
                annotation = annotation + ",multimapped"
            else:
                annotation = "multimapped"
        print(item[0]+"\t"+annotation)
