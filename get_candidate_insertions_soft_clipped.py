import pysam
import argparse
import sys
import re
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

def get_read_pos(record, min_insert, front):
    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)
    if not front:
        cigars = cigars[::-1]
    # Search until we find a softclip. Will happen at start or end. May have a hard clip before
    start = 0
    end = 0
    for cg in cigars:
        if cg.endswith('H'):
            step = int(cg[:cg.find("H")])
            end += step
        elif cg.endswith('S'):
            step = int(cg[:cg.find("S")])
            end += step
        else:
            break
    if end != 0:
        # Things are in the correct orientation
        end += start
        if abs(end - start) >= min_insert:
            return(start, end)
        else:
            return None
    else:
        return None

def get_ref_pos(record, front):
    cigars = re.findall('[0-9]*[A-Z]', record.cigarstring)
    if not front:
        cigars = cigars[::-1]
    # Search until we find a softclip. Will happen at start or end. May have a hard clip before
    start = 0
    for cg in cigars:
        if cg.endswith('H'):
            step = int(cg[:cg.find("H")])
            start += step
        elif cg.endswith('S'):
            step = int(cg[:cg.find("S")])
            start += step
        else:
            break
    if start != 0:
        if front:
            # Things are in the correct orientation
            # Softclips aren't on the reference. So we have to assume there is a perfect match and take the read position's base count
            return(record.reference_start - start, record.reference_start)
        else:
            return(record.reference_end, record.reference_end + start)
    else:
        return None

def get_read_length(cigarstring):
    # Gets the read length based on the Cigar String
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
    return count

parser = argparse.ArgumentParser( description='Extract reads with long insertions and soft clips')
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--original-bam', type=str, required=True)
parser.add_argument('--min-insertion-length', type=int, default=100)
parser.add_argument('--sc-length-filter', type=int, default=500)
parser.add_argument('--min-flank-size', required=False, default=100)
parser.add_argument('--merged', type=str, required=True)
args = parser.parse_args()

merged_reads = {}
with open(args.merged) as in_merge:
    for line in in_merge:
        line = line.strip()
        merged_reads[line] = 1

max_ref_gap_at_candidate = 5
sam_reader = pysam.AlignmentFile(args.original_bam)

mapped_count = {}
records_to_output = defaultdict(list)

for record in sam_reader.fetch():
    if record.is_unmapped:
        continue
    if record.query_name not in mapped_count:
        mapped_count[record.query_name] = 0
    mapped_count[record.query_name] += 1


print("\t".join(["chromosome", "reference_insertion_start", "reference_insertion_end", "read_name", "read_insertion_start", "read_insertion_end", "dust_score", "insertion_sequence", "pass_fail"]))
sam_reader = pysam.AlignmentFile(args.bam)
for record in sam_reader.fetch():
    if record.is_unmapped:
        continue

    read_annotation = "PASS"
    if record.mapq < 60:
        read_annotation = "mapq<60"        

    # Check if there a sufficently long soft clip on either end
    read_pos = get_read_pos(record, int(args.min_insertion_length), True)
    ref_pos = get_ref_pos(record, True)
    if read_pos is not None and ref_pos is not None:
        read_start = read_pos[0]
        read_end = read_pos[1]
        ref_start = ref_pos[1] - 1
        ref_end = ref_pos[1]
        insertion_sequence = record.query_sequence[read_start:read_end]
        sdust_score = calculate_sdust_score(insertion_sequence)
        annotation = read_annotation
        if not ((abs(ref_start - record.reference_start) > int(args.min_flank_size)) or 
            (abs(record.reference_end - ref_end) > int(args.min_flank_size))):
            if "PASS" in annotation:
                annotation = "flank_size"
            else:
                annotation = annotation + ",flank_size"
        if len(insertion_sequence) < args.sc_length_filter:
            if "PASS" in annotation:
                annotation = "sc_length"
            else:
                annotation = annotation + ",sc_length"
        records_to_output[record.query_name].append(("%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, sdust_score, insertion_sequence), annotation))

    read_pos = get_read_pos(record, int(args.min_insertion_length), False)
    ref_pos = get_ref_pos(record, False)
    if read_pos is not None and ref_pos is not None:
        read_len = get_read_length(record.cigarstring)
        read_start = read_len - read_pos[1]
        read_end = read_len - read_pos[0]
        ref_start = ref_pos[0]
        ref_end = ref_pos[0] + 1
        insertion_sequence = record.query_sequence[read_start:read_end]
        sdust_score = calculate_sdust_score(insertion_sequence)
        annotation = read_annotation
        if not ((abs(ref_start - record.reference_start) > int(args.min_flank_size)) or 
            (abs(record.reference_end - ref_end) > int(args.min_flank_size))):
            if "PASS" in annotation:
                annotation = "flank_size"
            else:
                annotation = annotation + ",flank_size"
        if len(insertion_sequence) < args.sc_length_filter:
            if "PASS" in annotation:
                annotation = "sc_length"
            else:
                annotation = annotation + ",sc_length"
        records_to_output[record.query_name].append(("%s\t%d\t%d\t%s\t%d\t%d\t%.1f\t%s" % (record.reference_name, ref_start, ref_end, record.query_name, read_start, read_end, sdust_score, insertion_sequence), annotation))

for read in records_to_output:
    for item in records_to_output[read]:
        annotation = item[1]
        #if mapped_count[read] > 1 and read not in merged_reads:
        #    if "PASS" not in annotation:
        #        annotation = annotation + ",multimapped"
        #    else:
        #        annotation = "multimapped"
        print(item[0]+"\t"+annotation)
        
