import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict

def get_overlap(start, end, int_start, int_end):
    if int_start <= start and int_end >= end:
        return end-start
    elif int_start > start and int_end < end:
        return int_end - int_start
    elif int_start <= start and int_end > start:
        return int_end - start
    elif int_start < end and int_end > end:
        return end-int_start
    return 0


def check_nearby(chrom, start, end, seen):
    ret = {}
    if chrom not in seen:
        return False
    i = start
    while i < end:
        if i in seen[chrom]["start"]:
            ret[seen[chrom]["start"][i]] = get_overlap(start,end,i,seen[chrom]["start"][i])
        if i in seen[chrom]["end"]:
            ret[seen[chrom]["end"][i]] = get_overlap(start,end,seen[chrom]["end"][i],i)
        i+=1
    # May be spanning within a block
    for item in seen[chrom]["start"]:
        if item <= start and seen[chrom]["start"][item] >= end:
            ret[seen[chrom]["start"][item]] = get_overlap(start,end,item,seen[chrom]["start"][item])
        if item >= start and seen[chrom]["start"][item] <= end:
            ret[seen[chrom]["start"][item]] = get_overlap(start,end,seen[chrom]["end"][item],item)
    return ret

def mapped_end_to_end(cigarstring):
    # Count up the position on the read until we get to a deletion
    count = 0
    read_start = 0
    read_end = 0
    match_count = 0
    large_indel = []
    clipped_bases = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
            match_count += int(cg[:cg.find("M")])
        if cg.endswith('='):
            count += int(cg[:cg.find("=")])
            match_count += int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
            match_count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
            if int(cg[:cg.find("I")]) > 50:
                large_indel.append([count, "I", int(cg[:cg.find("I")])])
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) > 50:
                large_indel.append([count, "D", int(cg[:cg.find("D")])])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            if match_count > 0 and read_end == 0:
                read_end = count
            count += int(cg[:cg.find("S")])
            clipped_bases += int(cg[:cg.find("S")])
            if match_count == 0:
                read_start = count
        elif cg.endswith('H'):
            if match_count > 0 and read_end == 0:
                read_end = count
            count += int(cg[:cg.find("H")])
            clipped_bases += int(cg[:cg.find("H")])
            if match_count == 0:
                read_start = count
    mapped = True
    if len(large_indel) > 0 or clipped_bases/count > 0.1:
        mapped = False
    ret = []
    pos = read_start
    for i in range(len(large_indel)):
        ret.append([pos, large_indel[i][0]])
        if large_indel[i][1] == "I":
            pos = large_indel[i][2]+large_indel[i][0]
        else:
            pos = large_indel[i][0]
    ret.append([pos,read_end])
    return (mapped, ret)
    
def get_family(annotation):
    if "L1" in annotation or "LINE" in annotation:
        return "LINE"
    elif "Alu" in annotation:
        return "Alu"
    elif "SVA" in annotation:
        return "SVA"
    else:
        return "Ambiguous"

def get_neb(avg_size, read_lengths, min_read_len, read_positions):
    total = 0
    for read in read_lengths:
        if read in read_positions:
            total += max(0, read_lengths[read]-1000-avg_size)
    return total/1000000000

def get_without_flank(cigarstring, read_start, read_end):
    read_count = 0
    ref_count = 0
    ref_start = -1
    ref_end = -1
    clipped_bases = 0
    last_op = 0
    for cg in re.findall('[0-9]*[A-Z=]', cigarstring):
        if cg.endswith('M'):
            read_count += int(cg[:cg.find("M")])
            ref_count += int(cg[:cg.find("M")])
            last_op = int(cg[:cg.find("M")])
        if cg.endswith('='):
            read_count += int(cg[:cg.find("=")])
            ref_count += int(cg[:cg.find("=")])
            last_op = int(cg[:cg.find("=")])
        elif cg.endswith('X'):
            read_count += int(cg[:cg.find("X")])
            ref_count += int(cg[:cg.find("X")])
            last_op = int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            read_count += int(cg[:cg.find("I")])
        elif cg.endswith('D'):
            ref_count += int(cg[:cg.find("D")])
        elif cg.endswith('S'):
            read_count += int(cg[:cg.find("S")])
            clipped_bases += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            read_count += int(cg[:cg.find("H")])
            clipped_bases += int(cg[:cg.find("H")])
        if read_count >= read_start and ref_start < 0:
            if ref_count > 0:
                if cg.endswith('I'):
                    ref_start = ref_count
                else:
                    ref_start = ref_count - (read_count - read_start)
            else:
                ref_start = 0
        if read_count >= read_end and ref_end < 0:
            if cg.endswith("H") or cg.endswith("S") or cg.endswith("I"):
                ref_end = ref_count
            else:
                ref_end = ref_count - (read_count - read_end)
    return ref_start, ref_end


parser = argparse.ArgumentParser( description='Identify polymorphic insertions using the reads aligned to the assembly haplotypes')
parser.add_argument('--tsv', required=True)
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--min-read-length', type=int, default=7000)
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument('--read-lengths', type=str, required=True)
parser.add_argument('--chrom', required=True)
args = parser.parse_args()

#giab_regions = {}
#with open(args.bed, 'r') as in_bed:
#    for line in in_bed:
#        row = line.strip().split('\t')
#        if row[0] == args.chrom:
#            if row[0] not in giab_regions:
#                giab_regions[row[0]] = {}
#                giab_regions[row[0]]["start"] = defaultdict(int)
#                giab_regions[row[0]]["end"] = defaultdict(int)
#            #print(row, file=sys.stderr)
#            giab_regions[row[0]]["start"][int(row[1])] = int(row[2])
#            giab_regions[row[0]]["end"][int(row[2])] = int(row[1])

#print("BED Input", file=sys.stderr)

read_lengths = {}
read_total_hc_bases = 0
read_flank_hc_bases = 0
bam_reader = pysam.AlignmentFile(args.bam)
count = 0
read_query_aln = {}
for record in bam_reader.fetch(args.chrom):
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.mapping_quality < 20 or record.infer_read_length() < args.min_read_length:
        continue
    if record.query_name in read_query_aln:
        read_query_aln[record.query_name] = max(read_lengths[record.query_name], record.query_alignment_length)
    else:
        read_query_aln[record.query_name] = record.query_alignment_length
    read_lengths[record.query_name] = record.infer_read_length()

with open(args.read_lengths, 'w') as out_reads:
    for read in read_query_aln:
        out_reads.write(read+"\t"+str(read_query_aln[read])+"\n")

