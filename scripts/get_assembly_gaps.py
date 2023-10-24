import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict
from intervaltree import Interval, IntervalTree



def mapped_end_to_end(cigarstring):
    # Count up the position on the read until we get to a deletion
    count = 0
    read_start = 0
    read_end = 0
    match_count = 0
    large_indel = False
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
                large_indel = True
        elif cg.endswith('D'):
            if int(cg[:cg.find("D")]) > 50:
                large_indel = True
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
    if large_indel or clipped_bases/count > 0.1:
        return (False, read_start, read_end)
    return (True, read_start, read_end)

def get_nearby(chrom, start, end, seen):
    i = start
    count = 0
    nearby = {}
    if chrom not in seen:
        return nearby
    while i < end:
        if i in seen[chrom]:
            nearby[count] = seen[chrom][i]
            count += 1
        i += 1
    return nearby

def get_nearby_with_position(chrom, start, end, seen):
    i = start
    nearby = {}
    count = 0
    if chrom not in seen:
        return nearby
    while i < end:
        if i in seen[chrom]:
            nearby[count] = [seen[chrom][i], i]
            count += 1
        i += 1
    return nearby
    
def get_polymorphic_in_window(chrom, start, end, seen, softclips, bam_ref_positions, window, cutoff, total_read_count, nearby):
    insert_reads = {}
    for i in nearby:
        for j in nearby[i][0]:
            if abs(nearby[i][1]-start) < window or abs(nearby[i][1]-end) < window:
                read = j.split(':')[1]
                insert_reads[read] = 1
    insert_count = len(insert_reads)
    softclip_count = 0
    if chrom in softclips:
        softclip_count = len(softclips[chrom][str(start)+":"+str(end)])
    if total_read_count == 0:
        return "NA"
    if (insert_count+softclip_count)/total_read_count > cutoff:
        return "Polymorphic_"+str(window)+"_"+str((insert_count+softclip_count)/total_read_count)
    return "NA"

def check_polymorphic(chrom, start, end, seen, softclips, bam_ref_positions):
    total_reads = {}
    for pos in bam_ref_positions[chrom]:
        pos_start = pos.split('-')[0]
        pos_end = pos.split('-')[1]
        if len(pos_start) == 0:
            continue
        if len(pos_end) == 0:
            continue
        pos_start = int(pos_start)
        pos_end = int(pos_end)
        if start > pos_start and end < pos_end:
            for item in bam_ref_positions[chrom][pos]:
                if item[2] < start and item[3] > end:
                    total_reads[item[0]] = 1
    total_read_count = len(total_reads)
    nearby = get_nearby_with_position(chrom, start-500, end+500, seen)
    # 500, 0.8
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, bam_ref_positions, 500, 0.8, total_read_count, nearby)
    if ret != "NA":
        return ret
    # 200, 0.5
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, bam_ref_positions, 200, 0.5, total_read_count, nearby)
    if ret != "NA":
        return ret
    # 100, 0.3
    ret = get_polymorphic_in_window(chrom, start, end, seen, softclips, bam_ref_positions, 100, 0.3, total_read_count, nearby)
    if ret != "NA":
        return ret
    return "NA"

parser = argparse.ArgumentParser( description='Identify polymorphic insertions using the reads aligned to the assembly haplotypes')
parser.add_argument('--tsv', required=True)
parser.add_argument('--hap1-bam', required=True)
parser.add_argument('--hap2-bam', required=True)
parser.add_argument('--max-distance', type=int, default=1000)
parser.add_argument('--fastq-suffix', required=True)
parser.add_argument('--fastq-folder', required=True)
parser.add_argument('--min-read-length', type=int, default=7000)
args = parser.parse_args()

mapped_regions = defaultdict(IntervalTree)
reader = pysam.AlignmentFile(args.hap1_bam)
for record in reader.fetch():
    mapped_regions[record.reference_name][record.reference_start:record.reference_end] = record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end)+"_"+record.query_name+":"+str(record.query_alignment_start)+"-"+str(record.query_alignment_end)
reader = pysam.AlignmentFile(args.hap2_bam)
for record in reader.fetch():
    mapped_regions[record.reference_name][record.reference_start:record.reference_end] = record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end)+"_"+record.query_name+":"+str(record.query_alignment_start)+"-"+str(record.query_alignment_end)

count = 0
inserts = defaultdict(IntervalTree)
read_lens = {}
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0 or row[2] == "reference_insertion_start" or "PASS" not in line:
            count = 1
            continue
        chrom_start = int(row[2]) - 500
        if chrom_start < 0:
            chrom_start = 0
        chrom_end = int(row[3]) + 500
        inserts[row[1]][chrom_start:chrom_end] = row[0]
        if row[0] not in read_lens:
            read_lens[row[0]] = {}
            # open up fastq
            fh = pysam.FastxFile(filename=row[0]+"/"+args.fastq_folder+"/"+row[0]+args.fastq_suffix)
            for entry in fh:
                read_lens[row[0]][entry.name] = len(entry.sequence)

    


count = 0
with open(args.tsv, 'r') as in_tsv:
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0 or row[2] == "reference_insertion_start":
            row.append("AssemblyCoverage\tSingleSample\tReadLenFilter")
            count = 1
            print('\t'.join(row))
            continue
        if "PASS" not in line:
            row.append("NA")
            row.append("NA")
            row.append("NA")
            continue
        chrom_start = int(row[2]) - args.max_distance
        if chrom_start < 0:
            chrom_start = 0
        chrom_end = int(row[3]) + args.max_distance
        nearby = mapped_regions[row[1]][chrom_start:chrom_end]
        mapped_count = 0
        mapped_str = ""
        for item in nearby:
            if mapped_count != 0:
                mapped_str += ";"
            mapped_str += item.data
            mapped_count += 1
        if mapped_count == 0:
            mapped_str = "NA"
        row.append(mapped_str)
        nearby = inserts[row[1]][int(row[2]):int(row[3])]
        samples = {}
        for item in nearby:
            samples[item.data] = 1
        if len(samples) > 1:
            row.append("MultiSample")
        else:
            row.append("SingleSample")
        if read_lens[row[0]][row[4]] < args.min_read_length:
            row.append("MinReadLength")
        else:
            row.append("PassReadLength")
        print("\t".join(row))


