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
parser.add_argument('--tsv-suffix', required=True)
parser.add_argument('--tsv-folder', required=True)
parser.add_argument('--bam-suffix', required=True)
parser.add_argument('--bam-folder', required=True)
parser.add_argument('--sample-list', required=True)
parser.add_argument('--hap-folder', required=True)
parser.add_argument('--hap1-suffix', type=str, required=True)
parser.add_argument('--hap2-suffix', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=1000)
args = parser.parse_args()

samples = args.sample_list.split(',')
sample_to_bam = {}
sample_to_tsv = {}
for sample in samples:
    sample_to_bam[sample] = sample+"/"+args.bam_folder+"/"+sample+args.bam_suffix
    sample_to_tsv[sample] = args.tsv_folder+"/"+sample+args.tsv_suffix

#seen = {}
#all_positions_list = defaultdict(list)
#for sample in sample_to_tsv:
#    tsv = sample_to_tsv[sample]
#    with open(tsv, 'r') as in_tsv:
#        count = 0
#        for line in in_tsv:
#            if count == 0:
#                count = 1
#                continue
#            row = line.strip().split('\t')
#            chrom = row[0]
#            start = int(row[1])
#            end =  int(row[2])
#            all_positions_list[row[0]].append((int(row[1])-500, int(row[2])+500))
#            read_insert = sample+":"+row[3]+":+:"+row[4]+"-"+row[5]
#            if chrom not in seen:
#                seen[chrom] = defaultdict(list)
#            updated_insert = read_insert+":"+chrom+":"+row[1]+"-"+row[2]
#            seen[chrom][start].append(updated_insert)
#
#print("Tsvs read in", file=sys.stderr)
#
#for chrom in all_positions_list:
#    sorted_by_lower_bound = sorted(all_positions_list[chrom], key=lambda tup: tup[0])
#    merged = []
#    for higher in sorted_by_lower_bound:
#        if not merged:
#            merged.append(higher)
#        else:
#            lower = merged[-1]
#            # test for intersection between lower and higher:
#            # we know via sorting that lower[0] <= higher[0]
#            if higher[0] <= lower[1]:
#                upper_bound = max(lower[1], higher[1])
#                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
#            else:
#                merged.append(higher)
#    all_positions_list[chrom] = merged
#
#all_positions = defaultdict(IntervalTree)
#for chrom in all_positions_list:
#    for pos in all_positions_list[chrom]:
#        start = pos[0]
#        end = pos[1]
#        all_positions[chrom][start:end] = 1
#print("All positions", file=sys.stderr)
#
#bam_ref_positions = {}
#softclips = {}
#count = 0
#for sample in sample_to_bam:
#    count = 0
#    reader = pysam.AlignmentFile(sample_to_bam[sample])
#    print(sample, file=sys.stderr)
#    for record in reader.fetch():
#        count += 1
#        if count % 500000 == 0:
#            print(count, file=sys.stderr)
#        cg_tuples = record.cigartuples
#        if cg_tuples[0][0] == 4 or cg_tuples[0][0] == 5:
#            # Record starts with a hard or softclip
#            if cg_tuples[0][1] >= 100:
#                # Have a clip that is larger than 100bp, check to see if it intersects a position we care about
#                nearby = get_nearby(record.reference_name, record.reference_start-50, record.reference_start+50, seen)
#                if len(nearby) > 0:
#                    # Have a nearby softclip, add it to the list for the position
#                    if record.reference_name not in softclips:
#                        softclips[record.reference_name] = defaultdict(set)
#                    for item in nearby:
#                        for j in nearby[item]:
#                            pos = j.split(':')[5]
#                            softclips[record.reference_name][pos].add(record.query_name)
#        elif cg_tuples[len(cg_tuples)-1][0] == 4 or cg_tuples[len(cg_tuples)-1][0] == 5:
#            # Record starts with a hard or softclip
#            if cg_tuples[len(cg_tuples)-1][1] >= 100:
#                # Have a clip that is larger than 100bp, check to see if it intersects a position we care about
#                nearby = get_nearby(record.reference_name, record.reference_end-50, record.reference_end+50, seen)
#                if len(nearby) > 0:
#                    # Have a nearby softclip, add it to the list for the position
#                    if record.reference_name not in softclips:
#                        softclips[record.reference_name] = defaultdict(set)
#                    for item in nearby:
#                        for j in nearby[item]:
#                            pos = j.split(':')[5]
#                            softclips[record.reference_name][pos].add(record.query_name)
#        nearby = all_positions[record.reference_name][record.reference_start:record.reference_end]
#        if len(nearby) > 0:
#            # Have a record in a region we care about
#            if record.reference_name not in bam_ref_positions:
#                bam_ref_positions[record.reference_name] = defaultdict(list)
#            for item in nearby:
#                #print(record.query_name+" "+record.reference_name+":"+str(item.begin)+"-"+str(item.end))
#                bam_ref_positions[record.reference_name][str(item.begin)+"-"+str(item.end)].append([record.query_name, record.mapping_quality, record.reference_start, record.reference_end])
#
#print("Bam Lists", file=sys.stderr)

sample_to_hap1 = {}
sample_to_hap2 = {}
for sample in samples:
    sample_to_hap1[sample] = args.hap_folder+"/"+sample+args.hap1_suffix
    sample_to_hap2[sample] = args.hap_folder+"/"+sample+args.hap2_suffix

count = 0
for sample in sample_to_hap1:
    polymorphic_reads = {}
    read_positions = defaultdict(list)
    # For each hap read in the alignments of the reads to the hap
    # If there is an end to end alignment of the read to either haplotype add it to the set
    reader = pysam.AlignmentFile(sample_to_hap1[sample])
    for record in reader.fetch():
        mapped, start, end = mapped_end_to_end(record.cigarstring)
        if mapped:
            polymorphic_reads[record.query_name] = 1
        read_positions[record.query_name].append([start,end])
    reader = pysam.AlignmentFile(sample_to_hap2[sample])
    for record in reader.fetch():
        mapped, start, end = mapped_end_to_end(record.cigarstring)
        if mapped:
            polymorphic_reads[record.query_name] = 1
        read_positions[record.query_name].append([start,end])
    # Read in tsv of inserts, If read mapped end to end with no insertion call made it is polymorphic
    # Add a polymorphic tag here and output
    with open(sample_to_tsv[sample], 'r') as in_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            if count == 0:
                row.append("AssemblyAligned\tSampleSimilarity")
                count = 1
                print("Sample\t"+'\t'.join(row))
                continue
            if "PASS" not in line:
                row.append("NA")
                row.append("NA")
                print(sample+"\t"+"\t".join(row))
                continue
            polymorphic = False
            count += 1
            reads = [row[3]]
            for read in reads:
                if read in polymorphic_reads:
                    polymorphic = True
                read_insert_start = int(row[4])
                read_insert_end = int(row[5])
                for position in read_positions[record.query_name]:
                    if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                        polymorphic = True
            if polymorphic:
                row.append("AssemblyPolymorphic")
            else:
                row.append("NotInAssembly")
            #row.append(check_polymorphic(row[0], int(row[1]), int(row[2]), seen, softclips, bam_ref_positions))
            print(sample+"\t"+"\t".join(row))


