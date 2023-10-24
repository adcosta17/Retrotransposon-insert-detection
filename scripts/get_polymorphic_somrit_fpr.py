import pysam
import argparse
import sys
import csv
import re
from collections import defaultdict



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

parser = argparse.ArgumentParser( description='Identify polymorphic insertions using the reads aligned to the assembly haplotypes')
parser.add_argument('--tsv', required=True)
parser.add_argument('--hap1', type=str, required=True)
parser.add_argument('--hap2', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--min-read-length', type=int, default=7000)
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument('--normalized', type=str, required=True)
args = parser.parse_args()

read_lengths = {}
bam_reader = pysam.AlignmentFile(args.bam)
count = 0
for record in bam_reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name in read_lengths and record.mapping_quality >= 20 and record.infer_read_length() >= args.min_read_length:
        read_lengths[record.query_name] = max(read_lengths[record.query_name], record.query_alignment_length)
    elif record.mapping_quality >= 20 and record.infer_read_length() >= args.min_read_length:
        read_lengths[record.query_name] = record.query_alignment_length


polymorphic_reads = {}
read_positions = defaultdict(list)
# For each hap read in the alignments of the reads to the hap
# If there is an end to end alignment of the read to either haplotype add it to the set

reader = pysam.AlignmentFile(args.hap1)
for record in reader.fetch():
    mapped, positions = mapped_end_to_end(record.cigarstring)
    if mapped:
        polymorphic_reads[record.query_name] = 1
    read_positions[record.query_name].extend(positions)

reader = pysam.AlignmentFile(args.hap2)
for record in reader.fetch():
    mapped, positions = mapped_end_to_end(record.cigarstring)
    if mapped:
        polymorphic_reads[record.query_name] = 1
    read_positions[record.query_name].extend(positions)

# Read in tsv of inserts, If read mapped end to end with no insertion call made it is polymorphic
# Add a polymorphic tag here and output

counts = {}
single_sample_counts = {}
sizes = {}
sizes["All"] = 0
sizes["LINE"] = 0
sizes["Alu"] = 0
sizes["SVA"] = 0
sizes["Ambiguous"] = 0
counts["LINE"] = 0
counts["All"] = 0
counts["Alu"] = 0
counts["SVA"] = 0
counts["Ambiguous"] = 0
single_sample_counts["LINE"] = 0
single_sample_counts["All"] = 0
single_sample_counts["Alu"] = 0
single_sample_counts["SVA"] = 0
single_sample_counts["Ambiguous"] = 0

with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            row.append("AssemblyAligned\tSamples")
            count = 1
            print("Sample\t"+'\t'.join(row))
            continue
        if args.sample not in row[5]:
            continue
        if row[6] == "PASS" and "Polymorphic" not in line:
            polymorphic = "PossiblyNovel_NotInAssembly"
            reads = row[5].split(',')
            read_samples = {}
            for read in reads:
                read_info = read.split(':')
                read_name = read_info[1]+":"+read_info[2]
                read_samples[read_info[0]] = 1
                if "Polymorphic" in polymorphic:
                    continue
                if read_name in polymorphic_reads:
                    polymorphic = "AssemblyPolymorphic_E2E"
                elif read_name in read_positions:
                    read_insert_start = int(read_info[4].split('-')[0])
                    read_insert_end = int(read_info[4].split('-')[1])
                    for position in read_positions[read_name]:
                        if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                            polymorphic = "AssemblyPolymorphic_Insert"
                else:
                    polymorphic = "NotInAssembly"
            row.append(polymorphic)
            single_sample = False
            if len(read_samples) == 1:
                row.append("SingleSample")
                single_sample = True
            else:
                row.append("MultiSample")
            if polymorphic == "PossiblyNovel_NotInAssembly":
                family = get_family(row[7])
                # Have a passing insert
                counts[family] += 1
                sizes[family] += len(row[4])
                counts["All"] += 1
                sizes["All"] += len(row[4])
                if single_sample:
                    single_sample_counts[family] += 1
                    single_sample_counts["All"] += 1
        else:
            row.append("NA\tNA")
        print(args.sample+"\t"+"\t".join(row))

avg_size = {} 
avg_size["All"] = 0
if counts["All"] > 0:
    avg_size["All"] = (sizes["All"])/(counts["All"])
avg_size["LINE"] = 0
if counts["LINE"] > 0:
    avg_size["LINE"] = (sizes["LINE"])/(counts["LINE"])
avg_size["Alu"] = 0
if counts["Alu"] > 0:
    avg_size["Alu"] = (sizes["Alu"])/(counts["Alu"])
avg_size["SVA"] = 0
if counts["SVA"] > 0:
    avg_size["SVA"] = (sizes["SVA"])/(counts["SVA"])
avg_size["Ambiguous"] = 0
if counts["Ambiguous"] > 0:
    avg_size["Ambiguous"] = (sizes["Ambiguous"])/(counts["Ambiguous"])
with  open(args.normalized, 'w') as out_norm:
    out_norm.write("Sample\tAll_Count\tAll_Single_Sample\tAll_NEB\tAll_Single_Sample_NEB\tLINE_Count\tLINE_Single_Sample\tLINE_NEB\tLINE_Single_Sample_NEB\tAlu_Count\tAlu_Single_Sample\tAlu_NEB\tAlu_Single_Sample_NEB\tSVA_Count\tSVA_Single_Sample\tSVA_NEB\tSVA_Single_Sample_NEB\tAmbiguous_Count\tAmbiguous_Single_Sample\tAmbiguous_NEB\tAmbiguous_Single_Sample_NEB\n")
    out_norm.write(args.sample)
    for item in avg_size:
        normalized_denom = get_neb(avg_size[item], read_lengths, args.min_read_length, read_positions)
        out_norm.write("\t"+str(counts[item])+"\t"+str(single_sample_counts[item])+"\t"+str(counts[item]/normalized_denom)+"\t"+str(single_sample_counts[item]/normalized_denom))
    out_norm.write("\n")
