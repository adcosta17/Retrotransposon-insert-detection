import pysam
import argparse
import sys
import csv
import re
from intervaltree import Interval, IntervalTree
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
        return ret
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


def mapped_end_to_end(cigarstring, read_length, reverse):
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
    if not reverse:
        return (mapped, ret)
    new_ret = []
    for item in ret:
        new_ret.append([read_length-item[1],read_length-item[0]])
    return(mapped,new_ret)
    
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
        total += max(0, read_lengths[read]-1000-avg_size)
    return total/1000000000


def get_type(sniffles_id):
    if "INS" in sniffles_id:
        return "INS"
    if "DEL" in sniffles_id:
        return "DEL"
    if "DUP" in sniffles_id:
        return "DUP"
    if "INV" in sniffles_id:
        return "INV"
    if "BND" in sniffles_id:
        return "BND"
    return "NA"


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


parser = argparse.ArgumentParser( description='Get Stats for DDTS Sniffles2 and Insertions HC regions')
parser.add_argument('--samples', required=True)
parser.add_argument('--hap1', type=str, required=True)
parser.add_argument('--hap2', type=str, required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--chroms', required=True)
parser.add_argument('--folder', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--tsv', required=True)
parser.add_argument('--normalized', type=str, required=True)
parser.add_argument('--sniffles', type=str, required=False)
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--min-read-length', type=int, default=7000)
args = parser.parse_args()

total_hc = 1
total_hc_flank = 1
read_query_aln = {}
for chrom in args.chroms.split(','):
    # Get the values from the files
    with open(args.sample+"/"+args.folder+"/"+args.sample+"_reads_"+chrom+".txt", 'r') as in_reads:
        for line in in_reads:
            row = line.strip().split('\t')
            if row[0] in read_query_aln:
                read_query_aln[row[0]] = max(read_query_aln[row[0]], int(row[1]))
            else:
                read_query_aln[row[0]] = int(row[1])

print("hc Input", file=sys.stderr)
print(total_hc, file=sys.stderr)
print(total_hc_flank, file=sys.stderr)

reads_to_use = {}
with open(args.tsv, 'r') as in_tsv:
    count = 0
    for line in in_tsv:
        row = line.strip().split('\t')
        if count == 0:
            count = 1
            continue
        if row[6] == "PASS" and "Polymorphic" not in line:
            reads = row[5].split(',')
            for read in reads:
                read_info = read.split(':')
                read_name = read_info[1]+":"+read_info[2]
                if args.sample == read_info[0]:
                    reads_to_use[read_name] = 1
                    #print(read_name)

print("reads Input", file=sys.stderr)

polymorphic_reads = {}
read_positions = defaultdict(list)
# For each hap read in the alignments of the reads to the hap
# If there is an end to end alignment of the read to either haplotype add it to the set

read_lengths = defaultdict(int)
seen_count = 0
reader = pysam.AlignmentFile(args.hap1)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use:
        continue
    mapped, positions = mapped_end_to_end(record.cigarstring, record.infer_read_length(), record.is_reverse)
    if mapped:
        polymorphic_reads[record.query_name] = 1
        seen_count += 1
    read_positions[record.query_name].extend(positions)
    read_lengths[record.query_name] = record.infer_read_length()
    #print(record.query_name +"\t"+str(record.query_alignment_start) + "\t"+ str(record.query_alignment_end) + "\t"+str(record.is_reverse)+"\t"+str(positions))

print("HAP1 Input", file=sys.stderr)

reader = pysam.AlignmentFile(args.hap2)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use:
        continue
    mapped, positions = mapped_end_to_end(record.cigarstring, record.infer_read_length(), record.is_reverse)
    if mapped:
        polymorphic_reads[record.query_name] = 1
    read_positions[record.query_name].extend(positions)
    read_lengths[record.query_name] = record.infer_read_length()
    #print(record.query_name +"\t"+str(record.query_alignment_start) + "\t"+ str(record.query_alignment_end) + "\t"+str(record.is_reverse)+"\t"+str(positions))

print("HAP2 Input", file=sys.stderr)

counts = {}
single_sample_counts = {}
sizes = defaultdict(list)
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
                    #print(read+"\te2e")
                elif read_name in read_positions:
                    read_insert_start = int(read_info[4].split('-')[0])
                    read_insert_end = int(read_info[4].split('-')[1])
                    if read_info[3] == '-':
                        read_insert_end = read_lengths[read_name] - int(read_info[4].split('-')[0])
                        read_insert_start = read_lengths[read_name] - int(read_info[4].split('-')[1])
                    for position in read_positions[read_name]:
                        if (read_insert_start - position[0]) > args.max_distance and  (position[1] - read_insert_end) > args.max_distance:
                            polymorphic = "AssemblyPolymorphic_Insert"
                            #print(read+"\tPolymorphic\t"+str(position))
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
                sizes[family].append(len(row[4]))
                sizes["All"].append(len(row[4]))
                counts[family] += 1
                counts["All"] += 1
                if single_sample:
                    single_sample_counts[family] += 1
        else:
            row.append("NA\tNA")
        print(args.sample+"\t"+"\t".join(row))

print("reads norm Input", file=sys.stderr)

avg_size = {} 
avg_size["All"] = 0
if counts["All"] > 0:
    avg_size["All"] = sum(sizes["All"])/len(sizes["All"])
avg_size["LINE"] = 0
if counts["LINE"] > 0:
    avg_size["LINE"] = sum(sizes["LINE"])/len(sizes["LINE"])
avg_size["Alu"] = 0
if counts["Alu"] > 0:
    avg_size["Alu"] = sum(sizes["Alu"])/len(sizes["Alu"])
avg_size["SVA"] = 0
if counts["SVA"] > 0:
    avg_size["SVA"] = sum(sizes["SVA"])/len(sizes["SVA"])
avg_size["Ambiguous"] = 0
if counts["Ambiguous"] > 0:
    avg_size["Ambiguous"] = sum(sizes["Ambiguous"])/len(sizes["Ambiguous"])


with open(args.normalized, 'w') as out_norm:
    out_norm.write("Sample\tAll_Count\tAll_Single_Sample\tAll_NEB\tAll_Single_Sample_NEB\tLINE_Count\tLINE_Single_Sample\tLINE_NEB\tLINE_Single_Sample_NEB\tAlu_Count\tAlu_Single_Sample\tAlu_NEB\tAlu_Single_Sample_NEB\tSVA_Count\tSVA_Single_Sample\tSVA_NEB\tSVA_Single_Sample_NEB\tAmbiguous_Count\tAmbiguous_Single_Sample\tAmbiguous_NEB\tAmbiguous_Single_Sample_NEB\n")
    out_norm.write(args.sample)
    for item in avg_size:
        normalized_denom = get_neb(avg_size[item], read_query_aln, args.min_read_length, read_positions)
        #normalized_denom = 1
        if normalized_denom > 0:
            out_norm.write("\t"+str(counts[item])+"\t"+str(single_sample_counts[item])+"\t"+str(counts[item]/normalized_denom)+"\t"+str(single_sample_counts[item]/normalized_denom))
        else:
            out_norm.write("\t"+str(counts[item])+"\t"+str(single_sample_counts[item])+"\t"+str(0)+"\t"+str(0))
    out_norm.write("\n")

print("Norm inserts output", file=sys.stderr)

