import pysam
import argparse
import sys
import re
import csv
from intervaltree import Interval, IntervalTree
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

parser = argparse.ArgumentParser( description='Get Stats for DDTS Sniffles2 pop VCF')
parser.add_argument('--samples', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--vcf', required=True)
parser.add_argument('--bam', required=True)
parser.add_argument('--bed', required=True)
parser.add_argument('--hap1', required=True)
parser.add_argument('--hap2', required=True)
args = parser.parse_args()

giab_regions = {}
with open(args.bed, 'r') as in_bed:
    for line in in_bed:
        row = line.strip().split('\t')
        if row[0] not in giab_regions:
            giab_regions[row[0]] = {}
            giab_regions[row[0]]["start"] = defaultdict(int)
            giab_regions[row[0]]["end"] = defaultdict(int)
        #print(row, file=sys.stderr)
        giab_regions[row[0]]["start"][int(row[1])] = int(row[2])
        giab_regions[row[0]]["end"][int(row[2])] = int(row[1])


read_lens = {}
samples_list = args.samples.split(',')
counts = {}
counts["All"] = defaultdict(int)
counts["INS"] = defaultdict(int)
counts["DEL"] = defaultdict(int)
counts["DUP"] = defaultdict(int)
counts["INV"] = defaultdict(int)
counts["BND"] = defaultdict(int)
counts_uniq = {}
counts_uniq["All"] = defaultdict(int)
counts_uniq["INS"] = defaultdict(int)
counts_uniq["DEL"] = defaultdict(int)
counts_uniq["DUP"] = defaultdict(int)
counts_uniq["INV"] = defaultdict(int)
counts_uniq["BND"] = defaultdict(int)

reads_to_use = {}
with open(args.vcf, 'r') as in_vcf:
    for line in in_vcf:
        if '#' in line:
            continue
        row = line.strip().split('\t')
        if not check_nearby(row[0], int(row[1])-1, int(row[1])+1, giab_regions):
            continue
        for item in row[7].split(';'):
            if "RNAMES" in item:
                name_list = item.split("RNAMES=")[1].split(',')
                for name in name_list:
                    reads_to_use[name] = 1


count = 0
polymorphic_reads = {}
reader = pysam.AlignmentFile(args.hap1)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use:
        continue
    mapped, positions = mapped_end_to_end(record.cigarstring)
    if mapped:
        polymorphic_reads[record.query_name] = 1

print("HAP1 Input", file=sys.stderr)

reader = pysam.AlignmentFile(args.hap2)
for record in reader.fetch():
    count += 1
    if count % 100000 == 0:
        print(count, flush=True, file=sys.stderr)
    if record.query_name not in reads_to_use:
        continue
    mapped, positions = mapped_end_to_end(record.cigarstring)
    if mapped:
        polymorphic_reads[record.query_name] = 1

print("HAP2 Input", file=sys.stderr)

i = 0
sample_bases = {}
sample_bases_hc = {}
read_total_hc_bases = 1
read_lens[args.sample] = defaultdict(int)
reader = pysam.AlignmentFile(args.bam)
i += 1
#for record in reader:
#    if record.mapping_quality > 0:
#        read_lens[args.sample][record.query_name] = record.infer_read_length()
#        ret = check_nearby(record.reference_name, record.reference_start, record.reference_end, giab_regions)
#        if len(ret) == 0:
#            continue
#        aln_total = 0
#        for item in ret:
#            aln_total += ret[item]
#        read_total_hc_bases += aln_total

read_sum = 1
#for read in read_lens[args.sample]:
#    read_sum += read_lens[args.sample][read]
sample_bases[args.sample] = read_sum/1000000000
sample_bases_hc[args.sample] = read_total_hc_bases/1000000000


with open(args.vcf, 'r') as in_vcf:
    for line in in_vcf:
        if '#' in line: 
            continue
        row = line.strip().split('\t')
        if not check_nearby(row[0], int(row[1])-1, int(row[1])+1, giab_regions):
            continue
        sv_type = get_type(row[2])
        if sv_type == "NA":
            continue
        polymorphic = False
        for item in row[7].split(';'):
            if "RNAMES" in item:
                name_list = item.split("RNAMES=")[1].split(',')
                for name in name_list:
                    if name in polymorphic_reads:
                        polymorphic = True
                        break
        if polymorphic:
            continue
        sample_count = defaultdict(int)
        for i in range(len(samples_list)):
            sample = samples_list[i]
            if "NULL" in row[9+i]:
                continue
            counts[sv_type][sample] += 1
            counts["All"][sample] += 1
            sample_count[sample] += 1
        print(line.strip())
        if len(sample_count) == 1:
            for sample in sample_count:
                counts_uniq["All"][sample] += 1
                counts_uniq[sv_type][sample] += 1
                #if sv_type == "INS" and sample == args.sample:
                #    print(line.strip())

print("Type\tSample\tAll\tAllNorm\tAllHC\tInsertions\tInsertionNorm\tInsertionHC\tDeletions\tDeletionsNorm\tDeletionsHC\tDuplications\tDuplicationsNorm\tDuplicationsHC\tInversions\tInversionsNorm\tInversionsHC\tTranslocations\tTranslocationsNorm\tTranslocationsHC")
for sample in samples_list:
    if sample != args.sample:
        continue
    print("Shared\t"+sample+"\t"+str(counts["All"][sample])+"\t"+str(counts["All"][sample]/sample_bases[sample])+"\t"+str(counts["All"][sample]/sample_bases_hc[sample])+"\t"+str(counts["INS"][sample])+"\t"+str(counts["INS"][sample]/sample_bases[sample])+"\t"+str(counts["INS"][sample]/sample_bases_hc[sample])+"\t"+str(counts["DEL"][sample])+"\t"+str(counts["DEL"][sample]/sample_bases[sample])+"\t"+str(counts["DEL"][sample]/sample_bases_hc[sample])+"\t"+str(counts["DUP"][sample])+"\t"+str(counts["DUP"][sample]/sample_bases[sample])+"\t"+str(counts["DUP"][sample]/sample_bases_hc[sample])+"\t"+str(counts["INV"][sample])+"\t"+str(counts["INV"][sample]/sample_bases[sample])+"\t"+str(counts["INV"][sample]/sample_bases_hc[sample])+"\t"+str(counts["BND"][sample])+"\t"+str(counts["BND"][sample]/sample_bases[sample])+"\t"+str(counts["BND"][sample]/sample_bases_hc[sample]))
    print("Unique\t"+sample+"\t"+str(counts_uniq["All"][sample])+"\t"+str(counts_uniq["All"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["All"][sample]/sample_bases_hc[sample])+"\t"+str(counts_uniq["INS"][sample])+"\t"+str(counts_uniq["INS"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["INS"][sample]/sample_bases_hc[sample])+"\t"+str(counts_uniq["DEL"][sample])+"\t"+str(counts_uniq["DEL"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["DEL"][sample]/sample_bases_hc[sample])+"\t"+str(counts_uniq["DUP"][sample])+"\t"+str(counts_uniq["DUP"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["DUP"][sample]/sample_bases_hc[sample])+"\t"+str(counts_uniq["INV"][sample])+"\t"+str(counts_uniq["INV"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["INV"][sample]/sample_bases_hc[sample])+"\t"+str(counts_uniq["BND"][sample])+"\t"+str(counts_uniq["BND"][sample]/sample_bases[sample])+"\t"+str(counts_uniq["BND"][sample]/sample_bases_hc[sample]))


