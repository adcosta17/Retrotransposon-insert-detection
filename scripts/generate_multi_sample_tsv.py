import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import threading
import copy
import pickle

def wrapper(func, args, res):
    res.append(func(*args))

def is_fail(pf_annotation):
    # Base case if pass
    if pf_annotation == "PASS":
        return False
    pf_annotation_split = pf_annotation.split(",")
    count = 0
    filters = ["in_centromere", "mapq_fraction"]
    for item in pf_annotation_split:
        if item in filters:
            count += 1
    if count == len(pf_annotation_split):
        return False
    return True

def get_soft_clip(cigartuples, reverse=False):
    # Gets the read length based on the Cigar String
    count = 0
    if reverse:
        cigartuples.reverse()
    for cg in cigartuples:
        if cg[0] == 4:
            count += cg[1]
        elif cg[0] == 5:
            count += cg[1]
        else:
            break
    return count

def same_family(annotation1, annotation2):
    if "ambiguous" in annotation1 or "ambiguous" in annotation2:
        return True
    if "LINE" in annotation1 and "LINE" in annotation2:
        return True
    if "SVA" in annotation1 and ("SVA" in annotation2 or "SINE" in annotation2):
        return True
    if "SINE" in annotation1 and ("SVA" in annotation2 or "SINE" in annotation2):
        return True
    if "ERV" in annotation1 and "ERV" in annotation2:
        return True
    return False

#from collections import namedtuple
#Arguments = namedtuple('Arguments',['sample','suffix','folder', 'bam_suffix', 'bam_folder','max_distance', 'threads'])
#args = Arguments('DDTS_0001_nn_n,DDTS_0002_nn_n,DDTS_0003_nn_n,DDTS_0004_nn_n,DDTS_0005_nn_n,DDTS_0006_nn_n,DDTS_0007_nn_n,DDTS_0008_nn_n,DDTS_0009_nn_n,DDTS_0010_nn_n,DDTS_0011_nn_n', '.all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv', 'winnow_read_analysis', '.sorted.phased.bam', 'winnow_phased', 20, 12)

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--sample', required=True)
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv")
parser.add_argument('--folder', default="read_analysis")
parser.add_argument('--bam-suffix', default=".sorted.phased.bam")
parser.add_argument('--bam-folder', default="phased")
parser.add_argument('--merged-suffix', default=".all.merged_reads.txt")
parser.add_argument('--merged-folder', default="filtered_mapped")
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--filters-to-exclude', default="NA")
parser.add_argument('--threads', type=int, default=1)
args = parser.parse_args()

threads_to_use = 1
if args.threads > 1:
    threads_to_use = args.threads

# Split string into list:
base_tree = {}

def get_base_tree(sample, folder, suffix):
    sample_tsv = "./"+sample+"/"+folder+"/"+sample+suffix
    ret_tree = {}
    print(sample_tsv,file=sys.stderr)
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            # get insertion coordinates, insert into intervaltree allowing for the maximum allowable distance
            if is_fail(row_args[8]):
                continue
            chrom = row_args[0]
            start = int(row_args[1])
            end = int(row_args[2])
            if chrom not in ret_tree:
                ret_tree[chrom] = {}
            if start not in ret_tree[chrom]:
                ret_tree[chrom][start] = defaultdict(list)
            ret_tree[chrom][start][end-start] = [row_args[9], sample+"_HP:"+str(row_args[19])]
    return ret_tree

results = []
thread_list = []
if threads_to_use > 1:
    for sample in args.sample.split(','):
        print(sample,file=sys.stderr)
        t = threading.Thread(target=wrapper, args=(get_base_tree, (sample, args.folder, args.suffix,), results))
        t.start()
        thread_list.append(t)
    for t in thread_list:
        t.join()
else:
    for sample in args.sample.split(','):
        results.append(get_base_tree(sample, args.folder, args.suffix))

chroms = {}
for res in results:
    for chrom in res:
        chroms[chrom] = 1

for res in results:
    for chrom in chroms:
        if chrom not in base_tree:
            base_tree[chrom] = {}
        if chrom in res:
            for start in res[chrom]:
                if start not in base_tree[chrom]:
                    base_tree[chrom][start] = res[chrom][start]
                else:
                    for pos in res[chrom][start]:
                        if pos in base_tree[chrom][start]:
                            base_tree[chrom][start][pos][1] = base_tree[chrom][start][pos][1]+","+res[chrom][start][pos][1]
                        else:
                            base_tree[chrom][start][pos] = res[chrom][start][pos]


print("Completed Pass 1",file=sys.stderr)
intervaltrees = {}
insert_sizes = defaultdict(IntervalTree)

def get_intervals_1(sample, folder, suffix, max_distance, base_tree):
    sample_tsv = "./"+sample+"/"+folder+"/"+sample+suffix
    interval_ret = {}
    print(sample_tsv,file=sys.stderr)
    insert_ret = defaultdict(IntervalTree)
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if is_fail(row_args[8]):
                continue
            chrom = row_args[0]
            start = int(row_args[1]) - max_distance
            end = int(row_args[2]) + max_distance
            # Check the base_tree to see if there are any inserts that fall within range
            nearby = []
            for i in range(start, end+1):
                if i in base_tree[chrom]:
                    # found a hit
                    for pos in base_tree[chrom][i]:
                        nearby.append(base_tree[chrom][i][pos])
            annotation = row_args[9]
            data_dict = defaultdict(int)
            for item in nearby:
                if same_family(annotation, item[0]):
                    data = item[1].split(',')
                    for d in data:
                        data_dict[d] += 1
            data_list = []
            for item in data_dict:
                data_list.append(item+"_"+str(data_dict[item]))
            start = int(row_args[1])
            end = int(row_args[2])
            if chrom not in interval_ret:
                interval_ret[chrom] = {}
            if start not in interval_ret[chrom]:
                interval_ret[chrom][start] = defaultdict(list)
            if end-start not in interval_ret[chrom][start]:
                interval_ret[chrom][start][end-start] = [annotation, ",".join(data_list)]
            insert_ret[chrom][start:end] = int(row_args[5]) - int(row_args[4])
    return (interval_ret, insert_ret)

results = []
thread_list = []
if threads_to_use > 1:
    for sample in args.sample.split(','):
        print(sample,file=sys.stderr)
        t = threading.Thread(target=wrapper, args=(get_intervals_1, (sample, args.folder, args.suffix, args.max_distance, base_tree,), results))
        t.start()
        thread_list.append(t)
    for t in thread_list:
        t.join()
else:
    for sample in args.sample.split(','):
        results.append(get_intervals_1(sample, args.folder, args.suffix, args.max_distance, base_tree))

for res in results:
    for chrom in chroms:
        insert_sizes[chrom] |= res[1][chrom]
        if chrom not in intervaltrees:
            intervaltrees[chrom] = {}
        if chrom in res[0]:
            for start in res[0][chrom]:
                if start not in intervaltrees[chrom]:
                    intervaltrees[chrom][start] = res[0][chrom][start]
                else:
                    for pos in res[0][chrom][start]:
                        if pos not in intervaltrees[chrom][start]:
                            intervaltrees[chrom][start][pos] = res[0][chrom][start][pos]

print("Completed Pass 2",file=sys.stderr)

def get_intervals_2(sample, folder, suffix, intervaltrees, max_distance):
    sample_tsv = "./"+sample+"/"+folder+"/"+sample+suffix
    local_tree = {}
    print(sample_tsv,file=sys.stderr)
    for chrom in intervaltrees:
        local_tree[chrom] = {}
        for start in intervaltrees[chrom]:
            local_tree[chrom][start] = defaultdict(str)
            for pos in intervaltrees[chrom][start]:
                local_tree[chrom][start][pos] = ""
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            #if count % 1000000 == 0:
            #    print(sample + " " + str(count), file=sys.stderr)
            #count += 1
            # get insertion coordinates, insert into intervaltree allowing for the maximum allowable distance
            if is_fail(row_args[8]):
                chrom = row_args[0]
                start = int(row_args[1]) - max_distance
                end = int(row_args[2]) + max_distance
                if chrom not in local_tree:
                    continue
                # Check the base_tree to see if there are any inserts that fall within range
                for i in range(start, end+1):
                    if i in local_tree[chrom]:
                        # found a hit
                        for pos in local_tree[chrom][i]:
                            if sample in local_tree[chrom][i][pos]:
                                vals = local_tree[chrom][i][pos].split('_')
                                n = int(vals[len(vals)-1])+1
                                local_tree[chrom][i][pos] = sample+"_fail_"+str(n)
                            else:
                                local_tree[chrom][i][pos] = sample+"_fail_1"
    return local_tree

results = []
thread_list = []
if threads_to_use > 1:
    for sample in args.sample.split(','):
        t = threading.Thread(target=wrapper, args=(get_intervals_2, (sample, args.folder, args.suffix, intervaltrees, args.max_distance,), results))
        t.start()
        thread_list.append(t)
    for t in thread_list:
        t.join()
else:
    for sample in args.sample.split(','):
        results.append(get_intervals_2(sample, args.folder, args.suffix, intervaltrees, args.max_distance))

for res in results:
    for chrom in chroms:
        if chrom in res:
            for start in res[chrom]:
                for pos in res[chrom][start]:
                    if pos in intervaltrees[chrom][start]:
                        if res[chrom][start][pos] != "":
                            intervaltrees[chrom][start][pos][1] = intervaltrees[chrom][start][pos][1]+','+res[chrom][start][pos]

print("Completed Pass 3",file=sys.stderr)

## Get merged reads
merged_reads = {}
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    sample_reads = "./"+sample+"/"+args.merged_folder+"/"+sample+args.merged_suffix
    merged_reads[sample] = defaultdict(str)
    with open(sample_reads) as input_reads:
        for line in input_reads:
            read = line.strip()
            merged_reads[sample][read] = 1

# Go through bam to get the haplotype counts for each insert position for each sample
# Count how many are haplotype 0, 1 or 2 for each sample and include that into the output tsv
haplotypes = {}
softclips = {}
regions = {}
for chrom in intervaltrees:
    if chrom not in haplotypes:
        haplotypes[chrom] = {}
    if chrom not in regions:
        regions[chrom] = {}
    if chrom not in softclips:
        softclips[chrom] = {}
    for start in intervaltrees[chrom]:
        for pos in intervaltrees[chrom][start]:
            end = start + pos
            if str(start)+"_"+str(end) not in haplotypes[chrom]:
                haplotypes[chrom][str(start)+"_"+str(end)] = defaultdict(int)
            if str(start)+"_"+str(end) not in softclips[chrom]:
                softclips[chrom][str(start)+"_"+str(end)] = defaultdict(int)
            r_start = int(round(start, -4))
            # Want to round down not up. If we've rounded up round it down
            if r_start > start:
                r_start -= 10000
            r_end = r_start + 10000
            regions[chrom][str(r_start)+'_'+str(r_end)] = 1

print("Haplotypes Setup",file=sys.stderr)

def get_softclip_haplotype(sample, bam_folder, bam_suffix, max_distance, intervaltrees, haplotypes, softclips, insert_sizes, regions, merged_reads):
    local_hap = copy.deepcopy(haplotypes)
    local_soft = copy.deepcopy(softclips)
    local_merged = copy.deepcopy(merged_reads)
    seen_reads = {}
    header = pysam.AlignmentFile("./"+sample+"/"+bam_folder+"/"+sample+bam_suffix).header
    sam_reader = pysam.AlignmentFile("./"+sample+"/"+bam_folder+"/"+sample+bam_suffix)
    for sq in header['SQ']:
        #print(sq['SN'],file=sys.stderr)
        if sq['SN'] not in intervaltrees:
            continue
        #print(sample + " " + sq['SN'], file=sys.stderr)
        for region in regions[sq['SN']]:
            r_start = int(region.split('_')[0])
            r_end = int(region.split('_')[1])
            for record in sam_reader.fetch(sq['SN'], r_start, r_end):
                if record.is_unmapped:
                    continue
                tags = record.get_tags()
                hap = "HP:0"
                for t in tags:
                    if t[0] == "HP":
                        hap = "HP:"+str(t[1])
                # now check each insert in the haplotype dict to see if its in the range of this record
                for i in range(record.reference_start-max_distance, record.reference_end+max_distance+1):
                    start = i
                    if start in intervaltrees[sq['SN']]:
                        for pos in intervaltrees[sq['SN']][start]:
                            end = start + pos
                            if record.query_name in local_merged and record.query_name not in seen_reads:
                                local_hap[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
                                seen_reads[record.query_name] = 1
                                continue
                            elif record.query_name in seen_reads:
                                continue
                            local_hap[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
                            # Check a slightly larger window for consistent soft clip position
                            if (record.reference_start - start <= max_distance and start < record.reference_start) or (start - record.reference_end <= max_distance and start > record.reference_end):
                                insert_size = 0
                                n_items = 0
                                if len(insert_sizes[sq['SN']][start:end]) > 0:
                                    for item in insert_sizes[sq['SN']][start:end]:
                                        insert_size += item.data
                                        n_items += 1
                                    insert_size = insert_size/n_items
                                if insert_size < max_distance:
                                	insert_size = max_distance
                                insert_size = min(insert_size, 100)
                                # Alignment starts or ends nearby. Check if there is a soft clip > insert size
                                if abs(start - record.reference_start) <= max_distance:
                                    size = get_soft_clip(record.cigartuples)
                                    if size > 0.25*insert_size:
                                        # Add softclip for this inserts
                                        local_soft[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
                                elif abs(end - record.reference_end) <= max_distance:
                                    size = get_soft_clip(record.cigartuples, True)
                                    if size > 0.25*insert_size:
                                        # Add softclip for this inserts
                                        local_soft[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
    return (local_soft, local_hap)

results = []
thread_list = []
if threads_to_use > 1:
    for sample in args.sample.split(','):
        print(sample,file=sys.stderr)
        t = threading.Thread(target=wrapper, args=(get_softclip_haplotype, (sample, args.bam_folder, args.bam_suffix, args.max_distance, intervaltrees, haplotypes, softclips, insert_sizes,regions,merged_reads[sample],), results))
        t.start()
        thread_list.append(t)
    for t in thread_list:
        t.join()
else:
    for sample in args.sample.split(','):
        results.append(get_softclip_haplotype(sample, args.bam_folder, args.bam_suffix, args.max_distance, intervaltrees, haplotypes, softclips, insert_sizes,regions,merged_reads[sample]))

for res in results:
    softclips_ret = res[0]
    haplotypes_ret =  res[1]
    for chrom in softclips_ret:
        for region in softclips_ret[chrom]:
            for s in softclips_ret[chrom][region]:
                softclips[chrom][region][s] += softclips_ret[chrom][region][s]
    for chrom in haplotypes_ret:
        for region in haplotypes_ret[chrom]:
            for s in haplotypes_ret[chrom][region]:
                haplotypes[chrom][region][s] += haplotypes_ret[chrom][region][s]


print("Haplotypes Added",file=sys.stderr)

print("chrom\tstart\tend\tannotation\tsamples\thaplotypes\tsoftclips")
for chrom in intervaltrees:
    for start in intervaltrees[chrom]:
        for pos in intervaltrees[chrom][start]:
            end = start + pos
            hap_list = []
            for h in haplotypes[chrom][str(start)+"_"+str(end)]:
                hap_list.append(h+"_"+str(haplotypes[chrom][str(start)+"_"+str(end)][h]))
            soft_list = []
            for s in softclips[chrom][str(start)+"_"+str(end)]:
                soft_list.append(s+"_"+str(softclips[chrom][str(start)+"_"+str(end)][s]))
            out_str = chrom+"\t"+str(start)+"\t"+str(end)+"\t"+intervaltrees[chrom][start][pos][0]+"\t"+intervaltrees[chrom][start][pos][1]+"\t"+",".join(hap_list)+"\t"+",".join(soft_list)
            print(out_str)

