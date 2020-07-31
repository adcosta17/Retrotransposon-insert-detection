import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict
from multiprocessing.pool import ThreadPool


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


def get_read_length(cigartuples):
    # Gets the read length based on the Cigar String
    count = 0
    for cg in cigartuples:
        if cg[0] == 0:
            count += cg[1]
        elif cg[0] == 5:
            count += cg[1]
        elif cg[0] == 1:
            count += cg[1]
        elif cg[0] == 6:
            count += cg[1]
        elif cg[0] == 4:
            count += cg[1]
        elif cg[0] == 8:
            count += cg[1]
    return count

def get_read_length_no_sc(cigartuples):
    count = 0
    for cg in cigartuples:
        if cg[0] == 0:
            count += cg[1]
        elif cg[0] == 1:
            count += cg[1]
        elif cg[0] == 6:
            count += cg[1]
        elif cg[0] == 8:
            count += cg[1]
    return count
    

def get_start_end_string(record_1, record_2):
    # Assume the cigars are oriented such the alignment is cigar1 - indel - cigar2
    # Remove any hard or soft clipped bases from the start of the read. This is our starting position
    # Next iterate over a cleaned cigar string to get how many bases in first half on read, add the indel size and then iterate over second cleaned cigar string
    # now we know length of read to traverse. Gen ends pos by adding length to start pos
    cg1_start = 0
    for cg in record_1.cigartuples:
        if cg[0] == 4:
            cg1_start += cg[1]
            break
        elif cg[0] == 5:
            cg1_start += cg[1]
    cg1_len = get_read_length_no_sc(record_1.cigartuples)
    cg1_end = cg1_start + cg1_len
    cg2_start = 0
    for cg in record_2.cigartuples:
        if cg[0] == 4:
            cg2_start += cg[1]
            break
        elif cg[0] == 5:
            cg2_start += cg[1]
    cg2_len = get_read_length_no_sc(record_2.cigartuples)
    cg2_end = cg2_start + cg2_len
    if record_1.is_reverse:
        tmp = get_read_length(record_1.cigartuples) - cg1_end
        cg1_end = get_read_length(record_1.cigartuples) - cg2_start
        cg2_start = tmp
    return [cg1_end, cg2_start]
    

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--sample', required=True)
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv")
parser.add_argument('--folder', default="read_analysis")
parser.add_argument('--bam-suffix', default=".sorted.phased.bam")
parser.add_argument('--bam-folder', default="phased")
parser.add_argument('--merged-reads', default=".all.merged_reads.txt")
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--filters-to-exclude', default="NA")
parser.add_argument('--insert-location-output', required=True)
parser.add_argument('--low-mapq-output', required=True)
parser.add_argument('--regions-list', required=True)
parser.add_argument('--threads', type=int, default=1)
parser.add_argument('--min-mapq', type=int, default=20)
args = parser.parse_args()


threads_to_use = 1
if args.threads > 1:
    threads_to_use = args.threads - 1

regions_split = {}
count = 0
with open(args.regions_list) as in_rl:
    for line in in_rl:
        if count % (threads_to_use/2) not in regions_split:
            regions_split[count % (threads_to_use/2)] = []
        regions_split[count % (threads_to_use/2)].append(line.strip())

# read in the list of merged reads, We want to exclude these from our check of missed split mappings
merged_reads = {}
for sample in args.sample.split(','):
    with open("./"+sample+"/filtered_mapped/"+sample+args.merged_reads) as in_mr:
        for line in in_mr:
            merged_reads[line.strip()] = 1

# Split string into list:
seen_count = 0
seen_reads = {}
low_mapq = {}
for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    # read reference insertions and set up interval trees
    sample_tsv = "./"+sample+"/"+args.folder+"/"+sample+args.suffix
    with open(sample_tsv) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            # get insertion coordinates, insert into intervaltree allowing for the maximum allowable distance
            if row_args[8] != "PASS":
                if "mapq<20" in row_args[8]:
                    if row_args[3] not in low_mapq:
                        low_mapq[row_args[3]] = {}
                    low_mapq[row_args[3]][row_args[1]+":"+row_args[2]] = row_args
            else:
                if row_args[3] not in seen_reads:
                    seen_reads[row_args[3]] = {}
                seen_reads[row_args[3]][row_args[4]+":"+row_args[5]] = row_args

def get_seen_reads(samples, regions, seen_reads):
    records_per_read = defaultdict(list)
    for sample in samples.split(','):
        sam_reader = pysam.AlignmentFile("./"+sample+"/"+args.bam_folder+"/"+sample+args.bam_suffix)
        for region in regions:
            for record in sam_reader.fetch(region=region):
                if record.is_unmapped or record.query_name not in seen_reads:
                    continue
                # Add record to the dict
                records_per_read[record.query_name].append(record)
    return records_per_read

def get_mapq_reads(samples, regions, mapq_reads):
    records_per_read = defaultdict(list)
    for sample in samples.split(','):
        sam_reader = pysam.AlignmentFile("./"+sample+"/"+args.bam_folder+"/"+sample+args.bam_suffix)
        for region in regions:
            for record in sam_reader.fetch(region=region):
                if record.is_unmapped or record.query_name not in mapq_reads:
                    continue
                if record.mapping_quality == 0:
                    records_per_read[record.query_name].append(record)
    return records_per_read

pool = ThreadPool(processes=int(threads_to_use))

results = []
for reg in regions_split:
    async_result = pool.apply_async(get_seen_reads, (args.sample, regions_split[reg], seen_reads))
    results.append(async_result)

results_mapq = []
for reg in regions_split:
    async_result = pool.apply_async(get_mapq_reads, (args.sample, regions_split[reg], low_mapq))
    results_mapq.append(async_result)

#combine results
records_per_read = defaultdict(list)
for i in results:
    res = i.get()
    for read in res:
        records_per_read[read].extend(res[read])

# have all of our reads, now check each one
with open(args.insert_location_output,'w') as out_insert_loc:
    for read in seen_reads:
        if read in merged_reads:
            continue
        if len(records_per_read[read]) > 1:
            # Have more than one alignment for this read. See if there are two that are near each other
            record_1 =  records_per_read[read][0]
            i = 1
            while i < len(records_per_read[read]):
                record_2 = records_per_read[read][i]
                if (record_1.reference_name == record_2.reference_name and 
                    record_1.is_reverse == record_2.is_reverse and
                    (record_1.reference_start - record_2.reference_end < 250 or
                     record_2.reference_start - record_1.reference_end < 250) and
                    (record_1.mapping_quality > args.min_mapq or record_2.mapping_quality > args.min_mapq)):
                    # These records line up enough to possibly support a split mapping insert
                    if record_1.reference_start > record_2.reference_start:
                        insert_pos = get_start_end_string(record_2, record_1)
                    else:
                        insert_pos = get_start_end_string(record_1, record_2)
                    # Have insert position for the read and where it should be. Intersect this with the tsv records I have seen
                    # If there is an insert called that is fully contained by these positions, flag it
                    for pos in seen_reads[read]:
                        start = int(pos.split(':')[0])
                        if insert_pos[0] < start and insert_pos[1] > start:
                            # Flag this insert
                            out_insert_loc.write("\t".join(seen_reads[read][pos])+"\n")
                record_1 = record_2
                i += 1

records_per_read = defaultdict(list)
for i in results_mapq:
    res = i.get()
    for read in res:
        records_per_read[read].extend(res[read])

# have all of our reads, now check each one
with open(args.low_mapq_output,'w') as out_low_mapq:
    for read in low_mapq:
        if len(records_per_read[read]) > 0:
            # Have more than one alignment for this read. See if there are two that are near each other
            i = 0
            while i < len(records_per_read[read]):
                record = records_per_read[read][i]
                for pos in low_mapq[read]:
                    if record.reference_start < int(low_mapq[read][pos][1]) and record.reference_end > int(low_mapq[read][pos][2]):
                        out_low_mapq.write("\t".join(low_mapq[read][pos])+"\n")
                i += 1 