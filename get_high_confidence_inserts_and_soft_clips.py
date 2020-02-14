import argparse
import pysam
import sys
import re
import gzip
from collections import defaultdict
from os import listdir
from os.path import isfile, join
from intervaltree import Interval, IntervalTree
from multiprocessing.pool import ThreadPool

def is_hc(intervaltrees, record):
    chrom = record[0]
    start = int(record[1])
    end = int(record[2])
    nearby = intervaltrees[chrom][start:end]
    if len(nearby) == 0:
        return False
    return True

parser = argparse.ArgumentParser( description='Filter inserts and soft clips based on their mapped fraction of the best repbase hit and if they fall into a high confidence region and are flanked')
parser.add_argument('--tsv', required=True)
parser.add_argument('--read-to-reference-bam', required=True)
parser.add_argument('--min-mapq', type=int, default=20)
parser.add_argument('--window-size', type=int, default=50)
parser.add_argument('--min-mapq-fraction', type=float, default=0)
parser.add_argument('--dust-fraction', type=float, default=0.5)
parser.add_argument('--threads', required=False, type=int, default=1)
parser.add_argument('--sc', type=bool, default=False)
args = parser.parse_args()

threads_to_use = 1
if args.threads > 1:
    threads_to_use = args.threads - 1

mapped_count = {}
annotated_count = {}
sam_reader = pysam.AlignmentFile(args.read_to_reference_bam)
for record in sam_reader.fetch():
    if record.is_unmapped:
        continue
    if record.query_name not in mapped_count:
        mapped_count[record.query_name] = 0
    if record.mapping_quality >= 20:
        mapped_count[record.query_name] += 1

# regions is a dict[chrom][region] = [(record, lowqual, total)]
# Each region for a chromsosome is 10kb long, tsv records where the insert falls in that region are placed in it
# During loop go through each chromosome's region list. Get bam records for region
# Add to the tsv record counts for each bam that intersects +/- 10 bp of insert
tsv_records = {}
regions = {}
pos = 0
# Read in tsv first
# get a list of reads we care about
with open(args.tsv) as in_tsv:
    count = 0
    for line in in_tsv:
        if count == 0:
            count = 1
            print(line.strip()+"\t"+"low_confidence_fraction")
            continue
        line = line.strip().split('\t')
        # Get bin start and end position, based on start position
        start = int(round(int(line[1]), -4))
        # Want to round down not up. If we've rounded up round it down
        if start > int(line[1]):
            start -= 10000
        start -= int(args.window_size)
        end = start + 10000 + 2*int(args.window_size)
        if args.sc:
            if mapped_count[line[3]] > 1:
                if line[8] == "PASS":
                    line[8] = "multimapped"
                else:
                    line[8] = line[8]+"multimapped"
        if line[0] not in regions:
            regions[line[0]] = defaultdict(list)
        regions[line[0]][str(start)+"-"+str(end)].append([line,0,0])

for chrom in regions:
    for region in regions[chrom]:
        if pos % int(threads_to_use) not in tsv_records:
            tsv_records[pos % int(threads_to_use)] = {}
        if chrom not in tsv_records[pos % int(threads_to_use)]:
            tsv_records[pos % int(threads_to_use)][chrom] = defaultdict(list)
        tsv_records[pos % int(threads_to_use)][chrom][region] = regions[chrom][region]
        pos += 1

contig_lengths = {}
header = pysam.AlignmentFile(args.read_to_reference_bam).header
for sq in header['SQ']:
    contig_lengths[sq['SN']] = int(sq['LN'])


# read in alignments for those reads only
def filter_hc(file_name, tsv_records, contig_lengths):
    sam_reader = pysam.AlignmentFile(file_name)
    hc_inserts_and_softclips = tsv_records
    for chrom in hc_inserts_and_softclips:
        for region in hc_inserts_and_softclips[chrom]:
            # Get the reads that align to this region
            try:
                start = region.split("-")[0]
                end = region.split("-")[1]
                if int(start) < 0:
                    start = "0"
                if int(end) > contig_lengths[chrom]:
                	end = str(contig_lengths[chrom])
                tmp_sam_reader = sam_reader.fetch(region=chrom+":"+start+"-"+end)
                for record in tmp_sam_reader:
                    # Go through all the reads in the region
                    for i in range(len(hc_inserts_and_softclips[chrom][region])):
                        line_arr = hc_inserts_and_softclips[chrom][region][i][0]
                        #if record.query_name == line_arr[3]:
                        #	print(str(record.reference_start) + " " + str(line_arr[1]) + " " + str(record.reference_end) + " " + str(line_arr[2]))
                        if ((int(record.reference_start) <= int(line_arr[1]) and int(record.reference_end) >= int(line_arr[1])) or
                           (int(record.reference_start) <= int(line_arr[1]) and int(record.reference_end) < int(line_arr[1]) and (int(line_arr[1]) - int(record.reference_end)) <= int(args.window_size)) or
                           (int(record.reference_start) <= int(line_arr[2]) and int(record.reference_end) >= int(line_arr[2])) or
                           (int(record.reference_start) > int(line_arr[2]) and int(record.reference_end) >= int(line_arr[2]) and (int(record.reference_start) - int(line_arr[2])) <= int(args.window_size))):
                            hc_inserts_and_softclips[chrom][region][i][2] += 1
                            if record.mapq < args.min_mapq:
                                hc_inserts_and_softclips[chrom][region][i][1] += 1
            except ValueError:
                pass
    return hc_inserts_and_softclips

pool = ThreadPool(processes=int(threads_to_use))

results = []
for reg in tsv_records:
    async_result = pool.apply_async(filter_hc, (args.read_to_reference_bam, tsv_records[reg], contig_lengths))
    results.append(async_result)


#combine results
for i in results:
    res = i.get()
    for chrom in res:
        for region in res[chrom]:
            for i in range(len(res[chrom][region])):
                dust = False
                if float(res[chrom][region][i][0][6]) < args.dust_fraction:
                    dust = True
                if res[chrom][region][i][2] == 0:
                    updated_line = res[chrom][region][i][0]
                    if updated_line[8] == "PASS":
                        if dust:
                            updated_line[8] = "mapq_fraction,dust"
                        else:
                            updated_line[8] = "mapq_fraction"
                    else:
                        if dust:
                            updated_line[8] = updated_line[8]+",mapq_fraction,dust"
                        else:
                            updated_line[8] = updated_line[8]+",mapq_fraction"
                    print("\t".join(updated_line) + "\t1")
                elif float(res[chrom][region][i][1])/res[chrom][region][i][2] < args.min_mapq_fraction:
                    if dust:
                        updated_line = res[chrom][region][i][0]
                        if updated_line[8] == "PASS":
                            updated_line[8] = "dust"
                        else:
                            updated_line[8] = updated_line[8]+",dust"
                        print("\t".join(updated_line)+ "\t"+str(float(res[chrom][region][i][1])/res[chrom][region][i][2]))
                    else:
                        print("\t".join(res[chrom][region][i][0])+ "\t"+str(float(res[chrom][region][i][1])/res[chrom][region][i][2]))
                else:
                    updated_line = res[chrom][region][i][0]
                    if updated_line[8] == "PASS":
                        if dust:
                            updated_line[8] = "mapq_fraction,dust"
                        else:
                            updated_line[8] = "mapq_fraction"
                    else:
                        if dust:
                            updated_line[8] = updated_line[8]+",mapq_fraction,dust"
                        else:
                            updated_line[8] = updated_line[8]+",mapq_fraction"
                    print("\t".join(updated_line) +"\t"+str(float(res[chrom][region][i][1])/res[chrom][region][i][2]))
