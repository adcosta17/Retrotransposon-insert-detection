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


def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation

parser = argparse.ArgumentParser( description='Filter inserts and soft clips based on their mapped fraction of the best repbase hit and if they fall into a high confidence region and are flanked')
parser.add_argument('--vcf', required=True)
parser.add_argument('--read-to-reference-bam', required=True)
parser.add_argument('--telomere-filter', required=True)
parser.add_argument('--centromere-filter', required=True)
parser.add_argument('--min-mapq', type=int, default=20)
parser.add_argument('--window-size', type=int, default=50)
parser.add_argument('--min-mapq-fraction', type=float, default=0)
parser.add_argument('--threads', required=False, type=int, default=1)
parser.add_argument('--sc', type=bool, default=False)
args = parser.parse_args()

threads_to_use = 1
if args.threads > 1:
    threads_to_use = args.threads - 1

telomeres = defaultdict(IntervalTree)
ceontromeres = defaultdict(IntervalTree)

with open(args.telomere_filter) as in_tf:
    count = 0
    for line in in_tf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        telomeres[chrom][start:end] = key

with open(args.centromere_filter) as in_cf:
    count = 0
    for line in in_cf:
        if count == 0:
            count = 1
            continue
        line_args = line.strip().split('\t')
        chrom = line_args[1]
        start = int(line_args[2])
        end = int(line_args[3])
        key = chrom + ":" + str(start) + "-" + str(end)
        ceontromeres[chrom][start:end] = key

# regions is a dict[chrom][region] = [(record, lowqual, total)]
# Each region for a chromsosome is 10kb + window_size long, tsv records where the insert falls in that region are placed in it
# During loop go through each chromosome's region list. Get bam records for region
# Add to the tsv record counts for each bam that intersects +/- window_size bp of insert
vcf_records = {}
regions = {}
pos = 0
# Read in tsv first
# get a list of reads we care about
vcf_in = pysam.VariantFile(args.vcf)
header_lines = str(vcf_in.header).strip().split('\n')
for line in header_lines:
    if "##" in line:
        print(line)
    else:
        print('##INFO=<ID=MAPQFRAC,Number=1,Type=String,Description="Mapping Quality Fraction">')
        print('##FILTER=<ID=in_centromere,Number=1,Type=String,Description="In centromere_filter">')
        print('##FILTER=<ID=in_telomere,Number=1,Type=String,Description="In Telomere filter">')
        print('##FILTER=<ID=mapq_fraction,Number=1,Type=String,Description="mapq_fraction filter">')
        print(line)
for rec in vcf_in.fetch():
    rec_row = str(rec).strip().split('\t')
    start = int(round(int(rec_row[1]), -4))
    # Want to round down not up. If we've rounded up round it down
    if start > int(rec_row[1]):
        start -= 10000
    start -= int(args.window_size)
    end = start + 10000 + 2*int(args.window_size)
    if len(ceontromeres[rec_row[0]][int(rec_row[1]):int(rec_row[1])+1]) > 0:
        rec_row[6] = "in_centromere"
    if len(telomeres[rec_row[0]][int(rec_row[1]):int(rec_row[1])+1]) > 0:
        rec_row[6] = "in_telomere"
    if rec_row[0] not in regions:
        regions[rec_row[0]] = defaultdict(list)
    regions[rec_row[0]][str(start)+"-"+str(end)].append([rec_row,0,0])

for chrom in regions:
    for region in regions[chrom]:
        if pos % int(threads_to_use) not in vcf_records:
            vcf_records[pos % int(threads_to_use)] = {}
        if chrom not in vcf_records[pos % int(threads_to_use)]:
            vcf_records[pos % int(threads_to_use)][chrom] = defaultdict(list)
        vcf_records[pos % int(threads_to_use)][chrom][region] = regions[chrom][region]
        pos += 1

contig_lengths = {}
header = pysam.AlignmentFile(args.read_to_reference_bam).header
for sq in header['SQ']:
    contig_lengths[sq['SN']] = int(sq['LN'])


# read in alignments for those reads only
def filter_hc(file_name, vcf_records, contig_lengths):
    sam_reader = pysam.AlignmentFile(file_name)
    hc_inserts_and_softclips = vcf_records
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
                        if ((int(record.reference_start) <= int(line_arr[1]) and int(record.reference_end) >= int(line_arr[1])) or
                           (int(record.reference_start) <= int(line_arr[1]) and int(record.reference_end) < int(line_arr[1]) and 
                            (int(line_arr[1]) - int(record.reference_end)) <= int(args.window_size)) or
                           (int(record.reference_start) <= int(line_arr[1])+1 and int(record.reference_end) >= int(line_arr[1])+1) or
                           (int(record.reference_start) > int(line_arr[1])+1 and int(record.reference_end) >= int(line_arr[1])+1 and 
                            (int(record.reference_start) - int(line_arr[1])+1) <= int(args.window_size))):
                            hc_inserts_and_softclips[chrom][region][i][2] += 1
                            if record.mapq < args.min_mapq:
                                hc_inserts_and_softclips[chrom][region][i][1] += 1
            except ValueError:
                pass
    return hc_inserts_and_softclips

pool = ThreadPool(processes=int(threads_to_use))

results = []
for reg in vcf_records:
    async_result = pool.apply_async(filter_hc, (args.read_to_reference_bam, vcf_records[reg], contig_lengths))
    results.append(async_result)


#combine results
for i in results:
    res = i.get()
    for chrom in res:
        for region in res[chrom]:
            for i in range(len(res[chrom][region])):
                if res[chrom][region][i][2] == 0:
                    updated_line = res[chrom][region][i][0]
                    updated_line[6] = "mapq_fraction"
                    updated_line[7] = updated_line[7] + ";MAPQFRAC=1"
                    print('\t'.join(updated_line))
                elif float(res[chrom][region][i][1])/res[chrom][region][i][2] < args.min_mapq_fraction:
                    updated_line = res[chrom][region][i][0]
                    updated_line[7] = updated_line[7] + ";MAPQFRAC=" + str(float(res[chrom][region][i][1])/res[chrom][region][i][2])
                    print('\t'.join(updated_line))
                else:
                    updated_line = res[chrom][region][i][0]
                    updated_line[6] = "mapq_fraction"
                    updated_line[7] = updated_line[7] + ";MAPQFRAC=" + str(float(res[chrom][region][i][1])/res[chrom][region][i][2])
                    print('\t'.join(updated_line))
