import argparse
import pysam
import sys
import re
import gzip
import pybedtools

from collections import defaultdict
from os import listdir
from os.path import isfile, join
from multiprocessing.pool import ThreadPool


# Gets the total read length based on the cigar string
def get_read_length(cigarstring):
    count = 0
    for cg in re.findall('[0-9]*[A-Z]', cigarstring):
        if cg.endswith('M'):
            count += int(cg[:cg.find("M")])
        elif cg.endswith('X'):
            count += int(cg[:cg.find("X")])
        elif cg.endswith('I'):
            count += int(cg[:cg.find("I")])
        elif cg.endswith('P'):
            count += int(cg[:cg.find("P")])
        elif cg.endswith('S'):
            count += int(cg[:cg.find("S")])
        elif cg.endswith('H'):
            count += int(cg[:cg.find("H")])
    return count

def full_read_ref_position(record):
    step = 0
    cg = re.findall('[0-9]*[A-Z]', record.cigarstring)
    if cg[0].endswith('H') and cg[1].endswith('S'):
        step = int(cg[0][:cg[0].find("H")])
        step += int(cg[1][:cg[1].find("S")])
    elif cg[0].endswith('H') or cg[0].endswith('S'):
        if cg[0].endswith('H'):
            step = int(cg[0][:cg[0].find("H")])
        else:
            step = int(cg[0][:cg[0].find("S")])
    ref_start = record.reference_start - step
    cg = cg[::-1]
    if cg[0].endswith('H') and cg[1].endswith('S'):
        step = int(cg[0][:cg[0].find("H")])
        step += int(cg[1][:cg[1].find("S")])
    elif cg[0].endswith('H') or cg[0].endswith('S'):
        if cg[0].endswith('H'):
            step = int(cg[0][:cg[0].find("H")])
        else:
            step = int(cg[0][:cg[0].find("S")])
    ref_end = record.reference_end + step
    return (ref_start, ref_end)


def get_read_pos(record, region):
    reg = region.split(":")[1].split("-")
    region_start = int(reg[0])
    region_end = int(reg[1])
    ref_pos = full_read_ref_position(record)
    # See if region start or region end are outside the read's start or end
    region_starts_before = False
    region_ends_after = False
    start_pos = 0
    end_pos = 0
    read_len = get_read_length(record.cigarstring)
    if region_start < ref_pos[0]:
        # Know this read clips the region at the start
        region_starts_before= True
        start_pos = 0
    elif ref_pos[1] < region_end:
        region_ends_after = True
        end_pos = read_len
    ref_count = record.reference_start
    read_count = 0
    for cg in re.findall('[0-9]*[A-Z]', record.cigarstring):
        step = 0
        if cg.endswith('M'):
            step = int(cg[:cg.find("M")])
            ref_count += step
            read_count += step
        elif cg.endswith('I'):
            step = int(cg[:cg.find("I")])
            read_count += step
        elif cg.endswith('D'):
            step = int(cg[:cg.find("D")])
            ref_count += step
        elif cg.endswith('N'):
            step = int(cg[:cg.find("N")])
            ref_count += step
        elif cg.endswith('S'):
            step = int(cg[:cg.find("S")])
            read_count += step
        elif cg.endswith('H'):
            step = int(cg[:cg.find("H")])
            read_count += step
        elif cg.endswith('='):
            step = int(cg[:cg.find("=")])
            ref_count += step
            read_count += step
        elif cg.endswith('X'):
            step = int(cg[:cg.find("X")])
            ref_count += step
            read_count += step
        # At each iteration need to check if we have reached the start of the region
        if ref_count <= region_start and not region_starts_before:
            start_pos = read_count
        if ref_count >= region_end and not region_ends_after:
            end_pos = read_count
            break
    if start_pos == end_pos:
        # if start and end are the same it means the region is deleted on the read
        # The read doesn't actually span the region
        return [0,0,0]
    return [start_pos, end_pos, read_len]

def get_min_flank(span_list):
    vals = []
    for item in span_list:
        vals.append(item[1])
        vals.append(item[2])
    return min(vals)

def get_max_div(span_list):
    vals = []
    for item in span_list:
        vals.append(item[4])
    return max(vals)

parser = argparse.ArgumentParser( description='Get all reads that flank any deleted regions based on alignments to original genome')
parser.add_argument('--bed', required=True)
parser.add_argument('--read-to-reference-bam', required=True)
parser.add_argument('--output-insert-list', required=True)
parser.add_argument('--output-sc-list', required=False)
parser.add_argument('--flank-size', required=False, type=int, default=100)
parser.add_argument('--threads', required=False, type=int, default=1)
args = parser.parse_args()

region_bed = pybedtools.BedTool(args.bed)
read_bam = pybedtools.BedTool(args.read_to_reference_bam)
intersection_res = region_bed.intersect(read_bam,bed=True,loj=True)

regions_of_interest = {}
pos = 0
for line in intersection_res:
    lines = str(line).split("\n")
    for l in lines:
        line_args = l.strip().split("\t")
        if len(line_args) > 7:
            if int(line_args[7]) >= 20:
                if (pos % int(args.threads)) not in regions_of_interest:
                    regions_of_interest[pos % int(args.threads)] = defaultdict(list)
                regions_of_interest[pos % int(args.threads)][line_args[0]+":"+str(line_args[1])+"-"+str(line_args[2])].append(line_args[6])
                pos += 1
mapped_count = {}

print("Computed Mapped Counts")
# Define the function to process what we want
def get_regions(file_name, flank_size, region_list, mapped_count):
    sam_reader = pysam.AlignmentFile(file_name)
    reads_that_span = {}
    for region in region_list:
        try:
            tmp_sam_reader = sam_reader.fetch(region=region)
            for record in tmp_sam_reader:
                # Check to see if the record intersects the region at all
                # Specifically if any hard or softclips at start or end of read have more than 100 bp in region
                if record.query_name in region_list[region]:
                    # Know read has a good alignment to this region
                    # Get the read position of where the region is based on the cigar string
                    pos_arr = get_read_pos(record, region)
                    if pos_arr != [0,0,0]:
                        if record.query_name not in reads_that_span:
                            reads_that_span[record.query_name] = {}
                        reads_that_span[record.query_name][region] = pos_arr
        except ValueError:
            pass
    return reads_that_span

pool = ThreadPool(processes=int(args.threads))

results = []
for reg in regions_of_interest:
    async_result = pool.apply_async(get_regions, (args.read_to_reference_bam, int(args.flank_size), regions_of_interest[reg], mapped_count))
    results.append(async_result)


#combine results
reads_that_span = {}
for i in results:
    res = i.get()
    for read in res:
        if read not in reads_that_span:
            reads_that_span[read] = {}
        for region in res[read]:
            reads_that_span[read][region] = res[read][region]

print("Got Reads that span")

sc_reads = []
with open(args.output_insert_list, 'w') as out_txt:
    for read in reads_that_span:
        for region in reads_that_span[read]:
            if (reads_that_span[read][region][0] >= args.flank_size and 
            	(reads_that_span[read][region][2] - reads_that_span[read][region][1]) >= args.flank_size):
                out_txt.write(read + "\t" + region + "\t"+ str(reads_that_span[read][region][0])+"\t"+ str(reads_that_span[read][region][1]) + "\n")
            elif (args.output_sc_list and 
            	 (reads_that_span[read][region][0] >= args.flank_size or 
            	 	(reads_that_span[read][region][2] - reads_that_span[read][region][1]) >= args.flank_size)):
                sc_reads.append(read + "\t" + region + "\t"+ str(reads_that_span[read][region][0])+"\t"+ str(reads_that_span[read][region][1]) + "\n")

if args.output_sc_list:
    with open(args.output_sc_list, 'w') as out_txt:
        for line in sc_reads:
            out_txt.write(line)
