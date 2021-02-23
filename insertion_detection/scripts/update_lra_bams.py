import pysam
import argparse
import sys
import re
from collections import defaultdict

def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0

parser = argparse.ArgumentParser( description='Update BAMs generated with LRA to have secondary mapping flags if more than 1 alignment')
parser.add_argument('--bam', type=str, required=True)
parser.add_argument('--output-bam', type=str, required=True)
args = parser.parse_args()

# Store the positions, qulaity and lengths of all read alignments
# Sort all of these to get the highest quality and then longest alignments first
# Select the first one as primary and flag the rest as secondary alignments
read_alignments = defaultdict(list)
header = pysam.AlignmentFile(args.bam).header
for sq in header['SQ']:
    sam_reader = pysam.AlignmentFile(args.bam)
    tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
    print(sq['SN'])
    for record in tmp_sam_reader:
        read_alignments[record.query_name].append((record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end), record.mapping_quality, record.query_alignment_length))

print("Got all Positions")

# For each read sort the list based on the quality and then the length
for read in read_alignments:
    records = sorted(read_alignments[read], key = lambda x: (-x[1], -x[2]))
    read_alignments[read] = records

print("Sorted Positions")

with pysam.AlignmentFile(args.output_bam, "wb", header=header) as out_bam:
    for sq in header['SQ']:
        sam_reader = pysam.AlignmentFile(args.bam)
        tmp_sam_reader = sam_reader.fetch(contig=sq['SN'])
        print(sq['SN'])
        for record in tmp_sam_reader:
            if record.reference_name+":"+str(record.reference_start)+"-"+str(record.reference_end) != read_alignments[record.query_name][0][0] and len(read_alignments[record.query_name][0]) > 1:
                # Update the flag
                a = pysam.AlignedSegment(header=header)
                a.query_name = record.query_name
                a.flag = record.flag + 256
                a.reference_id = get_id(header, record.reference_name)
                a.reference_start = record.reference_start
                a.mapping_quality = record.mapping_quality
                a.cigarstring = record.cigarstring
                a.query_sequence = record.query_sequence
                out_bam.write(a)
            else:
                out_bam.write(record)