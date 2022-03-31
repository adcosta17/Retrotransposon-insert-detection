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

parser = argparse.ArgumentParser( description='Remove mapping_artifacts and add haplotype tags to insert records')
parser.add_argument('--vcf', required=True)
parser.add_argument('--full-bam', required=True)
parser.add_argument('--insert-bam', required=True)
parser.add_argument('--min-mapq', type=int, default=20)
args = parser.parse_args()

read_set = defaultdict(int)

vcf_in = pysam.VariantFile(args.vcf)
for rec in vcf_in.fetch():
    rec_row = str(rec).strip().split('\t')
    for read in rec.info["RNAMES"]:
        read_set[read] = 1

# Read through the Full bam. Look for reads in the reads to use set
read_locations = {}
read_haplotypes = {}
sam_reader = pysam.AlignmentFile(args.full_bam)
for record in sam_reader.fetch():
    if record.is_unmapped or record.query_name not in read_set:
        continue
    if record.mapping_quality >= args.min_mapq:
        # Keep this mapping
        if record.query_name not in read_locations:
            read_locations[record.query_name] = []
        read_locations[record.query_name].append(record)

# Once we have the locations of where each read we care about aligns in the full bam, parse the insert bam
insert_locations = {}
sam_reader = pysam.AlignmentFile(args.insert_bam)
for record in sam_reader.fetch():
    if record.is_unmapped or record.mapping_quality < args.min_mapq:
        # Insert has not mapped well, Don't consider it
        continue
    if record.mapping_quality >= args.min_mapq:
        # Keep this mapping
        rec_name = record.query_name
        insert_locations[rec_name] = record

# Now go through the tsv again. Change the annotation as needed
vcf_in = pysam.VariantFile(args.vcf)
header_lines = str(vcf_in.header).strip().split('\n')
for line in header_lines:
    if "##" in line:
        print(line)
    else:
        print('##FILTER=<ID=possible_chimera,Number=1,Type=String,Description="possible_chimera filter">')
        print(line)
for rec in vcf_in.fetch():
    rec_row = str(rec).strip().split('\t')
    reads_to_use = []
    for read in rec.info["RNAMES"]:
        chimeric_fail = False
        if rec.id in insert_locations and read in read_locations:
            # have an insert for this read that mapped well enough
            if len(read_locations[read]) > 1:
                # Have a mapping for this insert, compare it to all of the reads mappings. If any intersect Flag this insert as a failure
                for record in read_locations[read]:
                    if (record.reference_name == rec.chrom and record.reference_start < rec.pos and record.reference_end > rec.pos +1):
                        # Have a read with multiple mappings that supports this insert, remove it
                        chimeric_fail = True
                        break
        if chimeric_fail:
            continue
        reads_to_use.append(read)
    if len(reads_to_use) == 0:
        rec_row[6] = "possible_chimera"
    print('\t'.join(rec_row))


