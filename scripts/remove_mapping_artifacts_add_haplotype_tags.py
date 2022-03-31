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
parser.add_argument('--tsv', required=True)
parser.add_argument('--full-bam', required=True)
parser.add_argument('--insert-bam', required=True)
parser.add_argument('--min-mapq', type=int, default=20)
parser.add_argument('--insert-filter', required=True)
args = parser.parse_args()

# Get all the inserts that will have mapped inserts
# Exclude those that have mapq below the required 20 and min_insert length

insert_filter = {}
with open(args.insert_filter) as in_filter:
    for row in in_filter:
        row_args = row.strip().split("\t")
        if row_args[0] not in insert_filter:
            insert_filter[row_args[0]] = {}
        insert_filter[row_args[0]][row_args[1]+":"+row_args[2]] = row_args

reads_to_use = defaultdict(int)
# Read in tsv first
# get a list of reads we care about
with open(args.tsv) as in_tsv:
    count = 0
    for line in in_tsv:
        if count == 0:
            count = 1
            continue
        line = line.strip().split('\t')
        if "mapq<20" in line[8] or "min_insertion_length" in line[8]:
            # Don't need to check these inserts. Add them to the output as is
            continue
        # Otherwise add the read to the list of reads we care about. Will need to parse through the full BAM
        # see where all the alignments of this read that are mapq20 or greater are and then compare them to the insert mapping locations
        reads_to_use[line[3]] = 1


# Read through the Full bam. Look for reads in the reads to use set
read_locations = {}
read_haplotypes = {}
sam_reader = pysam.AlignmentFile(args.full_bam, check_sq=False)
for record in sam_reader.fetch():
    if record.is_unmapped:
        continue
    tags = record.get_tags()
    for t in tags:
        if t[0] == "HP":
            if record.query_name not in read_haplotypes:
                read_haplotypes[record.query_name] = defaultdict(int)
            read_haplotypes[record.query_name][str(record.reference_start)+"-"+str(record.reference_end)] = int(t[1])
    if record.query_name not in reads_to_use:
        continue
    if record.mapping_quality >= args.min_mapq:
        # Keep this mapping
        if record.query_name not in read_locations:
            read_locations[record.query_name] = []
        read_locations[record.query_name].append(record)

# Once we have the locations of where each read we care about aligns in the full bam, parse the insert bam
insert_locations = {}
sam_reader = pysam.AlignmentFile(args.insert_bam, check_sq=False)
for record in sam_reader.fetch():
    if record.is_unmapped or record.mapping_quality < args.min_mapq:
        # Insert has not mapped well, Don't consider it
        continue
    if record.mapping_quality >= args.min_mapq:
        # Keep this mapping
        rec_name = record.query_name.split(":")[0]
        insert_pos = record.query_name.split(":")[1]
        if rec_name not in insert_locations:
            insert_locations[rec_name] = {}
        insert_locations[rec_name][insert_pos] = record

# Now go through the tsv again. Change the annotation as needed
with open(args.tsv) as in_tsv:
    count = 0
    for line in in_tsv:
        if count == 0:
            count = 1
            print(line.strip()+"\thaplotype")
            continue
        line_args = line.strip().split('\t')
        # Add read haplotype tag here
        if line_args[3] in read_haplotypes:
            found = False
            for pos in read_haplotypes[line_args[3]]:
                if int(pos.split("-")[0]) <= int(line_args[1]) and int(pos.split("-")[1]) >= int(line_args[2]):
                    line_args.append(str(read_haplotypes[line_args[3]][pos]))
                    found = True
            if not found:
                line_args.append("0")
        else:
            line_args.append("0")
        if line_args[0] in insert_filter:
            # Check to see if insert is in the insert filter set
            for item in insert_filter[line_args[0]]:
                if line_args[3] == insert_filter[line_args[0]][item][3] and abs(int(line_args[1]) - int(insert_filter[line_args[0]][item][1])) <= 20 and abs(int(line_args[2]) - int(insert_filter[line_args[0]][item][2])) <= 20:
                    # Failure
                    line_args[8] = update_annotation(line_args[8], "insert_filter")
                    break
        if "mapq<20" in line_args[8] or "min_insertion_length" in line_args[8]:
            # Don't need to check these inserts. Add them to the output as is
            print("\t".join(line_args))
            continue
        chimeric_fail = False
        if line_args[3] in insert_locations and line_args[3] in read_locations:
            # have an insert for this read that mapped well enough
            if line_args[4]+"-"+line_args[5] in insert_locations[line_args[3]] and len(read_locations[line_args[3]]) > 1:
                # Have a mapping for this insert, compare it to all of the reads mappings. If any intersect Flag this insert as a failure
                for record in read_locations[line_args[3]]:
                    if (record.reference_name == insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_name and
                        record.reference_start == insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_start and
                        record.reference_end == insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_end):
                        continue
                    if (record.reference_name == insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_name and
                        record.reference_start < insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_start and
                        record.reference_end > insert_locations[line_args[3]][line_args[4]+"-"+line_args[5]].reference_end):
                        chimeric_fail = True
                        break
        if chimeric_fail:
            line_args[8] = update_annotation(line_args[8], "possible_chimera")
        print("\t".join(line_args))


