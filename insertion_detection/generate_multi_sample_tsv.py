import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

def get_filters_to_exclude(pf_annotation, filters, repbase_annotation):
    # Base case if pass
    if pf_annotation == "PASS":
        return True
    # pf_annotation may have one or more filters
    # Need to check and elimiate each one found in filters
    # If they are all found in filters, then we can return True.
    # Otherwise there are things in
    if filters == "" or filters == "NA":
        return False
    pf_annotation_split = pf_annotation.split(",")
    count = 0
    for item in pf_annotation_split:
        if item in filters:
            # explicit check to ensure we don't include things that aren't mapped to repbase
            if item == "mapped_fraction":
                if repbase_annotation != "no_repbase_mapping":
                    count += 1
            else:
                count += 1
    if count == len(pf_annotation_split):
        return True
    return False

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


parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in a reference sample')
parser.add_argument('--sample', required=True)
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.tsv")
parser.add_argument('--folder', default="read_analysis")
parser.add_argument('--bam-suffix', default=".sorted.phased.bam")
parser.add_argument('--bam-folder', default="phased")
parser.add_argument('--max-distance', type=int, default=500)
parser.add_argument('--filters-to-exclude', default="NA")
args = parser.parse_args()

# Split string into list:
base_tree = defaultdict(IntervalTree)
seen_count = 0
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
                continue
            chrom = row_args[0]
            start = int(row_args[1]) - args.max_distance
            end = int(row_args[2]) + args.max_distance
            key = chrom + ":" + str(start) + "-" + str(end)
            base_tree[chrom][start:end] = [row_args[9], sample+"_HP:"+str(row_args[19])]

print("Completed Pass 1",file=sys.stderr)
seen_count = 0
intervaltrees = defaultdict(IntervalTree)
insert_sizes = defaultdict(IntervalTree)
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
                continue
            nearby = base_tree[row_args[0]][int(row_args[1]):int(row_args[2])]
            annotation = row_args[9]
            data_dict = defaultdict(int)
            if len(nearby) > 0:
                # have an insert already, update the data
                for item in nearby:
                    if same_family(annotation, item.data[0]):
                        data_dict[item.data[1]] += 1
            chrom = row_args[0]
            start = int(row_args[1]) - args.max_distance
            end = int(row_args[2]) + args.max_distance
            key = chrom + ":" + str(start) + "-" + str(end)
            data_list = []
            for item in data_dict:
                data_list.append(item+"_"+str(data_dict[item]))
            intervaltrees[chrom][start:end] = [annotation, ",".join(data_list)]
            insert_sizes[chrom][int(row_args[1]):int(row_args[2])] = int(row_args[5]) - int(row_args[4])

print("Completed Pass 2",file=sys.stderr)
seen_count = 0
# Third pass to get anything that overlaps the passing inserts in other samples but failed here 
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
                nearby = intervaltrees[row_args[0]][int(row_args[1]):int(row_args[2])]
                if len(nearby) > 0:
                    # have an insert already, update the data
                    for item in nearby:
                        intervaltrees[row_args[0]].remove(item)
                        data_list = []
                        found = False
                        for k in item.data[1].split(","):
                            if sample+"_fail" in k:
                                vals = k.split('_')
                                n = int(vals[len(vals)-1])+1
                                data_list.append(sample+"_fail_"+str(n))
                                found = True
                            else:
                                data_list.append(k)
                        if not found:
                            data_list.append(sample+"_fail_1")
                        intervaltrees[row_args[0]][item.begin:item.end] = [item.data[0], ",".join(data_list)]

print("Completed Pass 3",file=sys.stderr)

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
    for item in intervaltrees[chrom]:
        start = item.begin + args.max_distance
        end = item.end - args.max_distance
        if str(start)+"_"+str(end) not in haplotypes[chrom]:
            haplotypes[chrom][str(start)+"_"+str(end)] = defaultdict(int)
        if str(start)+"_"+str(end) not in softclips[chrom]:
            softclips[chrom][str(start)+"_"+str(end)] = defaultdict(int)
        r_start = int(round(start, -4))
        # Want to round down not up. If we've rounded up round it down
        if r_start > start:
            r_start -= 10000
        r_start -= 10
        r_end = r_start + 10000 + 20
        regions[chrom][str(r_start)+'_'+str(r_end)] = 1

print("Haplotypes Setup",file=sys.stderr)

for sample in args.sample.split(','):
    print(sample,file=sys.stderr)
    header = pysam.AlignmentFile("./"+sample+"/"+args.bam_folder+"/"+sample+args.bam_suffix).header
    sam_reader = pysam.AlignmentFile("./"+sample+"/"+args.bam_folder+"/"+sample+args.bam_suffix)
    for sq in header['SQ']:
        print(sq['SN'],file=sys.stderr)
        if sq['SN'] not in intervaltrees:
            continue
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
                for item in intervaltrees[sq['SN']][r_start:r_end]:
                    start = item.begin + args.max_distance
                    end = item.end - args.max_distance
                    if start >= record.reference_start and start <= record.reference_end:
                        haplotypes[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
                    # Check a slightly larger window for consistent soft clip position
                    if abs(start - record.reference_start) <= args.max_distance or abs(start - record.reference_end) <= args.max_distance:
                        insert_size = 0
                        if len(insert_sizes[sq['SN']][start:end]) > 0:
                            is_item = next(iter(insert_sizes[sq['SN']][start:end]))
                            insert_size = is_item.data
                        # Alignment starts or ends nearby. Check if there is a soft clip > insert size
                        if abs(start - record.reference_start) <= args.max_distance:
                            size = get_soft_clip(record.cigartuples)
                            if size > insert_size:
                                # Add softclip for this inserts
                                softclips[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1
                        elif abs(start - record.reference_end) <= args.max_distance:
                            size = get_soft_clip(record.cigartuples, True)
                            if size > insert_size:
                                # Add softclip for this inserts
                                softclips[sq['SN']][str(start)+"_"+str(end)][sample+"_"+hap] += 1


print("Haplotypes Added",file=sys.stderr)

print("chrom\tstart\tend\tannotation\tsamples\thaplotypes\tsoftclips")
for chrom in intervaltrees:
    for item in intervaltrees[chrom]:
        hap_list = []
        start = item.begin + args.max_distance
        end = item.end - args.max_distance
        for h in haplotypes[chrom][str(start)+"_"+str(end)]:
            hap_list.append(h+"_"+str(haplotypes[chrom][str(start)+"_"+str(end)][h]))
        soft_list = []
        for s in softclips[chrom][str(start)+"_"+str(end)]:
            soft_list.append(s+"_"+str(softclips[chrom][str(start)+"_"+str(end)][s]))
        out_str = chrom+"\t"+str(start)+"\t"+str(end)+"\t"+item.data[0]+"\t"+item.data[1]+"\t"+",".join(hap_list)+"\t"+",".join(soft_list)
        print(out_str)

