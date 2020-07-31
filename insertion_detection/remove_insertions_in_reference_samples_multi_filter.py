import pysam
import argparse
import sys
import csv
from intervaltree import Interval, IntervalTree
from collections import defaultdict

def get_inserts(all_inserts, chrom, start, end):
    mcount = 0
    for i in range(end-start):
        mcount = max(mcount, all_inserts[chrom][start+i])
    return mcount

def get_string(items):
    item_dict = defaultdict(int)
    for item in items:
        item_dict["in_sample_"+item.replace("_fail",'')] = 1
    new_items = []
    for item in item_dict:
        new_items.append(item)
    return ",".join(new_items)

def get_full_string(items, include_fail=False):
    item_dict = defaultdict(int)
    for item in items.split(','):
        if not include_fail and "fail" in item:
            continue
        item_dict["in_sample_"+item.replace("_fail",'')] += 1
    new_items = []
    for item in item_dict:
        new_items.append(item+"__"+str(item_dict[item]))
    return ",".join(new_items)

def get_sample_dict(items):
    new_items = defaultdict(int)
    for item in items.split(","):
        new_items[item] += 1
    return new_items

def get_sample_dict_no_fail(items):
    new_items = defaultdict(int)
    for item in items.split(","):
        new_item = item.replace("_fail",'')
        new_items[new_item] += 1
    return new_items

def is_only_sample(sample_dict, sample):
    count = 0
    for item in sample_dict:
        count += sample_dict[item]
    if count == sample_dict[sample]:
        return True
    return False

def in_multiple_samples(sample_dict, sample_count, sample):
    count = 0
    for item in sample_dict:
        # If there is a single insert in another sample we can allow it
        # Parse the item to remove the failure if there is one
        new_item = item.replace("_fail",'')
        if new_item == sample:
            continue
        if "fail" in item and new_item in sample_dict:
            if (sample_dict[new_item] + sample_dict[item]) >= 1:
                count += 1
        elif "fail" in item:
            if sample_dict[item] >= 1:
                count += 1
        else:
            if sample_dict[new_item] >= 1:
                count += 1
    if count >= sample_count:
        return True
    return False


def in_multiple_samples_no_fail(sample_dict, sample_count, sample, max_common):
    count = 0
    for item in sample_dict:
        # If there is a single insert in another sample we can allow it
        if item == sample or "fail" in item:
            continue
        if sample_dict[item] >= max_common:
            count += 1
    if count >= sample_count:
        return True
    return False


def get_insert_dict(insert_string):
    # Parse through the insert string
    insert_dict = {}
    for insert in insert_string.split(','):
        if "_fail_" in insert:
            sample = insert.split("_fail_")[0]
            count = int(insert.split("_fail_")[1])
            if sample not in insert_dict:
                insert_dict[sample] = defaultdict(int)
            insert_dict[sample]["fail"] = count
        elif "HP:" in insert:
            sample = insert.split("HP:")[0][:-1]
            hp = insert.split("HP:")[1].split('_')[0]
            count = int(insert.split("HP:")[1].split('_')[1])
            if sample not in insert_dict:
                insert_dict[sample] = defaultdict(int)
            insert_dict[sample][hp] = count
    return insert_dict

def get_haplotype_dict(haplotype_string):
    haplotype_dict = {}
    for insert in haplotype_string.split(','):
        if "HP:" in insert:
            sample = insert.split("HP:")[0][:-1]
            hp = insert.split("HP:")[1].split('_')[0]
            count = int(insert.split("HP:")[1].split('_')[1])
            if sample not in haplotype_dict:
                haplotype_dict[sample] = defaultdict(int)
            haplotype_dict[sample][hp] = count
    return haplotype_dict

def update_annotation(annotation, update):
    if annotation == "PASS":
        annotation = update
    else:
        annotation = annotation+","+update
    return annotation

def format_and_print_line_insert_dict(insert_dict, row_args, reason):
    ret = []
    for sample in insert_dict:
        ret.append("in_sample_"+sample)
    row_args[8] = update_annotation(row_args[8], ','.join(ret))
    row_args.append("insert\t"+reason+"_"+str(len(ret)))
    print("\t".join(row_args))

def get_filters(file):
    filter_dict = {}
    int_tree = defaultdict(IntervalTree)
    with open(file) as csvfile:
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            # Get the count per sample
            if row_args[0] not in filter_dict:
                filter_dict[row_args[0]] = {}
            if row_args[1]+"-"+row_args[2] in filter_dict[row_args[0]]:
                # Duplicate entry with different annotation, keep first one we see, ignore rest
                # Insert count should be identical between duplicates
                continue
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]] = {}
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["annotation"] = row_args[3]
            # parse insert list to get set of inserts per haplotype plus failures
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["inserts"] = get_insert_dict(row_args[4])
            # parse all haplotype counts to get total breakdown of haplotype stats per sample
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["haplotypes"] = get_haplotype_dict(row_args[5])
            if len(row_args) == 6:
                row_args.append('')
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["softclips"] = get_haplotype_dict(row_args[6])
            start = int(row_args[1]) - args.max_distance
            end = int(row_args[2]) + args.max_distance
            key = row_args[0] + ":" + str(start) + "-" + str(end)
            # setup string
            ret = []
            for sample in filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["inserts"]:
                ret.append("in_sample_"+sample)
            int_tree[row_args[0]][start:end] = ','.join(ret)
    return (filter_dict, int_tree)

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in any other sample')
parser.add_argument('--input', required=True)
#parser.add_argument('--sc', required=True)
parser.add_argument('--multi-reference-filter-xlarge-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-large-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-medium-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-small-tsv', type=str, required=True)
parser.add_argument('--telomere-filter', required=True)
parser.add_argument('--centromere-filter', required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument('--ref-sample', type=str, default="DDTS_0001_nn_n")
parser.add_argument('--reference-repeat-bed', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=100)
parser.add_argument('--max-common', type=int, default=1)
parser.add_argument('--sample-count', type=int, default=1)
parser.add_argument('--num-samples-needed-to-remove', type=int, default=1)
parser.add_argument('--insert-filter', required=True)
parser.add_argument('--low-mapq-filter', required=True)
args = parser.parse_args()

insert_filter = {}
with open(args.insert_filter) as in_filter:
    for row in in_filter:
        row_args = row.strip().split("\t")
        if row_args[0] not in insert_filter:
            insert_filter[row_args[0]] = {}
        insert_filter[row_args[0]][row_args[1]+":"+row_args[2]] = row_args

mapq_0_filter = {}
with open(args.low_mapq_filter) as in_filter:
    for row in in_filter:
        row_args = row.strip().split("\t")
        if row_args[0] not in mapq_0_filter:
            mapq_0_filter[row_args[0]] = {}
        mapq_0_filter[row_args[0]][row_args[1]+":"+row_args[2]] = row_args

# Read each multi-ref-filter-tsv and set up a dict for each
# Use exact positions, but store the number of inserts found nearby for each haplotype + fails
x_large_filter, intervaltree_xl = get_filters(args.multi_reference_filter_xlarge_tsv)
large_filter, intervaltree_large = get_filters(args.multi_reference_filter_large_tsv)
medium_filter, intervaltree_medium = get_filters(args.multi_reference_filter_medium_tsv)
small_filter, intervaltree_small = get_filters(args.multi_reference_filter_small_tsv)

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

# read input file and only emit records not near a record in the interval tree 
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            print(row.strip()+"\tsource\tsample_count")
            continue
        if float(row_args[18]) > 0.5:
            row_args[8] = update_annotation(row_args[8], "mapq_fraction")
        if len(ceontromeres[row_args[0]][int(row_args[1]):int(row_args[2])]) > 0 and "in_centromere" not in row_args[8]:
            row_args[8] = update_annotation(row_args[8], "in_centromere")
        if len(telomeres[row_args[0]][int(row_args[1]):int(row_args[2])]) > 0 and "in_telomere" not in row_args[8]:
            row_args[8] = update_annotation(row_args[8], "in_telomere")
        chrom = row_args[0]
        key = row_args[1]+"-"+row_args[2]
        # Check to see if there are any inserts nearby
        if chrom not in x_large_filter:
            # Alternative chromsosme that there were no passing inserts for, same as below
            row_args.append("insert\tpassed_nearby_check")
            print("\t".join(row_args))
            continue
        if key not in x_large_filter[chrom]:
            # This insert is already a fail for some other reason. Check the interval tree to get the nearest insert
            nearby = intervaltree_xl[chrom][int(row_args[1]):int(row_args[2])]
            if len(nearby) > 0:
                # If there is something that shows up in more than half the samples remove it
                item = next(iter(nearby))
                row_args[8] = update_annotation(row_args[8], item.data)
                row_args.append("insert\tother_failure_"+str(len(item.data.split(','))))
                print("\t".join(row_args))
                continue
            else:
                row_args.append("insert\tpassed_nearby_check")
                print("\t".join(row_args))
                continue
        else:
            # now need to check the haplotypes
            insert_dict = x_large_filter[chrom][key]["inserts"]
            # First check the referece control sample. If an insert found here flag it as in other samples
            ref_count = 0
            if args.ref_sample in insert_dict:
                ref_count = insert_dict[args.ref_sample]["0"] + insert_dict[args.ref_sample]["1"] + insert_dict[args.ref_sample]["2"]
            if ref_count > 0:
                # have an insert in ref_control sample, filter is out
                format_and_print_line_insert_dict(insert_dict, row_args, "ref_control_sample")
                continue
            # Remaining inserts should only be ones for which no passing insert is found in the ref control sample
            insert_dict = large_filter[chrom][key]["inserts"]
            total_inserts = 0
            for insert in insert_dict[args.sample]:
                if insert == "fail":
                    continue
                total_inserts += insert_dict[args.sample][insert]
            total_haplotyped = 0
            haplotype_dict = large_filter[chrom][key]["haplotypes"]
            for insert in haplotype_dict[args.sample]:
                total_haplotyped += haplotype_dict[args.sample][insert]
            if len(row_args) == 19:
                hp = row_args[18]
            else:
                hp = row_args[19]
            if hp == "0":
                # Insert is on a read that isn't haplotyped. Check to see if there are inserts in the region that are haplotyped
                if total_inserts == insert_dict[args.sample][hp] and total_haplotyped == haplotype_dict[args.sample][hp]:
                    # All inserts and reads are from same haplotype on sample and no haplotype has been called
                    # If at least 30% of the reads have inserts, then its likley that this is more polymorphic than a hotspot
                    if float(total_inserts)/total_haplotyped > 0.3:
                        format_and_print_line_insert_dict(insert_dict, row_args, "No_haplotype_likely_polymorphic")
                        continue
                    # Otherwise we need to check every other sample and see if they have the same trend at this position.
                    # If all samples have no haplotyped reads at this position, and more than 30% of reads have an insert:
                    # Flag as polymorphic rather than hotspot
                    all_sample_insert_count = total_inserts
                    all_sample_hap_count = total_haplotyped
                    other_inserts_hap = False
                    other_samples_hap = False
                    for sample in insert_dict:
                        if sample == args.sample or sample == args.ref_sample:
                            continue
                        all_sample_insert_count += insert_dict[sample]["0"]
                        if insert_dict[sample]["1"] != 0 or insert_dict[sample]["2"] != 0:
                            other_inserts_hap = True
                        all_sample_hap_count += haplotype_dict[sample]["0"]
                        if haplotype_dict[sample]["1"] != 0 or haplotype_dict[sample]["2"] != 0:
                            other_samples_hap = True
                    if other_inserts_hap or other_samples_hap:
                        # There are haplotyped reads in the other samples at this position, default to a failure case as below
                        format_and_print_line_insert_dict(insert_dict, row_args, "likely_alignment_mistake_multi_sample")
                        continue
                    elif float(all_sample_insert_count)/all_sample_hap_count > 0.3:
                        format_and_print_line_insert_dict(insert_dict, row_args, "No_haplotype_likely_polymorphic_multi_sample")
                        continue
                    else:
                        # Need to now check a tight window for failure inserts on the reference sample.
                        ref_count = 0
                        insert_dict_small = small_filter[chrom][key]["inserts"]
                        if args.ref_sample in insert_dict_small:
                            ref_count = insert_dict_small[args.ref_sample]["0"] + insert_dict_small[args.ref_sample]["1"] + insert_dict_small[args.ref_sample]["2"] + insert_dict_small[args.ref_sample]["fail"]
                        if ref_count > 0:
                            # have an insert in ref_control sample, filter is out
                            format_and_print_line_insert_dict(insert_dict_small, row_args, "ref_control_sample_fail")
                            continue
                        found = False
                        if row_args[0] in insert_filter:
                            # Check to see if insert is in the insert filter set
                            for item in insert_filter[row_args[0]]:
                                if row_args[3] == insert_filter[row_args[0]][item][3] and abs(int(row_args[1]) - int(insert_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(insert_filter[row_args[0]][item][2])) <= 20:
                                    # Failure
                                    found = True
                                    format_and_print_line_insert_dict(insert_dict_small, row_args, "insert_filter_fail")
                                    break
                        if row_args[0] in mapq_0_filter and not found:
                            for item in mapq_0_filter[row_args[0]]:
                                if abs(int(row_args[1]) - int(mapq_0_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(mapq_0_filter[row_args[0]][item][2])) <= 20:
                                    found = True
                                    format_and_print_line_insert_dict(insert_dict_small, row_args, "mapq_0_filter_fail")
                                    break
                        if not found:
                            row_args.append("insert\thap_0_pass")
                            print("\t".join(row_args))
                        continue
                else:
                    # We have inserts found within the window on other reads in the region and sample that are haployped
                    # Likely an error. Flag and remove
                    if total_inserts != insert_dict[args.sample][hp]:
                        format_and_print_line_insert_dict(insert_dict, row_args, "likely_alignment_mistake_single")
                        continue
                    else:
                        ref_count = 0
                        insert_dict_small = small_filter[chrom][key]["inserts"]
                        if args.ref_sample in insert_dict_small:
                            ref_count = insert_dict_small[args.ref_sample]["0"] + insert_dict_small[args.ref_sample]["1"] + insert_dict_small[args.ref_sample]["2"] + insert_dict_small[args.ref_sample]["fail"]
                        if ref_count > 0:
                            # have an insert in ref_control sample, filter is out
                            format_and_print_line_insert_dict(insert_dict_small, row_args, "ref_control_sample_fail")
                            continue
                        found = False
                        if row_args[0] in insert_filter:
                            # Check to see if insert is in the insert filter set
                            for item in insert_filter[row_args[0]]:
                                if row_args[3] == insert_filter[row_args[0]][item][3] and abs(int(row_args[1]) - int(insert_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(insert_filter[row_args[0]][item][2])) <= 20:
                                    # Failure
                                    found = True
                                    format_and_print_line_insert_dict(insert_dict_small, row_args, "insert_filter_fail")
                                    break
                        if row_args[0] in mapq_0_filter and not found:
                            for item in mapq_0_filter[row_args[0]]:
                                if abs(int(row_args[1]) - int(mapq_0_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(mapq_0_filter[row_args[0]][item][2])) <= 20:
                                    found = True
                                    format_and_print_line_insert_dict(insert_dict_small, row_args, "mapq_0_filter_fail")
                                    break
                        if not found:
                            row_args.append("insert\thap_0_pass_2")
                            print("\t".join(row_args))
                        continue
            # Insert here is on a read that is haplotyped
            # See if a majority of reads that are of this haplotype contain an insert and if the insert is shared between both haplotypes
            # If we see > 50% of reads of the haplotype and only this haplotype with an insert than polymorphic
            # If inserts found here on both haplotypes, then dig deeper
            hp_insert_fraction = float(insert_dict[args.sample][hp])/haplotype_dict[args.sample][hp]
            other_hp_insert_fraction = 0
            if hp == "1" and haplotype_dict[args.sample]["2"] > 0:
                other_hp_insert_fraction = float(insert_dict[args.sample]["2"])/haplotype_dict[args.sample]["2"]
            elif hp == "2" and haplotype_dict[args.sample]["1"] > 0:
                other_hp_insert_fraction = float(insert_dict[args.sample]["1"])/haplotype_dict[args.sample]["1"]
            # Check to see if insert is haplotype specific, and at what fraction
            # make sure that at least half the reads are haplotyped at this position
            total_called_hp = haplotype_dict[args.sample]["1"] + haplotype_dict[args.sample]["2"]
            hp_total_fraction = float(total_called_hp)/(total_called_hp + haplotype_dict[args.sample]["0"])
            if hp_insert_fraction > 0.5 and other_hp_insert_fraction < 0.1:# and hp_total_fraction > 0.5:
                # insert is likely haplotype specific within sample. Call polymporphic
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_haplotype_specific_single_hp")
                continue
            elif hp_insert_fraction > 0.5 and other_hp_insert_fraction > 0.5:# and hphp_total_fraction > 0.5:
                # Insert appears in both haplotypes
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_polymorphic_both_hp")
                continue
            # If we get here we have an insert that is haplotyped, but not enough inserts on a haplotype to call it polymorphic
            # Check the other samples and repeat checks. If we do get more than 50% of inserts on a haplotype over all samples at region, call polymorphic
            hp_insert_count = insert_dict[args.sample][hp]
            hp_total_count = haplotype_dict[args.sample][hp]
            other_hp_insert_count = 0
            other_hp_total_count = 0
            all_read_count = total_called_hp + haplotype_dict[args.sample]["0"]
            if hp == "1" and haplotype_dict[args.sample]["2"] > 0:
                other_hp_insert_count = insert_dict[args.sample]["2"]
                other_hp_total_count = haplotype_dict[args.sample]["2"]
            elif hp == "2" and haplotype_dict[args.sample]["1"] > 0:
                other_hp_insert_count = insert_dict[args.sample]["1"]
                other_hp_total_count = haplotype_dict[args.sample]["1"]

            for sample in insert_dict:
                if sample == args.sample or sample == args.ref_sample:
                    continue
                hp_insert_count += insert_dict[sample][hp]
                hp_total_count += haplotype_dict[sample][hp]
                if hp == "1":
                    other_hp_insert_count += insert_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["2"]
                elif hp == "2":
                    other_hp_insert_count += insert_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["1"]
                all_read_count += haplotype_dict[sample]["0"] + haplotype_dict[sample]["1"] + haplotype_dict[sample]["2"]
                total_called_hp += haplotype_dict[sample]["1"] + haplotype_dict[sample]["2"]
            hp_insert_fraction = 0
            other_hp_insert_fraction = 0
            total_called_hp_fraction = 0
            if hp_total_count > 0:
                hp_insert_fraction = float(hp_insert_count)/hp_total_count
            if other_hp_total_count > 0:
                other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
            if all_read_count > 0:
                total_called_hp_fraction = float(total_called_hp)/all_read_count
            if hp_insert_fraction > 0.3 and other_hp_insert_fraction < 0.1:# and total_called_hp_fraction > 0.5:
                # Likely polymorphic when looking at all non reference control samples
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_polymorphic_single_hp_multi_sample")
                continue
            elif hp_insert_fraction > 0.3 and other_hp_insert_fraction > 0.3:# and total_called_hp_fraction > 0.5:
                # Likely polymorphic when looking at all non reference control samples, both hp
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_polymorphic_both_hp_multi_sample")
                continue
            # If here we don't have enough haplotyped reads to make the call, or we don't see enough inserts specific to a haplotype
            # Check to see if there are many haplotype specific softclips
            insert_dict = large_filter[chrom][key]["inserts"]
            haplotype_dict = large_filter[chrom][key]["haplotypes"]
            softclip_dict = large_filter[chrom][key]["softclips"]
            hp_insert_count = insert_dict[args.sample][hp]
            if args.sample in softclip_dict:
                hp_insert_count += softclip_dict[args.sample][hp]
            hp_total_count = haplotype_dict[args.sample][hp]
            other_hp_insert_count = 0
            other_hp_total_count = 0
            if hp == "1" and haplotype_dict[args.sample]["2"] > 0:
                other_hp_insert_count = insert_dict[args.sample]["2"]
                if args.sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[args.sample]["2"]
                other_hp_total_count = haplotype_dict[args.sample]["2"]
            elif hp == "2" and haplotype_dict[args.sample]["1"] > 0:
                other_hp_insert_count = insert_dict[args.sample]["1"]
                if args.sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[args.sample]["1"]
                other_hp_total_count = haplotype_dict[args.sample]["1"]
            for sample in insert_dict:
                if sample == args.sample or sample == args.ref_sample:
                    continue
                hp_insert_count += insert_dict[sample][hp]
                if sample in softclip_dict:
                    hp_insert_count += softclip_dict[sample][hp]
                hp_total_count += haplotype_dict[sample][hp]
                if hp == "1":
                    other_hp_insert_count += insert_dict[sample]["2"]
                    if sample in softclip_dict:
                        other_hp_insert_count += softclip_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["2"]
                elif hp == "2":
                    other_hp_insert_count += insert_dict[sample]["1"]
                    if sample in softclip_dict:
                        other_hp_insert_count += softclip_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["1"]
            hp_insert_fraction = 0
            other_hp_insert_fraction = 0
            total_called_hp_fraction = 0
            if hp_total_count > 0:
                hp_insert_fraction = float(hp_insert_count)/hp_total_count
            if other_hp_total_count > 0:
                other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
            if hp_insert_fraction > 0.75 and other_hp_insert_fraction < 0.2:# and total_called_hp_fraction > 0.5:
                # Likely polymorphic when looking at all non reference control samples
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_polymorphic_single_hp_multi_sample_softclip")
                continue
            elif hp_insert_fraction > 0.75 and other_hp_insert_fraction > 0.75:# and total_called_hp_fraction > 0.5:
                # Likely polymorphic when looking at all non reference control samples, both hp
                format_and_print_line_insert_dict(insert_dict, row_args, "liklely_polymorphic_both_hp_multi_sample_softclip")
                continue
            # These should be real inserts
            # Last check to see if small distance failures in reference
            ref_count = 0
            insert_dict_small = small_filter[chrom][key]["inserts"]
            if args.ref_sample in insert_dict_small:
                ref_count = insert_dict_small[args.ref_sample]["0"] + insert_dict_small[args.ref_sample]["1"] + insert_dict_small[args.ref_sample]["2"] + insert_dict_small[args.ref_sample]["fail"]
            if ref_count > 0:
                # have an insert in ref_control sample, filter is out
                format_and_print_line_insert_dict(insert_dict_small, row_args, "ref_control_sample_fail")
                continue
            found = False
            if row_args[0] in insert_filter:
                # Check to see if insert is in the insert filter set
                for item in insert_filter[row_args[0]]:
                    if row_args[3] == insert_filter[row_args[0]][item][3] and abs(int(row_args[1]) - int(insert_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(insert_filter[row_args[0]][item][2])) <= 20:
                        # Failure
                        found = True
                        format_and_print_line_insert_dict(insert_dict_small, row_args, "insert_filter_fail")
                        break
            if row_args[0] in mapq_0_filter and not found:
                for item in mapq_0_filter[row_args[0]]:
                    if abs(int(row_args[1]) - int(mapq_0_filter[row_args[0]][item][1])) <= 20 and abs(int(row_args[2]) - int(mapq_0_filter[row_args[0]][item][2])) <= 20:
                        found = True
                        format_and_print_line_insert_dict(insert_dict_small, row_args, "mapq_0_filter_fail")
                        break
            if not found:
                row_args.append("insert\tpassed_hp_check")
                print("\t".join(row_args))
