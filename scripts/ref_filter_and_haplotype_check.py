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
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["haplotypes"]= {}
            if len(row_args) > 5:
                filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["haplotypes"] = get_haplotype_dict(row_args[5])
            filter_dict[row_args[0]][row_args[1]+"-"+row_args[2]]["softclips"] = {}
            if len(row_args) > 6:
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


def in_ref_sample(insert_dict, ref_sample, row_args, message):
    ref_count = 0
    if ref_sample in insert_dict:
        ref_count = insert_dict[ref_sample]["0"] + insert_dict[ref_sample]["1"] + insert_dict[ref_sample]["2"]
    if ref_count > 0:
        # have an insert in ref_control sample, filter is out
        #format_and_print_line_insert_dict(insert_dict, row_args, message)
        return True
    return False

def in_ref_sample_small(insert_dict, ref_sample, row_args, message):
    ref_count = 0
    if ref_sample in insert_dict:
        ref_count = insert_dict[ref_sample]["0"] + insert_dict[ref_sample]["1"] + insert_dict[ref_sample]["2"]
        if "fail" in insert_dict[ref_sample]:
            ref_count += insert_dict[ref_sample]["fail"]
    if ref_count > 0:
        # have an insert in ref_control sample, filter is out
        format_and_print_line_insert_dict(insert_dict, row_args, message)
        return True
    return False


def get_total_inserts(insert_dict, passed_sample):
    total_inserts = 0
    for sample in insert_dict:
        if sample == passed_sample:
            continue
        total_inserts += insert_dict[sample]["0"]
        total_inserts += insert_dict[sample]["1"]
        total_inserts += insert_dict[sample]["2"]
    return total_inserts


def filter_non_hp(insert_dict, ref_sample, haplotype_dict, softclip_dict, row_args, annotation, cutoff):
    to_print = False
    if int(row_args[1]) == 74727222 and row_args[0] == "chr1":
        to_print = True
    all_sample_insert_count = 0
    all_sample_hap_count = 0
    for sample in insert_dict:
        if sample == ref_sample:
            continue
        all_sample_insert_count += insert_dict[sample]["0"]
        all_sample_insert_count += insert_dict[sample]["1"]
        all_sample_insert_count += insert_dict[sample]["2"]
        if sample in haplotype_dict:
            all_sample_hap_count += haplotype_dict[sample]["0"]
            all_sample_hap_count += haplotype_dict[sample]["1"]
            all_sample_hap_count += haplotype_dict[sample]["2"]
    if all_sample_hap_count > 0:
        if to_print:
            print("Base " + row_args[3] + " " + str(all_sample_insert_count) + " " + str(all_sample_hap_count), file=sys.stderr)
        if float(all_sample_insert_count)/all_sample_hap_count >= cutoff:
            format_and_print_line_insert_dict(insert_dict, row_args, annotation)
            return True
    # now check softclips
    all_sample_insert_count = 0
    all_sample_hap_count = 0
    for sample in insert_dict:
        if sample == ref_sample:
            continue
        all_sample_insert_count += insert_dict[sample]["0"]
        all_sample_insert_count += insert_dict[sample]["1"]
        all_sample_insert_count += insert_dict[sample]["2"]
        if sample in softclip_dict:
            all_sample_insert_count += softclip_dict[sample]["0"]
            all_sample_insert_count += softclip_dict[sample]["1"]
            all_sample_insert_count += softclip_dict[sample]["2"]
        if sample in haplotype_dict:
            all_sample_hap_count += haplotype_dict[sample]["0"]
            all_sample_hap_count += haplotype_dict[sample]["1"]
            all_sample_hap_count += haplotype_dict[sample]["2"]
    if all_sample_hap_count > 0:
        if to_print:
            print("Softclips " + row_args[3] + " " + str(all_sample_insert_count) + " " + str(all_sample_hap_count), file=sys.stderr)
        if float(all_sample_insert_count)/all_sample_hap_count >= cutoff:
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_softclip")
            return True
    # Need to now check if it is following a haplotyped pattern
    # If its a hap0 it may be hap specific and being missed. Meaning it may fail a larger check overall 
    # Think 50% of reads are one hap and less than 100% of those show/support an insert
    # Need to account for this
    for hp in ["1", "2"]:
        all_sample_insert_count = 0
        all_sample_hap_count = 0
        all_sample_other_hp_insert_count = 0
        all_sample_other_hp_count = 0
        for sample in insert_dict:
            if sample == ref_sample:
                continue
            all_sample_insert_count += insert_dict[sample]["0"]
            all_sample_insert_count += insert_dict[sample][hp]
            if hp == "1":
                all_sample_other_hp_insert_count += insert_dict[sample]["2"]
            else:
                all_sample_other_hp_insert_count += insert_dict[sample]["1"]
            if sample in softclip_dict:
                all_sample_insert_count += softclip_dict[sample]["0"]
                all_sample_insert_count += softclip_dict[sample][hp]
                if hp == "1":
                    all_sample_other_hp_insert_count += softclip_dict[sample]["2"]
                else:
                    all_sample_other_hp_insert_count += softclip_dict[sample]["1"]
            if sample in haplotype_dict:
                all_sample_hap_count += haplotype_dict[sample]["0"]
                all_sample_hap_count += haplotype_dict[sample][hp]
                if hp == "1":
                    all_sample_other_hp_count += haplotype_dict[sample]["2"]
                else:
                    all_sample_other_hp_count += haplotype_dict[sample]["1"]
        if all_sample_hap_count > 0:
            if to_print:
                print("Softclips hp " + hp + " " + row_args[3] + " " + str(all_sample_insert_count) + " " + str(all_sample_hap_count), file=sys.stderr)
            if float(all_sample_insert_count)/all_sample_hap_count >= cutoff:
                format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_softclip_hp"+hp)
                return True
    return False

parser = argparse.ArgumentParser( description='Remove insertions that are near insertions in any other sample')
parser.add_argument('--input', required=True)
#parser.add_argument('--sc', required=True)
parser.add_argument('--multi-reference-filter-xlarge-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-large-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-medium-tsv', type=str, required=True)
parser.add_argument('--multi-reference-filter-small-tsv', type=str, required=True)
parser.add_argument('--telomere-filter', required=True)
parser.add_argument('--centromere-filter', required=True)
parser.add_argument('--mapq0-filter', required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument('--ref-sample', type=str, default="DDTS_0001_nn_n")
parser.add_argument('--reference-repeat-bed', type=str, required=True)
parser.add_argument('--max-distance', type=int, default=100)
parser.add_argument('--max-common', type=int, default=1)
parser.add_argument('--sample-count', type=int, default=1)
parser.add_argument('--num-samples-needed-to-remove', type=int, default=1)
args = parser.parse_args()


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

mapq0_filter = {}
with open(args.mapq0_filter) as in_filter:
    for row in in_filter:
        row_args = row.strip().split("\t")
        if row_args[0] not in mapq0_filter:
            mapq0_filter[row_args[0]] = {}
        mapq0_filter[row_args[0]][row_args[1]+":"+row_args[2]] = row_args

# read input file and only emit records not near a record in the interval tree
with open(args.input) as csvfile:
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        failure = False
        if count == 0:
            count = 1
            print(row.strip()+"\tsource\tsample_count\ttype")
            continue
        to_print = False
        #if row_args[3] == "b8ce41be-5e49-4739-937c-4e67c17846c9" or row_args[3] == "6d04ec05-6b3d-4c51-ba3b-60ef22c4a859" or row_args[3] == "d717c2ef-06b8-4d6e-9b0a-a928d0c1ef77" or row_args[3] == "f7b9ad3d-f06d-4a9f-9d50-5205ae8186cf":
        #    to_print = True
        #    print(row_args[3], file=sys.stderr)
        chrom = row_args[0]
        key = row_args[1]+"-"+row_args[2]
        # Check to see if there are any inserts nearby
        if chrom not in x_large_filter:
            # Alternative chromsosme that there were no passing inserts for, same as below
            row_args.append("insert\ton_alt_chrom\tNA")
            print("\t".join(row_args))
            continue
        if key not in x_large_filter[chrom]:
            # This insert is already a fail for some other reason. Check the interval tree to get the nearest insert
            nearby = intervaltree_xl[chrom][int(row_args[1]):int(row_args[2])]
            if len(nearby) > 0:
                # If there is something that shows up in more than half the samples remove it
                item = next(iter(nearby))
                row_args[8] = update_annotation(row_args[8], item.data)
                row_args.append("insert\tother_failure_"+str(len(item.data.split(',')))+"\tNA")
                print("\t".join(row_args))
                continue
            else:
                row_args.append("insert\tpassed_nearby_check\tNA")
                print("\t".join(row_args))
                continue
        # now need to check the haplotypes
        insert_dict = x_large_filter[chrom][key]["inserts"]
        haplotype_dict = x_large_filter[chrom][key]["haplotypes"]
        softclip_dict = x_large_filter[chrom][key]["softclips"]
        # First check the referece control sample. If an insert found here flag it as in other samples
        annotation = ""
        if in_ref_sample(insert_dict, args.ref_sample, row_args,"ref_control_sample"):
            annotation = "ref_control_sample"
        # Remaining inserts should only be ones for which no passing insert is found in the ref control sample
        if len(row_args) == 19:
            hp = row_args[18]
        else:
            hp = row_args[19]
        if hp == "0":
            # No haplotype for this insert.
            if filter_non_hp(insert_dict, args.ref_sample, haplotype_dict, softclip_dict, row_args, annotation+"_"+"No_haplotype_likely_polymorphic_multi_sample_500", 0.8):
                continue
            insert_dict = large_filter[chrom][key]["inserts"]
            haplotype_dict = large_filter[chrom][key]["haplotypes"]
            softclip_dict = large_filter[chrom][key]["softclips"]
            if filter_non_hp(insert_dict, args.ref_sample, haplotype_dict, softclip_dict, row_args, annotation+"_"+"No_haplotype_likely_polymorphic_multi_sample_200", 0.5):
                continue
            insert_dict = medium_filter[chrom][key]["inserts"]
            haplotype_dict = medium_filter[chrom][key]["haplotypes"]
            softclip_dict = medium_filter[chrom][key]["softclips"]
            if filter_non_hp(insert_dict, args.ref_sample, haplotype_dict, softclip_dict, row_args, annotation+"_"+"No_haplotype_likely_polymorphic_multi_sample_100", 0.3):
                continue
            #insert_dict = small_filter[chrom][key]["inserts"]
            #haplotype_dict = small_filter[chrom][key]["haplotypes"]
            #softclip_dict = small_filter[chrom][key]["softclips"]
            #if filter_non_hp(insert_dict, args.ref_sample, haplotype_dict, softclip_dict, row_args, annotation+"_"+"No_haplotype_likely_polymorphic_multi_sample_50", 0.25):
            #    continue
            if annotation == "ref_control_sample":
                format_and_print_line_insert_dict(insert_dict, row_args, annotation)
                continue
            if in_ref_sample_small(small_filter[chrom][key]["inserts"], args.ref_sample, row_args, annotation+"_"+"ref_control_sample_fail"):
                continue
            if row_args[0] in mapq0_filter:
                # Check to see if insert is in the insert filter set
                mapq0_found = False
                for item in mapq0_filter[row_args[0]]:
                    if abs(int(row_args[1]) - int(mapq0_filter[row_args[0]][item][1])) <= 500 and abs(int(row_args[2]) - int(mapq0_filter[row_args[0]][item][2])) <= 500:
                        format_and_print_line_insert_dict(insert_dict, row_args, "mapq0_filter")
                        mapq0_found = True
                        break
                if mapq0_found:
                    continue
            total_inserts = get_total_inserts(insert_dict, args.sample)
            if total_inserts == 0:
                row_args.append("insert\thap_0_pass\tnovel")
            else:
                row_args.append("insert\thap_0_pass\thotspot")
            print("\t".join(row_args))
            continue
        # Inserts here have a haplotype. Check if they follow the expected pattern
        insert_dict = x_large_filter[chrom][key]["inserts"]
        haplotype_dict = x_large_filter[chrom][key]["haplotypes"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in haplotype_dict:
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if hp_insert_fraction >= 0.8 and other_hp_insert_fraction < 0.2:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_500")
            continue
        elif hp_insert_fraction >= 0.8 and other_hp_insert_fraction >= 0.8:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_500")
            continue
        # If here we don't have enough haplotyped reads to make the call, or we don't see enough inserts specific to a haplotype
        # Check to see if there are many haplotype specific softclips
        softclip_dict = x_large_filter[chrom][key]["softclips"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in softclip_dict:
                hp_insert_count += softclip_dict[sample][hp]
                hp_insert_count += softclip_dict[sample]["0"]
            if sample in haplotype_dict:    
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if hp_insert_fraction > 0.8 and other_hp_insert_fraction < 0.2:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_softclip_500")
            continue
        elif hp_insert_fraction > 0.8 and other_hp_insert_fraction > 0.8:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_softclip_500")
            continue
        ##
        ## Check using the large filter next
        ##
        insert_dict = large_filter[chrom][key]["inserts"]
        haplotype_dict = large_filter[chrom][key]["haplotypes"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in haplotype_dict:
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if hp_insert_fraction >= 0.5 and other_hp_insert_fraction < 0.5:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_200")
            continue
        elif hp_insert_fraction >= 0.5 and other_hp_insert_fraction >= 0.5:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_200")
            continue
        # If here we don't have enough haplotyped reads to make the call, or we don't see enough inserts specific to a haplotype
        # Check to see if there are many haplotype specific softclips
        softclip_dict = large_filter[chrom][key]["softclips"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in softclip_dict:
                hp_insert_count += softclip_dict[sample][hp]
                hp_insert_count += softclip_dict[sample]["0"]
            if sample in haplotype_dict:    
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["1"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if hp_insert_fraction > 0.5 and other_hp_insert_fraction < 0.5:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_softclip_200")
            continue
        elif hp_insert_fraction > 0.5 and other_hp_insert_fraction > 0.5:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_softclip_200")
            continue
        # These should be real inserts
        # Last check to see if small distance failures in reference
        if annotation == "ref_control_sample":
            format_and_print_line_insert_dict(insert_dict, row_args, annotation)
            continue
        if in_ref_sample_small(small_filter[chrom][key]["inserts"], args.ref_sample, row_args, annotation+"_"+"ref_control_sample_fail"):
            continue
        ##
        ## Check using the Medium filter next
        ##
        insert_dict = medium_filter[chrom][key]["inserts"]
        haplotype_dict = medium_filter[chrom][key]["haplotypes"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            if to_print:
                print("Sample: "+sample, file=sys.stderr)
                print(insert_dict[sample][hp], file=sys.stderr)
                print(insert_dict[sample]["0"], file=sys.stderr)
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in haplotype_dict:
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
                if to_print:
                    print(haplotype_dict[sample][hp], file=sys.stderr)
                    print(haplotype_dict[sample]["0"], file=sys.stderr)
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if to_print:
            print("After", file=sys.stderr)
            print(hp_total_count, file=sys.stderr)
            print(hp_insert_count, file=sys.stderr)
            print(hp_insert_fraction, file=sys.stderr)
            print(other_hp_insert_fraction, file=sys.stderr)
        if hp_insert_fraction >= 0.3 and other_hp_insert_fraction >= 0.3:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_100")
            continue
        elif hp_insert_fraction >= 0.3:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_100")
            continue
        #elif hp_insert_fraction >= 0.4 and other_hp_insert_fraction >= 0.4:
        #    # Likely polymorphic when looking at all non reference control samples, both hp
        #    format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_100")
        #    continue
        # If here we don't have enough haplotyped reads to make the call, or we don't see enough inserts specific to a haplotype
        # Check to see if there are many haplotype specific softclips
        softclip_dict = medium_filter[chrom][key]["softclips"]
        hp_insert_count = 0
        hp_total_count = 0
        other_hp_insert_count = 0
        other_hp_total_count = 0
        for sample in insert_dict:
            if sample == args.ref_sample:
                continue
            hp_insert_count += insert_dict[sample][hp]
            hp_insert_count += insert_dict[sample]["0"]
            if sample in softclip_dict:
                hp_insert_count += softclip_dict[sample][hp]
                hp_insert_count += softclip_dict[sample]["0"]
            if sample in haplotype_dict:    
                hp_total_count += haplotype_dict[sample][hp]
                hp_total_count += haplotype_dict[sample]["0"]
            if hp == "1":
                other_hp_insert_count += insert_dict[sample]["2"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["2"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["2"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
            elif hp == "2":
                other_hp_insert_count += insert_dict[sample]["1"]
                if sample in softclip_dict:
                    other_hp_insert_count += softclip_dict[sample]["1"]
                if sample in haplotype_dict:
                    other_hp_total_count += haplotype_dict[sample]["1"]
                    other_hp_total_count += haplotype_dict[sample]["0"]
        hp_insert_fraction = 0
        other_hp_insert_fraction = 0
        if hp_total_count > 0 and hp_insert_count > 1:
            hp_insert_fraction = float(hp_insert_count)/hp_total_count
        if other_hp_total_count > 0:
            other_hp_insert_fraction = float(other_hp_insert_count)/other_hp_total_count
        if hp_insert_fraction > 0.3 and other_hp_insert_fraction > 0.3:
            # Likely polymorphic when looking at all non reference control samples, both hp
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_softclip_100")
            continue
        elif hp_insert_fraction > 0.3:
            # Likely polymorphic when looking at all non reference control samples
            format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_single_hp_multi_sample_softclip_100")
            continue
        #elif hp_insert_fraction > 0.3 and other_hp_insert_fraction > 0.3:
        #    # Likely polymorphic when looking at all non reference control samples, both hp
        #    format_and_print_line_insert_dict(insert_dict, row_args, annotation+"_"+"liklely_polymorphic_both_hp_multi_sample_softclip_100")
        #    continue
        # These should be real inserts
        # Last check to see if small distance failures in reference
        if annotation == "ref_control_sample":
            format_and_print_line_insert_dict(insert_dict, row_args, annotation)
            continue
        if in_ref_sample_small(small_filter[chrom][key]["inserts"], args.ref_sample, row_args, annotation+"_"+"ref_control_sample_fail"):
            continue
        ## Final check to see if we've missed anything
        insert_dict = x_large_filter[chrom][key]["inserts"]
        haplotype_dict = x_large_filter[chrom][key]["haplotypes"]
        if row_args[0] in mapq0_filter:
                # Check to see if insert is in the insert filter set
                mapq0_found = False
                for item in mapq0_filter[row_args[0]]:
                    if abs(int(row_args[1]) - int(mapq0_filter[row_args[0]][item][1])) <= 500 and abs(int(row_args[2]) - int(mapq0_filter[row_args[0]][item][2])) <= 500:
                        format_and_print_line_insert_dict(insert_dict, row_args, "mapq0_filter")
                        mapq0_found = True
                        break
                if mapq0_found:
                    continue
        total_inserts = get_total_inserts(insert_dict, args.sample)
        if total_inserts == 0:
            row_args.append("insert\tpassed_hp_check\tnovel")
        else:
            row_args.append("insert\tpassed_hp_check\thotspot")
        print("\t".join(row_args))
