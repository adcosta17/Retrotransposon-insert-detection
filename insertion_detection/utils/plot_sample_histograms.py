import argparse
import sys
import datetime
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

parser = argparse.ArgumentParser( description='Generate plot of insert lengths and samples found near')
parser.add_argument('--output-dir', required=True)
parser.add_argument('--suffix', default=".all.read_insertions.repbase_annotated.high_confidence.chimeric_filtered.tsv")
parser.add_argument('--folder', default="read_analysis")
parser.add_argument('--samples', required=True)
args = parser.parse_args()

insert_sizes = {}
insert_sizes["LINE"] = [0] * 101
insert_sizes["SINE"] = [0] * 101
insert_sizes["SVA"] = [0] * 101
insert_sizes["ERV"] = [0] * 101

for sample in args.samples.split(","):
    print(sample,file=sys.stderr)
    sample_tsv = "./"+sample+"/"+args.folder+"/"+sample+args.suffix
    with open(sample_tsv) as csvfile:
        #reader = csv.DictReader(csvfile, delimiter="\t")
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            if row_args[8] == "PASS":
                i = int(len(row_args[7])/100)
                if i > 100:
                    i = 100
                if "LINE" in row_args[9]:
                    insert_sizes["LINE"][i] += 1
                elif "SINE" in row_args[9]:
                    insert_sizes["SINE"][i] += 1
                elif "SVA" in row_args[9]:
                    insert_sizes["SVA"][i] += 1
                elif "ERV" in row_args[9]:
                    insert_sizes["ERV"][i] += 1

# Generate histograms based on counts
x = np.array(range(101))
y = np.array(insert_sizes["LINE"])
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/insert_sizes_LINE.png')
plt.clf()
y = np.array(insert_sizes["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_SINE.png')
plt.clf()
y = np.array(insert_sizes["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ERV.png')
plt.clf()
y = np.array(insert_sizes["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_SVA.png')
plt.clf()

insert_sizes = {}
insert_sizes["LINE"] = [0] * 101
insert_sizes["SINE"] = [0] * 101
insert_sizes["SVA"] = [0] * 101
insert_sizes["ERV"] = [0] * 101
insert_sizes_no_sample_1 = {}
insert_sizes_no_sample_1["LINE"] = [0] * 101
insert_sizes_no_sample_1["SINE"] = [0] * 101
insert_sizes_no_sample_1["SVA"] = [0] * 101
insert_sizes_no_sample_1["ERV"] = [0] * 101

sample_counts = [0] * 12
sample_counts2 = [0] * 11

sample_counts_no_sample_1 = [0] * 12
sample_counts_no_sample_1_2 = [0] * 11

sample_counts_annotated = [0] * 12
sample_counts_annotated_2 = [0] * 11

sample_counts_annotated_no_sample_1 = [0] * 12
sample_counts_annotated_no_sample_1_2 = [0] * 11

sample_counts_families = {}
sample_counts_families["LINE"] = [0] * 12
sample_counts_families["SINE"] = [0] * 12
sample_counts_families["SVA"] = [0] * 12
sample_counts_families["ERV"] = [0] * 12

sample_counts_families_no_sample1 = {}
sample_counts_families_no_sample1["LINE"] = [0] * 12
sample_counts_families_no_sample1["SINE"] = [0] * 12
sample_counts_families_no_sample1["SVA"] = [0] * 12
sample_counts_families_no_sample1["ERV"] = [0] * 12

for sample in args.samples.split(","):
    print(sample, file=sys.stderr)
    sample_tsv = "./"+sample+"/"+args.folder+"/"+sample+".all.read_insertions_and_soft_clipped.repbase_annotated.high_confidence.chimeric_filtered.reference_ava_filtered.tsv"
    with open(sample_tsv) as csvfile:
        #reader = csv.DictReader(csvfile, delimiter="\t")
        count = 0
        for row in csvfile:
            row_args = row.strip().split("\t")
            if count == 0:
                count = 1
                continue
            i = int(len(row_args[7])/100)
            if i > 100:
                i = 100
            if row_args[8] == "PASS":
                if "LINE" in row_args[9]:
                    insert_sizes["LINE"][i] += 1
                    insert_sizes_no_sample_1["LINE"][i] += 1
                    sample_counts_families["LINE"][0] += 1
                    sample_counts_families_no_sample1["LINE"][0] += 1
                elif "SINE" in row_args[9]:
                    insert_sizes["SINE"][i] += 1
                    insert_sizes_no_sample_1["SINE"][i] += 1
                    sample_counts_families["SINE"][0] += 1
                    sample_counts_families_no_sample1["SINE"][0] += 1
                elif "SVA" in row_args[9]:
                    insert_sizes["SVA"][i] += 1
                    insert_sizes_no_sample_1["SVA"][i] += 1
                    sample_counts_families["SVA"][0] += 1
                    sample_counts_families_no_sample1["SVA"][0] += 1
                elif "ERV" in row_args[9]:
                    insert_sizes["ERV"][i] += 1
                    insert_sizes_no_sample_1["ERV"][i] += 1
                    sample_counts_families["ERV"][0] += 1
                    sample_counts_families_no_sample1["ERV"][0] += 1
                sample_counts[0] += 1
                sample_counts2[0] += 1
                sample_counts_annotated[0] += 1
                sample_counts_annotated_2[0] += 1
                sample_counts_no_sample_1[0] += 1
                sample_counts_no_sample_1_2[0] += 1
                sample_counts_annotated_no_sample_1[0] += 1
                sample_counts_annotated_no_sample_1_2[0] += 1
            else:
                only_in_sample = True
                sc = defaultdict(int)
                for item in row_args[8].split(','):
                    if "in_sample" in item:
                        sc[item.replace("_fail",'')] += 1
                    else:
                        # Failure other than in sample, need to skip
                        only_in_sample = False
                count = 0
                in_sample1 = False
                for item in sc:
                    count += 1
                    if "in_sample_DDTS_0001" in item:
                        in_sample1 = True
                if count > 0:
                    #print(str(count))
                    #print(sc)
                    sample_counts[count] += 1
                    if not in_sample1:
                        sample_counts_no_sample_1[count] += 1
                    if row_args[9] != "no_repbase_mapping":
                        sample_counts_annotated[count] += 1
                        if not in_sample1:
                            sample_counts_annotated_no_sample_1[count] += 1
                    if count < 11:
                        sample_counts2[count] += 1
                        if not in_sample1:
                            sample_counts_no_sample_1_2[count] += 1
                        if row_args[9] != "no_repbase_mapping":
                            sample_counts_annotated_2[count] += 1
                            if not in_sample1:
                                sample_counts_annotated_no_sample_1_2[count] += 1
                    if only_in_sample:
                        if "LINE" in row_args[9]:
                            sample_counts_families["LINE"][count] += 1
                            if not in_sample1:
                                insert_sizes_no_sample_1["LINE"][i] += 1
                                sample_counts_families_no_sample1["LINE"][count] += 1
                        elif "SINE" in row_args[9]:
                            sample_counts_families["SINE"][count] += 1
                            if not in_sample1:
                                insert_sizes_no_sample_1["SINE"][i] += 1
                                sample_counts_families_no_sample1["SINE"][count] += 1
                        elif "SVA" in row_args[9]:
                            sample_counts_families["SVA"][count] += 1
                            if not in_sample1:
                                insert_sizes_no_sample_1["SVA"][i] += 1
                                sample_counts_families_no_sample1["SVA"][count] += 1
                        elif "ERV" in row_args[9]:
                            sample_counts_families["ERV"][count] += 1
                            if not in_sample1:
                                insert_sizes_no_sample_1["ERV"][i] += 1
                                sample_counts_families_no_sample1["ERV"][count] += 1

# Generate histograms based on counts
x = np.array(range(101))
y = np.array(insert_sizes["LINE"])
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/insert_sizes_ref_filtered_LINE.png')
plt.clf()
y = np.array(insert_sizes["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_SINE.png')
plt.clf()
y = np.array(insert_sizes["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_ERV.png')
plt.clf()
y = np.array(insert_sizes["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_SVA.png')
plt.clf()

y = np.array(insert_sizes_no_sample_1["LINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_no_sample_1_LINE.png')
plt.clf()
y = np.array(insert_sizes_no_sample_1["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_no_sample_1_SINE.png')
plt.clf()
y = np.array(insert_sizes_no_sample_1["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_no_sample_1_ERV.png')
plt.clf()
y = np.array(insert_sizes_no_sample_1["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/insert_sizes_ref_filtered_no_sample_1_SVA.png')
plt.clf()

# Generate histograms based on counts
x = np.array(range(12))
y = np.array(sample_counts)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts.png')
plt.clf()

x = np.array(range(12))
y = np.array(sample_counts_annotated)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_annotated.png')
plt.clf()

x = np.array(range(12))
y = np.array(sample_counts_no_sample_1)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_no_sample_1.png')
plt.clf()

x = np.array(range(12))
y = np.array(sample_counts_annotated_no_sample_1)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_annotated_no_sample_1.png')
plt.clf()

y = np.array(sample_counts_families["LINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_all_LINE.png')
plt.clf()
y = np.array(sample_counts_families["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_all_SINE.png')
plt.clf()
y = np.array(sample_counts_families["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_all_SVA.png')
plt.clf()
y = np.array(sample_counts_families["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_all_ERV.png')
plt.clf()


y = np.array(sample_counts_families_no_sample1["LINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_no_sample_1_LINE.png')
plt.clf()
y = np.array(sample_counts_families_no_sample1["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_no_sample_1_SINE.png')
plt.clf()
y = np.array(sample_counts_families_no_sample1["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_no_sample_1_SVA.png')
plt.clf()
y = np.array(sample_counts_families_no_sample1["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'/sample_counts_no_sample_1_ERV.png')
plt.clf()

x = np.array(range(11))
y = np.array(sample_counts2)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_0_10.png')
plt.clf()

x = np.array(range(11))
y = np.array(sample_counts_annotated_2)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_0_10_annotated.png')
plt.clf()


x = np.array(range(11))
y = np.array(sample_counts_no_sample_1_2)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_0_10_no_sample_1.png')
plt.clf()

x = np.array(range(11))
y = np.array(sample_counts_annotated_no_sample_1_2)
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'/sample_counts_0_10_annotated_no_sample_1.png')
plt.clf()
