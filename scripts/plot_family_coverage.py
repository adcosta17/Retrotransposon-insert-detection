import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import mappy as mp

plt.rcParams.update({'font.size': 22})

parser = argparse.ArgumentParser( description='Generate plot of coverage across repeat family by percentile')
parser.add_argument('--samples', required=True)
parser.add_argument('--seqs', required=True)
parser.add_argument('--output-dir', required=True)
args = parser.parse_args()

ctg_sizes = {}
with pysam.FastxFile(args.seqs) as fh:
    for entry in fh:
        ctg_sizes[entry.name] = len(entry.sequence)
        #print(entry.name+"\t"+str(len(entry.sequence)))

family_sizes = defaultdict(list)
all_family_sizes = defaultdict(list)
positions = {}
positions["LINE"] = [0]*101
positions["SINE"] = [0]*101
positions["SVA"] = [0]*101
positions["All_LINE"] = [0]*101
positions["All_SINE"] = [0]*101
positions["All_SVA"] = [0]*101

a = mp.Aligner(args.seqs)  # load or build index
if not a: raise Exception("ERROR: failed to load/build index")

for sample in args.samples.split(','):
    print(sample)
    sample_family_sizes = defaultdict(list)
    sample_all_family_sizes = defaultdict(list)
    sample_positions = {}
    sample_positions["LINE"] = [0]*101
    sample_positions["SINE"] = [0]*101
    sample_positions["SVA"] = [0]*101
    sample_positions["All_LINE"] = [0]*101
    sample_positions["All_SINE"] = [0]*101
    sample_positions["All_SVA"] = [0]*101
    with open(sample+"/winnow_realign_read_analysis/"+sample+".read_insertions.repbase_annotated.mapq_ct_filtered.ma_filtered.ref_filtered_haplotype_checked.tsv") as csvfile:
        count = 0
        for row in csvfile:
            if count == 0:
                count += 1
                continue
            row_args = row.strip().split("\t")
            if len(row_args[7]) > 7000 or len(row_args[7]) < 100:
                continue
            if row_args[8] == "PASS":
                if "LINE" in row_args[9]:
                    family_sizes["LINE"].append(len(row_args[7]))
                    sample_family_sizes["LINE"].append(len(row_args[7]))
                    for hit in a.map(row_args[7]):
                        if "LINE" in hit.ctg:
                            start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                            end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                            i = start
                            while i <= end:
                                positions["LINE"][i] += 1
                                sample_positions["LINE"][i] += 1
                                i += 1
                            break
                elif "SINE" in row_args[9]:
                    family_sizes["SINE"].append(len(row_args[7]))
                    sample_family_sizes["SINE"].append(len(row_args[7]))
                    for hit in a.map(row_args[7]):
                        if "SINE" in hit.ctg:
                            start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                            end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                            i = start
                            while i <= end:
                                positions["SINE"][i] += 1
                                sample_positions["SINE"][i] += 1
                                i += 1
                            break
                elif "SVA" in row_args[9]:
                    family_sizes["SVA"].append(len(row_args[7]))
                    sample_family_sizes["SVA"].append(len(row_args[7]))
                    for hit in a.map(row_args[7]):
                        if "SVA" in hit.ctg:
                            start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                            end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                            i = start
                            while i <= end:
                                positions["SVA"][i] += 1
                                sample_positions["SVA"][i] += 1
                                i += 1
                            break
                elif "ERV" in row_args[9]:
                    family_sizes["ERV"].append(len(row_args[7]))
                    sample_family_sizes["ERV"].append(len(row_args[7]))
            if "LINE" in row_args[9]:
                all_family_sizes["LINE"].append(len(row_args[7]))
                sample_all_family_sizes["LINE"].append(len(row_args[7]))
                for hit in a.map(row_args[7]):
                    if "LINE" in hit.ctg:
                        start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                        end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                        i = start
                        while i <= end:
                            positions["All_LINE"][i] += 1
                            sample_positions["All_LINE"][i] += 1
                            i += 1
                        break
            elif "SINE" in row_args[9]:
                all_family_sizes["SINE"].append(len(row_args[7]))
                sample_all_family_sizes["SINE"].append(len(row_args[7]))
                for hit in a.map(row_args[7]):
                    if "SINE" in hit.ctg:
                        start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                        end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                        i = start
                        while i <= end:
                            positions["All_SINE"][i] += 1
                            sample_positions["All_SINE"][i] += 1
                            i += 1
                        break
            elif "SVA" in row_args[9]:
                all_family_sizes["SVA"].append(len(row_args[7]))
                sample_all_family_sizes["SVA"].append(len(row_args[7]))
                for hit in a.map(row_args[7]):
                    if "SVA" in hit.ctg:
                        start = int(hit.r_st/ctg_sizes[hit.ctg]*100)
                        end = int(hit.r_en/ctg_sizes[hit.ctg]*100)
                        i = start
                        while i <= end:
                            positions["All_SVA"][i] += 1
                            sample_positions["All_SVA"][i] += 1
                            i += 1
                        break
            elif "ERV" in row_args[9]:
                all_family_sizes["ERV"].append(len(row_args[7]))
                sample_all_family_sizes["ERV"].append(len(row_args[7]))
    os.mkdir(args.output_dir+"/"+sample)
    plt.hist(sample_family_sizes["LINE"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("LINE Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Passing_Inserts_LINE.png')
    plt.clf()
    plt.hist(sample_family_sizes["SINE"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("Alu Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Passing_Inserts_Alu.png')
    plt.clf()
    plt.hist(sample_family_sizes["SVA"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("SVA Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Passing_Inserts_SVA.png')
    plt.clf()
    plt.hist(sample_family_sizes["ERV"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("ERV Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Passing_Inserts_ERV.png')
    plt.clf()
    plt.hist(sample_all_family_sizes["LINE"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("LINE Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_All_Inserts_LINE.png')
    plt.clf()
    plt.hist(sample_all_family_sizes["SINE"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("Alu Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_All_Inserts_Alu.png')
    plt.clf()
    plt.hist(sample_all_family_sizes["SVA"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("SVA Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_All_Inserts_SVA.png')
    plt.clf()
    plt.hist(sample_all_family_sizes["ERV"], 250, facecolor='blue', alpha=0.5)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("ERV Insert Size Distribution")
    fig = plt.gcf()
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_All_Inserts_ERV.png')
    plt.clf()
    x = np.array(range(101))
    y = np.array(sample_positions["LINE"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "1500", "3000", "4500", "6000"])
    plt.xlabel("Location")
    plt.ylabel("Coverage")
    plt.title("LINE Element Coverage Breakdown")
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_Passing_Inserts_LINE.png')
    plt.clf()
    y = np.array(sample_positions["SINE"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "75", "150", "225", "300"])
    plt.xlabel("Location")
    plt.ylabel("Coverage")
    plt.title("Alu Element Coverage Breakdown")
    plt.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_Passing_Inserts_Alu.png')
    plt.clf()
    y = np.array(sample_positions["SVA"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "400", "800", "1200", "1600"])
    plt.xlabel("Location")
    plt.ylabel("Coverage")
    plt.title("SVA Element Coverage Breakdown")
    plt.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_Passing_Inserts_SVA.png')
    plt.clf()
    y = np.array(sample_positions["All_LINE"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "1500", "3000", "4500", "6000"])
    plt.xlabel("Location")
    plt.title("LINE Element Coverage Breakdown")
    plt.ylabel("Coverage")
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5)
    fig.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_All_Inserts_LINE.png')
    plt.clf()
    y = np.array(sample_positions["All_SINE"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "75", "150", "225", "300"])
    plt.xlabel("Location")
    plt.ylabel("Coverage")
    plt.title("Alu Element Coverage Breakdown")
    plt.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_All_Inserts_Alu.png')
    plt.clf()
    y = np.array(sample_positions["All_SVA"])
    plt.bar(x, height= y)
    plt.xticks([0, 25, 50, 75, 100], ["0", "400", "800", "1200", "1600"])
    plt.xlabel("Location")
    plt.ylabel("Coverage")
    plt.title("SVA Element Coverage Breakdown")
    plt.savefig(args.output_dir+"/"+sample+"/"+sample+'_Coverage_All_Inserts_SVA.png')
    plt.clf()


# Generate histograms based on counts
#print(family_sizes["LINE"])
#print(max(family_sizes["LINE"]))
plt.hist(family_sizes["LINE"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("LINE Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'Passing_Inserts_LINE.png')
plt.clf()
plt.hist(family_sizes["SINE"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("Alu Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'Passing_Inserts_Alu.png')
plt.clf()
plt.hist(family_sizes["SVA"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("SVA Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'Passing_Inserts_SVA.png')
plt.clf()
plt.hist(family_sizes["ERV"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("ERV Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'Passing_Inserts_ERV.png')
plt.clf()

plt.hist(all_family_sizes["LINE"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("LINE Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'All_Inserts_LINE.png')
plt.clf()
plt.hist(all_family_sizes["SINE"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("Alu Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'All_Inserts_Alu.png')
plt.clf()
plt.hist(all_family_sizes["SVA"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("SVA Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'All_Inserts_SVA.png')
plt.clf()
plt.hist(all_family_sizes["ERV"], 250, facecolor='blue', alpha=0.5)
plt.xlabel("Size")
plt.ylabel("Count")
plt.title("ERV Insert Size Distribution")
fig = plt.gcf()
fig.savefig(args.output_dir+'All_Inserts_ERV.png')
plt.clf()


# Generate histograms based on counts
x = np.array(range(101))
y = np.array(positions["LINE"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "1500", "3000", "4500", "6000"])
plt.xlabel("Location")
plt.ylabel("Coverage")
plt.title("LINE Element Coverage Breakdown")
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'Coverage_Passing_Inserts_LINE.png')
plt.clf()
y = np.array(positions["SINE"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "75", "150", "225", "300"])
plt.xlabel("Location")
plt.ylabel("Coverage")
plt.title("Alu Element Coverage Breakdown")
plt.savefig(args.output_dir+'Coverage_Passing_Inserts_Alu.png')
plt.clf()
y = np.array(positions["SVA"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "400", "800", "1200", "1600"])
plt.xlabel("Location")
plt.ylabel("Coverage")
plt.title("SVA Element Coverage Breakdown")
plt.savefig(args.output_dir+'Coverage_Passing_Inserts_SVA.png')
plt.clf()
y = np.array(positions["All_LINE"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "1500", "3000", "4500", "6000"])
plt.xlabel("Location")
plt.title("LINE Element Coverage Breakdown")
plt.ylabel("Coverage")
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'Coverage_All_Inserts_LINE.png')
plt.clf()
y = np.array(positions["All_SINE"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "75", "150", "225", "300"])
plt.xlabel("Location")
plt.ylabel("Coverage")
plt.title("Alu Element Coverage Breakdown")
plt.savefig(args.output_dir+'Coverage_All_Inserts_Alu.png')
plt.clf()
y = np.array(positions["All_SVA"])
plt.bar(x, height= y)
plt.xticks([0, 25, 50, 75, 100], ["0", "400", "800", "1200", "1600"])
plt.xlabel("Location")
plt.ylabel("Coverage")
plt.title("SVA Element Coverage Breakdown")
plt.savefig(args.output_dir+'Coverage_All_Inserts_SVA.png')
plt.clf()