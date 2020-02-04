import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser( description='Generate plot of coverage across repeat family by percentile')
parser.add_argument('--input', required=True)
parser.add_argument('--output-dir', required=True)
parser.add_argument('--repbase-fasta', type=str, required=True)
args = parser.parse_args()

repbase_lengths = defaultdict(int)
# Read in repbase fasta and get lengths for each sequence
with open(args.repbase_fasta) as rep_in:
    line = rep_in.readline()
    while line != '':
        name = line.strip().split(">")[1]
        line = rep_in.readline()
        if line == '':
            break
        seq_len = len(line.strip())
        repbase_lengths[name] = seq_len
        line = rep_in.readline()


family_bins = {}
family_bins["LINE"] = [0] * 101
family_bins["SINE"] = [0] * 101
family_bins["SVA"] = [0] * 101
family_bins["ERV"] = [0] * 101

family_bins_insert_pos = {}
family_bins_insert_pos["LINE"] = [0] * 101
family_bins_insert_pos["SINE"] = [0] * 101
family_bins_insert_pos["SVA"] = [0] * 101
family_bins_insert_pos["ERV"] = [0] * 101

with open(args.input) as csvfile:
    #reader = csv.DictReader(csvfile, delimiter="\t")
    count = 0
    for row in csvfile:
        row_args = row.strip().split("\t")
        if count == 0:
            count = 1
            continue
        if row_args[8] == "PASS" and row_args[19] == "insert":
            start = int((int(row_args[16])/repbase_lengths[row_args[9]])*100)
            end = int((int(row_args[17])/repbase_lengths[row_args[9]])*100)
            i = start
            while i <= end:
                if "LINE" in row_args[9]:
                    family_bins["LINE"][i] += 1
                elif "SINE" in row_args[9]:
                    family_bins["SINE"][i] += 1
                elif "SVA" in row_args[9]:
                    family_bins["SVA"][i] += 1
                elif "ERV" in row_args[9]:
                    family_bins["ERV"][i] += 1
                i += 1

            l_insert = int(row_args[5]) - int(row_args[4])
            start = int(((int(row_args[13])-int(row_args[4]))/l_insert)*100)
            end = int(((int(row_args[14])-int(row_args[4]))/l_insert)*100)
            i = start
            while i <= end:
                if "LINE" in row_args[9]:
                    family_bins_insert_pos["LINE"][i] += 1
                elif "SINE" in row_args[9]:
                    family_bins_insert_pos["SINE"][i] += 1
                elif "SVA" in row_args[9]:
                    family_bins_insert_pos["SVA"][i] += 1
                elif "ERV" in row_args[9]:
                    family_bins_insert_pos["ERV"][i] += 1
                i += 1

# Generate histograms based on counts
x = np.array(range(101))
y = np.array(family_bins["LINE"])
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'LINE.png')
plt.clf()
y = np.array(family_bins["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'SINE.png')
plt.clf()
y = np.array(family_bins["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'ERV.png')
plt.clf()
y = np.array(family_bins["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'SVA.png')
plt.clf()

y = np.array(family_bins_insert_pos["LINE"])
plt.bar(x, height= y)
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(args.output_dir+'LINE_inserts.png')
plt.clf()
y = np.array(family_bins_insert_pos["SINE"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'SINE_inserts.png')
plt.clf()
y = np.array(family_bins_insert_pos["ERV"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'ERV_inserts.png')
plt.clf()
y = np.array(family_bins_insert_pos["SVA"])
plt.bar(x, height= y)
plt.savefig(args.output_dir+'SVA_inserts.png')
plt.clf()