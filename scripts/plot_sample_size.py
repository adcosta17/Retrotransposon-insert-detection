import argparse
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pysam
from collections import defaultdict
import mappy as mp
import matplotlib.font_manager as font_manager
import matplotlib


parser = argparse.ArgumentParser( description='Generate plots of sample size simulations')
parser.add_argument('--input-dir', required=True)
parser.add_argument('--output-dir', required=True)
args = parser.parse_args()

means = np.arange(0.1, 5.1, 0.1)
effects = np.arange(1.1, 5.1, 0.1)

for m in means:
    # For each mean plot the min samples needed to get a significant result 95% of the time for each effect size
    # Line plot, X axis is effect size, Y Axis is min samples. Capped at 5000
    mean = round(m,2)
    print(mean)
    x = []
    x_label = []
    y = []
    for e in effects:
        base_mean = 0
        val = 5000
        effect = round(e,2)
        with open(args.input_dir+"/Sample_"+str(mean)+"_"+str(effect)+"_.txt", 'r') as in_file:
            for line in in_file:
                row = line.strip().split("\t")
                coverage = row[0]
                base_mean = round(float(row[2]),2)
                counts = np.array(list(map(int, row[4:])))
                #print(counts)
                val = np.percentile(counts, 95)
                if val == 5000:
                    val = 0
        x.append(base_mean)
        x_label.append(str(base_mean))
        y.append(val)
    Y = np.array(y)
    X = np.array(x)
    plt.plot(X, Y)
    #plt.xticks([0, 25, 50, 75, 100], ["0", "400", "800", "1200", "1600"])
    plt.xlabel("Single Treatment Mean")
    plt.ylabel("Min Samples for 95% Significant")
    plt.title("Mean "+str(mean)+" for RT Treatment vs Single Treatment Min Samples")
    plt.savefig(args.output_dir+'/Samples_'+str(mean)+".png")
    plt.clf()
