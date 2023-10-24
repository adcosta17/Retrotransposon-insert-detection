import pysam
import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser( description='Evaluate the number of samples needed to have a significant difference ')
parser.add_argument('--iterations', required=False, type=int, default=5000)
parser.add_argument('--max-samples', required=False, type=int, default=5000)
parser.add_argument('--depth', required=False, type=int, default=30)
parser.add_argument('--effect', required=False, type=float, default=2.0)
parser.add_argument('--mean', required=False, type=float, default=1.0)
args = parser.parse_args()

e = args.effect
base_mean = args.mean * e
treated_mean = args.mean
min_sizes = defaultdict(int)
for i in range(args.iterations):
    all_base = np.random.poisson(int(base_mean), args.max_samples)
    all_treated = np.random.poisson(int(treated_mean), args.max_samples)
    min_size = args.max_samples
    for j in range(args.max_samples):
        base_set = all_base[0:j+1]
        treated_set = all_treated[0:j+1]
        pval = stats.ttest_ind(base_set, treated_set,equal_var=False).pvalue
        #print(pval)
        if min_size == args.max_samples and pval < 0.05:
            min_size = j+1
            break
    min_sizes[i] = min_size
    #print(min_sizes[i])
print(str(args.depth)+"\t"+str(treated_mean)+"\t"+str(base_mean)+"\t"+str(e)+"\t"+"\t".join([str(x) for x in list(min_sizes.values())]))
