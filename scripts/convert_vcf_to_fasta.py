import pysam
import argparse
import sys
import csv
from collections import defaultdict

parser = argparse.ArgumentParser( description='Convert a .read_insertions.tsv file to fasta')
parser.add_argument('--input', required=True)
args = parser.parse_args()

vcf = pysam.VariantFile(args.input)
for rec in vcf.fetch():
    if rec.info["SVTYPE"] == "INS" and int(rec.info["SVLEN"]) > 50:
        if "INS" in rec.alts[0]:
            continue
        print(">%s\n%s" % (rec.id, rec.alts[0]))
