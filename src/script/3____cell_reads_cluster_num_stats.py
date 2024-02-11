# Get cell reads count, cell cluster count from bam file

# Date: 06042022
# log: 06162022, changed cluster regex, add -lib* before there is no -lib info in the MB barcode
#      06162022, add number of reads in each cluster


import pandas as pd
import openpyxl
from collections import defaultdict
import numpy as np
import gzip
import pysam
import re
import sys
import argparse
import csv
from pathlib import Path
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-inbam', type=str, help="sorted bam file needed to be dedupped")
parser.add_argument('-outdir', type=str, help="directory to write cell reads and cell cluster")

args = parser.parse_args()

INBAMFILE = args.inbam

# name inherit from bam file
file_name = Path(INBAMFILE).stem

# Path function will take care of black slash, default will remove it
OUT_CELLREADS = Path(args.outdir) / (file_name + ".cell_reads.csv")
OUT_CELLCLU = Path(args.outdir) / (file_name + ".cell_clusters.csv")
OUT_CLUREADS = Path(args.outdir) / (file_name + ".clusters_size.csv")


# ---- stream in bam file and count -------

cell_reads_dict = defaultdict(int)
cell_clusters_dict = defaultdict(set)
cell_clusters_size = defaultdict(int)

inbam = pysam.AlignmentFile(INBAMFILE, "rb")

with inbam as bam:
  for read in bam.fetch(until_eof = True):
    groups = re.search('^.*\|(BC.*)-([ATCGN]{16}-lib.*)#[ATCGN]{12}$', read.qname)
    CB = groups[1]
    MB = groups[2]
    cell_reads_dict[CB]+=1 # number of reads from each cell
    cell_clusters_dict[CB].add(MB) # unique MB in each cell
    cell_clusters_size[CB+","+CB+"-"+MB]+=1 # number of reads in each cluster
inbam.close()

# get number of clusters per cell
clusters_per_cell = {k: len(v) for k, v in cell_clusters_dict.items()}


# ---- write dictionary to csv --------
pd.DataFrame.from_dict(data=cell_reads_dict, orient='index').set_axis(['num_reads'],axis=1).sort_values('num_reads',ascending=False).to_csv(OUT_CELLREADS)

pd.DataFrame.from_dict(data=cell_clusters_size, orient='index').set_axis(['cluster_num_reads'],axis=1).sort_values('cluster_num_reads',ascending=False).to_csv(OUT_CLUREADS)

pd.DataFrame.from_dict(data=clusters_per_cell, orient='index').set_axis(['num_clusters'],axis=1).sort_values('num_clusters',ascending=False).to_csv(OUT_CELLCLU)


