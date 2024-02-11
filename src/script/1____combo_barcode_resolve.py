# 1____DNA_end_barcode_resolve.py

import pandas as pd
import openpyxl
from collections import defaultdict
import numpy as np
import gzip 
import pybktree
import re
import os
import h5py # /dataOS/wenxingzhao/software/anaconda/install/bin/pip install xxx
import time
import psutil
import logging
from Bio import SeqIO
from Levenshtein import distance as levenshtein_distance
import datetime
import sys
import argparse
import logging

# ---- parse arguments --------

# Create the parser and add arguments
parser = argparse.ArgumentParser()
parser.add_argument('-R1', type=str, help="filepath to R1, include 10xBC and UMI")
parser.add_argument('-R2', type=str, help="filepath to R2, include cell barcode and RNA/DNA insert")
parser.add_argument('-outRNA', type=str, help="filepath to the barcode resolved RNA insertion fastq file")
parser.add_argument('-outDNA', type=str, help="filepath to the barcode resolved DNA insertion fastq file")
parser.add_argument('-lib', type=str, help="lib/cluster name when used in sequencing")
parser.add_argument('-log', type=str, help="filepath to the log file")

args = parser.parse_args()

R1_file = args.R1
R2_file = args.R2
outRNA_file = args.outRNA
outDNA_file = args.outDNA
lib_name = args.lib
logfile = args.log


##### tmp modification ######
## add index as complex barcode ## 

#basename = os.path.basename(R1_file)
sample_idx = lib_name # lib5


logging.basicConfig(level=logging.DEBUG, filename=logfile, filemode="a+",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")

logging.info('---- start ----')




# ---- dictionary of reference cell barcode ----
BCD_DIR = "/snakemake/barcode"

BC3 = pd.read_csv(BCD_DIR+"/dna_barcode_shuf_3.txt",lineterminator="\n", header=None)[0].tolist()
BC3_dict = dict(zip(BC3, ["BC3_" + str(x) for x in np.arange(1, len(BC3))]))

BC2 = pd.read_csv(BCD_DIR+"/dna_barcode_shuf_2.txt",lineterminator="\n", header=None)[0].tolist()
BC2_dict = dict(zip(BC2, ["BC2_" + str(x) for x in np.arange(1, len(BC2))]))

BC1 = pd.read_csv(BCD_DIR+"/dna_barcode_shuf_1.txt",lineterminator="\n", header=None)[0].tolist()
BC1_dict = dict(zip(BC1, ["BC1_" + str(x) for x in np.arange(1, len(BC1))]))


# ---- 10x barcode whitelist -----
WL = "/snakemake/barcode/3M-february-2018.txt.gz"
my_file = gzip.open(WL, "rb")
content = my_file.read().split(b'\n')
content = [x.decode() for x in content]
whitelist_10x = set(content)

logging.info("---- read in 10X BC ----")

# -------- add helper functions -------------
def hamming(str1, str2):
    count = 0 
    assert (len(str1)==len(str2)), "length doesn't match"
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            count += 1
    return(count)

def SUFFIX(pattern):
    return(pattern[1:])

def neighbors(pattern, d):
    if d == 0:
        return(pattern)
    if len(pattern) == 1:
        return(['A','T','C','G'])
    nlist = []
    suffixn = neighbors(SUFFIX(pattern),d)
    for string in suffixn:
        if hamming(string, SUFFIX(pattern)) < d:
            for x in ['A','T','C','G']:
                nlist.append(x+string)
        else:
            nlist.append(pattern[0]+string)
    return(nlist)


def n_mis_neighbors(words_dict, d):
  '''
  words dict with keys the reference barcode and values the barcode name
  '''
  hash_ = {}
  words = list(words_dict.keys())
  for i in words:
    word = i
    name = words_dict[i]
    # remaining bcd string
    bcd_nei = neighbors(word, d)
    for b in bcd_nei:
      hash_[b] = name
  return hash_

BC3_hash = n_mis_neighbors(BC3_dict, 2) # allow two mismatches in barcode
BC2_hash = n_mis_neighbors(BC2_dict, 2)
BC1_hash = n_mis_neighbors(BC1_dict, 2)

# -------- RNA/DNA adaptor sequences -------------
UMI5_all = neighbors("AAAAA",5) # correspond to NNNNN

DNA_adp = set(neighbors("ACAACGCACAGTGTCTAG",2))
RNA_adp_1 = [y + x for y in neighbors("CGCTT", 1) for x in UMI5_all] 
RNA_adp = set([h+t for h in RNA_adp_1 for t in neighbors("ATAGCATT", 1)])


logging.info("---- finish constructing barcode ----")





# -------- Function: decode barcode sequences -------------

# BC3 is variable length 
# we will have a function take in full barcode sequences 
# return total barcode length, and barcode combination

def debarcode(seq, qual):
  # BC3: 14-17bp,  BC2:14 bp,  BC1:14 bp
  # input seq length: 17+7+14+7+14 = 59
  BC3_14bp = seq[0:14]
  BC3_15bp = seq[0:15]
  BC3_16bp = seq[0:16]
  BC3_17bp = seq[0:17]
  
  BC3_idx = np.argwhere([x in BC3_hash for x in [BC3_14bp, BC3_15bp, BC3_16bp, BC3_17bp]]).tolist()
  
  if( len(BC3_idx) == 0):
    return("BCFAIL") # insert index is 0
    
  else:
    true_BC3_seq = [BC3_14bp, BC3_15bp, BC3_16bp, BC3_17bp][BC3_idx[0][0]]
    BC3_len = len(true_BC3_seq)
    BC3 = BC3_hash.get(true_BC3_seq, "NA")
    
    BC2_seq = seq[(BC3_len+7): (BC3_len+7+14)]
    BC2 = BC2_hash.get(BC2_seq, "NA")
    
    BC1_seq = seq[(BC3_len+7+14+7): (BC3_len+7+14+7+14)]
    BC1 = BC1_hash.get(BC1_seq, "NA")
    
    BCD = BC3 + "-" + BC2 + "-" + BC1

    if "NA" in BCD:
      return("BCFAIL")
    
    else:
      DNA_adaptor = seq[(BC3_len+54) : (BC3_len+72)] # 54 = 7+14+7+14+7+5, adaptor_len=18bp
      RNA_adaptor = seq[(BC3_len+49) : (BC3_len+67)] # 49 = 7+14+7+14+7, adaptor_len=18bp
      
      if DNA_adaptor in DNA_adp:
        TYPE = "DNA"
        insert = seq[(BC3_len+72) : ] # 72 = 7+14+7+14+7+5+18
        quality = qual[(BC3_len+72) : ]
        return((TYPE, BCD, insert, quality))
      
      elif RNA_adaptor in RNA_adp:
        TYPE = "RNA"
        insert = seq[(BC3_len+70) : ] # keep consistent with DNA end so easy for fastq quality
        quality = qual[(BC3_len+70) : ]
        return((TYPE, BCD, insert, quality))
      
      else:
        return("NOADP")

      


# -------- read in files -------------

total_reads = 0
wrong_BC = 0
no_adp = 0
wrong_10x = 0
final_DNA_reads = 0
final_RNA_reads = 0

try:
  with gzip.open(R1_file,'r') as R1: #after phiX removal, change gzip.open to open
    with gzip.open(R2_file,'r') as R2:
      with open(outRNA_file,'w') as RNAout:
        with open(outDNA_file,'w') as DNAout:
          while True:
            R1_fq = []
            R2_fq = []
            total_reads += 1
            for i in range(4):
              R1_line = R1.readline()
              R2_line = R2.readline()
              if not R1_line and not R2_line:
                  break
              R1_fq.append(R1_line.decode().strip('\n'))
              R2_fq.append(R2_line.decode().strip('\n'))
              
            if len(R1_fq) == 0: # break at EOF 
              break
              
            R1_seq = R1_fq[1].strip('\n')
            R2_seq = R2_fq[1].strip('\n')
            R2_qual = R2_fq[3].strip('\n')
            R2_ID = R2_fq[0].split(" ")[0]
            
            
            if ("N" in R1_seq) or ("N" in R2_seq):
              pass
            
            bc10x = R1_seq[0:16]
            
            if bc10x in whitelist_10x:
              
              R1_UMI_seq = R1_seq[16:28]
              
              R2_resolve = debarcode(R2_seq, R2_qual)
              
              if R2_resolve == "BCFAIL":
                wrong_BC += 1
                continue
              
              if R2_resolve == "NOADP":
                no_adp += 1
                continue

              out_str = R2_ID + "|" + R2_resolve[1] + "-" + bc10x + "-" + sample_idx, "#" + R1_UMI_seq + "\n" + R2_resolve[2] + "\n+\n" + R2_resolve[3] + "\n"
              
              if R2_resolve[0]=="DNA":
                final_DNA_reads += 1
                DNAout.writelines(out_str)
              else:
                final_RNA_reads += 1
                RNAout.writelines(out_str)
                
            else:
              wrong_10x +=1
              continue

        DNAout.close()        
      RNAout.close()
        
except IOError as e:
    print('Operation failed: %s' % e.strerror)
    

# print total reads in raw fastq file
logging.info("total lines in file "+ R1_file +" : " + str(total_reads))
logging.info("incorrect 10x barcode (not perfect match) :"+ str(wrong_10x))
logging.info("reads incomplete CB :"+ str(wrong_BC))
logging.info("Final DNA reads :"+ str(final_DNA_reads))
logging.info("Final RNA reads :"+ str(final_RNA_reads))

logging.info("---- finished ----")
