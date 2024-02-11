# TEMP_____remove_duplicates.py

import pandas as pd
import openpyxl
from collections import defaultdict
import numpy as np
import gzip
import pysam
import re
import sys
import argparse
from Levenshtein import distance as levenshtein_distance

def hamming(str1, str2):
    count = 0
    assert (len(str1)==len(str2)), "length doesn't match"
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            count += 1
    return(count)



parser = argparse.ArgumentParser()
parser.add_argument('-inbam', type=str, help="sorted bam file needed to be dedupped")
parser.add_argument('-outbam', type=str, help="dedupped bam file")

args = parser.parse_args()

INBAMFILE = args.inbam
OUTBAMFILE = args.outbam

DIST = 8 # two reads distance no less than 8 to be regard as unique 
UMIDIST = 2 # maxium distance of UMI difference to be two umis


inbam = pysam.AlignmentFile(INBAMFILE, "rb")
outbam = pysam.AlignmentFile(OUTBAMFILE, "wb", template=inbam)


prev_bcd = ""
prev_achr = ""
prev_umi = ""
prev_aend = 0
bcd_umi_dict = defaultdict(list)  # to avoid the orginal bam file have different bcd lines not consecutively ordered, in that case it's not enough to only check for the previous line

with inbam as bam:
  with outbam as out:

    for read in bam.fetch(until_eof = True):
      bcd = re.search('^.*\|(BC.*)#[ATCGN]{12}$', read.qname).group(1)
      umi = re.search('^.*#([ATCGN]{12})$', read.qname).group(1)
      aend = read.aend - read.alen # read aligned end minus aligned length => aligned start
      achr = read.reference_name
      
      if( achr == prev_achr and (aend - prev_aend)<DIST ):
        if( bcd in bcd_umi_dict.keys() and any([levenshtein_distance(umi,x)<=UMIDIST for x in bcd_umi_dict[bcd]]) ):
            pass
        else:
          out.write(read)
          bcd_umi_dict[bcd].append(umi)
      else:
        out.write(read)
        bcd_umi_dict = defaultdict(list)
        bcd_umi_dict[bcd]=[umi]
      
      prev_bcd =  bcd 
      prev_umi =  umi 
      prev_aend = aend
      prev_achr = achr
      




      
