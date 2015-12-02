# get the light and dark allele from a OR file
# reads file and prints the major and minor allele, as well as the light allele for each population

import sys, os, re
#import numpy as np
import gzip

#from scipy import stats
#from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description="""
reads file and prints the major and minor allele, as well as the light allele for each population
"""
                                 )

parser.add_argument("--in","-i", dest="infile", help="file with OR", required=True)

args = parser.parse_args()
infile = vars(args)['infile']
if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
for line in inf:
    if re.match("\#",line):
        continue
    line=line.rstrip()
    fields=line.split()
    if fields[0] == "":
        continue
    pops = (len(fields)-3)/4
    lights=[]
    for i in range(0,pops):
        if fields[4+i*4] == 'na':
            lights.append("-")
        elif  float(fields[4+i*4]) <= 1:
            lights.append(fields[2][1])
        elif float(fields[4+i*4]) > 1:
            lights.append(fields[2][0])
        else:
            lights.append("-")
    print "\t".join(fields[0:3])+"\t"+"\t".join(lights)
