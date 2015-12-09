def v_num_of_null_smaller (x,y):
    return float(np.sum( y < x))

#!/usr/bin/env python
import sys, os, glob, re
#import math
import numpy as np
from datetime import datetime
from numpy import array, zeros
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
#from scipy.stats import chisquare
import argparse
import gzip
parser = argparse.ArgumentParser(description='calculate and add q fdr values to a sync file with the p value at the end of each line.') 
parser.add_argument("--sync","-s", dest="obsfile", help="sync file with observed P values as last column", required=True)
parser.add_argument("--null","-n", dest="nullfile", help="file with null P values",required=True)
parser.add_argument("--out","-o", dest="outfile", help="outfile",required=True)
parser.add_argument("--verbose","-v",  action="store_true", dest="verbose", help="print operation performed to standard out (default: False)",default=False)
args = parser.parse_args()
infile = vars(args)['obsfile']
nullfile = vars(args)['nullfile']
outfile =  vars(args)['outfile']
verbose =  vars(args)['verbose']
if re.search("\.b?gz$",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")	
obs_pvals = []
chr_pos_ref = []
if verbose:
    start_time = datetime.now()
    print "Reading cmhout file at "+str(start_time)
for line in inf:
    entries = line.rstrip().split("\t")
    chr_pos_ref.append("\t".join(entries[0:3]))
    try:
        obs_pvals.append(float(entries[-1]))
    except:
        print "Problem with last field of:\n"+line
inf.close()

obs_pvals = np.array(obs_pvals)
if verbose:
    print str(len(obs_pvals))+" lines read"
    print "Reading null Pvalue file at "+str(datetime.now())
counter = 0
if re.search("\.b?gz$",nullfile):
    inf = gzip.open(nullfile,'rb')
else:
    inf = open(nullfile,"r")	
null_pvals = [float(x.rstrip()) for x in inf]
inf.close()
null_pvals=np.sort(np.array(null_pvals))
if verbose:
    print str(len(null_pvals))+" lines read"
    print "calculating ranks for observed Pvalues at "+str(datetime.now())
r_obs_P = np.argsort(obs_pvals.argsort())
r_obs_P += 1
n_ratio = float(len(obs_pvals))/len(null_pvals)
null_smaller = np.searchsorted(null_pvals,obs_pvals)
q_value = null_smaller * n_ratio / r_obs_P
#v_num_of_null_smaller = np.vectorize(lambda ( x , null_pvals) : float(np.sum( null_pvals < x)))
if verbose:
    print "Calculating q values and  writing output files at "+str(datetime.now())
if outfile == "stdout":
    outf=sys.stdout
else:
    outf = open(outfile,"w")
#for i,chrpos in enumerate(chr_pos_ref):
for i in np.argsort(obs_pvals):
    #q_value = null_smaller * n_ratio / r_obs_P[i]
    print >>outf, "{0}\t{1}\t{2}".format(chr_pos_ref[i],obs_pvals[i],q_value[i])
    if verbose and (i+1)%10000 == 0:
        print str(i+1)+" lines written"
outf.close()
if verbose:
    print "Done after "+str(datetime.now() - start_time)

    
