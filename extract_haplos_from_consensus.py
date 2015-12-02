
import argparse
from argparse import RawTextHelpFormatter
import sys
import re
import os

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Extracts the information on the provided haplotypes (--samples) from a haplotype conensus format and prints the resulting haplotype consensus file into the stdout.

""")
parser.add_argument('--consensus',"-c", required=True, dest='consensus', type=str, help="Haplotype file in consensus format.")
parser.add_argument('--samples',"-s", required=True, dest='samples', type=str, help="Number of haplotypes that should be in the output consensus file comma separated, e.g. 3-14,17,20-25,27-31,32-43,46,50-60 ")

args = parser.parse_args()


#-------------------
# function
#-------------------


def translate_samples(samples):
	a=samples.split(",")
	s=[]
	for i in a:
		x=i.split("-")
		if len(x) ==2:
			s.extend(range(int(x[0])-1,int(x[1])))
		else:
			s.append(int(x[0])-1)
	#s=map(int,samples.split(","))
	return s
	

####################
# main
####################

samples=translate_samples(args.samples)
#print samples

f=open(args.consensus,"r")

for l in f:
	chr,pos,ref,ref_r,haplo,qual,info=l.split()
	h_out=[];i_out=[]
#	print haplo
	for i in samples:
		h_out.append(haplo[i-1])
		i_out.append(qual[i-1])
	h_out="".join(h_out)
	i_out="".join(i_out)
	print "\t".join([chr,pos,ref,ref_r,h_out,i_out,info])
	

