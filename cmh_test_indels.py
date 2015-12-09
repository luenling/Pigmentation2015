def cmh(replicate_tables,confint=True):
	"""
	gets replicate table an optional boolean confint and returns a string with either pValue, odds ratio and upper and lower conf. interval odds ratio, or just pV and odds ratio, depending on the value of confint
	"""
	# import rpy2.robjects as robjects
	# from rpy2.robjects.packages import importr
	#from rpy2.robjects import r
	#thanks to robert or martin for example code
	#r=robjects.r
	response_robj = robjects.IntVector(replicate_tables)
	replicatecount = int(float(len(replicate_tables))/4.0)
	dim_robj=robjects.IntVector([2,2,replicatecount])
	response_rar=robjects.r['array'](response_robj,dim=dim_robj)
	#testres=r['mantelhaen.test'](response_rar)
	testres=r['try_CMHtest'](response_rar)
	if len(testres) > 1: 
		pvalue=testres[2][0]
		odds_int=[testres[3][0],testres[3][1]]
		odds_ratio=testres[4][0]
	else:
		pvalue=np.nan
		odds_int=[np.nan,np.nan]
		odds_ratio=np.nan
	#assert(pvalue<=1)
	#assert(pvalue>=0)
	if confint:
		return "\t".join([str(pvalue),str(odds_ratio),str(odds_int[0]),str(odds_int[1])])
	else:
		return "\t".join([str(pvalue),str(odds_ratio)]) 


def get_replicate_tables(rep_dict):
	replicate_tables = []
	for r in rep_dict:
		for t in rep_dict[r]:
			replicate_tables.extend(t)
	return replicate_tables


if __name__ == '__main__':
	#get one random table per replicate
	import sys, os, re
	import numpy as np
	#import Sync_parser
	#from numpy import array, zeros
	#import math
	import gzip
	import argparse
	import rpy2.robjects as robjects
	from rpy2.robjects.packages import importr
	parser = argparse.ArgumentParser(description="""create P values, odds ratio and conf intervals for a CMH test on an indel count sync file. 
 The file should have a format akin to a tab delimited sync file:
Chromosome\tposition\tallele1\tallele2,3,4...\t(count_a1:count_a2:count:...)*n
Only the first two colon separated count values are taken for each pools, so 100:20:0:0 is treated the same as 100:20 and the resulting OR is always for the the test a1:a2.  
 """) 

	parser.add_argument("--pops","-p", dest="pops", help="the population pairs in the CMH file to test, mutual exclusive with reps, if set, reps is ignored, also do multiple cmh tests with the same using \"|\" (eg. \"1:2|1:3,1:4\")",required=True)
	parser.add_argument("--in","-i", dest="infile", help="indel file to be read", required=True)
	parser.add_argument("--confint","-c", dest="confint",  action='store_false', help="confidence intervals", default=True)
	args = parser.parse_args()
	infile = vars(args)['infile']
	pops = vars(args)['pops']
	confint = vars(args)['confint']
	#thanks to robert or martin for example code
	# create R function for try
	r=robjects.r
	robjects.r('''
        try_CMHtest <- function(count_table) {
        try( mantelhaen.test(count_table),
        silent=TRUE
        )
        }
        ''')
	mult=False
	if pops:
		if '|' in pops:
			pops = pops.split('|')
		else:
			pops = [ pops ]
		pops = [  x.split(',')  for x in pops]
		pops = [ [  int(z) for y in x for z in y.split(":") ]  for x in pops ]
		for i in pops:
			assert(len(i)%2 == 0) or sys.exit('some populations not paired up :'+str(i))
		pops = [ np.array(x) for x in pops ]
		pops = [ x -1 for x in pops]
	else:
		sys.exit("populations not set")
	if re.search("\.b?gz",infile):
		inf = gzip.open(infile,'rb')
	else:
		inf = open(infile,"r")

	for line in inf:
		if re.match("\s*\#",line):
			# comment lines
			continue
		line=line.rstrip()
		entries=line.split()
		try:
			counts=[ [ int(z) for z in x.split(",")] for x in entries[4:]]
		except:
			continue
		thiscmhp=[]
		for i in pops:
			replicate_tables=[]
			for idx in i: 
				replicate_tables.extend(counts[idx])
			thiscmhp.append(cmh(replicate_tables,confint))
		thiscmhp="\t".join(thiscmhp)
		print "\t".join(entries[:4])+"\t"+thiscmhp
	inf.close()
