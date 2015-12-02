import sys, os, re
import numpy as np
import argparse
import gzip

#from scipy import stats
#Author: Lukas Endler

def split_syncline(syncline):
    if syncline:
        syncline = syncline.rstrip()
        entries=syncline.split()
    else:
        entries = False
    return entries

def get_next_line(inf):
    try:
        l=inf.next()
    except:
        l=False
    return l

def catchup_chrom(inf,chrom):
    entries=[]
    while(True):
        l=get_next_line(inf)
        if (not l):
            entries=l
            break
        entries=split_syncline(l)		
        if (entries[0] >= chrom):
            break
    return entries

def catchup_line_num(inf,line_num,chrom):
    entries=[]
    while(True):
        l=get_next_line(inf)
        if (not l):
            entries=l
            break
        entries=split_syncline(l)
        if (int(entries[1]) >= line_num) or (entries[0] != chrom):
            break
    return entries

parser = argparse.ArgumentParser(description='read an allele frequency file in Sync Format (CHR BPS Maj/Minor Alleles) and one or more files containing the mapping of Dmel positions to Simulans as provided by Ram\'s Mauve scripts, and tell for each position whether the simulans allele corresponds to the major (1), or minor (2) allele, whether it is not matching at all (0), or whether the position does not map between the species genomes (-1)') 
parser.add_argument("-a","--af", dest="affile", help="allele frequency file, needs to be sorted the same way as mapfile", required=True)
parser.add_argument("-m","--maps", dest="mapfiles", help="comma seperated list of files with D. mel - D. sim. genome position mappings", required=True)
parser.add_argument("--oa", action="store_true", dest="outall", help="output the allele, if mapped, not a digit code", default=False)


args = parser.parse_args()
affile = vars(args)['affile']
mapfiles = vars(args)['mapfiles'].split(",")
outall = vars(args)['outall']

if re.search("\.b?gz",affile):
    af = gzip.open(affile,"rb")
else:
    af = open(affile,"r")
mfh =[ gzip.open(fn,"rb") if re.search("\.b?gz",fn) else open(fn,"r") for fn in mapfiles ]
# read first lines of all entries
mf_entries = [ split_syncline(get_next_line(fh)) for fh in mfh]
for line in af:
    (chrom, bps, alleles) = line.split()[0:3]
    bps=int(bps)
    #print  "{}\t{}\t{}\t{}".format(chrom, bps, alleles,mf_entries)
    all_cons = [] # allele conservation (-1: not mapping, 0: no allele, 1: major is ancestral, 2: minor is ancestral)
    for ind in range(0,len(mfh)):
        # catch up if necessary
        if (mf_entries[ind] and mf_entries[ind][0] < chrom):
            # chrom lagging behind
            mf_entries[ind]=catchup_chrom(mfh[ind],chrom)
        if (mf_entries[ind] and (mf_entries[ind][0] == chrom) and  (int(mf_entries[ind][1]) < bps) ):
            # bps lagging behind
            mf_entries[ind]=catchup_line_num(mfh[ind],bps,chrom) 
        if not(mf_entries[ind]) or mf_entries[ind][-1] == "-":
            all_cons.append("-1")
        else:
            all = alleles.find(mf_entries[ind][-1])+1
            if (all and outall):
                all=mf_entries[ind][-1]
            all_cons.append(str(all))
    print "{}\t{}\t{}\t{}".format(chrom, bps, alleles,"\t".join(all_cons))

af.close()
for fh in mfh:
    fh.close()
