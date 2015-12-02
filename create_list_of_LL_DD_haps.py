# reads a light allele sync file and a haplotype file and gives teh proportion of light light and dark dark haplotypes 

def read_snp_file(infile,pop):
    """
    reads a light allele file and returns a dict with chrom->bps->[ light, dark]
    """
    snp_dict=defaultdict(defaultdict)
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line=line.rstrip()
        fields=line.split()
        light= fields[pop]
        if light != "-":
            snp_dict[fields[0]][fields[1]]=[ light, fields[2][1-fields[2].index(light)]]
    inf.close()
    return snp_dict


import sys, os, re
import numpy as np
import gzip

#from scipy import stats
from collections import defaultdict
#np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
import argparse
parser = argparse.ArgumentParser(description="""
 reads a light allele sync file and a haplotype file and gives the proportion of light light, dark dark, other haplotypes and r2
output:
chrom snp1 snp2 [popn: tot_count frctLL fDD fOther r2] total_count tot_fLL tot_fDD tot_fOther nan 
"""
                                 )

parser.add_argument("--in","-i", dest="infile", help="file with haplotypes", required=True)
parser.add_argument("--ld","-l", dest="ldfile", help="file with light alleles", required=True)
parser.add_argument("--pos","-p", dest="pop", help="population in light allele file", required=True)
parser.add_argument("--mc","-m", dest="min_count", help="minimal count for a snp pair default=25", default=25)

args = parser.parse_args()
infile = vars(args)['infile']
ldfile = vars(args)['ldfile']
pop = int(vars(args)['pop'])
min_count = int(vars(args)['min_count'])

pop=pop+2
snp_dict=read_snp_file(ldfile,pop)

if re.search("\.b?gz$",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
snp=False
chrom=False
snp2=False
snp_key=False
hap_dict=defaultdict(lambda: defaultdict(list)) # dictionary of chrom->snp1-snp2_>[l1l2,d1d2,[ [ cLL1,cll2,cll3 ... ],...,[ cdd1, cdd2, .. ],[cO1,cO2,..] ] ]
r2_dict=defaultdict(lambda: defaultdict(list)) # dictionary of chrom->snp1-snp2_>[r21, ...]
for line in inf:
    # skip lines until start of snp found
    if re.match("\#\#\#SNPs:",line):
        try: # gives the right pattern to split
            (chrom,snp,snp2) = re.split("[:,]",line.split()[2])
        except: # wrong number of snps or something else fishy - skip bloody snp
            snp=False
        if snp and not(snp in snp_dict[chrom].keys() and snp2 in snp_dict[chrom].keys()):
            snp = False
        if snp: # it is a good old useable snp
            snp_key=snp+"-"+snp2
            ll=snp_dict[chrom][snp][0]+snp_dict[chrom][snp2][0]
            dd=snp_dict[chrom][snp][1]+snp_dict[chrom][snp2][1]
            hap_dict[chrom][snp_key].extend([ll,dd])
        continue
    if snp and re.match("^-*[ACTG]+[ACTG-]*\s+",line):
        fields=line.split()
        pops= (len(fields)-1)/4
        if not snp_key in r2_dict[chrom].keys():
            r2_dict[chrom][snp_key]=[ float(fields[3+x*4].replace("-","nan")) for x in range(0,pops) ]
        if len(hap_dict[chrom][snp_key]) == 2:
            # add fields with zeros
            hap_dict[chrom][snp_key].append(np.zeros((3,pops),dtype=int))
        haplo=fields[0].replace("-","")
        try:
            counts = np.array([ int(fields[1+x*4].replace("-","0")) for x in range(0,pops)],dtype=int)
        except:
            sys.exit("problem reading this:\n"+line)
        if haplo == hap_dict[chrom][snp_key][0]:
            # LL
            hap_dict[chrom][snp_key][2][0] += counts
        elif haplo == hap_dict[chrom][snp_key][1]:
            # DD
            hap_dict[chrom][snp_key][2][1] += counts
        else:
            hap_dict[chrom][snp_key][2][2] += counts
    elif snp:
        snp=False
    # continue

for chrom in sorted(hap_dict.keys()):
    for snp_key in sorted(hap_dict[chrom].keys(),key=lambda x: int(x.split("-")[0])):
        snps=snp_key.split("-")
        # get totals for ll, dd, o   
        tots=np.sum(hap_dict[chrom][snp_key][2],axis=1)
        # get total sum
        total=np.sum(hap_dict[chrom][snp_key][2])
        # get totals for each pop
        pop_tots=np.sum(hap_dict[chrom][snp_key][2],axis=0)
        snp_str=chrom+"\t"+snps[0]+"\t"+snps[1]
        if total <= min_count*2:
            continue
        for i in range(0,len(hap_dict[chrom][snp_key][2][0])):
            if pop_tots[i] <= min_count:
                snp_str += "\tnan"*5
            else:
                fracts=hap_dict[chrom][snp_key][2][:,i]/float(pop_tots[i])
                snp_str += "\t"+str(pop_tots[i])+"\t"+"\t".join([ "{:.2f}".format(x) for x in fracts])+"\t"+ "{:.2f}".format(r2_dict[chrom][snp_key][i])
        fracts=tots/float(total)
        snp_str +=  "\t"+str(total)+"\t"+"\t".join([ "{:.2f}".format(x) for x in fracts])+"\t"+"{:.2f}".format(np.median(r2_dict[chrom][snp_key][i]))
        print snp_str
                
                
                       
