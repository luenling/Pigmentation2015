# takes a file with SNP pairs
# format:
# chrom pos1 pos2 freq1 freq2 r2
# takes a SNP file with Chrom pos freq recrate
# sample n snps -> look for pairs (matched in distance & freq & rec rate) until n reached
# get r2 for EU and SA for all these pairs

def calc_r2(c_AB,tot,fA,fB):
    if fA >= 0.999 or fB >= 0.999:
        return 1
    return (float(c_AB)/tot - fA*fB)**2/(fA*fB*(1-fA)*(1-fB))


def get_r2(span_nucs,snp1,snp2,min_int=10):
    nucs= ["A","C","G","T"]
    maj1=nucs[ np.argmax(snp1)]
    maj2=nucs[ np.argmax(snp2)]
    tot=float(sum(span_nucs.values()))
    r2=0.0
    if tot >= min_int:
        fAx=sum([ y for x,y in span_nucs.items() if x[0] == maj1])/tot
        fxB=sum([ y for x,y in span_nucs.items() if x[1] == maj2])/tot
        cAB=span_nucs[(maj1,maj2)]
        r2=calc_r2(cAB,tot,fAx,fxB)
    else:
        r2=np.nan
    return r2

def get_readspan_and_nucs(chrom,bps,bamfh,snp_span,snp1,snp2,map_qual=20,base_qual=20):
    """
    gets a chrom and a sorted SNP pair, and a bam file handle returns the counts of the each nucleotide-pair as a dict and the ACTG counts for each SNP
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101)
    for pile in bamfh.pileup(region=chrom,start=(bps[0]-1),end=bps[1]):
        if pile.pos == bps[0]-1:
            pile_col1=pile.pileups
        elif pile.pos == bps[1]-1:
            pile_col2=pile.pileups
            break
    #b=[x.alignment.seq[x.qpos] for x in pile_col ]
    #counts=[ b.count(x) for x in ["A","C","G","T"]]
    #read_dict=defaultdict(list)
    snp_nucs=[]
    read_dict=defaultdict()
    for pile_read in pile_col1:
        # check whether matched, mapping qual and base qual alright
        if  pile_read.is_del or (pile_read.alignment.mapq < map_qual) or (ord(pile_read.alignment.qual[pile_read.qpos])- 33) < base_qual :
            continue
        # assign read_read
        read_dict[pile_read.alignment.qname]=pile_read.alignment.seq[pile_read.qpos]
        snp_nucs.append(pile_read.alignment.seq[pile_read.qpos])
    snp1 +=np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)
    snp_nucs=[]
    for pile_read in pile_col2:
        # check whether matched, mapping qual and base qual alright
        if  pile_read.is_del or (pile_read.alignment.mapq < map_qual) or (ord(pile_read.alignment.qual[pile_read.qpos])- 33) < base_qual :
            continue
        # assign read_read
        if pile_read.alignment.qname in read_dict.keys():
            span_nucs[(read_dict[pile_read.alignment.qname],pile_read.alignment.seq[pile_read.qpos])] += 1 
        snp_nucs.append(pile_read.alignment.seq[pile_read.qpos])
    snp2 +=np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)
    bamfh.reset()
    return np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)


def read_cmh_file(infile,popsA,popsB,rr):
    """
    reads a cmh file and returns a dictionary of chromosomes with a key pointing at a numpy array of bp positions and one at the AFs for each chromsome key and recombination rate:
    eg cmhdict['2L']['bps'] = array[BPs],  cmhdict['2L']['af1'] = array[afs] and cmhdict['2L']['af2']=array[af2] cmhdict['2L']['rr']=array[recrate]
    and pVs array : [ [pVs,CHR,Afs,index] ]
    for the afs the average of the floats at the populations given in the argument pops is taken  
    """
    # intermediate lists of values
    cmh_dict= defaultdict(lambda: defaultdict(list)) #dictionary of chroms and different subfields defaulting to lists
    if re.search("\.b?gz$",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    popsA += 2
    #popsB = popsB.copy()
    popsB += 2
    rr += 2
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line.rstrip()
        fields=line.split()
        cmh_dict[fields[0]]['bps'].append(int(fields[1]))
        if popsA.ndim: # not just one column
            af1=sum([ float(fields[x]) for x in popsA ])/len(popsA)
        else:
            af1=float(fields[popsA])
        if popsB.ndim: # not just one column
            af2=sum([ float(fields[x]) for x in popsB ])/len(popsB)
        else:
            af2=float(fields[popsB])
        cmh_dict[fields[0]]['af1'].append( af1 )        
        cmh_dict[fields[0]]['af2'].append( af2 )
        try:
            recrate=float(fields[rr])
        except:
            recrate=0.0
        cmh_dict[fields[0]]['rr'].append( recrate )        
        #pVs.append([ pV,fields[0], afs])
    # convert to np arrays
    for i in cmh_dict.keys():
        for j in cmh_dict[i].keys():
            cmh_dict[i][j] = np.array(cmh_dict[i][j])
    for i in cmh_dict.keys():
        sort_vector = cmh_dict[i]['bps'].argsort()
        for j in cmh_dict[i].keys():
            cmh_dict[i][j] = cmh_dict[i][j][sort_vector]
    inf.close()
    return cmh_dict


def create_bin_indeces(snp_dict,bins,key):
    """
    creates dictionary of the bin index positions in the AF histogram for each chromosome with np arrays of the indeces of all SNPs in that frequency bin for sampling.
    """    
    for i in snp_dict.keys():
        dig=np.digitize(snp_dict[i][key],bins)
        # bin indeces are shifted +1 against histogram count indeces
        dig -= 1
        indx_bins=defaultdict(list)
        for j,x in enumerate(dig):
            indx_bins[x].append(j)
        for j in indx_bins.keys():
            indx_bins[j]=np.array(indx_bins[j])
        snp_dict[i][key+'_idx']=indx_bins
    return True

def get_idxs(snp_dict,chrom,afbins,rrbin,pop):
    # get list of indeces of snps with the right afs, rr, pop
    af_key="af"+pop+"_idx"
    af1_idxs=np.intersect1d(snp_dict[chrom][af_key][afbins[0]],snp_dict[chrom]['rr_idx'][rrbin])
    af2_idxs=np.intersect1d(snp_dict[chrom][af_key][afbins[1]],snp_dict[chrom]['rr_idx'][rrbin])
    return [af1_idxs,af2_idxs]

def get_rand_idxs(snp_dict,chrom,idxs,dist,n,v=False):
    """
    gets the snp_dict, chromosome, indeces that for the afs and rr for snp1 and snp2, the distance and the sample size and gives a list of pairs of indeces back   
    """
    # first pick one of idxs[0]
    db=[dist*0.75,dist*1.25]
    # get BPS[0]
    # calc snp_dict[chrom]['bps']
    # maximum number of unique indeces to be gotten
    ml=len(idxs[0])
    tries = 0
    idx1=set()
    pairs=[]
    while (tries <= 100) or len(idx1) < n:
        # first condition just to make sure it ends somehow
        idx1_test=set()
        num_left= n - len(idx1)
        while len(idx1) < ml and len(idx1_test) < num_left:
            # create set to test for potential snps to add
            add_set=set(np.random.random_integers(0,ml-1,size=num_left))
            idx1_test.update(add_set.difference(idx1))
        for idx in idx1_test:
            bps1=snp_dict[chrom]['bps'][idxs[0][idx]]
            bps2=snp_dict[chrom]['bps'][idxs[1]]
            bps2 = abs(bps2 - bps1)
            # check which idxs[1] in distance dist +/- 0.25 * dist
            idx2 = np.where( ( bps2 >= db[0]) & (bps2 <= db[1] ) )[0]
            if len(idx2) > 0:
                # snp pair to add:
                # take random from snp2 in right distance
                idx2 = np.random.choice(idx2,1)
                pairs.append(np.array([idxs[0][idx],idxs[1][idx2]]))
                idx1.add(idx)
            else:
                continue
            # get SNPs in distance
        tries += 1
        if v and (tries%10 == 0):
            print "tries:{}\tidxs:{}\tpairs:{}\t".format(tries,len(idx1),len(pairs))
    return pairs
    
import pysam
import copy
import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import argparse
import gzip
parser = argparse.ArgumentParser(description="reads a file with SNPs with allele frequencies (AF) and pValues and creates random samples of the same size, chromosome and AF distribution as the SNPs to test for correlation of logORs given inthe same file. If given three logORs per population, also consideres the CI for checking consistency. Distance +/- 50%")

parser.add_argument("-i","--in", dest="inA", help="file with pos, base afs and rec rate", required=True)
parser.add_argument("-t","--tar", dest="inB", help="file with SNPs to sample (format CHR SNP1 SNP2)", required=True)
parser.add_argument("--pA", dest="popsA",  help="populations to use for the average AFs of pop1  (default: \"4,5,6\")",default="4,5,6")
parser.add_argument("--pB", dest="popsB",  help="populations to use for the average AFs of popB  (default: \"1,2,3\")",default="1,2,3")
parser.add_argument("--rr", dest="rr",  help="population postion of the recombination rate (int)", required=True)
parser.add_argument("--perc", dest="perc",  help="binning for AFs (default: \"0.0,0.1,0.25,0.5,0.75,0.9,1\")",default="0.0,0.1,0.25,0.5,0.75,0.9,1")
parser.add_argument("-n", dest="n",  help="number of random samples to create for each snp pair", default="1000")
#parser.add_argument("-o", dest="outfile",  help="outputfile", required=True)
parser.add_argument("-r","--rec", dest="rec", help="list with thresholds for very low,low, medium, high and so on recombination columns def. \"0.01,1.68,3.21\" (default: False)", default="0.01,1.68,3.21")

# for haplotypes
parser.add_argument("--bA", dest="bam_filesA", help="bam files, list comma separated for pop A", required=True)
parser.add_argument("--bB", dest="bam_filesB", help="bam files for pop B, list comma separated", required=True)
parser.add_argument("--mq", dest="map_qual", help="mapping quality threshold", default=20)
parser.add_argument("--min_int", dest="min_int", help="minimal intersection depths for SNP combinations (default=10)", default=10)
parser.add_argument("--min_cov", dest="min_cov", help="minimal coverage for SNPs (int)", default=10)
parser.add_argument("--bq", dest="base_qual", help="base quality threshold", default=20)


sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)



args = parser.parse_args()
# for haplotyping
bam_filesA = vars(args)['bam_filesA'].split(",")
bam_filesB = vars(args)['bam_filesB'].split(",")
map_qual = vars(args)['map_qual']
base_qual = vars(args)['base_qual']
min_int = vars(args)['min_int']
min_cov = vars(args)['min_cov']
# for random sampling
n = int(vars(args)['n'])
inA = vars(args)['inA']
inB = vars(args)['inB']
rr = int(vars(args)['rr'])
perc = np.array(vars(args)['perc'].split(","), dtype=float )
popsA = np.array( vars(args)['popsA'].split(","),dtype=int )
popsB = np.array( vars(args)['popsB'].split(","),dtype=int )
reclevs=vars(args)['rec'].split(",")
reclevs=[0.0]+reclevs+["inf"]
# array of afs bins
#abins=np.arange(0,1.0000000001,perc)
abins=perc
# array of rec rate bins
rbins=np.array(reclevs,dtype=float)
# read target SNP pair file
tar_snps=[]
with open(inB) as f:
    for line in f:
        if re.match("\#",line) or re.match("\s+",line):
            continue
        line=line.rstrip()
        entries=line.split("\t")
        entry=defaultdict()
        entry['CHR']=entries[0]
        entry['BP1']=int(entries[1])
        entry['BP2']=int(entries[2])
        entry['dist']=abs(entry['BP1']-entry['BP2'])
        tar_snps.append(entry)
# read cmh file
snp_dict= read_cmh_file(inA,popsA,popsB,rr)

# add bins for af1
create_bin_indeces(snp_dict,abins,"af1")
# add bins for af2
create_bin_indeces(snp_dict,abins,"af2")
# add bins for rr
create_bin_indeces(snp_dict,rbins,"rr")
#snp_dict['rbins']=rbins
#snp_dict['abins']=abins

# get frequencies and recrate for comparison SNPs
for snp_pair in tar_snps:
    chrom=snp_pair['CHR']
    bp1=snp_pair['BP1']
    bp2=snp_pair['BP2']
    idx1=np.searchsorted(snp_dict[chrom]['bps'],bp1)-1
    idx2=np.searchsorted(snp_dict[chrom]['bps'],bp2)-1
    snp_pair['af1']=[ snp_dict[snp_pair['CHR']]['af1'][idx1],snp_dict[snp_pair['CHR']]['af1'][idx2]]
    snp_pair['af2']=[ snp_dict[snp_pair['CHR']]['af2'][idx1],snp_dict[snp_pair['CHR']]['af2'][idx2]]
    snp_pair['af1bin']=np.digitize(snp_pair['af1'],abins)-1
    snp_pair['af2bin']=np.digitize(snp_pair['af2'],abins)-1
    #snp_pair['afbin']=[ np.array(list(set([snp_pair['af1bin'][x],snp_pair['af2bin'][x]]))) for x in [0,1] ]
    snp_pair['avg1']=np.mean(snp_pair['af1'])
    snp_pair['avg2']=np.mean(snp_pair['af2'])
    snp_pair['rr']=(snp_dict[snp_pair['CHR']]['rr'][idx1] + snp_dict[snp_pair['CHR']]['rr'][idx2])/2.0
    snp_pair['rrbin']=np.digitize([ snp_pair['rr'] ],rbins)[0]-1

nucs= ["A","C","G","T"]
print "CHR\tBPS1\tBPS2\tdist\tr2A\tr2B\tnRand\tmedian_rA\tmean_rA\tstd_rA\tmedian_rB\tmean_rB\tstd_rB\twil_P\tq_rA\tq_rB\tq_rAB"
bam_fhsA = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_filesA]
bam_fhsB = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_filesB]

for snp_pair in tar_snps: 
    # get SNPs with the required chrom, afs and rr
    # create samples for af1
    idxs1=get_idxs(snp_dict,snp_pair['CHR'],snp_pair['af1bin'],snp_pair['rrbin'],"1")    
    idxs2=get_idxs(snp_dict,snp_pair['CHR'],snp_pair['af2bin'],snp_pair['rrbin'],"2")    
    idxs=[ np.unique(np.concatenate((idxs1[x],idxs2[x]))) for x in [0,1]]  
    rnd_idxs=get_rand_idxs(snp_dict,snp_pair['CHR'],idxs,snp_pair['dist'],n)
    rnd_idxs.sort(key=min)
    for x in rnd_idxs: x.sort()
    # calculate r2 values for
    # for each snp pair in popA:
    chrom=snp_pair['CHR']
    sbps=[ snp_pair['BP1'], snp_pair['BP2'] ]
    span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
    # span_nucsB=defaultdict(int); snp1B=np.zeros(4,dtype=int); snp2B=np.zeros(4,dtype=int)
    for fh in bam_fhsA:
        get_readspan_and_nucs(chrom,sbps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
        #fh.close()
    r2A=get_r2(span_nucs,snp1,snp2,min_int)
    span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
    for fh in bam_fhsB:
        get_readspan_and_nucs(chrom,sbps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
        #fh.close()
    r2B=get_r2(span_nucs,snp1,snp2,min_int)
    rand_rA=[]
    rand_rB=[]
    # bam_fhsA = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_listA]
    # bam_fhsB = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_listB]
    for bps in rnd_idxs:
        bps=snp_dict[chrom]['bps'][bps]
        span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
        for fh in bam_fhsA:
            get_readspan_and_nucs(chrom,bps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
        rand_rA.append(get_r2(span_nucs,snp1,snp2,min_int))
        span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
        for fh in bam_fhsB:
            get_readspan_and_nucs(chrom,bps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
        rand_rB.append(get_r2(span_nucs,snp1,snp2,min_int))
    rand_rA=np.array(rand_rA)
    rand_rB=np.array(rand_rB)
    inds = np.where(  ( np.isfinite(rand_rA)) & ( np.isfinite(rand_rB)) )
    rand_rA=np.array(rand_rA[inds])
    rand_rB=np.array(rand_rB[inds])
    wilcox=stats.wilcoxon(rand_rA,rand_rB)[1]
    m_rA=np.median(rand_rA)
    mean_rA=np.mean(rand_rA)
    std_rA=np.std(rand_rA)
    q_rA=float(len(rand_rA[rand_rA >= r2A]))/len(rand_rA)
    m_rB=np.median(rand_rB)
    mean_rB=np.mean(rand_rB)
    std_rB=np.std(rand_rB)
    q_rB=float(len(rand_rB[rand_rB >= r2B]))/len(rand_rB)
    q_rAB=float(len(rand_rB[ (rand_rA >= r2A)  | (rand_rB >= r2B) ]))/len(rand_rB)
    
    print "{}\t{}\t{}\t{}\t{:.4g}\t{:.4g}\t{}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}\t{:.4g}".format(chrom,sbps[0],sbps[1],snp_pair['dist'],r2A,r2B,len(rand_rA),m_rA,mean_rA,std_rA,m_rB,mean_rB,std_rB,wilcox,q_rA,q_rB,q_rAB)
    
    
for bamfh in  bam_fhsA:
    bamfh.close()
for bamfh in  bam_fhsB:
    bamfh.close()
    

# def crapola:
#     # create samples for af2
#     rnd_idxs2=get_rand_idxs(snp_dict,snp_pair['CHR'],idxs,dist,n)
#     rnd_idxs2.sort(key=min)
#     for x in rnd_idxs2: x.sort()
#     rand2_rA=[]
#     rand2_rB=[]
#     # bam_fhsA = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_listA]
#     # bam_fhsB = [ pysam.Samfile(bam_file,'rb') for bam_file in bam_listB]
#     for bps in rnd_idxs2:
#         bps=snp_dict[chrom]['bps'][bps]
#         span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
#         for fh in bam_fhsA:
#             get_readspan_and_nucs(chrom,bps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
#         rand2_rA.append(get_r2(span_nucs,snp1,snp2,min_int))
#         span_nucs=defaultdict(int); snp1=np.zeros(4,dtype=int); snp2=np.zeros(4,dtype=int)
#         for fh in bam_fhsB:
#             get_readspan_and_nucs(chrom,bps,fh,span_nucs,snp1,snp2,map_qual,base_qual)
#         rand2_rB.append(get_r2(span_nucs,snp1,snp2,min_int))
#     rand2_rA=np.array(rand2_rA)
#     rand2_rB=np.array(rand2_rB)
#     rand2_rA=np.array(rand2_rA[inds])
#     rand2_rB=np.array(rand2_rB[inds])
#     rand_all_A=np.array([ rand1_rA, rand2_rA])
#     rand_all_B=np.array([ rand1_rB, rand2_rB])
#     wil_2=stats.wilcoxon(rand1_rA,rand1_rB)[1]
#     wil_comb=stats.wilcoxon(rand_all_rA,rand_all_rB)[1]
    # inds = np.where(  ( np.isfinite(rand2_rA)) & ( np.isfinite(rand2_rB)) )
