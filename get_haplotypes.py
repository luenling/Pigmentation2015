import pysam
from distutils.version import LooseVersion
import numpy as np
import copy
import gzip
#from scipy import stats
import sys
import re
from collections import defaultdict
import argparse

def read_snps(infile):
    """
    reads a cmh or sync like file and creates a dict with chrom->[ [ bps, {} ] ]  
    """
    snp_dict=defaultdict(list)  #dictionary of chroms with positions and pVals of snps
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
        snp_dict[fields[0]].append([ int(fields[1]), {} ])
    inf.close()
    return snp_dict

def get_reads_and_nucs(chrom,bps,snp_indx,bamfh,read_dict,map_qual=20,base_qual=20):
    """
    gets a SNP position (chrom & bps), SNP indx, a bam file handle and a read_dict ( defaultdict(list)), adds the reads with snp index and nuc to the read_dict and returns an array with the ACTG counts
    read_dict[readname]->[ [snp indx,nucleotide] ]    
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101) 
    # have to go through loop, I think
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups
            break
    #b=[x.alignment.seq[x.qpos] for x in pile_col ]
    #counts=[ b.count(x) for x in ["A","C","G","T"]]
    #read_dict=defaultdict(list)
    snp_nucs=[]
    for pile_read in pile_col:
        # check whether matched, mapping qual and base qual alright
        if  pile_read.is_del or (pile_read.alignment.mapq < map_qual) or (ord(pile_read.alignment.qual[pile_read.qpos])- 33) < base_qual :
            continue
        # assign read_read
        read_dict[pile_read.alignment.qname].append([snp_indx,pile_read.alignment.seq[pile_read.qpos]])
        snp_nucs.append(pile_read.alignment.seq[pile_read.qpos])
    return np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)

def get_reads_and_nucs_new(chrom,bps,snp_indx,bamfh,read_dict,map_qual=20,base_qual=20):
    """
    function for the new pysam synthax (pysam 0.8.0 and up), thanks to the maintainers for breaking :( 
    gets a SNP position (chrom & bps), SNP indx, a bam file handle and a read_dict ( defaultdict(list)), adds the reads with snp index and nuc to the read_dict and returns an array with the ACTG counts
    read_dict[readname]->[ [snp indx,nucleotide] ]    
    """
    # create pileup and move iterator to the right position - problem due to pysam
    # if you use pysam 0.6, you will have to subtract a readlength and a bit from start (eg. start=bps-101) 
    # have to go through loop, I think
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups
            break
    #b=[x.alignment.seq[x.qpos] for x in pile_col ]
    #counts=[ b.count(x) for x in ["A","C","G","T"]]
    #read_dict=defaultdict(list)
    snp_nucs=[]
    for pile_read in pile_col:
        # check whether matched, mapping qual and base qual alright
        if  pile_read.is_del or (pile_read.alignment.mapping_quality < map_qual) or pile_read.alignment.query_qualities[pile_read.query_position] < base_qual :
            continue
        # assign read_read
        read_dict[pile_read.alignment.query_name].append([snp_indx,pile_read.alignment.query_sequence[pile_read.query_position]])
        snp_nucs.append(pile_read.alignment.query_sequence[pile_read.query_position])
    return np.array([ snp_nucs.count(x) for x in ["A","C","G","T"]],dtype=int)


def create_hap_dict(read_dict):
    """
    creates a dict of all haplotypes with the structure: (snp ids) -> ( alleles ) -> counts 
    """
    hap_dict=defaultdict(lambda: defaultdict(lambda: 0) )
    for read_id in read_dict:
        snps = sorted(read_dict[read_id],key=lambda x: x[0])
        snp_idxs=[x[0] for x in snps]
        snp_alleles=[x[1] for x in snps]
        if len(set(snp_idxs)) != len(snp_idxs):
            # if SNPs duplicated : read pair overlapping, check for consistency - if not make X
            keep_idx=[0]
            for i in range(1,len(snps)):
                if snp_idxs[i] == snp_idxs[i-1]:
                    if snp_alleles[i] != snp_alleles[i-1]:
                        snp_alleles[i-1] == 'N'
                else:
                    keep_idx.append(i)
            snp_idxs=[snp_idxs[x] for x in keep_idx]
            snp_alleles=[snp_alleles[x] for x in keep_idx]
        if len(snp_idxs) < 2:
            continue
        hap_dict[tuple(snp_idxs)][tuple(snp_alleles)] +=1
    return hap_dict
    
def add_read_counts(set_a,set_b,hap_dict):
    """
    takes 2 sorted snp idx tuples and halplotype dicts (snp_idx_tuple -> allele tuples -> count ) with tdx_a being a subset of b and adds the counts of b to the appropriate a haplotypes
    eg. : (0,1) -> { (A,C):1, (G,C):2, (A,T):1 }, (0,1,2) -> {(A,C,T):1, (A,C,G):2, (A,T,T):1} => (0,1) -> { (A,C):1+3, (G,C):2, (A,T):1+1 }
    """
    # get subset positions:
    isec=tuple(sorted(set(set_a).intersection(set_b)))    
    if isec != set_a:
        return False
    isec_in_b =  [set_b.index(x) for x in isec]
    # go through a and add counts to hap_dict[isec][alleles]
    for i in hap_dict[set_b]:
        hap_dict[set_a][tuple([i[x] for x in isec_in_b])] += hap_dict[set_b][i]

def calc_r2(c_AB,tot,fA,fB):
    return (float(c_AB)/tot - fA*fB)**2/(fA*fB*(1-fA)*(1-fB))

################
###   MAIN   ###
################

### for testing try this bamfile
### bam_file="/Volumes/vetgrid10/Data/Vienna_2010/Realigned/pop1_earlylate_merged_RG_sanger_real.bam"
### chrom="2L"
### bps=40829


parser = argparse.ArgumentParser(description="""
gets a list of bam files, names for the pools in the bam files, and a sync-like file with snps (only chromosome and basepos needed) and prints a file with the following for all pairs (or even ) of SNPs that share reads:
...
""") 
parser.add_argument("-b", dest="bam_files", help="bam files, list comma separated, multiple files per pool seperated by colon (eg. \"dark1bam:dark2.bam:dark3,light1.bam:light2.bam,base.bam\")", required=True)
parser.add_argument("-p", dest="pool_names", help="names for pools, list comma separated", required=True)
parser.add_argument("-s", dest="snp_file", help="snp file", required=True)
parser.add_argument("--mq", dest="map_qual", help="mapping quality threshold", default=20)
parser.add_argument("--min_int", dest="min_int", help="minimal intersection depths for SNP combinations (int)", default=10)
parser.add_argument("--min_cov", dest="min_cov", help="minimal coverage for SNPs (int)", default=10)
parser.add_argument("--pairs", dest="pairs", action="store_true", help="print pairs of SNPs only (default: print all, only works with full out)", default=False)
parser.add_argument("--bq", dest="base_qual", help="base quality threshold", default=20)
parser.add_argument("-f", dest="full_out", action="store_true", help="create hard to parse, long winded output format", default=False)
#parser.add_argument("-s", dest="snp_file", help="snp file", required=True)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbose state on stderr", default=False)


args = parser.parse_args()
snp_file = vars(args)['snp_file']
bam_files = vars(args)['bam_files'].split(",")
map_qual = int(vars(args)['map_qual'])
base_qual = int(vars(args)['base_qual'])
min_int = int(vars(args)['min_int'])
min_cov = int(vars(args)['min_cov'])
pairs = vars(args)['pairs']
full_out = vars(args)['full_out']
pool_names = vars(args)['pool_names'].split(",")
pool_haps = {}
#snp_file = '/Volumes/vetgrid10/Data/Vienna_2010/Realigned/Haplos/test.vcf'
#bam_files = ['/Volumes/vetgrid10/Data/Vienna_2010/Realigned/pop1_earlylate_merged_RG_sanger_real.bam','/Volumes/vetgrid10/Data/Vienna_2010/Realigned/pop2_earlylate_merged_RG_sanger_real.bam']
#map_qual = 20
#base_qual = 20
#min_int = 10
#min_cov = 10
#pairs = True
#full_out = False
#pool_names = ["t1","t2"]
#pool_haps = {}

if not full_out:
    # set pairs to True for special output
    pairs = True


#bam_files="/Volumes/Temp/Lukas/Data/SA_A7/BGI_100_101/BGI101_sanger_A7_VL_r5.picnodupl.filtered.mq20_chrfilt_RG_real.bam,/Volumes/Temp/Lukas/Data/SA_A7/BGI_100_101/BGI101_sanger_A7_VVD_r4.picnodupl.filtered.mq20_chrfilt_RG_real.bam,/Volumes/Temp/Lukas/Data/SA_A7/BGI_100_101/BGI101_sanger_A7_L_r5.picnodupl.filtered.mq20_chrfilt_RG_real.bam"
#pool_names="A7_VL_r5,A7_VVD_r4,A7_L_r5"
# check pysam version as synthax changed
if LooseVersion(pysam.__version__) >= "0.8.1":
    pysam_new = True
else:
    pysam_new=False
if len(pool_names) != len(bam_files):
    sys.exit("Length pool names and bam files not equal")
if not full_out:
    # print header line
    fields=["af1","af2","hfr","r2","N"]
    print "chr\tbps1\tbps2\ta1\ta2\t"+"\t".join([ x+"_"+y for x in pool_names for y in fields ])
snp_dict=read_snps(snp_file)
for chrom in sorted(snp_dict.keys()):
    for bam_idx,bam_list in enumerate(bam_files):
        read_dict=defaultdict(list)
        for bam_file in bam_list.split(":"):
            if pysam_new:
                bamfh=pysam.AlignmentFile(bam_file,'rb')
            else:
                bamfh=pysam.Samfile(bam_file,'rb')
            # get dictionary of reads with snp idx and alleles
            for snp_idx in range(0,len(snp_dict[chrom])):
                if not(pool_names[bam_idx] in snp_dict[chrom][snp_idx][1].keys()):
                    snp_dict[chrom][snp_idx][1][pool_names[bam_idx]] = np.zeros(4,dtype=int)
                if pysam_new:
                    snp_dict[chrom][snp_idx][1][pool_names[bam_idx]]+=get_reads_and_nucs_new(chrom,snp_dict[chrom][snp_idx][0],snp_idx,bamfh,read_dict,map_qual=map_qual,base_qual=base_qual)
                else:
                    snp_dict[chrom][snp_idx][1][pool_names[bam_idx]]+=get_reads_and_nucs(chrom,snp_dict[chrom][snp_idx][0],snp_idx,bamfh,read_dict,map_qual=map_qual,base_qual=base_qual)
            bamfh.close()
        # transform into a dictionary of different snp combinations occuring
        hap_dict=create_hap_dict(read_dict)
        # get all snp combinations found on reads/readpairs
        snp_sets=sorted(hap_dict.keys(),key=lambda x: (len(x),x))
        isec_sets=set([])
        # run through to see which intersection sets need to be added to get all subsets
        for i in range(0,len(snp_sets)-1):
            if len(snp_sets[i]) > 2: # add all tuples that are subsets, just to have all snp pairs for r2 calculations
                for sub_pair in [ (snp_sets[i][x],snp_sets[i][y]) for x in range(0,len(snp_sets[i])-1) for y in range(x+1,len(snp_sets[i]))]:
                    isec_sets.add(sub_pair)
                    # look at all combinations with later starting/longer sets
            for j in range(i+1,len(snp_sets)):
                # get intersections between sets with len > 2 that are not already in the set
                if len(snp_sets[j]) < 3 or len(set(snp_sets[i]).intersection(snp_sets[j])) < 2:
                    # no possible overlap > 1 or if overlap is < 2
                    continue
                if len(set(snp_sets[i]).intersection(snp_sets[j])) != len(snp_sets[i]):
                    if tuple(sorted(set(snp_sets[i]).intersection(snp_sets[j]))) in snp_sets:
                        # set is already in hap_dict
                        continue
                    else:
                        isec_sets.add(tuple(sorted(set(snp_sets[i]).intersection(snp_sets[j]))))
                        continue
        # add intersection sets to set of snp index sets, unique and sort
        snp_sets.extend(isec_sets)
        snp_sets=sorted(set(snp_sets),key=lambda x: (len(x),x))
        # go through all snp idx sets (sorted by length) and add counts for longer supersets
        for i in range(0,len(snp_sets)-1):
            # look at all combinations with later starting/longer sets
            for j in range(i+1,len(snp_sets)):
                if not ( snp_sets[j] in hap_dict.keys() ) or len(set(snp_sets[i]).intersection(snp_sets[j])) != len(snp_sets[i]):
                    # empty second set, or overlap != intersection set
                    continue
                add_read_counts(snp_sets[i],snp_sets[j],hap_dict)
        pool_haps[pool_names[bam_idx]]=copy.deepcopy(hap_dict)
    #############
    # output part
    #############
    if full_out:
        # print parameters
        print "Parameters used"
        print "Mapping quality:"+str(map_qual)+"\tBase Qual:"+str(base_qual)
        print "Minimal Intersection count:"+str(min_int)+"\tMin. coverage:"+str(min_cov)
        # print bam files:
        print "Bam files used for pools:"
        for bam_idx,bam_list in enumerate(bam_files):
            print pool_names[bam_idx]+":\t"+",".join([ bam_file for bam_file in bam_list.split(":") ])+"\n"
            # print SNPs (number, Chrom, pos, ....)
        print "SNP File:"
        if re.search("\.b?gz",snp_file):
            snp_inf = gzip.open(snp_file,'rb')
        else:
            snp_inf = open(snp_file,"r")
        line_num=1
        for line in snp_inf:        
            print str(line_num)+"\t"+line.rstrip()
            line_num += 1
        snp_inf.close()
        print "SNPs considered"
        print " \t \t \t \t \t "+" \t \t \t \t".join(pool_names)
        print "Number\tChrom\tPos"+"\tDepth\tA\tC\tG\tT"*len(pool_names)
    snp_depth=defaultdict(list)
    for i,j in enumerate(snp_dict[chrom]):
        for pool in pool_names:
            snp_depth[pool].append(sum(j[1][pool]))
            if full_out:
                print "{}\t{}\t{}\t".format(i+1,chrom,j[0])+"\t".join([ "{}\t{}".format(snp_depth[pool][i],"\t".join([str(x) for x in j[1][pool]])) for pool in pool_names ] )
    if full_out:
        # print Haplotypes
        print " "
        print "Haplotypes contig: "+chrom
        print "Cnt: count of intersecting reads supporting haplotype,Tot: total count of reads stretching SNPs, Frc: fraction Cnt/Tot, r2: r2 calculated using Method 1 from LDx article (Feder et al. 2012),afs: allele frequencies for each SNP, mDp: mean coverage of SNPs in haplotype"
        # print " "*(len(snp_dict[chrom])+ 2)+" \t  \t  \t     \t    \t   \t".join(pool_names)
        # print "".join(["-" for x in range(0,len(snp_dict[chrom]))])+"\tCnt\tTot\tFrc\tr2\tafs\tmDp"*len(pool_names)
        print " "*(len(snp_dict[chrom])+ 2)+" \t  \t     \t   \t".join(pool_names)
        print "".join(["-" for x in range(0,len(snp_dict[chrom]))])+"\tCnt\tFrc\tr2\tafs"*len(pool_names)
    snp_sets=[]
    for pool in pool_names:
        snp_sets.extend(pool_haps[pool].keys())
    snp_sets=sorted(set(snp_sets),key=lambda x: x)
    nucs = dict(zip(["A","C","G","T"],[0,1,2,3]))
    for snp_set in snp_sets:
        if pairs and len(snp_set) > 2:
            # only output pairs
            continue
        # all haplotypes
        haplo_sets=[]
        for pool in pool_names:
            haplo_sets.extend(pool_haps[pool][snp_set].keys())
        haplo_sets=sorted(set(haplo_sets),key=lambda x: x)
        # now a call to any combination of snp_set and haplo should create a new entry with count=0, I love defaultdict ;)
        tot={}
        for pool in pool_names:
            tot[pool]=sum([ pool_haps[pool][snp_set][x] for x in haplo_sets])
        #tot=sum([ hap_dict[snp_set][x] for x in hap_dict[snp_set].keys() ])
        if min_int > max([ tot[x] for x in pool_names ]):
            # intersection count too low in all pools
            continue
        if max([ snp_depth[y][x] for x in snp_set for y in pool_names]) < min_cov:
            # snp depths too low in all pools
            continue
        if full_out:
            print "###SNPs: "+",".join([str(x+1) for x in snp_set])+"\t"+chrom+":"+",".join([ str(snp_dict[chrom][x][0]) for x in snp_set])+"\tdistance:"+str(snp_dict[chrom][snp_set[-1]][0]-snp_dict[chrom][snp_set[0]][0])

        # get depth of reads averaged over snps in set
        avg_depth={ pool:sum([ snp_depth[pool][x] for x in snp_set])/len(snp_set) for pool in pool_names }
        # sort haplo sets descending by average number of reads found for each combination 
        haplo_sets.sort(key= lambda x: -1*np.mean([ pool_haps[pool][snp_set][x] for pool in pool_names]) )
        for haplo in haplo_sets:
            if full_out:
                # create string of dashes (one dash per snp)
                hap_str=["-"]*len(snp_dict[chrom])
            else:
                # create first head of string: chrom, bps1, bps 2, allele 1, allele 2
                hap_str=chrom+"\t"+"\t".join([ str(snp_dict[chrom][x][0]) for x in snp_set]) +"\t"+ "\t".join(haplo)
            for i in range(0,len(snp_set)):
                # go through every haplotype to create a string
                if full_out:
                    hap_str[snp_set[i]]=haplo[i]
                    hap_str="".join(hap_str)
            for pool in pool_names:
                if tot[pool] < min_int or min([ snp_depth[pool][x] for x in snp_set]) < min_cov:
                    if full_out:
                        hap_str+="\tNA"*4
                    else:
                        hap_str+="\tNA"*5
                    continue
                if len(snp_set) == 2 and pool_haps[pool][snp_set][haplo] > 0:
                    # calculate r2
                    c_AB=pool_haps[pool][snp_set][haplo]
                    fA=float(sum([ pool_haps[pool][snp_set][x] for x in pool_haps[pool][snp_set] if x[0] == haplo[0] ]))/tot[pool]
                    fB=float(sum([ pool_haps[pool][snp_set][x] for x in pool_haps[pool][snp_set] if x[1] == haplo[1] ]))/tot[pool]
                    r2= 0 if (fA >= 0.99999 or fB >= 0.99999 ) else calc_r2(c_AB,tot[pool],fA,fB)
                    r2="{:.3f}".format(r2)
                else:
                    r2="NA"
                allelefs=[]
                for i in range(0,len(snp_set)):
                    allelefs.append( float(snp_dict[chrom][snp_set[i]][1][pool][nucs[haplo[i]]])/snp_depth[pool][snp_set[i]])
                #hap_str += "\t{}\t{}\t{:.2f}\t{}\t{}\t{:.1f}".format(pool_haps[pool][snp_set][haplo],tot[pool],float(pool_haps[pool][snp_set][haplo])/tot[pool],r2,",".join(["{:.2f}".format(x) for x in allelefs]),avg_depth[pool])
                if full_out:
                    hap_str += "\t{}\t{:.2f}\t{}\t{}".format(pool_haps[pool][snp_set][haplo],float(pool_haps[pool][snp_set][haplo])/tot[pool],r2,",".join(["{:.2f}".format(x) for x in allelefs]))
                else:
                    hap_str += "\t{}\t{:.2f}\t{}\t{}".format("\t".join(["{:.2f}".format(x) for x in allelefs]),float(pool_haps[pool][snp_set][haplo])/tot[pool],r2,tot[pool])
            print hap_str
        #spacer bewtween SNP sets



    
              
                



