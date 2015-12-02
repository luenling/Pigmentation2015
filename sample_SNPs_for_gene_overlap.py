def read_cmh_file(infile,pops=[ 1 ]):
    """
    reads a cmh file and returns a dictionary of chromosomes with a key pointing at a numpy array of bp positions and one at the AFs for each chromsome key:
    eg cmhdict['2L']['bps'] = array[BPs],  cmhdict['2L']['afs'] = array[pV]
    for the afs the average of the floats at the populations given in the argument pops is taken  
    """
    # intermediate lists of values
    cmh_dict= defaultdict(lambda: defaultdict(list)) #dictionary of chroms and different subfields defaulting to lists
    inf = open(infile,"r")
    pops = pops.copy()
    pops += 2
    n_pops = len(pops)
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line.rstrip()
        fields=line.split()
        cmh_dict[fields[0]]['bps'].append(int(fields[1]))
        cmh_dict[fields[0]]['afs'].append( sum([ float(fields[x]) for x in pops ])/n_pops  )
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

def add_genesets(snp_dict,gene_file):
    """
    reads a sync file with a space separated list of fbgn numbers for each position and adds them to a snp dictionary
    """
    inf = open(gene_file,"r")
    for i in snp_dict.keys():
        snp_dict[i]['genes']=np.empty(len(snp_dict[i]['bps']), dtype=set)
    for line in inf:
        if re.match("\#",line):
            continue
        line.rstrip()
        fields=line.split()
        if len(fields) < 3:
            continue
        bps=int(fields[1])
        if fields[0] in snp_dict.keys():
            idx = snp_dict[fields[0]]['bps'].searchsorted(bps)
            if (idx < len(snp_dict[fields[0]]['bps'])) and snp_dict[fields[0]]['bps'][idx] == bps:
                snp_dict[fields[0]]['genes'][idx]=set([ x for x in fields[2:] ])
    return True



def create_af_histo(snp_dict,perc):
    """
    creates a new entry for each chromosome holding an array with counts for each afs bin
    50.0/perc bins are created
    """
    for i in  snp_dict.keys():
        hist=np.histogram(snp_dict[i]['afs'], bins=np.ceil(100.0/perc), range=(0.0,1.0))
        # get rid of bins with afs == 0 in the spectrum
        vals=[]
        bins=[]
        for j,x in enumerate(hist[0]):
            if x > 0:
                if len(vals) == 0:
                    bins.append(hist[1][j])
                elif (bins[-1] != hist[1][j] ):
                    vals.append(0)
                    bins.append(hist[1][j])
                vals.append(x)
                bins.append(hist[1][j+1])
        snp_dict[i]['afs_hist'] = { 'counts':np.array(vals),'bins':np.array(bins) }
    return True

def create_bin_indeces(snp_dict,sig_snp):
    """
    creates dictionary of the bin index positions in the AF histogram for each chromosome with np arrays of the indeces of all SNPs in that frequency bin for sampling.
    """
    for i in sig_snp.keys():
        dig=np.digitize(snp_dict[i]['afs'],sig_snp[i]['afs_hist']['bins'])
        # bin indeces are shifted +1 against histogram count indeces
        dig -= 1
        indx_bins=defaultdict(list)
        for j,x in enumerate(dig):
            indx_bins[x].append(j)
        for j in indx_bins.keys():
            indx_bins[j]=np.array(indx_bins[j])
        snp_dict[i]['bin_idx']=indx_bins
    return True
            
def create_random_sample(idx_bins,count_array):
    """
    uses np.random.choice to sample SNP indeces from the bin index structure  according to a , around a 100 times slower than self created way. the replace=False creates the performance hit
    """
    idxs=[]
    for i,x in enumerate(count_array):
        if x > 0:
            idxs.extend(np.random.choice(idx_bins[i],size=x,replace=False))
    return idxs

def create_random_sample_alt(idx_bins,count_array):
    """
    uses unique_sample_of_int to sample SNP indeces from the bin index structure  according to a
    """
    idxs=[]
    for i,x in enumerate(count_array):
        if x > 0:
            idxs.extend([ idx_bins[i][ind] for ind in unique_sample_of_int(len(idx_bins[i])-1,x) ] )
    return idxs

def unique_sample_of_int(max,size):
    """
    samples a set of unique integers with length size from the interval \[0,max\]
    """
    idxs=set()
    num_left = size - len(idxs)
    while num_left > 0:
        idxs = idxs.union(set(np.random.random_integers(0,max,size=num_left)))
        num_left = size - len(idxs)
    return idxs

def get_random_geneIDs(snp_dict,sig_snp):
    geneIDs=[]
    for chrom in sig_snp.keys():
        idxs = create_random_sample_alt(snp_dict[chrom]['bin_idx'],sig_snp[chrom]['afs_hist']['counts'])
        for idx in idxs:
            if snp_dict[chrom]['genes'][idx]:
                geneIDs.extend(snp_dict[chrom]['genes'][idx])
    return set(geneIDs)

def get_gene_list(infile):
    inf = open(infile,"r")
    geneIDs=[]
    for line in inf:
        if re.match("\#",line):
            continue
        line.rstrip()
        fields=line.split()
        if len(fields) > 2:
            geneIDs.extend([ x for x in fields[2:] ])        
    inf.close()
    return set(geneIDs)

def get_sig_gene_set(snp_dict,sig_snp):
    """
    returns set of all geneIDs linked to the SNPs of sig_snp in snp_dict    
    """
    geneIDs=[]
    for chrom in sig_snp.keys():
        for bps in sig_snp[chrom]['bps']:
            idx = snp_dict[chrom]['bps'].searchsorted(bps)
            if (idx < len(snp_dict[chrom]['bps'])) and snp_dict[chrom]['bps'][idx] == bps and snp_dict[chrom]['genes'][idx]:
                geneIDs.extend(snp_dict[chrom]['genes'][idx])
    return set(geneIDs)

import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description="reads a file with SNPs with allele frequencies (AF) and annotations of gene names for each SNP, and creates random samples of the same size, chromosome and AF distribution as the SNPs from another file - eg. the significant SNPs. It gives a histogram of overlapping genes found in the samples.\n It can either be used with two SNP files A and B, or with one SNP file A and a list of target genes. To calculate a P Value the observed number of overlapping genes has to be given.")

parser.add_argument("--sA", dest="inA", help="syncfile A with all SNPs and AFs", required=True)
parser.add_argument("--aA", dest="aA", help="syncfile A SNPs and gene ID annotations", required=True)
parser.add_argument("--sigA", dest="sigA", help="syncfile of significant SNPs of A", required=True)
parser.add_argument("--sB", dest="inB", help="syncfile B with all SNPs and AFs", default=False)
parser.add_argument("--aB", dest="aB", help="syncfile B SNPs and gene ID annotations", default=False)
parser.add_argument("--sigB", dest="sigB", help="syncfile of significant SNPs of B", default=False)
parser.add_argument("-t", dest="target",  help="list of target genes for overlap, not used for A-B sampling",default=False)
parser.add_argument("--pops","-p", dest="pops",  help="populations in A and B to use for the average AFs (default: \"4,5,6\")",default="4,5,6")
parser.add_argument("--perc", dest="perc",  help="percent step for AF binning (default: 1.0)",default=1.0)
parser.add_argument("-n", dest="rnd_samp",  help="number of random samples to create", default="10000")
parser.add_argument("-o", dest="outfile",  help="outputfile", required=True)
parser.add_argument("-v", dest="verb",  action="store_true", help="outputfile", default=False)
parser.add_argument("--plot", dest="plot",  action="store_true", help="create a plot of the drawn values (default: False)", default=False)
parser.add_argument("--obs", dest="obs", help="observed number of overlapping genes (default: False)", default=False)

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


args = parser.parse_args()
inA = vars(args)['inA']
aA = vars(args)['aA']
sigA = vars(args)['sigA']
inB = vars(args)['inB']
aB = vars(args)['aB']
sigB = vars(args)['sigB']
target = vars(args)['target']
pops = np.array( vars(args)['pops'].split(","),dtype=int )
rnd_samp = int(vars(args)['rnd_samp'])
perc = float(vars(args)['perc'])
outfile = vars(args)['outfile']
verb= vars(args)['verb']
plot= vars(args)['plot']
obs= vars(args)['obs']
if obs:
    obs=int(obs)
# create SNP dict from file
if verb:
    print "create SNP dict from file "+ inA
snp_dictA= read_cmh_file(inA,pops)
if verb:
    print "adding gene sets from "+ aA
add_genesets(snp_dictA,aA)
# read signifcant snps
if verb:
    print "reading sign genes "+ sigA
sig_snpsA = read_cmh_file(sigA,pops)
# create histogram
if verb:
    print "creating AF histogram and index bins"
create_af_histo(sig_snpsA,perc)
# create index bins
create_bin_indeces(snp_dictA,sig_snpsA)
num_sig_genesA=get_sig_gene_set(snp_dictA,sig_snpsA)
# second snp file given:
if inB:
    if ( not aB ) or ( not sigB ):
        sys.exit("not all files required for B given")
    if verb:
        print "create SNP dict from file "+ inB
    snp_dictB= read_cmh_file(inB,pops)
    if verb:
        print "adding gene sets from "+ aB
    add_genesets(snp_dictB,aB)
    # read signifcant snps
    if verb:
        print "reading sign genes "+ sigB
    sig_snpsB = read_cmh_file(sigB,pops)
    # create histogram
    if verb:
        print "creating AF histogram and index bins"
    create_af_histo(sig_snpsB,perc)
    # create index bins
    create_bin_indeces(snp_dictB,sig_snpsB)
    num_sig_genesB=get_sig_gene_set(snp_dictB,sig_snpsB)

if target and (not inB):
    # read file with gene IDs
    if verb:
        print "reading gene set from " + target
    target_genes=get_gene_list(target)
    
# create dict with gene ID counts
gene_counter=defaultdict(int)
#create random samples
if verb:
    print "random sampling " 

for i in range(0,rnd_samp):
    if verb and i%1000 == 0:
        print "sample " + str(i)
    rnd_genes=get_random_geneIDs(snp_dictA,sig_snpsA)
    if inB:
        target_genes=get_random_geneIDs(snp_dictB,sig_snpsB)
    com_gen_num=len(target_genes.intersection(rnd_genes))
    gene_counter[com_gen_num] += 1
#unravel gene_count to np.array, adding 6 zeros at the end to get nicer histograms:
gene_count=np.zeros(max(gene_counter.keys())+6,dtype=int)
for i in gene_counter.keys():
    gene_count[i] = gene_counter[i]
# write out results as histogram file
with open(outfile,'w') as of:
    if verb:
        print "writing outputfile" + outfile
    pVstr = ""
    if obs != False:
        pVstr = "p value for overlap of " + str(obs) + " genes = " + str(float(gene_count[obs:].sum())/gene_count.sum())
    print >>of, "\# "+str(rnd_samp)+" random samples" + pVstr
    print >>of,"Overlap\tCount"
    for i,x in enumerate(gene_count):
        print >>of, "{}\t{}\t".format(i,x)

if plot :
    if verb:
        print "creating histogramm plot " + outfile + ".sampling_hist.pdf"    
    plt.figure()
    plt.plot(np.arange(0,len(gene_count)),gene_count,'bx',linestyle="None")
    plt.ylabel("Counts")
    plt.xlabel("Overlapping Genes")
    if obs != False:
        if verb:
            print "p value for overlap of " + str(obs) + " genes = " + str(float(gene_count[obs:].sum())/gene_count.sum())     
        plt.axvline(obs,color='r',linestyle='dashed')
        plt.text(obs+1,gene_count.max()/2,"P("+str(obs)+")="+str(float(gene_count[obs:].sum())/gene_count.sum()))
    plt.savefig(outfile+".sampling_hist.pdf")


