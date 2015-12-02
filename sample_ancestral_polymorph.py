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

def add_anc_polyms(snp_dict,ap_file):
    """
    reads a sync file with mapped polymorphisms in D. sim from a file and adds the info to the snp_dict as a list [ 0/1 , 00/01/.. ]
    """
    inf = open(ap_file,"r")
    for i in snp_dict.keys():
        snp_dict[i]['ap']=np.empty(len(snp_dict[i]['bps']), dtype=list)
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
                snp_dict[fields[0]]['ap'][idx]= [ int(fields[5]), fields[6]  ]
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
            

def create_random_sample_alt(idx_bins,count_array):
    """
    uses unique_sample_of_int to sample SNP indeces from the bin index structure  according to afs
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

def get_random_aps(snp_dict,sig_snp,chroms):
    """
    return a dictionary of chromosomes with a numpy array of three entries (zero,one and two shared alleles)
    """
    aps = { x : np.zeros(3,dtype=int) for x in chroms }
    aps['total']= np.zeros(3,dtype=int)
    for chrom in chroms:
        idxs = create_random_sample_alt(snp_dict[chrom]['bin_idx'],sig_snp[chrom]['afs_hist']['counts'])
        for idx in idxs:
            if snp_dict[chrom]['ap'][idx]:
                # add to the right field
                aps[chrom][ 2 -snp_dict[chrom]['ap'][idx][1].count("0") ] += 1
                aps['total'][ 2 -snp_dict[chrom]['ap'][idx][1].count("0") ] += 1 
    return aps


import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser(description="reads a file with SNPs with allele frequencies (AF), and creates random samples of the same size, chromosome and AF distribution as the SNPs from another file - eg. the significant SNPs. It gives a histogram of ancestral polymorphisms/ancestral allleles found in the samples.\n  To calculate a P Value the observed number of overlapping genes has to be given.")

parser.add_argument("--sA", dest="inA", help="syncfile A with all SNPs and AFs", required=True)
parser.add_argument("--sigA", dest="sigA", help="syncfile of significant SNPs of A", required=True)
parser.add_argument("--aA", dest="aA", help="syncfile A with mapping to D. sim.", required=True)
parser.add_argument("--pops","-p", dest="pops",  help="populations in A to use for the average AFs (default: \"4,5,6\")",default="4,5,6")
parser.add_argument("--perc", dest="perc",  help="percent step for AF binning (default: 1.0)",default=1.0)
parser.add_argument("-n", dest="rnd_samp",  help="number of random samples to create", default="10000")
parser.add_argument("-o", dest="outfile",  help="outputfile", required=True)
parser.add_argument("-v", dest="verb",  action="store_true", help="outputfile", default=False)
#parser.add_argument("--plot", dest="plot",  action="store_true", help="create a plot of the drawn values (default: False)", default=False)
parser.add_argument("-q","--quant", dest="quant", help="quantiles in percent to output apart from the median, the minimum and the maximum (default: \"0.01,2.5,97.5,99.99\")", default="0.01,2.5,97.5,99.99")
parser.add_argument("--chrs", dest="chroms", help="chromosomes to collect values for (default: \"X,2L,2R,3L,3R,4\" )", default="X,2L,2R,3L,3R,4")

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


args = parser.parse_args()
inA = vars(args)['inA']
aA = vars(args)['aA']
sigA = vars(args)['sigA']
pops = np.array( vars(args)['pops'].split(","),dtype=int )
rnd_samp = int(vars(args)['rnd_samp'])
perc = float(vars(args)['perc'])
outfile = vars(args)['outfile']
verb= vars(args)['verb']
#plot= vars(args)['plot']
chroms=  vars(args)['chroms'].split(",")
quant=  [ float(x) for x in vars(args)['quant'].split(",") ]
quant.insert(0,0.0)
quant.insert(0,50.0)
quant.append(100.0)
# create SNP dict from file
if verb:
    print "create SNP dict from file "+ inA
snp_dictA= read_cmh_file(inA,pops)
if verb:
    print "adding polymorphism data from "+ aA
add_anc_polyms(snp_dictA,aA)
# read signifcant snps
if verb:
    print "reading sign snps "+ sigA
sig_snpsA = read_cmh_file(sigA,pops)
# get number of signif SNPs
sig_snp_num={ x:len(sig_snpsA[x]['bps']) for x in sig_snpsA.keys() }
sig_snp_num['all_chroms']=np.sum([ sig_snp_num[x] for x in sig_snpsA.keys() ])
sig_snp_num['total']=np.sum([ sig_snp_num[x] for x in chroms])
# create histogram
if verb:
    print "creating AF histogram and index bins"
create_af_histo(sig_snpsA,perc)
# create index bins
create_bin_indeces(snp_dictA,sig_snpsA)

# create dict with gene ID counts
ap_counter= { x: np.zeros((rnd_samp,3),dtype=int) for x in chroms }
ap_counter['total']= np.zeros((rnd_samp,3),dtype=int)
#create random samples
if verb:
    print "random sampling " 

for i in range(0,rnd_samp):
    if verb and i%10000 == 0:
        print "sample " + str(i)
    rnd_aps=get_random_aps(snp_dictA,sig_snpsA,chroms)
    for chrom in rnd_aps.keys():
        ap_counter[chrom][i] = rnd_aps[chrom]
# write out results as histogram file
out_str="\t".join([ x+str(y) for x in ["Z_","N_","O_","T_","GEO_"]  for y in quant ])
out_str = "Chrom\t"+out_str
with open(outfile,'w') as of:
    if verb:
        print "writing outputfile" + outfile
    print >>of, "\# "+str(rnd_samp)+" random samples, sample size: "+str(sig_snp_num['all_chroms'])+"; snps on sampled chroms:"+ str(sig_snp_num['total'])
    print >>of, out_str
    # total = False
    for chrom in sorted(ap_counter.keys()):
        # if type(total) == np.ndarray:
        #     total += ap_counter[chrom]
        # else:
        #     total =  ap_counter[chrom].copy()
        output="\t".join([ str(stats.scoreatpercentile(sig_snp_num[chrom]-(ap_counter[chrom][:,2]+ap_counter[chrom][:,1]),y)) for y in quant ])
        output+="\t"+"\t".join([ str(stats.scoreatpercentile(ap_counter[chrom][:,x],y)) for x in range(0,3) for y in quant ])
        output+="\t"+"\t".join([ str(stats.scoreatpercentile(ap_counter[chrom][:,2]+ap_counter[chrom][:,1],y)) for y in quant ])
        #print "{}\t{}".format(chrom,output)
        print >>of, "{}\t{}".format(chrom,output)    
    # output="\t".join([ str(stats.scoreatpercentile(sig_snp_num['tot_chrom']-(total[:,2]+total[:,1]),y)) for y in quant ])
    # output+="\t"+"\t".join([ str(stats.scoreatpercentile(total[:,x],y)) for x in range(0,3) for y in quant ])
    # output+="\t"+"\t".join([ str(stats.scoreatpercentile(total[:,2]+total[:,1],y)) for y in quant ])
    # print >>of, "total\t"+output
                         

# if plot :
#     if verb:
#         print "creating histogramm plot " + outfile + ".sampling_hist.pdf"    
#     plt.figure()
#     plt.plot(np.arange(0,len(gene_count)),gene_count,'bx',linestyle="None")
#     plt.ylabel("Counts")
#     plt.xlabel("Overlapping Genes")
#     if obs != False:
#         if verb:
#             print "p value for overlap of " + str(obs) + " genes = " + str(float(gene_count[obs:].sum())/gene_count.sum())     
#         plt.axvline(obs,color='r',linestyle='dashed')
#         plt.text(obs+1,gene_count.max()/2,"P("+str(obs)+")="+str(float(gene_count[obs:].sum())/gene_count.sum()))
#     plt.savefig(outfile+".sampling_hist.pdf")


