def read_cmh_file(infile,pops,popsB,pop_pV,pVs,cI):
    """
    reads a cmh file and returns a dictionary of chromosomes with a key pointing at a numpy array of bp positions and one at the AFs for each chromsome key:
    eg cmhdict['2L']['bps'] = array[BPs],  cmhdict['2L']['afs'] = array[afs] and cmhdict['2L']['lOR1']=array[lOR1] cmhdict['2L']['lOR2']=array[lor2] cmhdict['2L']['pV'] = array[pV]
    and pVs array : [ [pVs,CHR,Afs,index] ]
    for the afs the average of the floats at the populations given in the argument pops is taken  
    """
    # intermediate lists of values
    cmh_dict= defaultdict(lambda: defaultdict(list)) #dictionary of chroms and different subfields defaulting to lists
    if re.search("\.b?gz$",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    pops += 2
    #popsB = popsB.copy()
    popsB = [ x+2 for x in popsB ]
    pop_pV += 2
    n_pops = len(pops)
    #load dictionary
    for line in inf:
        if re.match("\#",line):
            continue
        line.rstrip()
        fields=line.split()
        try:
            pV=float(fields[pop_pV])
        except:
            continue
        if not( pV >= 0 and pV <= 1.0):
            if not ( np.isnan(pV)):
                print "Dodgy Pv in: " + line 
            continue
        cmh_dict[fields[0]]['bps'].append(int(fields[1]))
        afs=sum([ float(fields[x]) for x in pops ])/n_pops
        cmh_dict[fields[0]]['afs'].append( afs )        
        cmh_dict[fields[0]]['pV'].append( pV )        
        if cI:
            cmh_dict[fields[0]]['lOR1'].append( float(fields[popsB[0][0]]) )
            cmh_dict[fields[0]]['lOR2'].append( float(fields[popsB[1][0]]) )
            cmh_dict[fields[0]]['lORCI1a'].append(float(fields[popsB[0][1]]))
            cmh_dict[fields[0]]['lORCI1b'].append(float(fields[popsB[0][2]]))
            cmh_dict[fields[0]]['lORCI2a'].append(float(fields[popsB[1][1]]))
            cmh_dict[fields[0]]['lORCI2b'].append(float(fields[popsB[1][2]]))
        else:
            cmh_dict[fields[0]]['lOR1'].append( float(fields[popsB[0]]) )
            cmh_dict[fields[0]]['lOR2'].append( float(fields[popsB[1]]) )          
        #pVs.append([ pV,fields[0], afs])
    # convert to np arrays
    for i in cmh_dict.keys():
        for j in cmh_dict[i].keys():
            cmh_dict[i][j] = np.array(cmh_dict[i][j])
    for i in cmh_dict.keys():
        sort_vector = cmh_dict[i]['bps'].argsort()
        for j in cmh_dict[i].keys():
            cmh_dict[i][j] = cmh_dict[i][j][sort_vector]
        cmh_dict[i]['lOR1']=np.nan_to_num(np.log(cmh_dict[i]['lOR1']))
        cmh_dict[i]['lOR2']=np.nan_to_num(np.log(cmh_dict[i]['lOR2']))
        #cmh_dict[i]['scOR']=np.sign(cmh_dict[i]['lOR2'])*np.sign(cmh_dict[i]['lOR1'])
        if cI:
            cmh_dict[i]['lORCI1a']=np.nan_to_num(np.log(cmh_dict[i]['lORCI1a']))
            cmh_dict[i]['lORCI1b']=np.nan_to_num(np.log(cmh_dict[i]['lORCI1b']))
            cmh_dict[i]['lORCI2a']=np.nan_to_num(np.log(cmh_dict[i]['lORCI2a']))
            cmh_dict[i]['lORCI2b']=np.nan_to_num(np.log(cmh_dict[i]['lORCI2b']))
            cmh_dict[i]['lORCI1']=np.sign(cmh_dict[i]['lORCI1a'])*(np.sign(cmh_dict[i]['lORCI1a']) == np.sign(cmh_dict[i]['lORCI1b']))
            cmh_dict[i]['lORCI2']=np.sign(cmh_dict[i]['lORCI2a'])*(np.sign(cmh_dict[i]['lORCI2a']) == np.sign(cmh_dict[i]['lORCI2b']))
            cmh_dict[i]['sORCI']=cmh_dict[i]['lORCI1']*cmh_dict[i]['lORCI2']
    inf.close()
    for i in cmh_dict.keys():
        for j,pV in enumerate(cmh_dict[i]['pV']):
            pVs.append([ pV, i, cmh_dict[i]['afs'][j], j])
    # sort pVs array
    pVs.sort(key=lambda x: x[0])
    return cmh_dict

def create_af_histo(snp_dict,perc):
    """
    creates a new entry for each chromosome holding an array with counts for each afs bin
    100.0/perc bins are created
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
        #print str(max) +" "+str(size) +" "+ str(num_left) +" "+str(idxs) 
        num_left = size - len(idxs)
    return idxs

def get_random_snps(snp_dict,sig_snp):
    snps=[]
    for chrom in sig_snp.keys():
        idxs = create_random_sample_alt(snp_dict[chrom]['bin_idx'],sig_snp[chrom]['afs_hist']['counts'])
        for idx in idxs:
            snps.append(chrom+str(snp_dict[chrom]['bps'][idx]))
    return set(snps)

def get_random_lOR(snp_dict,sig_snp,cI):
    lORs=[]
    for chrom in sig_snp.keys():
        idxs = create_random_sample_alt(snp_dict[chrom]['bin_idx'],sig_snp[chrom]['afs_hist']['counts'])
        #print str(snp_dict[chrom]['lOR1'][idxs])+"\n"+str(snp_dict[chrom]['lOR2'][idxs])
        for idx in idxs:
            if cI:
                lORs.append([snp_dict[chrom]['lOR1'][idx],snp_dict[chrom]['lOR2'][idx],snp_dict[chrom]['sORCI'][idx] ])
            else:
                lORs.append([snp_dict[chrom]['lOR1'][idx],snp_dict[chrom]['lOR2'][idx]])
    # print str(lORs)
    return np.array(lORs,dtype=float)

def get_snp_dict_by_ranks(pVs,rank):
    """
    gets a sorted pVs array with chrom and afs and returns a dict with chrom and afs 
    """
    snp_dict= defaultdict(lambda: defaultdict(list))
    for entry in pVs[0:rank]:
        snp_dict[entry[1]]['afs'].append( entry[2] )
    for chrom in snp_dict.keys():
        snp_dict[chrom]['afs'] = np.array(snp_dict[chrom]['afs'])
    return snp_dict

def get_lORs_by_ranks(snp_dict,pVs,rank,cI):
    """
    gets a sorted pVs array with chrom and afs and idxs and the snp_dict and returns a logOR array and the highest pV 
    """
    lORs=[]
    for entry in pVs[0:rank]:
        chrom=entry[1]
        idx=entry[3]
        if cI:
            lORs.append([snp_dict[chrom]['lOR1'][idx],snp_dict[chrom]['lOR2'][idx],snp_dict[chrom]['sORCI'][idx] ])
        else:
            lORs.append([snp_dict[chrom]['lOR1'][idx],snp_dict[chrom]['lOR2'][idx]])
    return [ np.array(lORs,dtype=float), pVs[rank-1][0] ]

import sys, os, re
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import argparse
import gzip
parser = argparse.ArgumentParser(description="reads a file with SNPs with allele frequencies (AF) and pValues and creates random samples of the same size, chromosome and AF distribution as the SNPs to test for correlation of logORs given inthe same file. If given three logORs per population, also consideres the CI for checking consistency.")

parser.add_argument("-i","--in", dest="inA", help="file with pos, base af and pV", required=True)
parser.add_argument("--ranks", dest="ranks", help="list of ranks for cutoff values for sampling (default:\"10,50,100,500,1000,5000,10000,50000,100000\")", default="10,50,100,500,1000,5000,10000,50000,100000")
#parser.add_argument("--sB", dest="inB", help="file with lOR to compare", default=False)
parser.add_argument("--pV", dest="pop_pV", help="position of pV", required=True)
parser.add_argument("--pA", dest="popsA",  help="populations to use for the average AFs (default: \"4,5,6\")",default="4,5,6")
parser.add_argument("--lOR", dest="popsB",  help="position in of the logORs to compare, if 2x3 given assumes avg(logOR) and confidence interval (default: \"1,2\"; for CI:\"4:5:6,1:2:3\")",default="1,2")
parser.add_argument("--detRanks", dest="detRanks",  help="list of ranks for cutoff values without random sampling (if given creates outfile.det)",default=False)
#parser.add_argument("--pB", dest="popsB",  help="populations in B to compare lOR, if 2x3 given assumes avg(logOR) and confidence interval (default: \"1,2\"; for CI:\"4:5:6,1:2:3\")",default="1,2")
parser.add_argument("--perc", dest="perc",  help="percent step for AF binning (default: 5.0)",default=5.0)
parser.add_argument("-n", dest="rnd_samp",  help="number of random samples to create for each rank cutoff", default="10000")
parser.add_argument("-o", dest="outfile",  help="outputfile", required=True)
parser.add_argument("-v", dest="verb",  action="store_true", help="verbouse output", default=False)

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


args = parser.parse_args()
inA = vars(args)['inA']
ranks =np.array( vars(args)['ranks'].split(","),dtype=int )
detRanks =vars(args)['detRanks']
#inB = vars(args)['inB']
popsA = np.array( vars(args)['popsA'].split(","),dtype=int )
pop_pV = int( vars(args)['pop_pV'])
popsB = vars(args)['popsB'].split(",")
cI=False
if ":" in popsB[0]:
    cI=True
    popsB = [np.array(x.split(":"),dtype=int ) for x in popsB   ]
    assert(len(popsB[0])==3 and len(popsB[1])==3)
else:
    popsB=np.array(popsB,dtype=int)
rnd_samp = int(vars(args)['rnd_samp'])
perc = float(vars(args)['perc'])
outfile = vars(args)['outfile']
verb= vars(args)['verb']
# array for pVs
pVs=[]
# create SNP dict from file
if verb:
    print "create SNP dict from file "+ inA
snp_dictA= read_cmh_file(inA,popsA,popsB,pop_pV,pVs,cI)

if detRanks != False:
    of=open(outfile+".det",'w')
    detRanks=np.array( detRanks.split(","),dtype=int )
    rho=[]
    max_pV=[]
    fr_id=[]
    fr_ci_id=[]
    for rank in detRanks:
        [ lORs, mpV ]=get_lORs_by_ranks(snp_dictA,pVs,rank,cI)
        max_pV.append(np.log10(mpV))
        rho.append(stats.spearmanr(np.array(lORs[:,0:2]))[0])
        sign_lORs=np.sign(lORs[:,0])*np.sign(lORs[:,1])
        fr_id.append(len(sign_lORs[sign_lORs > 0])/float(len(sign_lORs)))
        if cI:
            sign_lORCI=lORs[:,2]
            fr_ci_id.append(float(len(sign_lORCI[sign_lORCI > 0]))/len(sign_lORCI))
    for j,rank in enumerate(detRanks):
        if cI:
            if j == 0:
                print >>of,"rank\tpV\tfid\trho\tfid_ci"
            print >>of,"{}\t{}\t{}\t{}\t{}".format(rank,max_pV[j],fr_id[j],rho[j],fr_ci_id[j])
        else:
            if j == 0:
                print >>of,"rank\tfid\trho"
            print >>of,"{}\t{}\t{}\t{}".format(rank,max_pV[j],fr_id[j],rho[j])
    of.close()
    
of=open(outfile,'w')
# go through rank one by one
percentiles=[0.05,0.5,2.5,97.5,99.5,99.95]
fr_per=["fr_id_"+str(x) for x in percentiles]
rho_per=["rho_"+str(x) for x in percentiles]
fr_ci_per=["fr_ci_id_"+str(x) for x in percentiles]
if cI:
    print >>of,"rank\tfr_id_mean\t"+"\t".join(fr_per)+"\trho_id_mean\t"+"\t".join(rho_per)+"\tfr_ci_id_mean\t"+"\t".join(fr_ci_per)
else:
    print >>of,"rank\tfr_id_mean\t"+"\t".join(fr_per)+"\trho_id_mean\t"+"\t".join(rho_per)
#create_af_histo(snp_dictA,perc)
# go through rank one by one
for rank in ranks:
    # create histogram
    sig_snpsA=get_snp_dict_by_ranks(pVs,rank)
    if verb:
        print "creating AF histogram and index bins"
    create_af_histo(sig_snpsA,perc)
    # create index bins
    create_bin_indeces(snp_dictA,sig_snpsA)
    # create arrays with rhos and fractions of identical effect
    rho=[]
    fr_id=[]
    fr_ci_id=[]
    # random sampling
    if verb:
        print "random sampling rank "+str(rank)
        #print str(sig_snpsA)
    for i in range(0,rnd_samp):
        if verb and i%1000 == 0:
            print "rank " +str(rank) + " sample " + str(i)
        lORs=get_random_lOR(snp_dictA,sig_snpsA,cI)
        rho.append(stats.spearmanr(np.array(lORs[:,0:2]))[0])
        sign_lORs=np.sign(lORs[:,0])*np.sign(lORs[:,1])
        fr_id.append(len(sign_lORs[sign_lORs > 0])/float(len(sign_lORs)))
        if cI:
            sign_lORCI=lORs[:,2]
            fr_ci_id.append(float(len(sign_lORCI[sign_lORCI > 0]))/len(sign_lORCI))
    rho_perc=[ "{0:.4f}".format(x) for x in  np.percentile(rho,percentiles) ]
    fr_id_perc=[ "{0:.4f}".format(x) for x in np.percentile(fr_id,percentiles)]
    if cI:
        fr_ci_id_perc=[ "{0:.4f}".format(x) for x in np.percentile(fr_ci_id,percentiles)]
        print >>of,"{}\t{}\t{}\t{}\t{}\t{}\t{}".format(rank,np.mean(fr_id),"\t".join(fr_id_perc),np.mean(rho),"\t".join(rho_perc),np.mean(fr_ci_id),"\t".join(fr_ci_id_perc))
    else:
        print >>of,"{}\t{}\t{}\t{}\t{}".format(rank,np.mean(fr_id),"\t".join(fr_id_perc),np.mean(rho),"\t".join(rho_perc))
    of.flush()

of.close()


# #unravel gene_count to np.array, adding 6 zeros at the end to get nicer histograms:
# snp_count=np.zeros(max(snp_counter.keys())+6,dtype=int)
# for i in snp_counter.keys():
#     snp_count[i] = snp_counter[i]
# # write out results as histogram file
# with open(outfile,'w') as of:
#     if verb:
#         print "writing outputfile" + outfile
#     pVstr = ""
#     if obs != False:
#         pVstr = "p value for overlap of " + str(obs) + " snps = " + str(float(snp_count[obs:].sum())/snp_count.sum())
#     print >>of, "\# "+str(rnd_samp)+" random samples" + pVstr
#     print >>of,"Overlap\tCount"
#     for i,x in enumerate(snp_count):
#         print >>of, "{}\t{}\t".format(i,x)
