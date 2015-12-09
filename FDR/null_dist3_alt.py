def get_indeces(nlpVals,inds=100):
    """
    gets a numpy array of sorted -log10 transformed pValues and returns an array of indeces of pValues randomly sampled from inds intervals spaced evenly on from 0 to max(nlpVals) 
    """
    # gets an array of pValues and returns an array of bins with counts
    max_val=nlpVals.max()
    max_val=nlpVals[len(nlpVals)*5e-5]
    indices=[]
    intervals=np.linspace(0,max_val,inds)
    for i in range(0,len(intervals) -1):
        interval=np.where((nlpVals > intervals[i]) & (nlpVals <= intervals[i+1]))[0]
        if len(interval):
            indices.append(interval[ np.random.randint(0,len(interval)) ])
    indices=np.array(indices)
    return indices

def dist_med(obs,null):
    """
    takes arrays and calculates the mean square distance from the first median 
    returns a float. 
    """
    distances = ((obs-null)**2)/len(obs)
    dist = np.sum(distances)
    return dist

def generate_beta_null_table(colsums, a=1, b=1):
    from numpy.random import beta, binomial #,seed
    #colsums = [x==0 and x+1 or x for x in colsums]
    #allele1
    # for testing purposes only
    # seed(100)
    if colsums[0] == 0:
        n11=0
    else:
        p11 = beta(a, b)
        p11=float(p11)
        n11 = binomial(colsums[0],p11)
    n21= colsums[0] -n11
    #allele2    
    if colsums[1] ==0:
        n12=0
    else:
        p12 = beta(a,b)
        p12=float(p12)
        n12 = binomial(colsums[1],p12)
    n22 = colsums[1]-n12
    return n11,n12,n21,n22
    
def generate_random_table(rowsums, colsums, dim=[2,2]):
    from numpy.random import hypergeometric
    #from numpy import array
    n11= hypergeometric(rowsums[0], rowsums[1],  colsums[0], size=1)[0]
    n21= colsums[0] -n11
    n12 = rowsums[0]-n11
    n22 = rowsums[1]-n21
    #table =array([[n11, n21],[ n12, n22]])
    return n11,n12,n21,n22
    
def cmh(replicate_tables):
    #from rpy2.robjects.packages import importr
    #    from rpy2.robjects import r
    #thanks to robert or martin for example code
    #r=robjects.r
    response_robj = robjects.IntVector(replicate_tables)
    replicatecount = int(float(len(replicate_tables))/4.0)
    dim_robj=robjects.IntVector([2,2,replicatecount])
    response_rar=robjects.r['array'](response_robj,dim=dim_robj)
    testres=r['try_CMHtest'](response_rar)
    if len(testres) > 1:
        pvalue=testres[2][0]
    else:
        pvalue=1.0
    if np.isnan(pvalue):
        pvalue=1.0
    #testres=r['mantelhaen.test'](response_rar)
    #pvalue=testres[2][0]
    #assert(pvalue<=1)
    #assert(pvalue>=0)
    return pvalue

def create_hist(pVals,bins=1000,log=False,max_val=False):
    """
    gets a numpy array of pValues and returns an array of bins with counts. if log=True, p values are first transformed to their negative decadic logarithm and then binned. If max_val provided binning performed from 0 to max_val.
    """
    # gets an array of pValues and returns an array of bins with counts
    hist_count = zeros(bins)
    if (log):
        nlp = -1.0*np.log10(pVals)
        if (max_val):
            nlp_max=max_val
        else:
            nlp_max = nlp.max()
        pVals = pVals/nlp_max
    for i in pVals:
        j=int(bins*i)
        if j >= bins: 
        # brute force setting (not great)
            j = bins - 1
        hist_count[j] += 1
    return hist_count

#need sync parser
class SyncLineParser(object):
    def __init__(self,syncline, replicates=None, includedel=False):
        '''
        SyncLineParser
        pass it a single line from a sync file
        optional options:
            which (default=None) will indicate which columns have the phenotype of interest
            includedel  (default=False) set to True indicates that the last two columns of sync should be included
        other attributes set on init:
                self.whichcols #
                self.chr 
                self.pos 
                self.ref 
                self.seqs #sync columns with bp information
                self.cmhp
                self.replicates #dictionary with "which"
        functions:
            self.get_two_major_alleles()
                 sets self.majminor, gets the two major alleles, ignoring anything higher than a di-allellic site
            self.get_pop_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            self.get_overall_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            
        '''
        self.includedel = includedel
        if self.includedel:
            self.ncol = 6
        else:
            self.ncol = 4
        self.whichcols=replicates #labels for replicates
        #include last two columns if True
        #parse syncline
        sline = syncline.split()
        self.chr =sline.pop(0)
        self.pos =sline.pop(0)
        self.ref =sline.pop(0)
        self.seqs = [[int (x) for x in y.split(':')[:self.ncol] ] for y in sline if ':' in y]
        if ':' not in sline[-1]: #CMH output
            self.cmhp =sline.pop()
        else:
            self.cmhp =None
        #make dictionary with information for phenotype or replicate
        self.majminor = None
        self.pop_allele_freqs =None
        self.overall_allele_freqs=None
        self.seq_dict =None
        self.reduced_seq = None
        self.reduced_seq_dict = None


    def get_seq_dict(self):
        if not self.seq_dict:
            self.seq_dict = {}.fromkeys(set(self.whichcols))
            for r in range(0, len(self.seqs)):
                rep = self.whichcols[r]
                self.seq_dict.setdefault(rep,[]).append(self.seqs)
        return self.seq_dict

    def get_two_major_alleles(self):
        if not self.majminor:
            whichcols = range(0, self.ncol)
            allele_totals = [sum([y[x] for y in self.seqs]) for x in whichcols]
            allele_totals =zip(allele_totals,whichcols )
            allele_totals.sort(reverse=True)
            self.majminor = [allele_totals[0][1], allele_totals[1][1]]
        return self.majminor

    def get_pop_allele_freqs(self, two =True):
        if not self.pop_allele_freqs:
            if two and not self.majminor:
                self.get_two_major_alleles()
            if two:
                whichcols = self.majminor
            else:
                whichcols = range(0, self.ncol)
            reduced_seq = [[float(x[y]) for y in whichcols] for x in self.seqs] #reduce 
            pop_totals = [float(sum(x)) for x in reduced_seq]
            self.pop_allele_freqs =[]
            for i in range(len(self.seqs)):
                self.pop_allele_freqs.append([x/pop_totals[i] for x in reduced_seq])
        return self.pop_allele_freqs
    
    def get_overall_allele_freqs(self, two=True):    
        if not self.overall_allele_freqs:
            if not self.pop_allele_freqs:
                self.get_pop_allele_freqs(two)
            whichcols = range(0, self.pop_allele_freqs[0])
            self.overall_allele_freqs =[sum([y[x] for y in self.pop_allele_freqs]) for x in whichcols]
        return self.overall_allele_freqs
    
    def get_reduced_seq (self):
        if not self.reduced_seq:
            if not self.majminor:
                self.get_two_major_alleles()
            self.reduced_seq = [[x[y] for y in self.majminor] for x in self.seqs] #reduce 
        return self.reduced_seq

    def get_reduced_seq_dict(self):
        if not self.reduced_seq_dict:
            if not self.reduced_seq:
                self.get_reduced_seq()
            self.reduced_seq_dict = {}
            for r,seq in enumerate(self.reduced_seq):
                rep = self.whichcols[r]
                self.reduced_seq_dict.setdefault(rep,[]).append(seq)
        return self.reduced_seq_dict

def get_replicate_tables(rep_dict):
    replicate_tables = []
    for r in rep_dict:
        for t in rep_dict[r]:
            replicate_tables.extend(t)
    return replicate_tables



def chi2distance (obs_pval_hist,null_pval_hist):
    """
    takes two histograms with identical bin numbers and calculates the chisquare distance between them.
    returns a float. 
    """
    #sum1 = (obs_pval_hist+null_pval_hist)*2
    chi2 = np.divide((obs_pval_hist-null_pval_hist)**2,(obs_pval_hist+null_pval_hist)*2)
    chisum = np.sum(np.ma.masked_invalid(chi2))
    return chisum

def kull_leib_distance (obs_pval_hist,null_pval_hist):
    """
    takes two histograms (numpy arrays) with identical bin numbers and calculates the average of the Kullback_Leibler divergences in both directions KL(hist1-hist2) and KL(hist2,hist1) between them. returns a float
    """
    kl_on = 0
    kl_no = 0
    #calculate total number of points for normalisation
    tot_o=np.sum(obs_pval_hist)
    tot_n=np.sum(null_pval_hist)
    # normalise
    obs_p_norm = obs_pval_hist / tot_o
    null_p_norm = null_pval_hist / tot_n
    # use masked arrays to exclude inf and NaN
    division = np.divide(obs_p_norm,null_p_norm)
    kl_on = np.sum( np.ma.masked_invalid(obs_p_norm * np.log(division)))
    division = np.divide(null_p_norm,obs_p_norm)
    kl_no = np.sum(np.ma.masked_invalid(null_p_norm * np.log(division)))
    # go two ways as not symmetric and take the mean
    return (kl_on+kl_no)/2

        
if __name__ == '__main__':
    #get one random table per replicate
    import sys, os,re
    import gzip
    import numpy as np
    #import Sync_parser
    from numpy import array, zeros
    np.seterr(all=None, divide='ignore', over=None, under=None, invalid='ignore')
    import math
    import rpy2.robjects as robjects
    from numpy import array, zeros
    #from scipy.stats import chisquare
    import argparse
    parser = argparse.ArgumentParser(description='create P values for a CMH test on a sync file with the alleles shuffled between control and case according to a beta distribution') 

    parser.add_argument("--alpha","-a", dest="alpha", help="the value of alpha to be tested", required=True)
    parser.add_argument("--reps","-r", dest="replicates", help="the positions of replicates in the CMH file (def. \'1,2,3,1,2,3\')",default='1,2,3,1,2,3')
    parser.add_argument("--in","-i", dest="infile", help="CMH file to be read", required=True)
    parser.add_argument("--app", dest="appendix", help="outfile appendix (optional)",default="")
    parser.add_argument("--b2a", dest="b2a", help="value for beta to alpha ratio (default beta=alpha)",default="1")
    parser.add_argument("--log", action="store_true", dest="log", help="also calculate logarithmic binning (default off)",default=False)
    parser.add_argument("--coverages", dest="coverages", help="coverages of populations for adapting beta (default none)",default=False)
    parser.add_argument("--chi2", action="store_true", dest="calc_dist", help="calculate chi square distance to obs pVal histogram (default False)",default=False)
    parser.add_argument("--500k", action="store_true", dest="five_hun", help="only do the first 500000 snps (default False)",default=False)
    parser.add_argument("--10x", action="store_true", dest="tenx", help="calculate ten times as many pValues and do not output the distances (sets chi2 to false)(default False)",default=False)
    parser.add_argument("--dontreuse", action="store_false", dest="reuse", help="do not use existing p values",default=True)
#    parser.add_argument("--adjust_beta", dest="adjust_beta", help="adjust beta to fit mean coverage (give mean coverages of populations eg, [100,])",default="0")

    args = parser.parse_args()
    infile = vars(args)['infile']
    log = vars(args)['log']
    reuse = vars(args)['reuse']
    fh_K = vars(args)['five_hun']
    calc_dist=False
    if (vars(args)['calc_dist']):
        calc_dist=True
    # if tenx set do 10 pValues for each site and not calculate any distances
    tenx=vars(args)['tenx']
    rep_pV=1 # repeat random pV for each locus rep_pV times
    if (tenx):
        #calc_dist=False
        rep_pV=10
    alph = float(vars(args)['alpha'])
    b = alph
    if (float(vars(args)['b2a']) > 0.0):
        #setting b to beta value
        b = float(vars(args)['b2a'])*alph
    #replicates is specific to the way this data is organized, 0 is a placeholder and removed 
    #replicates =[1,2,1,2]
    replicates = [ int(x) for x in vars(args)['replicates'].split(',')]
    reps = list(set(replicates))
    # calculate directory of betas
    if (vars(args)['coverages']):
        covs = [ float(x) for x in vars(args)['coverages'].split(',')]
    else:
        # set all coverages to 1
        covs = [1.0 for x in range(len(replicates))]
    # each population gets its own beta
    if (len(covs) != len(replicates)):
        sys.exit("coverages not same length as replicates")
    betas={}
    for i in reps:
        # get indeces of replicates
        indx = [x for x in range(len(replicates)) if replicates[x] == i]
        if (len(indx) != 2):
            sys.exit("some populations not in pairs")
        # calculate coverage ratio
        cov_ratio = float(covs[indx[0]])/(covs[indx[0]]+covs[indx[1]])
        # calculate beta for each replicate
        betas[i]= alph*(1.0-cov_ratio)/cov_ratio
    beta_string = "beta"
    # create string representation for beta values used  
    for i in sorted(list(betas)):
        beta_string += "_r{0}_{1:.2f}".format(i,betas[i])    
    file_appendix =vars(args)['appendix']
    if not file_appendix and tenx:
        file_appendix = "_10x"
    if  (float(vars(args)['b2a']) != 1):
        outfile_name = os.path.split(infile)[-1]+"_a_"+str(alph)+"_b_"+str(b)+"_nullP.out"+file_appendix
    elif (vars(args)['coverages']):
        outfile_name = os.path.split(infile)[-1]+"_a_"+str(alph)+"_"+beta_string+"_nullP.out"+file_appendix
    else:
        outfile_name = os.path.split(infile)[-1]+"_"+str(alph)+"_nullP.out"+file_appendix
    out = open(outfile_name, 'w')
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")
    #begin parsing sync file
    count = 0
    # for calcualting chi2 distance
    obs_pvals= []
    null_pvals = []
    r=robjects.r
    robjects.r('''
        try_CMHtest <- function(count_table) {
        try( mantelhaen.test(count_table),
        silent=TRUE
        )
        }        
        ''')
    for line in inf:
        if (re.match('\s*#',line)):
            continue
        syncline = SyncLineParser(line, replicates)
        reducedseq = syncline.get_reduced_seq_dict()
        #check a few
        if syncline.cmhp and reuse:
            if calc_dist:
                obs_pvals.append(float(syncline.cmhp))
            if count< 10:
                replicate_tables = get_replicate_tables(reducedseq)
                thiscmhp = cmh(replicate_tables)
                assert abs(thiscmhp-float(syncline.cmhp)) < 0.005, "CMH discrepancy.  Check replicate order"+str(thiscmhp)+" "+str(syncline.cmhp)
        else:
            if calc_dist:
                replicate_tables = get_replicate_tables(reducedseq)
                #print line
                thiscmhp = cmh(replicate_tables)
                obs_pvals.append(thiscmhp)
        count +=1
        if count%100000 == 0:
            print "Done with "+ str(count) +" sites"
        
        for _ in range(0,rep_pV):
            random_tables = []
            for rep in reps: #get random 2x2 table for each replicate
                # rowsums= reducedseq[rep].sum(axis=1)
                # colsums= reducedseq[rep].sum(axis=0)
                #rep = reps[r] #use this replicate
                rowsums = [sum(y) for y in reducedseq[rep]] #sum of 2 alleles within replicate
                colsums=[]    #sum of alleles across replicates
                for i in [0,1]: 
                    sumcol =0
                    for j in reducedseq[rep]:
                        sumcol+=j[i]
                    colsums.append(sumcol)
                t= generate_beta_null_table(colsums, a=alph, b=betas[rep])
                random_tables.extend(t)
            random_cmhp = cmh(random_tables)    
            if calc_dist:
                null_pvals.append(random_cmhp)
            print >> out, random_cmhp
        if fh_K and count == 500000:
            break
    if calc_dist:
        obs_pvals=np.array(obs_pvals)
        null_pvals=np.array(null_pvals) 
        obs_pval_hist = create_hist(obs_pvals)
        null_pval_hist = create_hist(null_pvals)
        out_file= os.path.split(infile)[-1]+"chisquare.out"
        chisum = chi2distance(obs_pval_hist,null_pval_hist)
        #kuld = kull_leib_distance(obs_pval_hist,null_pval_hist)
        if len(obs_pvals) == len(null_pvals): # only get distance from 1st median if length of vectors equal  
            obs_pvals.sort()
            nl_obs=-1*np.log10(obs_pvals)
            null_pvals.sort()
            nl_null=-1*np.log10(null_pvals)
            indeces=get_indeces(nl_null)
            null_sam=nl_null[indeces]
            obs_sam=nl_obs[indeces]
            dist=dist_med(obs_sam,null_sam)
            outstring = "{0}\t{1}\t{2:.2f}\t{3:.2f}".format(alph,beta_string,chisum,dist)
        else:
            outstring = "{0}\t{1}\t{2:.2f}\tNA".format(alph,beta_string,chisum)
        with open(out_file,"a") as f:
            beta_string=",".join([str(x)+":"+str(betas[x]) for x in sorted(list(betas)) ])
            if (log):
                max_val = min(min(obs_pvals),min(null_pvals))
                max_val=math.log10(max_val)
                obs_pval_hist = create_hist(obs_pvals,bins=1000,log=True,max_val=max_val)
                null_pval_hist = create_hist(null_pvals,bins=1000,log=True,max_val=max_val)
                chisum_log = chi2distance(obs_pval_hist,null_pval_hist)
                #kuld_log = kull_leib_distance(obs_pval_hist,null_pval_hist)
                outstring += "\t{0:.2f}".format(chisum_log)
            print >>f, outstring
    inf.close()
