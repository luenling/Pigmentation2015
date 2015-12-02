def cmh(replicate_tables,confint=True):
    """
    gets replicate table an optional boolean confint and returns a string with either pValue, odds ratio and upper and lower conf. interval odds ratio, or just pV and odds ratio, depending on the value of confint
    if only one replicate given, uses fishers exact test
    """
    # import rpy2.robjects as robjects
    # from rpy2.robjects.packages import importr
    #from rpy2.robjects import r
    #thanks to robert or martin for example code
    #r=robjects.r
    response_robj = robjects.IntVector(replicate_tables)
    replicatecount = int(float(len(replicate_tables))/4.0)
    if replicatecount == 1:
        dim_robj=robjects.IntVector([2,2])
    else:
        dim_robj=robjects.IntVector([2,2,replicatecount])
    response_rar=robjects.r['array'](response_robj,dim=dim_robj)
        #testres=r['mantelhaen.test'](response_rar)
    if replicatecount == 1:
        testres=r['try_fishertest'](response_rar)
    else:
        testres=r['try_CMHtest'](response_rar)
    if len(testres) > 1:
        results=dict(zip(testres.names,list(testres)))
        pvalue=results['p.value'][0]
        odds_int=[results['conf.int'][0],results['conf.int'][1]]
        odds_ratio=results['estimate'][0]
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
            if twoonly and not self.majminor:
                self.get_two_major_alleles()
            if twoonly:
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


if __name__ == '__main__':
    #get one random table per replicate
    import sys, os, re
    import numpy as np
    #import Sync_parser
    from numpy import array, zeros
    import math
    import gzip
    import argparse
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    parser = argparse.ArgumentParser(description='create P values, ods ratio and conf intervals for a CMH test on a sync file') 

    parser.add_argument("--reps","-r", dest="replicates", help="the positions of replicates in the CMH file (def. \'1,2,3,1,2,3\')",default='1,2,3,1,2,3')
    parser.add_argument("--pops","-p", dest="pops", help="the population pairs in the CMH file to test, mutual exclusive with reps, if set, reps is ignored, also do multiple cmh tests with the same using \"|\" (eg. \"1:2|1:3,1:4\")",default=False)
    parser.add_argument("--in","-i", dest="infile", help="CMH file to be read", required=True)
    parser.add_argument("--confint","-c", dest="confint",  action='store_false', help="confidence intervals", default=True)
    nucs = [x for x in "ATCG"]
    args = parser.parse_args()
    infile = vars(args)['infile']
    replicates = [ int(x) for x in vars(args)['replicates'].split(',')]
    reps = list(set(replicates))
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
        try_fishertest <- function(count_table) {
        try( fisher.test(count_table),
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
        for i in reps:
            # get indeces of replicates
            indx = [x for x in range(len(replicates)) if replicates[x] == i]
            if (len(indx) != 2):
                sys.exit("some populations not in pairs")
    if re.search("\.b?gz",infile):
        inf = gzip.open(infile,'rb')
    else:
        inf = open(infile,"r")

    for line in inf:
        if re.match("\s*\#",line):
            # comment lines
            continue
        line=line.rstrip()
        if type(pops) != bool:
            syncline = SyncLineParser(line)
            reducedseq = syncline.get_reduced_seq()
            thiscmhp=[]
            for i in pops:
                replicate_tables=[]
                for idx in i: 
                    replicate_tables.extend(reducedseq[idx])
                thiscmhp.append(cmh(replicate_tables,confint))
            thiscmhp="\t".join(thiscmhp)
        else:
            syncline = SyncLineParser(line, replicates)
            reducedseq = syncline.get_reduced_seq_dict()
            replicate_tables = get_replicate_tables(reducedseq)
            thiscmhp = cmh(replicate_tables,confint)
        alleles = "".join([nucs[x] for x in syncline.majminor])
        print syncline.chr+"\t"+str(syncline.pos)+"\t"+alleles+"\t"+thiscmhp
    inf.close()
