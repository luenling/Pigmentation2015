import sys, os, re
from collections import defaultdict
import argparse
import gzip
import collections
import numpy as np
from scipy.stats import fisher_exact

#need sync parser
class SyncLineParser(object):
    def __init__(self,syncline, replicates=None, includedel=False, strands=False, min_count = 0, min_fract = 0):
        '''
        SyncLineParser
        pass it a single line from a sync file
        optional options:
            which (default=None) will indicate which columns have the phenotype of interest
            includedel  (default=False) set to True indicates that the last two columns of sync should be included
        other attributes set on init:
                self.whichcols #
                self.min_count : all alleles having a count < mincount over all populations are set to zero (default:0) 
                self.min_fract : all alleles having a frequency below min_fract in a population are set to zero in that population (default: 0)
                self.chr 
                self.pos 
                self.ref 
                self.seqs #sync columns with bp information
                self.cmhp
        functions:
            self.get_two_major_alleles()
                sets self.majminor, gets the two major alleles, ignoring anything higher than a di-allellic site
            self.get_pop_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            self.get_overall_allele_freqs(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            self.get_reduced_allele_freq_dict(two=True)
                gets overall allele freqs from seq information across all replicates
                if two is set to True, major and minor alleles only
            
        '''
        self.includedel = includedel
        self.strands = strands
        if self.includedel:
            self.ncol = 6
        else:
            self.ncol = 4
        self.whichcols=replicates #labels for replicates
        #include last two columns if True
        #parse syncline
        sline = syncline.split()
        self.min_count = min_count
        self.min_fract = min_fract
        self.chr =sline.pop(0)
        self.pos =sline.pop(0)
        self.ref =sline.pop(0)
        self.seqs = [ np.array([int (x) for x in y.split(':')[:self.ncol] ]) for y in sline if ':' in y]
        if self.strands:
            self.fw_strand = [ np.array([int (x) for x in y.split(':')[6:6+4] ]) for y in sline if ':' in y]
            self.strand_bias=None
        # get rid of alleles with count less than min_count
        if min_count:
            self.set_seqs_mc()
        if min_fract:
            self.set_seqs_mf()
        if ':' not in sline[-1]: #CMH output
            self.cmhp =sline.pop()
        else:
            self.cmhp =None
        #make dictionary with information for phenotype or replicate
        self.majminor = None
        self.coverages = None
        self.pop_allele_freqs =None
        self.overall_allele_freqs=None
        self.seq_dict =None
        self.reduced_seq = None
        self.reduced_seq_dict = None
        self.reduced_af_dict = None
        

    def set_seqs_mc(self):
        seq_counts = [ np.array([ x[y] for x in self.seqs ]).sum() for y in range(0, self.ncol) ]
        for seq in self.seqs:
            for i,count in enumerate(seq_counts):
                if count > 0 and count < self.min_count:
                    seq[i] = 0
        return self.seqs

    def set_seqs_mf(self):
        fracts = [ x[:self.ncol]/float(x[:self.ncol].sum()) for x in self.seqs ]
        for i,seq in enumerate(self.seqs):
            for j in range(0, self.ncol) :
                if seq[j] > 0 and ( fracts[i][j] < self.min_fract):
                    seq[j] = 0
        return self.seqs

        
    def get_seq_dict(self):
        if self.seq_dict == None:
            self.seq_dict = collections.defaultdict(list)
            for r in range(0, len(self.seqs)):
                self.seq_dict[self.whichcols[r]].append(self.seqs[r])
        return self.seq_dict

    def get_two_major_alleles(self):
        if not self.majminor:
            whichcols = range(0, self.ncol)
            allele_totals = np.array([sum([y[x] for y in self.seqs]) for x in whichcols])
            # get the highest ranked columns in reverse order
            if sorted(allele_totals)[-2] > 0:
                self.majminor =list(allele_totals.argsort()[-1:-3:-1])
            else:
                self.majminor = [ allele_totals.argsort()[-1] ]
        return self.majminor
    def get_pop_coverages(self, two = True):
        if not self.coverages:
            if two and not self.majminor:
                self.get_two_major_alleles()
            if two:
                whichcols = self.majminor
            else:
                whichcols = range(0, self.ncol)
            reduced_seq = np.array([[float(x[y]) for y in whichcols] for x in self.seqs]) #reduce
            self.coverages =  [ x.sum() for x in reduced_seq]
        return self.coverages

    def get_pop_allele_freqs(self, two = True):            
        if not self.pop_allele_freqs:
            if two and not self.majminor:
                self.get_two_major_alleles()
            if two:
                whichcols = self.majminor
            else:
                whichcols = range(0, self.ncol)
                #print whichcols
            reduced_seq = np.array([[float(x[y]) for y in whichcols] for x in self.seqs]) #reduce
            #pop_totals =  [ x.sum() for x in reduced_seq]
            pop_totals = self.get_pop_coverages(two)
            self.pop_allele_freqs =[]
            for i in range(len(self.seqs)):
                if (pop_totals[i] > 0):
                    freq=reduced_seq[i]/pop_totals[i]
                else:
                    freq=np.zeros(len(reduced_seq[i]))+np.nan
                self.pop_allele_freqs.append(freq)
        return self.pop_allele_freqs
    
    def get_overall_allele_freqs(self, two=True):   
        if not self.overall_allele_freqs:
            self.overall_allele_freqs=[]
            if not self.pop_allele_freqs:
                self.get_pop_allele_freqs(two)
            num_pop=len(self.pop_allele_freqs)
            self.overall_allele_freqs = [sum([y[x] for y in self.pop_allele_freqs])/num_pop for x in range(0, len(self.pop_allele_freqs[0]))]
        return self.overall_allele_freqs
    
    def get_reduced_seq (self, two = True):
        if self.reduced_seq == None:
            if not self.majminor:
                self.get_two_major_alleles()
            self.reduced_seq = np.array([x[self.majminor] for x in self.seqs]) #reduce 
        return self.reduced_seq

    def get_reduced_seq_dict(self):
        if self.reduced_seq_dict == None:
            if self.reduced_seq == None:
                self.get_reduced_seq()
            self.reduced_seq_dict = {}
            for r,seq in enumerate(self.reduced_seq):
                rep = self.whichcols[r]
                if rep in self.reduced_seq_dict:
                    self.reduced_seq_dict[rep]=np.vstack([self.reduced_seq_dict[rep],seq])
                else:
                    self.reduced_seq_dict[rep]=seq
        return self.reduced_seq_dict
            
    def get_reduced_allele_freq_dict(self,two=True):
        if self.reduced_af_dict == None:
            if not self.pop_allele_freqs:
                self.get_pop_allele_freqs(two)
            self.reduced_af_dict = collections.defaultdict(list)
            for r,af in enumerate(self.pop_allele_freqs):
                rep = self.whichcols[r]
                self.reduced_af_dict[rep].append(af)
        return self.reduced_af_dict
    
    def get_strand_bias(self, two=True):
        """
        returns an array with three entries: first and array with all individual strandbias fisher pV, second the total  strandbias fisher pV, and last an array with the ratios of the supporting reads on the less common strand to the more common strand for the major and minor alleles. If there is anything going wrong with the first fisher test - maybe not a true snp or something messy, it returns False.
        """
        if not self.strands or not two:
            # no stranded information read
            return False
        if not self.majminor:
            self.get_two_major_alleles()
        whichcols = self.majminor
        num_pop=len(self.seqs)
        fw_rev = np.array([ [ self.fw_strand[indx][whichcols], self.seqs[indx][whichcols] - self.fw_strand[indx][whichcols] ] for indx in range(0,num_pop)],dtype=int)
        fw_rev_tot = fw_rev.sum(axis=0)
        try:
            pop_pv = [ fisher_exact(x)[1] for x in fw_rev ]
        except:
            self.strand_bias=False
            return self.strand_bias
        tot_pv = fisher_exact(fw_rev_tot)[1]
        min_supp_read = np.divide(fw_rev_tot.min(axis=0),fw_rev_tot.max(axis=0),dtype=float)
        self.strand_bias = [ pop_pv, tot_pv, min_supp_read ]
        return self.strand_bias

        # whichcols = syncline.majminor
        # num_pop=len(syncline.seqs)
        # cont_table=np.zeros((2,2),dtype=int)
        # for indx in range(0,num_pop):
        #     cont_table += np.array([ syncline.fw_strand[indx][whichcols], syncline.seqs[indx][whichcols] - syncline.fw_strand[indx][whichcols] ])      
# for line in inf:
#     syncline = Sync_parser.SyncLineParser(line,strands=True)
#     syncline.get_strand_bias()
#     print str(syncline.strand_bias)
#     print str(syncline.seqs)
#     print str(syncline.fw_strand)


def read_inversion_file(infile):
    """
    reads martins inversion file (format: CHR INV BPS ALLELE) and creates a dict with chrom->bps->Inv->Allele]
    """
    snp_dict=defaultdict(lambda: defaultdict(defaultdict))
    inf = open(infile,"r")
    #load dictionary
    for line in inf:
        line=line.rstrip()
        fields=line.split()
        snp_dict[fields[0]][fields[2]][fields[1]]=fields[3]
    inf.close()
    return snp_dict

parser = argparse.ArgumentParser(description='reads martin\'s inversion file and a sync/cmh file and gives the fraction of the alleles fixed in the inversion.') 

parser.add_argument("--inv", dest="inv", help="file with alleles segregating with inversions", required=True)
parser.add_argument("--infile", dest="infile", help="sync or cmh file to get positions from, can be gzipped", required=True)
parser.add_argument("--header", dest="head",  action="store_true", help="write header (default=False)", default=False)
parser.add_argument("--count", dest="count",  action="store_true", help="add overall count column (default=False)", default=False)
args = parser.parse_args()
infile = vars(args)['infile']
inv=vars(args)['inv']
head=vars(args)['head']
count=vars(args)['count']
nucs={'A':0,'T':1,'C':2,'G':3}
snp_dict=read_inversion_file(inv)

if re.search("\.b?gz",infile):
    inf = gzip.open(infile,'rb')
else:
    inf = open(infile,"r")
for line in inf:
    line=line.rstrip()
    fields=line.split()
    if ( fields[0] in snp_dict and fields[1] in snp_dict[ fields[0] ]):
        #print line
        entry = SyncLineParser(line)
        afs = entry.get_pop_allele_freqs(two=False)
        if head:            
            if count:
                pops = [ "pop"+str(x)+"\tcov"+str(x) for x in range(1,len(afs)+1) ]
            else:
                pops = [ "pop"+str(x) for x in range(1,len(afs)+1) ]
                
            # if entry.cmhp:
            #     pops.append("P")
            print "CHR\tINV\tBPS\tREF\tALT\t"+"\t".join(pops)
            head=False
        for x in snp_dict[ fields[0] ][ fields[1] ].keys():
            allele=snp_dict[ fields[0] ][ fields[1] ][x]
            afs_x=[str(i[nucs[allele]]) for i in afs]
            if count:
                cts_x = [str(i) for i in entry.get_pop_coverages() ]
                afs_x = ["\t".join((i,j)) for (i,j) in zip(afs_x,cts_x)]
            out_str=fields[0]+"\t"+x+"\t"+fields[1]+"\t"+entry.ref+"\t"+allele+"\t"+"\t".join(afs_x)
            print out_str
inf.close()




