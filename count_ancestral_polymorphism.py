import sys, os, re
import numpy as np
import argparse
from collections import defaultdict

#Author: Lukas Endler
parser = argparse.ArgumentParser(description="""
reads a ancestral polymorphism file and outputs counts per chromosome
2L 1 00 xxx
2L 1 10 yyy
2L 1 01 zzz
2L 2 00 ...
no allele overlap: 00
combinations of the major (1), or minor (2) allele of simulans matching the major and minor of mel (position):
eg. major mel = major sim, minor mel != minor sim : 10
 minor mel = minor sim, major mel != major sim : 2
 minor mel = major sim: 1
  major mel = minor sim: 20
""") 
parser.add_argument("-i","--inf", dest="infile", help="ancestral polym comp file", required=True)
parser.add_argument("-f","--header", dest="header",  action="store_true", help="print header default:False", default=False)
parser.add_argument("-w","--win_size", dest="winsize", help="window size (int, default:0 (no windows))", default=0)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbouse state on stderr", default=False)
parser.add_argument("-r","--rec", dest="rec", help="list with two thresholds for high,low and medium recombination columns eg. \"1.68,3.21\" (default: False)", default=False)

args = parser.parse_args()
infile = vars(args)['infile']
winsize = int(vars(args)['winsize'])
verb = vars(args)['verb']
header = vars(args)['header']
rec=vars(args)['rec']
if rec:
    rec=[ float(x) for x in rec.split(",") ]

inf=open(infile,"r")
    
pol_dict=defaultdict(lambda: defaultdict( lambda: defaultdict( lambda: (defaultdict(int)))))

for line in inf:
    fields=line.split()
    chr=fields[0]
    if rec:
        rec_rate= float(fields[-1])
        if rec_rate == 0 :
            win="zero"
        elif rec_rate < rec[0]:
            win="low"
        elif rec_rate > rec[1]:
            win = "high"
        elif  rec_rate >= rec[0] and rec_rate <= rec[1]:
            win = "medium"
        else:
            print "rec rate is fishy :"+str(rec_rate)
            sys.exit()
    else:
        bps=int(fields[1])
        if winsize:
            win=(np.floor(bps/winsize)+0.5)*winsize
        else:
            win=0
    pol_sim=fields[5]
    a_state=fields[6]
    pol_dict[chr][win][pol_sim][a_state] += 1

val_dict={ '1': {'00': 0, '10':1, '01': 2} , '2':  {'00': 3, '10':4, '01': 5, '20':4, '02': 5, '12':6, '21': 6}   }
if header:
    if rec:
        print "CHR\tRecRate\tm0\tm10\tm01\tp00\tp10\tp01\tp12"
    else:
        print "CHR\tWin\tm0\tm10\tm01\tp00\tp10\tp01\tp12"
for chr in sorted(pol_dict.keys()):
    for win in sorted(pol_dict[chr].keys()):
        vals = np.zeros(7,dtype=int)
        for pol_sim in sorted(pol_dict[chr][win].keys()):
            for a_state in sorted(pol_dict[chr][win][pol_sim].keys()):
                vals[ val_dict[pol_sim][a_state]  ] += pol_dict[chr][win][pol_sim][a_state]             
        print "{}\t{}\t{}".format(chr,win,"\t".join([ str(x) for x in vals]))
