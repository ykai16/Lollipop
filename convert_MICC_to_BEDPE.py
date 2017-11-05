# !/bash
# Author: Yan Kai
# Date: May 29th, 2017
# This script is to prepare the loops in BEDPE format from loop file from MICC.
# MICC format: chr	start	end	chr	start	end	peakA	peakB	cA	cB	cAB	-log10(1-PostProb)	fdr
# BEDPE format: chrom1    start1    end1    chrom2    start2    end2    IAB    fdr    strand1    strand2

import sys, re, os
import numpy as np
import math
import pandas as pd
import HTSeq
import GenomeData
from optparse import OptionParser

def main(argv):
    parser = OptionParser()
    parser.add_option("-i", "--loops", action="store", type="string", dest="loops", metavar="<file>", help="loops identified from MICC")
    parser.add_option("-c", "--count", action="store", type="int", dest='cutoff', metavar="<int>", help="the cutoff to select loops: only loops with PETs above the given cutoff would be chosen")
    parser.add_option("-f", "--fdr", action="store", type="float", dest='fdr', metavar="<float>", help="the fdr to select loops: only loops with fdr below the given threshould would be chosen")
    parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="loops in BEDPE format following the criteria")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)

    micc = pd.read_table(opt.loops)
    PETcut = opt.cutoff
    FDRcut = opt.fdr
    hg19_chroms = GenomeData.hg19_chroms
    outfile = open(opt.output, 'w')
    outline = 'chrom1'+'\t'+'start1'+'\t'+'end1'+'\t'+'chrom2'+'\t'+'start2'+'\t'+'end2'+'\t'+'IAB'+'\t'+'FDR'+'\t'+'strand1'+'\t'+'strand2'+'\n'
    outfile.write(outline)


    for index, row in micc.iterrows():
        chrom1 = row[0]
        chrom2 = row[3]
        Iab = row[10]
        FDR = row[12]
        if chrom1 == chrom2 and chrom1 in hg19_chroms and Iab >= PETcut and FDR <= FDRcut:
            outline = chrom1+'\t'+str(row[1])+'\t'+str(row[2])+'\t'+chrom2+'\t'+str(row[4])+'\t'+str(row[5])+'\t'+str(Iab)+'\t'+str(FDR)+'\t'+'.'+'\t'+'.'+'\n'
            outfile.write(outline)
    outfile.close()

if __name__ == "__main__":
	main(sys.argv)
