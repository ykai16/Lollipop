# Author: Yan Kai
# Date: May 29th, 2017
# This script is to prepare the high-confidence positive and negative long-range interactions among CTCF binding sites, as revealed from ChIP-Seq or ChIA-PET.
# The input to generate positive data is the loops determined by CTCF ChIA-PET experiment.
# The negative data are randomly sampled from CTCF binding sites.

import sys, re, os
import GenomeData
import numpy as np
import bisect
import math
import pandas as pd
import HTSeq
from optparse import OptionParser

def find_summits_in_anchors(anchor, chrom_summits):
    """
    anchor = HTSeq.iv
    chrom_summits = []  # list of summit position on one chrom
    """
    chrom = anchor.chrom
    overlapped = 0
    for summit in chrom_summits:
        pos = HTSeq.GenomicPosition(chrom, summit,'.')
        if pos.overlaps(anchor):
            ans = summit
            overlapped = 1
            break
    if overlapped == 0:
        ans = (anchor.start + anchor.end)/2
    return str(ans)

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--peak", action="store", type="string", dest="peak", metavar="<file>", help="the CTCF peaks or summits in BED format")
    parser.add_option("-a", "--chiapet", action="store", type="string", dest="chiapet", metavar="<file>", help="the interaction file from ChIA-PET experiment")
    parser.add_option("-c", "--hic", action="store", type="string", dest="hic", metavar="<file>", help="the CTCF interactions identified by hic data")
    parser.add_option("-o", "--train", action="store", type="string", dest="training", metavar="<file>", help="the resulting file with positive and sampled negative interactions for training")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)


    peak = pd.read_table(opt.peak, header=None)
    chia = pd.read_table(opt.chiapet)
    hic = pd.read_table(opt.hic)
    chroms = GenomeData.hg19_chroms
    outfile = open(opt.training,'w')

    # The length constraint of loops
    minLength = 10000
    maxLength = 1000000

    true_loops = {} #true_interaction = {'chrXX':[(summit1,summit2),(summit1,summit2)...]}
    less_sig_loops = {} # To ensure a negative loops is not a less significant true loop.less_sig_loops = {'chrXX':[(summit1,summit2),(summit1,summit2)...]}
    bs_pool = {} # the binding sites pool for generating the negative loops. bs_pool = {'chrom':[summit1, summit2,...]}

    for index, row in peak.iterrows():
        chrom = row[0]
        if chrom in chroms and chrom != 'chrY':
            if chrom not in bs_pool.keys():
                bs_pool[chrom] = []
            summit = (row[1]+row[2])/2
            bs_pool[chrom].append(summit)



    # Sort the binding site list
    for chrom in bs_pool.keys():
        bs_pool[chrom].sort()

    # get the hic loops. HiC loops are used as a supplementary file to ensure that the randomly generated negative loops are not true loops identified in HiC.
    hic_loops = {}
    for index, row in hic.iterrows():
        chrom = row['chrom']
        anchor1 = HTSeq.GenomicInterval(chrom, row['start1']-1000, row['start1']+1000, '.')
        anchor2 = HTSeq.GenomicInterval(chrom, row['start2']-1000, row['start2']+1000,'.')
        if chrom not in hic_loops.keys():
            hic_loops[chrom] = []
        anchor1_summit = find_summits_in_anchors(anchor1, bs_pool[chrom])
        anchor2_summit = find_summits_in_anchors(anchor2, bs_pool[chrom])
        hic_loops[chrom].append((anchor1_summit,anchor2_summit))

    outline = 'chrom'+'\t'+'start1'+'\t'+'start2'+'\t'+'response'+'\t'+'length'+'\n'
    outfile.write(outline)
    loop_length = []

    for index, row in chia.iterrows():
        IAB = row['IAB']
        FDR = row['FDR']
        if (row['chrom1'] == row['chrom2'] and row['chrom1'] in chroms and row['chrom1'] != 'chrY'):
            chrom = row['chrom1']
            anchor1 = HTSeq.GenomicInterval(chrom, row['start1'], row['end1'],'.')
            anchor2 = HTSeq.GenomicInterval(chrom, row['start2'], row['end2'],'.')
            # Get the summit position of the anchors
            if chrom in bs_pool.keys():
                anchor1_summit = find_summits_in_anchors(anchor1, bs_pool[chrom])
                anchor2_summit = find_summits_in_anchors(anchor2, bs_pool[chrom])
            else:
                anchor1_summit = 'NaN'
                anchor2_summit = 'NaN'

        # distance is the genomic length between the two motifs. To focus on long-range interactions,
        # we required that distance >= 10kb and <= 1m
            if (anchor1_summit != 'NaN' and anchor2_summit != 'NaN'):
                if (int(anchor1_summit) > int(anchor2_summit)):
                    temp = anchor1_summit
                    anchor1_summit = anchor2_summit
                    anchor2_summit = temp
                distance = int(anchor2_summit) - int(anchor1_summit)

                if (distance >= minLength and distance <= maxLength):
                    if (IAB >= 2 and FDR <= 0.05):
                        loop_length.append(distance)
                        if chrom not in true_loops.keys():
                            true_loops[chrom] = []

                        true_loops[chrom].append((int(anchor1_summit), int(anchor2_summit)))
                        outline = chrom+'\t'+anchor1_summit+'\t'+anchor2_summit+'\t'+str(1)+'\t'+str(distance)+'\n'
                        outfile.write(outline)
                    else:
                        if chrom not in less_sig_loops.keys():
                            less_sig_loops[chrom] = []
                        less_sig_loops[chrom].append((int(anchor1_summit), int(anchor2_summit)))

    Ratio = 5 # Ratio = 5 means 5 negative interaction will be generated for each positive interaction.
    NumNeg = len(loop_length)*Ratio # NumNeg is the totoal number of negative interactions.


    negative_interactions = []
    selected_neg = []

    # Generate the negative interactions pool
    total = 0
    for chrom in true_loops.keys():
        for i_left in xrange(len(bs_pool[chrom])-1):
            m_left = bs_pool[chrom][i_left]
            for i_right in xrange(i_left+1, len(bs_pool[chrom])):
                m_right = bs_pool[chrom][i_right]
                length = m_right - m_left
                if length >= minLength and length <= maxLength and (m_left, m_right) not in true_loops[chrom]:
                    if chrom in less_sig_loops.keys():
                        if (m_left, m_right) not in less_sig_loops[chrom]:
                            iv = HTSeq.GenomicInterval(chrom, m_left, m_right, '.')
                            negative_interactions.append(iv)
                            total += 1
                    else:
                        iv = HTSeq.GenomicInterval(chrom, m_left, m_right, '.')
                        negative_interactions.append(iv)
                        total += 1
    print 'There are '+str(total)+' negative loops in total'



    selected_neg = np.random.choice(negative_interactions, NumNeg, replace=False)
    for iv in selected_neg:
        length = iv.end - iv.start
        outline = iv.chrom+'\t'+str(iv.start)+'\t'+str(iv.end)+'\t'+str(0)+'\t'+str(length)+'\n'
        outfile.write(outline)


    outfile.close()

if __name__ == "__main__":
	main(sys.argv)
