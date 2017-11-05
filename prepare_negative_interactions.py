# !/bash
# Author: Yan Kai
# This script is to prepare the positive and negative long-range interactions among CTCF motifs.
# The positive data is from CTCF ChIA-PET experiment
# The negative data is randomly sampled from the motifs involved in loops.
# The negative training data are sampled in a way that they have the same length distribution as the positive interactions.

import sys, re, os
import numpy as np
import bisect
import math
import pandas as pd
import HTSeq
from optparse import OptionParser

def find_motifs_in_loops(motif_list, chia_start, chia_stop, chip_list):

    starts = []
    anchor_motif = 'NaN'
    mleft = bisect.bisect_left(motif_list,chia_start)
    mright = bisect.bisect_right(motif_list, chia_stop)
    num_motif_on_anchor = mright - mleft

    if (num_motif_on_anchor == 1):
        anchor_motif = str(motif_list[mleft])
        starts.append(int(anchor_motif))
    elif (num_motif_on_anchor > 1):
        sleft = bisect.bisect_left(chip_list, chia_start)
        sright = bisect.bisect_right(chip_list, chia_stop)
        num_summit_on_anchor = sright - sleft
        for j in xrange(mleft,mright):
            starts.append(motif_list[j])

        if (num_summit_on_anchor == 1):
            summit = chip_list[sleft]
            dist = []
            for motif_index in xrange(mleft, mright):
                distance = abs(motif_list[motif_index] - summit) # The motif closest to the summit is chosen.
                dist.append(distance)
            if (min(dist) <= 1000):
                anchor_motif = str(motif_list[np.argmin(dist)])
    return anchor_motif





def main(argv):
    parser = OptionParser()
    parser.add_option("-m", "--motif", action="store", type="string", dest="motif", metavar="<file>", help="the CTCF motif file in BED format")
    parser.add_option("-p", "--chipseq", action="store", type="string", dest="chipseq", metavar="<file>", help="the CTCF ChIP-Seq peak file in BED format")
    parser.add_option("-a", "--chiapet", action="store", type="string", dest="chiapet", metavar="<file>", help="the interaction file from ChIA-PET experiment")
    parser.add_option("-c", "--hic", action="store", type="string", dest="hic", metavar="<file>", help="the CTCF motif interactions found by hic data")
    parser.add_option("-o", "--train", action="store", type="string", dest="training", metavar="<file>", help="the resulting file with positive and sampled negative interactions for training")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)


    motif = pd.read_table(opt.motif, header=None)
    chip = pd.read_table(opt.chipseq, header=None)
    chia = pd.read_table(opt.chiapet)
    hic = pd.read_table(opt.hic)
    outfile = open(opt.training,'w')

    # get the hic loops. HiC loops are used a supplementary file to ensure that the randomly generated negative loops are not true loops identified in HiC.
    hic_loops = {}
    for index, row in hic.iterrows():
        chrom = row['chrom']
        anchor1 = row['start1']
        anchor2 = row['start2']
        if chrom not in hic_loops.keys():
            hic_loops[chrom] = []
        hic_loops[chrom].append((anchor1,anchor2))

# The goal is to construct the long-range interactions among CTCF motifs, thus we need to assign CTCF motifs to the
# anchor regions of ChIA-PET. The strategy is the following: If the anchor is marked by only one motif, replace the
# anchor with the motif; If the anchor is marked by more than one motif, but this anchor is marked by only one CTCF
# ChIP-Seq summit and this summit overlaps with one motif, replace the anchor with the motif.
    motif_length = 18
    true_interaction = {} # true_interaction = {'chrXX':[(start1,start2),(start1,start2)...]}
    motif_dic = {} # motif_dic = {'chrXX':[motif_start1, motif_start2,...]}
    motif_pool = {} # the motif pool for generating the negative loops. motif_pool = {'chrom':[start1, start2,...]}
    for index, row in motif.iterrows():
        chrom = row[0]
        if chrom not in motif_dic.keys():
            motif_dic[chrom] = []
            true_interaction[chrom] = []
            motif_pool[chrom] = []
        motif_dic[chrom].append(row[1])

    chip_dic = {} # chip_dic = {'chrXX':[summit1, summit2, summit3...]}

    for index, row in chip.iterrows():
        chrom = row[0]
        if chrom not in chip_dic.keys():
            chip_dic[chrom] = []
        summit = (row[1]+row[2])/2
        chip_dic[chrom].append(summit)



    # Sort the motif list and chip list
    for chrom in motif_dic.keys():
        if chrom in chip_dic.keys():
            motif_dic[chrom].sort()
            chip_dic[chrom].sort()


    outline = 'chrom'+'\t'+'start1'+'\t'+'start2'+'\t'+'response'+'\t'+'length'+'\n'
    outfile.write(outline)
    loop_length = []

    for index, row in chia.iterrows():
        if (row['chrom1'] == row['chrom2'] and row['chrom1'] != 'chrM'):
            chrom = row['chrom1']
            starts1 = []
            starts2 = []

            # Get the motifs responsible for the anchors
            if chrom in chip_dic.keys():
                anchor1_motif = find_motifs_in_loops(motif_dic[chrom], row['start1'], row['stop1'], chip_dic[chrom])
                anchor2_motif = find_motifs_in_loops(motif_dic[chrom], row['start2'], row['stop2'], chip_dic[chrom])
            else:
                anchor1_motif = 'NaN'
                anchor2_motif = 'NaN'
            if anchor1_motif != 'NaN':
                motif_pool[chrom].append(int(anchor1_motif))
            if anchor2_motif != 'NaN':
                motif_pool[chrom].append(int(anchor2_motif))

        # distance is the distance between the two motifs. To focus on long-range interactions,
        # we required that distance >= 10kb and <= 3m
            if (anchor1_motif != 'NaN' and anchor2_motif != 'NaN'):
                if (int(anchor1_motif) > int(anchor2_motif)):
                    temp = anchor1_motif
                    anchor1_motif = anchor2_motif
                    anchor2_motif = temp
                distance = int(anchor2_motif) - int(anchor1_motif)

                if (distance >= 10000 and distance <= 3000000):
                    loop_length.append(distance)
                    true_interaction[chrom].append((int(anchor1_motif), int(anchor2_motif)))
                    outline = chrom+'\t'+anchor1_motif+'\t'+anchor2_motif+'\t'+str(1)+'\t'+str(distance)+'\n'
                    outfile.write(outline)


# Then construct the negative interactions
    N = 100
    bin_length = 3000000/N
    density, bin_edges = np.histogram(loop_length, bins = N, range = (0,3000000),density = True)
    density = density * bin_length

    Ratio = 5 # Ratio = 5 means 5 negative interaction will be generated for each positive interaction.
    NumNeg = len(loop_length)*Ratio # NumNeg is the totoal number of negative interactions.
    #neg_loop_length = np.random.choice(bin_edges[1:N+1],NumNeg, replace=True, p=density)




    num_samples_by_distance = {} # total number of samples in each distance category
    chroms = motif_dic.keys()

    distances = bin_edges[1:N+1]
    NumSample = NumNeg
    negative_interactions = {} # negative_interactions = {'distance':[loop1,loop2...]} loop = HTSeq.Genomic(chrom, start1, start2,'.')
    selected_neg = {} #same as negative_interactions
    for i in xrange((len(distances))):
        distance = distances[i]
        num_samples_by_distance[distance] = math.ceil(NumSample * density[i])
        negative_interactions[distance] = []

    # Generate the negative interactions pool
    total = 0
    for chrom in motif_pool.keys():
        for i_left in xrange(len(motif_pool[chrom])-1):
            m_left = motif_pool[chrom][i_left]
            for i_right in xrange(i_left+1, len(motif_pool[chrom])):
                m_right = motif_pool[chrom][i_right]
                length = m_right - m_left
                if length >= 10000 and length <= 3000000 and (m_left, m_right) not in true_interaction[chrom] and (m_left, m_right) not in hic_loops[chrom]:
                    distance = bin_edges[bisect.bisect_right(bin_edges[1:N+1],length)+1]
                    iv = HTSeq.GenomicInterval(chrom, m_left, m_right, '.')
                    negative_interactions[distance].append(iv)
                    total += 1
                    outline = chrom+'\t'+str(m_left)+'\t'+str(m_right)+'\t'+str(0)+'\t'+str(length)+'\n'
                    outfile.write(outline)
    print 'There are '+str(total)+' negative loops in total'




    outfile.close()

if __name__ == "__main__":
	main(sys.argv)
