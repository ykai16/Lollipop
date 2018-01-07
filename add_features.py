# !/bash
# Author: Yan Kai
# This script is to prepare features for the generated loops for training purpose.
# Two inputs: 1. The generated positive and negative loops in the format: chrom+start1+start2+length+response
#             2. The infomation table that contains the complete path for necessary raw BED file and peak file.
# Output is: The training data, i.e the interactions plus all the listed features.
import sys, re, os
import numpy as np
import bisect
import GenomeData
import collections
import time
from optparse import OptionParser
import pandas as pd
import HTSeq

def assign_motif_pattern(strand1, strand2):
    if (strand1 == '-' and strand2 == '+'):
        return '2'
    elif (strand1 == '+' and strand2 == '-'):
        return '4'
    else:
        return '3'


def find_motif_pattern(anchor1, anchor2, motif):
    """
     Input:
         anchor = HTSeq.GenomicInterval(chrom,summit-ext, summit+ext,'.')
         motif = {'chromXX':{start:(strand, score)}}
     Output:
         a tuple (pattern, avg_motif_strength, std_motif_strength)
     Rules to assign motif pattern:
     1. Both anchors have no motif, assign 0;
     2. One anchor has no motif, no matter how many motifs the other anchor may have, assign 1;
     3. Both anchors have 1 motif: no ambuguity, divergent=2;tandem=3; convergent=4
     4. Anchors have multiple motifs: in each anchor, choose the one with the highest motif strength
    """
    chrom = anchor1.chrom

    starts = list(motif[chrom].keys())

    num1 = 0
    starts1 = []
    scores1 = []
    for start in starts:
        if start >= anchor1.start and start <= anchor1.end:
            num1 += 1
            starts1.append(start)
            scores1.append(motif[chrom][start][1])
    num2 = 0
    starts2 = []
    scores2 = []
    for start in starts:
        if start >= anchor2.start and start <= anchor2.end:
            num2 += 1
            starts2.append(start)
            scores2.append(motif[chrom][start][1])

    if num1 ==0 and num2 == 0:
        return (0,0,0)
    else:
        if (num1*num2 == 0):
            if num1 == 0:
                return (1, max(scores2)/2.0, np.std([0, max(scores2)]))
            else:
                return (1, max(scores1)/2.0, np.std([0, max(scores1)]))
        else:
            if (num1 == 1 and num2 == 1):
                strand1 = motif[chrom][starts1[0]][0]
                strand2 = motif[chrom][starts2[0]][0]
                pattern = assign_motif_pattern(strand1, strand2)
                return (pattern, np.mean([max(scores1),max(scores2)]), np.std([max(scores1),max(scores2)]))
            else:
                index1 = scores1.index(max(scores1))
                strand1 = motif[chrom][starts1[index1]][0]
                index2 = scores2.index(max(scores2))
                strand2 = motif[chrom][starts2[index2]][0]
                pattern = assign_motif_pattern(strand1, strand2)
                return (pattern, np.mean([max(scores1),max(scores2)]), np.std([max(scores1),max(scores2)]))


def add_motif_pattern(train, Peak):
    """
    This function is to add the motif pattern feature for interacting anchors in training data.
    Peak is the complete path for the table containing all these information for all motifs. The format is:
    chrom + start + end + strand + pvalue + score + phastCon
    Return a training data with added motif pattern feature
    """
    print "Preparing the motif pattern feature..."
    info = pd.read_table(Peak)
    motif = {}  # motif = {'chromXX':{start:(strand, score)}}
    for index, row in info.iterrows():
        chrom = row['chrom']
        start = row['start']
        strand = row['strand']
        score = row['score']

        if chrom not in motif.keys():
            motif[chrom] = {}
        motif[chrom][start] = (strand, score)

    motif_pattern = []
    avg_motif_strength_list = []
    std_motif_strength_list = []

    for index, row in train.iterrows():
        chrom = row['chrom']
        start1 = row['start1']
        ext = 500
        anchor1 = HTSeq.GenomicInterval(chrom, start1-ext, start1+ext, '.')
        start2 = row['start2']
        anchor2 = HTSeq.GenomicInterval(chrom, start2-ext, start2+ext, '.')

        (pattern, avg_motif_strength, std_motif_strength) = find_motif_pattern(anchor1, anchor2, motif)
        motif_pattern.append(pattern)
        avg_motif_strength_list.append(avg_motif_strength)
        std_motif_strength_list.append(std_motif_strength)

    train['motif_pattern'] = pd.Series(motif_pattern, index = train.index)
    train['avg_motif_strength'] = pd.Series(avg_motif_strength_list, index = train.index)
    train['std_motif_strength'] = pd.Series(std_motif_strength_list, index = train.index)
    return train


def add_anchor_conservation(train, chroms, Peak):
    """
    To add the feature of sequence conservation on anchors.
    Peak is the folder that contains the phastCon.wig files of all chroms.
    The name format of PhastCon of each chrom is chrXX.phastCons100way.wigFix
    """
    print 'Preparing the sequence conservation of anchors...'
    starts = {}
    cvg = {}

    ext = 20
    training = [] # [chrom_train_DF, chrom2_train_DF, ...]
    for chrom in chroms:
        chrom_train = train[train['chrom'] == chrom].copy()
        chrom_train.reset_index(inplace=True)
        print 'Read in phastCon in '+chrom+'...'

        # Read in the phastCon track
        cvg = [0]*GenomeData.hg19_chrom_lengths[chrom]
        phastCon = Peak+'/'+chrom+'.phastCons100way.wigFix'
        wiggle = open(phastCon,'r')
        for line in wiggle:
            if line[0] == 'f':
                i = 0
                start = int(line.strip().split(' ')[2].split('=')[1])
            else:
                signal = line.strip().split(' ')[0]
                if signal == 'NA':
                    signal = 0
                else:
                    signal = float(signal)
                cvg[start + i] = signal
                i += 1
        wiggle.close()

        AvgCons = []
        DevCons = []
        for index, row in chrom_train.iterrows():
            con1 = sum(cvg[(row['start1']-ext): (row['start1']+ext)])
            con2 = sum(cvg[(row['start2']-ext): (row['start2']+ext)])
            AvgCons.append((con1+con2)/2.0)
            DevCons.append(np.std([con1, con2]))
        chrom_train['avg_conservation'] = pd.Series(AvgCons)
        chrom_train['std_conservation'] = pd.Series(DevCons)
        training.append(chrom_train)

    new_train = pd.concat(training, ignore_index=True)

    return new_train



def add_local_feature(signal, train, BED):
    """
    This function is to calculate the signal values on summit positions (summit +/- 2kb)
    """
    print "Preparing the local features of "+str(signal)+'...'
    extension = 2000
    fragment = 150 # This is the size of the ChIP fragment, usually it is 150.
    shift = fragment/2
    BED_reader = open(BED,'r')
    read_info = {}  # read_info = {'chrom':[start1, start2,...]}  actually 'start' here is the mid-point of a fragment

    read_number = 0
    for line in BED_reader:
        read_number += 1
        pline = line.strip()
        sline = pline.split('\t')
        chrom = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        strand = sline[5]

        if chrom not in read_info.keys():
            read_info[chrom] = []
        if (strand == '+'):
            start = start + shift
        else:
            start = end - shift

        read_info[chrom].append(start)
    BED_reader.close()

    # make sure reads are sorted
    for chrom in read_info.keys():
        read_info[chrom] = sorted(read_info[chrom])
        
    RPKMs1 = []
    RPKMs2 = []
    for index, row in train.iterrows():
        chrom = row['chrom']
        start1 = row['start1']
        start2 = row['start2']

        count1 = bisect.bisect_right(read_info[chrom], start1+extension) - bisect.bisect_left(read_info[chrom], start1-extension)
        count2 = bisect.bisect_right(read_info[chrom], start2+extension) - bisect.bisect_left(read_info[chrom], start2-extension)

        RPKM1 = float(count1)/(float(read_number)*2*extension)*1000000000
        RPKM2 = float(count2)/(float(read_number)*2*extension)*1000000000

        RPKMs1.append((RPKM1+RPKM2)/2.0)
        RPKMs2.append(np.std([RPKM1, RPKM2]))

    signal1 = 'avg_'+str(signal)
    signal2 = 'std_'+str(signal)
    train[signal1] = pd.Series(RPKMs1, index = train.index)
    train[signal2] = pd.Series(RPKMs2, index = train.index)
    return train

def add_gene_expression(train, Peak):
    """
    This function is to add the gene expression value of the looped region as a feature.The gene expression file's format is:
    gene_id   locus   value
    A1BG    chr19:coordinate1-coordiate2   1.31

    """
    exp_file = pd.read_table(Peak)
    gene_exp = {}  # {'chrom':{iv1:fpkm1,iv2:fpkm2...}}
    for index, row in exp_file.iterrows():
        gene = row['gene_id']
        region = row['locus']
        fpkm = row['value']
        chrom = region.split(':')[0]
        start = int(region.split(':')[1].split('-')[0])
        end = int(region.split(':')[1].split('-')[1])
        iv = HTSeq.GenomicInterval(chrom, start, end, '.')
        if chrom not in gene_exp.keys():
            gene_exp[chrom] = {}
        gene_exp[chrom][iv] = fpkm

    loop_expressions = []
    for index, row in train.iterrows():
        chrom = row['chrom']
        start1 = row['start1']
        start2 = row['start2']
        iv = HTSeq.GenomicInterval(chrom, start1, start2)
        loop_expression = 0
        for gene in gene_exp[chrom].keys():
            if gene.overlaps(iv):
                loop_expression += gene_exp[chrom][gene]
        loop_expressions.append(loop_expression)
    train['expression'] = pd.Series(loop_expressions, index = train.index)
    return train


def add_regional_feature_by_reads(signal, train, anchors, BED):
    """
    This function is to add the in-between and loop-flanking features for an interacting pair from raw reads BED files, whose format is:
    chrom+start+end.

    """
    print "Preparing the in-between and loop-flanking features of "+str(signal)+'...'
    BED = open(BED, 'r')
    shift = 75
    signal_dic = {}  # signal_dic = {"chrXX":[start1, start2, ...]} 'start' here are mid-point of one fragment
    read_number = 0
    for line in BED:
        read_number += 1
        pline = line.strip()
        sline = pline.split('\t')
        chrom = sline[0]
        start = int(sline[1])
        end = int(sline[2])
        strand = sline[5]

        if chrom not in signal_dic.keys():
            signal_dic[chrom] = []
        if strand == '+':
            start = start + shift
        else:
            start = end - shift
        signal_dic[chrom].append(start)
    BED.close()
    
    # make sure reads are sorted
    for chrom in signal_dic.keys():
        signal_dic[chrom] = sorted(signal_dic[chrom])

    in_between = []
    upstream = []
    downstream = []
    for index, row in train.iterrows():
        chrom = row['chrom']
        start1 = row['start1']
        start2 = row['start2']

        index1 = anchors[chrom].index(start1)
        index2 = anchors[chrom].index(start2)
        if index1 != 0:
            up_motif = anchors[chrom][index1 - 1]
            up_count = bisect.bisect_right(signal_dic[chrom], start1) - bisect.bisect_left(signal_dic[chrom], up_motif)
            up_strength = float(up_count)/float(abs(up_motif-start1)*read_number)*1e+9
        else:
            up_strength = 0
        upstream.append(up_strength)
        if index2 != (len(anchors[chrom])-1):
            down_motif = anchors[chrom][index2 + 1]
            down_count = bisect.bisect_right(signal_dic[chrom], down_motif) - bisect.bisect_left(signal_dic[chrom], start2)
            down_strength = float(down_count)/float(abs(down_motif-start2)*read_number)*1e+9
        else:
            down_strength = 0
        downstream.append(down_strength)

        strength = 0
        count = bisect.bisect_right(signal_dic[chrom], start2) - bisect.bisect_left(signal_dic[chrom], start1)

        strength = float(count)/float(abs(start2-start1)*read_number)*1e+9
        in_between.append(strength)

    in_between_signal = signal+'_in-between'
    train[in_between_signal] = pd.Series(in_between, index=train.index)

    upstream_signal = signal+'_left'
    train[upstream_signal] = pd.Series(upstream, index=train.index)

    downstream_signal = signal+'_right'
    train[downstream_signal] = pd.Series(downstream, index=train.index)

    return train




def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--training", action="store", type="string", dest="training", metavar="<file>", help="the training interactions generated in previous step")
	parser.add_option("-t", "--table", action="store", type="string", dest="info_table", metavar="<file>", help="the infomation table contains the paths for necessary files.")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="the output file with all the interactions and calculated features.")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
        	sys.exit(1)

	start_time = time.time()
	chroms = GenomeData.hg19_chroms
	#chroms = ['chr1']
	train = pd.read_table(opt.training)
	train = train.sort_values(by=['chrom','start1','start2'], axis = 0, ascending=[1,1,1])
	info_table = pd.read_table(opt.info_table)
	anchors = {} # anchors = {'chr':set(summit1, summit2,)}
	
    # Generate the anchors pool from the anchors of positive loops
	for index,row in train.iterrows():
	        chrom = row['chrom']
	        if chrom not in anchors.keys():
	            anchors[chrom] = set()
	        anchors[chrom].add(row['start1'])
	        anchors[chrom].add(row['start2'])
	for chrom in anchors.keys():
		anchors[chrom] = list(anchors[chrom])
		anchors[chrom].sort()

	for index, row in info_table.iterrows():
		signal = row['Signal']
		BED = row['BED']
		Peak = row['Peak']

		if (signal == 'Motif'):
			train = add_motif_pattern(train, Peak)
		
		
 	for index, row in info_table.iterrows():
		signal = row['Signal']
		BED = row['BED']
		Peak = row['Peak']

		if (signal == 'PhastCon'):
		    train = add_anchor_conservation(train, chroms, Peak)
		elif (signal == 'Gene expression'):
			train = add_gene_expression(train, Peak)
		else:
			if (BED != 'No'):
				train = add_local_feature(signal, train, BED)
				train = add_regional_feature_by_reads(signal, train, anchors, BED)



	train.to_csv(opt.outfile,sep='\t', index=False)
	end_time = time.time()
	elapsed = end_time-start_time
	print "Time elapased: "+str(elapsed)+'seconds'


if __name__ == "__main__":
	main(sys.argv)
