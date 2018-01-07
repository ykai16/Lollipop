import sys, re, os
import numpy as np
import bisect
import collections
import pandas as pd
import HTSeq
import operator
import GenomeData

def prepare_interactions(data, chrom, i, start_list, distance_distal, distance_proximal):
    """
    This function is to prepare the potential interactions for one motif with all the downstream motifs within a certain
    range.
    i: the index of the particular summits in the motif list
    start_list: the list of all the summits in one chrom,it's sorted
    distance_distal/proximal: only the interactions whose lengh is in the range of [proximal, distal] will be kept.

    Returned: a pandas data frame with 4 columns:
    chrom+start1+start2+length
    All the chroms and start1 will be identifcal with each other in this data frame.

    """

    data1 = pd.DataFrame()
    start1 = start_list[i]
    chromosome = []
    left = []
    right = []
    length = []
    for j in xrange(i+1,len(start_list)):
        start2 = start_list[j]
        interval = start2 - start1
        if (interval >= distance_proximal and interval <=  distance_distal):
            chromosome.append(chrom)
            left.append(int(start1))
            right.append(int(start2))
            length.append(int(interval))
        else:
            break

    data1['chrom'] = pd.Series(chromosome)
    data1['peak1'] = pd.Series(left)
    data1['peak2'] = pd.Series(right)
    data1['length'] = pd.Series(length)

    data = pd.concat([data,data1],ignore_index=True)
    return data



def prepare_reads_info(signal_table):
    """
    This function is to prepare the reads info from raw .bed files for the local features.

    Returned: read_info = {factor:{chrom:[start1, start2, start3]}}

    """
    read_info = {}  # read_info = {factor:{chrom:[start1, start2, start3]}}
    read_numbers = {} # read_numbers = {'H3K4me1':read_number, ...}
    shift = 75   # half of the fragment size
    for index, row in signal_table.iterrows():
        factor = row['Signal']
        BED = row['BED']
        if (factor  != 'Motif' and factor != 'Gene expression' and factor != 'PhastCon'):
            BED_reader = open(BED,'r')
            read_info[factor] = {}
            read_number = 0
            for line in BED_reader:
                read_number += 1
                pline = line.strip()
                sline = pline.split('\t')
                chrom = sline[0]
                start = int(sline[1])
                end = int(sline[2])
                strand = sline[5]
                if chrom not in read_info[factor].keys():
                    read_info[factor][chrom] = []

                if (strand == '+'):
                    start = start + shift
                else:
                    start = end - shift

                read_info[factor][chrom].append(start)
            read_numbers[factor] = read_number
            for chrom in read_info[factor].keys():
                read_info[factor][chrom] = sorted(read_info[factor][chrom])
            BED_reader.close()

    return (read_info, read_numbers)



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
    Return the data with added motif pattern feature and motif strength
    """
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
        start1 = row['peak1']
        ext = 500
        anchor1 = HTSeq.GenomicInterval(chrom, start1-ext, start1+ext, '.')
        start2 = row['peak2']
        anchor2 = HTSeq.GenomicInterval(chrom, start2-ext, start2+ext, '.')

        (pattern, avg_motif_strength, std_motif_strength) = find_motif_pattern(anchor1, anchor2, motif)
        motif_pattern.append(pattern)
        avg_motif_strength_list.append(avg_motif_strength)
        std_motif_strength_list.append(std_motif_strength)

    train['motif_pattern'] = pd.Series(motif_pattern, index = train.index)
    train['avg_motif_strength'] = pd.Series(avg_motif_strength_list, index = train.index)
    train['std_motif_strength'] = pd.Series(std_motif_strength_list, index = train.index)
    return train


def add_anchor_conservation(train, chrom, Peak):
    """
    To add the feature of sequence conservation on anchors.
    Peak is the folder that contains the phastCon.wig files of all chroms.
    The name format of PhastCon of each chrom is chrXX.phastCons100way.wigFix
    """
    starts = {}
    cvg = {}

    ext = 20
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
    for index, row in train.iterrows():
        con1 = sum(cvg[(int(row['peak1'])-ext): (int(row['peak1'])+ext)])
        con2 = sum(cvg[(int(row['peak2'])-ext): (int(row['peak2'])+ext)])
        AvgCons.append((con1+con2)/2.0)
        DevCons.append(np.std([con1, con2]))
    train['avg_conservation'] = pd.Series(AvgCons)
    train['std_conservation'] = pd.Series(DevCons)

    return train


def add_features(data, anchor_motifs, read_info, read_numbers, signals):
    """
    This function is to add both the local and inbetween features to the data.
    read_info = {factor:{chrom:[start1, start2, start3]}}
    inbetween_signals = {factor:{'ChrXX':{summit:peak_height}}}
    anchor_motifs = {'chrXX':[start1, start2...]}
    """
    extension = 2000

    for factor in signals:
        print "Preparing features for "+str(factor)+'...'
        avg_signal = 'avg_'+str(factor)
        std_signal = 'std_'+str(factor)
        inbetween_signal = str(factor)+'_in-between';upstream_signal = str(factor)+'_upstream';downstream_signal = str(factor)+'_downstream'

        reads = read_info[factor]
        read_number = read_numbers[factor]
        anchor1_RPKM = []
        anchor2_RPKM = []
        in_between = [];upstream = [];downstream = []

        for index, row in data.iterrows():
            chrom = row['chrom']
            start1 = int(row['peak1'])
            start2 = int(row['peak2'])

            # Get the RPKM read counts on anchors
            count1 = bisect.bisect_right(reads[chrom], start1+extension) - bisect.bisect_left(reads[chrom], start1-extension)
            count2 = bisect.bisect_right(reads[chrom], start2+extension) - bisect.bisect_left(reads[chrom], start2-extension)
            RPKM1 = float(count1)/(float(read_number)*2*extension)*1000000000
            RPKM2 = float(count2)/(float(read_number)*2*extension)*1000000000
            anchor1_RPKM.append(np.mean([RPKM1, RPKM2]))
            anchor2_RPKM.append(np.std([RPKM1, RPKM2]))

            # Get the RPKM values of the looped regions
            strength = 0
            count = bisect.bisect_right(reads[chrom], start2) - bisect.bisect_left(reads[chrom], start1)
            strength = float(count)/float(abs(start2 - start1)*read_number)*1e+9
            in_between.append(strength)
            index1 = anchor_motifs[chrom].index(start1)
            index2 = anchor_motifs[chrom].index(start2)
            if index1 != '0':
                up_motif = anchor_motifs[chrom][index1 - 1]
                up_count = bisect.bisect_right(reads[chrom],start1) - bisect.bisect_left(reads[chrom], up_motif)
                up_strength = float(up_count)/float(abs(up_motif-start1)*read_number)*1e+9
            else:
                up_strength = 0
            upstream.append(up_strength)
            if index2 != (len(anchor_motifs[chrom])-1):
                down_motif = anchor_motifs[chrom][index2 + 1]
                down_count = bisect.bisect_right(reads[chrom], down_motif) - bisect.bisect_left(reads[chrom], start2)
                down_strength = float(down_count)/float(abs(down_motif-start2)*read_number)*1e+9
            else:
                down_strength = 0
            downstream.append(down_strength)

        data[avg_signal] = pd.Series(anchor1_RPKM, index = data.index)
        data[std_signal] = pd.Series(anchor2_RPKM, index = data.index)
        data[inbetween_signal] = pd.Series(in_between, index = data.index)
        data[upstream_signal] = pd.Series(upstream, index = data.index)
        data[downstream_signal] = pd.Series(downstream, index = data.index)
    return data

def add_gene_expression(data, signal_table):
    """
    This function is to add the gene expression value of the looped region as a feature.The gene expression file's format is:
    gene_id   locus   value
    A1BG    chr19:coordinate1-coordiate2   1.31

    """
    for index,row in signal_table.iterrows():
        if row['Signal'] == 'Gene expression':
            print 'Preparing features for gene expression...'
            exp_file = pd.read_table(row['Peak'])
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
            for index, row in data.iterrows():
                chrom = row['chrom']
                start1 = row['peak1']
                start2 = row['peak2']
                iv = HTSeq.GenomicInterval(chrom, start1, start2)
                loop_expression = 0
                for gene in gene_exp[chrom].keys():
                    if gene.overlaps(iv):
                        loop_expression += gene_exp[chrom][gene]
                loop_expressions.append(loop_expression)
            data['expression'] = pd.Series(loop_expressions, index = data.index)
    return data

def prepare_features_for_interactions(data, summits, signal_table, raw_features):
    """

    data is a pandas dataframe with chrom+start1+start2+length.
    raw_features is a tuple of the motif information and coverage vectors for all factors.
    raw_features = (signals, read_info, read_numbers)
    motifs_info = {'chromXX':{start:(strand, score, phastCon)}}
    signals is a list of available factors with the same order as in the signal_table
    """

    signals = raw_features[0]
    read_info = raw_features[1]
    read_numbers = raw_features[2]

    data = add_features(data, summits, read_info, read_numbers, signals)
    data = add_gene_expression(data, signal_table)


    return data


