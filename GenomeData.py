#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng, Dustin E Schones and Keji Zhao
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

"""
This module contains classes of genome data, e.g. chromsomes
per species, the size of the chromosomes, etc.
"""

import re, os, sys, shutil

GenomeDataError = "Error in GenomeData class";


test_chroms = ['chr6']


mm9_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chrX', 'chrY']

mm10_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
             'chr18','chr19','chrX', 'chrY']


hg19_chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
	     'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
	     'chr18','chr19','chr20','chr21','chr22','chrX','chrY']


mm9_chrom_lengths = {'chr1':197195432, 'chr2':181748087, 'chr3':159599783, 'chr4':155630120, 'chr5':152537259, 'chr6':149517037,
                     'chr7':152524553, 'chr8':131738871, 'chr9':124076172,
                     'chr10':129993255, 'chr11':121843856, 'chr12':121257530,
                     'chr13':120284312, 'chr14':125194864, 'chr15':103494974,
                     'chr16':98319150, 'chr17':95272651, 'chr18':90772031,
                     'chr19':61342430, 'chrX':166650296, 'chrY':15902555}

mm10_chrom_lengths = {'chr1':195471971, 'chr2':182113224, 'chr3':160039680,
                     'chr4':156508116, 'chr5':151834684, 'chr6':149736546,
                     'chr7':145441459, 'chr8':129401213, 'chr9':124595110,
                     'chr10':130694993, 'chr11':122082543, 'chr12':120129022,
                     'chr13':120421639, 'chr14':124902244, 'chr15':104043685,
                     'chr16':98207768, 'chr17':94987271, 'chr18':90702639,
                     'chr19':61431566, 'chrX':171031299, 'chrY':91744698}


hg19_chrom_lengths = {'chr1':249250621,  'chr2':243199373, 'chr3':198022430,
		      'chr4':191154276, 'chr5':180915260, 'chr6':171115067,
		      'chr7':159138663,  'chr8':146364022, 'chr9':141213431,
		      'chr10':135534747, 'chr11':135006516, 'chr12':133851895,
		      'chr13':115169878, 'chr14':107349540, 'chr15':102531392,
		      'chr16':90354753,  'chr17':81195210,  'chr18':78077248,
		      'chr19':59128983,  'chr20':63025520,  'chr21':48129895,
		      'chr22':51304566,  'chrX':155270560,  'chrY':59373566}

