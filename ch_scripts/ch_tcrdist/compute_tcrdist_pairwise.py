"""
This file was exported from tcr-dist package with some modifications (https://github.com/phbradley/tcr-dist)

MIT License

Copyright (c) 2017 Philip Harlan Bradley and Jeremy Chase Crawford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
## the BLOSUM62 scoring matrix
##
## Amino acid substitution matrices from protein blocks.
## Henikoff S, Henikoff JG.
## Proc Natl Acad Sci U S A. 1992 Nov 15;89(22):10915-9.
## PMID: 1438297

#!/usr/local/bin python3

if __name__ == '__main__' and __package__ is None:
    import __init__ as ch_init
    __package__ = ch_init.get_package()

import os
import sys
print(sys.version)
import argparse
from functools import reduce

from ..ch_tcrdist.ch_functions_multiprocess import *
from ..ch_tcrdist.ch_base import *
from ..ch_tcrdist import ch_read_blosum as blosum
from ..ch_tcrdist import ch_read_files

#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will compute distances between all "
                                             "the tcrs and write out a distances matrix. The order of rows and columns"
                                             "in the distance matrix will match that in the input file."
                                             "It will also compute the rank scores for each tcr with respect to all the "
                                             "epitopes present")

parser.add_argument("-clones_file", type=str, default="sew_imgted.txt",
                    help="input file path.")
parser.add_argument("-distfile_prefix", type=str)
parser.add_argument("-outfile", type=str)
parser.add_argument("-organism", type=str, default='HomoSapiens')
parser.add_argument("-intrasubject_nbrdists", type=bool, default=2)
parser.add_argument("-clones_files", type=str)
parser.add_argument("-epitope_prefixes", type=str)
parser.add_argument("-nbrdist_percentiles", type=str, default='5;10;25',
                    help="nbrdist_percentiles.")
parser.add_argument("-chains", type=str, default='B',
                    help="chains.")
parser.add_argument("-ch_distance_type", type=str, default='default')
parser.add_argument("-matrix", type=str, help='path to matrix')
parser.add_argument("-distance_params", type=str,
                    default='gap_penalty_v_region:4,gap_penalty_cdr3_region:8,weight_cdr3_region:3,align_cdr3s:False,'
                            'trim_cdr3s:True,scale_factor:1.0439137134052388')
parser.add_argument("-epitope_col", type=str, default='antigen.epitope')   #####!!!!!!!
parser.add_argument("-multiproc", type=int, default=2)
parser.add_argument("-vmatrix", type=str, default=None)

#work with args
args = parser.parse_args()
clones_file = args.clones_file
organism = args.organism
intrasubject_nbrdists = args.intrasubject_nbrdists
new_nbrdists = not intrasubject_nbrdists
nbrdist_percentiles = list(map(int, args.nbrdist_percentiles.split(';')))
chains = args.chains
epitope_col = args.epitope_col
multiproc = args.multiproc
vmatrix = args.vmatrix


ch_distance_type = args.ch_distance_type
if 'default' in ch_distance_type:
    sign = 1
    import ch_tcr_distances_default as ch_tcr_distances
else:
    sign = -1 #just little hack for distance computation. Biopython is looking for alignments with higheest score. So I made score < 0 and dist = -score
    import ch_tcr_distances_custom as ch_tcr_distances

if not args.matrix:
    matrix = blosum.bsd4
else:
    if 'shugay' in ch_distance_type:
        matrix = blosum.compute_matrix_shugai(args.matrix)
    else:
        matrix = blosum.compute_cdrmatrix(blosum.read_blosum(args.matrix), sign=sign)

distance_params = ch_tcr_distances.DistanceParams(distance_matrix=matrix, config_string=args.distance_params)

if not args.clones_files:
    clones_files=[clones_file]
else:
    clones_files=args.clones_files

if not args.epitope_prefixes:
    epitope_prefixes = ['']*len(clones_files)
else:
    epitope_prefixes = args.epitope_prefixes

if not args.outfile:
    outfile = '{}_{}_{}_pairwise_distances.tsv'.format(make_dirprefix(clones_files[0], 'dists'), organism, chains)
else:
    outfile = args.outfile

if not args.distfile_prefix:
    distfile_prefix = distfile_prefix = make_dirprefix(clones_files[0], 'dists')
else:
    distfile_prefix = args.distfile_prefix
#

if 'shugay' in ch_distance_type:
    rep_dists = ch_tcr_distances.compute_all_v_region_distances_shugay(vmatrix, organism, distance_params)
else:
    rep_dists = ch_tcr_distances.compute_all_v_region_distances(organism, distance_params)


tcr_col = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha',
           'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta',
           'v.alpha', 'v.beta', epitope_col, 'species']

all_tcrs = ch_read_files.read_tcr(clones_file, organism, chains, epitope_col, tcr_col)

tcr_col_end = ['cdr3.1', 'v.1', '{}.1'.format(epitope_col),
               'cdr3.2', 'v.2', '{}.2'.format(epitope_col),
               'distance']


#sys.exit()

with open(outfile, 'w') as out:
    out.write('{}\n'.format('\t'.join(tcr_col_end)))
    for iti, tcr in all_tcrs.iterrows():
        print(iti)
        for _, ntcr in all_tcrs.iterrows():
            distance = ch_tcr_distances.compute_distance(tcr['tcr_info'], ntcr['tcr_info'],
                                                         chains, distance_params, rep_dists=rep_dists)
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(tcr['cdr3.beta'], tcr['v.beta'], tcr[epitope_col],
                                                            ntcr['cdr3.beta'], ntcr['v.beta'], ntcr[epitope_col],
                                                            distance))

print('Computation of distances is completed')
