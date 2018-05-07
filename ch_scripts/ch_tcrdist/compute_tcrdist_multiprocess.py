#!/usr/local/bin python3
import time
ts = time.time()
if __name__ == '__main__' and __package__ is None:
    import __init__ as ch_init
    __package__ = ch_init.get_package()

import os
import sys
print(sys.version)
import argparse
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

parser.add_argument("-clones_file", type=str, default="cmvtet/cmvtet_imgted_b.txt",
                    help="input file path.")
parser.add_argument("-db_file", type=str, default="cmvtet/vdjdb_imgted.txt",
                    help="input database file path.")
parser.add_argument("-distfile_prefix", type=str)
parser.add_argument("-outfile", type=str, default="custom_computation_cmvtet_b.txt")
parser.add_argument("-organism", type=str, default='HomoSapiens')
parser.add_argument("-intrasubject_nbrdists", type=bool, default=2)
parser.add_argument("-clones_files", type=str)
parser.add_argument("-epitope_prefixes", type=str)
parser.add_argument("-nbrdist_percentiles", type=str, default='10',#'5;10;25',
                    help="nbrdist_percentiles.")
parser.add_argument("-chains", type=str, default='B',
                    help="chains.")
parser.add_argument("-ch_distance_type", type=str, default='custom_bsd4neg')
parser.add_argument("-matrix", type=str, help='path to matrix', default='bsd4neg')
parser.add_argument("-distance_params", type=str,
                    default='gap_penalty_v_region:4,gap_penalty_cdr3_region:8,weight_cdr3_region:3,align_cdr3s:True,' #!!!
                            'trim_cdr3s:True,scale_factor:1.0439137134052388')
parser.add_argument("-epitope_col", type=str, default='epitope.seq')   #####!!!!!!!
parser.add_argument("-multiproc", type=int, default=2)
parser.add_argument("-vmatrix", type=str, default=None)

#work with args
args = parser.parse_args()
clones_file = args.clones_file
db_file = args.db_file
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
    elif 'bsd4neg' in ch_distance_type:
        matrix = blosum.bsd4neg
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
    outfile = '{}_{}_{}_nbrdists.tsv'.format(make_dirprefix(clones_files[0], 'dists'), organism, chains)
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

sample_tcrs = ch_read_files.read_tcr(clones_file, organism, chains, epitope_col)
db_tcrs = ch_read_files.read_tcr(db_file, organism, chains, epitope_col)

sample_tcrs.to_csv('{}_{}_{}_all_tcrs.tsv'.format(make_dirprefix(clones_files[0], 'all_tcrs'), organism, chains))
db_tcrs.to_csv('{}_{}_{}_all_tcrs.tsv'.format(make_dirprefix(db_file, 'all_tcrs'), organism, chains))


#all_tcrs, all_tcrs_col = get_groups(sample_tcrs, epitope_col, tcr_col, chains, nbrdist_percentiles, distance_params, rep_dists, multiproc, vdjdb_tcrs, ch_distance_type)

iti = 0

numepitopes = len(pd.unique(db_tcrs[epitope_col]))

def compute_dists(index):
    ntcr = sample_tcrs.loc[index,]  # tcrr
    edists = []
    for tindex, tcr in epigroup.iterrows():
        if tcr['tcr_info'] == ntcr['tcr_info']:
            print(tcr['tcr_info'], ntcr['tcr_info'], 'continued')
            continue
        edists.append(ch_tcr_distances.compute_distance(tcr['tcr_info'], ntcr['tcr_info'],
                                                        chains, distance_params, rep_dists=rep_dists))
        if edists[-1] == 0:
            print(tcr['tcr_info'], ntcr['tcr_info'], 'added to dists')
    edists.sort()
    for nbrdist_percentile in nbrdist_percentiles:
        #nbrdist = ch_tcr_distances.sort_and_compute_nbrdist_from_distances(edists, nbrdist_percentile,
        #                                                                   dont_sort=True)
        wtd_nbrdist = ch_tcr_distances.sort_and_compute_weighted_nbrdist_from_distances(edists,
                                                                                        nbrdist_percentile,
                                                                                        dont_sort=True)
        #ntcr['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
        #    nbrdist)
        ntcr['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
            wtd_nbrdist)
    return ntcr


for epitope, epigroup in db_tcrs.groupby(epitope_col):
    if len(epigroup) < 2:
        iti += 1
        continue
    iti += 1
    print('epitope {}, {}/{}'.format(epitope, iti, numepitopes))
    for nbrdist_percentile in nbrdist_percentiles:
        #sample_tcrs['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        sample_tcrs['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
        #tcr_col.append('{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile))
        tcr_col.append('{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile))

    pool = Pool(multiproc)
    results = pool.map(compute_dists, [i for i in range(len(sample_tcrs.index))])
    pool.close()
    pool.join()
    sample_tcrs = pd.concat(results, axis=1).T

print('Computation of distances is completed')
sample_tcrs.to_csv(outfile, sep='\t', index=None, columns=tcr_col)

te = time.time()
print(te-ts)