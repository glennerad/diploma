import argparse
from ch_base import *
import ch_read_blosum

#==============================

curr_version = 1.0
parser = argparse.ArgumentParser(description="Pipeline, that will compute only tcr_distances/custom_tcr_distances from dataset")

parser.add_argument("-distance_type", type=str, nargs="+", default="default_shugay".split(";"))
parser.add_argument("-epitope_col", type=str, default="epitope.seq")
parser.add_argument("-multiproc", type=int, default=60)
parser.add_argument("-nlim", type=int, default=30)
parser.add_argument("-datasets", type=str, nargs="+", default=("tcr_subst/vdjdb.beta_imgted.txt".split(";")),
                    help="datasets.")
parser.add_argument("-aamatrix", type=str, default=("tcr_subst/beta.lor.flat"),
                    help="aamatrix.")
parser.add_argument("-vmatrix", type=str, default=("tcr_subst/v.scores.flat"),
                    help="vmatrix.")
parser.add_argument("-organisms", type=str, nargs="+", default="HomoSapiens;MusMusculus".split(";")) #HomoSapiens;
parser.add_argument("-chains", type=str, nargs="+", default="B".split(";")) #A;B;AB
args = parser.parse_args()

distance_type = args.distance_type
epitope_col = args.epitope_col
multiproc = args.multiproc
nlim = args.nlim
datasets = args.datasets
aamatrix = args.aamatrix
vmatrix = args.vmatrix
chains = args.chains
organisms = []
for i in args.organisms:
    if i in org_to_vdjdb_org:
        organisms.append(org_to_vdjdb_org[i])
    else:
        organisms.append(i)


distmatrix = ch_read_blosum.compute_matrix_shugai(aamatrix)
maxvalue = distmatrix[max(distmatrix, key=lambda x: distmatrix[x])]
gapvalues = [int(maxvalue)+1, 2*(int(maxvalue)+1)]


def compute_dists(dataset, diststype, scriptdir, organism, chain, cdr3w, vw, gvw, gcw, dtype, aamatrix=None, vmatrix=None):
    outfile = "{}_{}_{}_{}_{}_{}_nbrdists.tsv".format(make_dirprefix(dataset, diststype, scriptdir), organism, chain, cdr3w, vw, gcw) #'dists_shugay' '../../imgt_work/'
    distance_params = "gap_penalty_v_region:{}," \
                      "gap_penalty_cdr3_region:{}," \
                      "weight_cdr3_region:{}," \
                      "weight_v_region:{}," \
                      "align_cdr3s:False," \
                      "trim_cdr3s:True," \
                      "scale_factor:1.0439137134052388".format(gvw, gcw, cdr3w, vw)
    print(
        '1 START ch_compute_distances_multiprocess.py; organism: {}; weight_cdr3_region: {}; weight_v_region: {}; gap_penalty_cdr3: {}; method: {}\n'.format(
            organism, cdr3w, vw, gcw, dtype))
    cmd = "python3 ch_compute_distances_multiprocess.py -clones_file {} " \
          "-outfile {} " \
          "-organism {} " \
          "-chains {} " \
          "-ch_distance_type {} " \
          "-distance_params {} " \
          "-epitope_col {} " \
          "-multiproc {} " \
          "-matrix {} " \
          "-vmatrix {}".format(
        dataset, outfile, organism, chain, dtype, distance_params, epitope_col, multiproc, aamatrix, vmatrix)
    os.system(cmd)
    print('END ch_compute_distances_multiprocess.py\n=================\n')


cdr3weights = [1]
vweights = [1, 5, 10, 15]
for organism in organisms:
    for dtype in distance_type:
        for dataset in datasets:
            for chain in chains:
                if dtype == 'default':
                    compute_dists(dataset, 'dists_shugay', '../../imgt_work/', organism, chain, 3, 1, 4, 8, dtype)
                else:
                    for cdr3w in cdr3weights:
                        for vw in vweights:
                            for gvalue in gapvalues:
                                compute_dists(dataset, 'dists_shugay', '../../imgt_work/', organism, chain, cdr3w, vw, gapvalues[0], gvalue, dtype, aamatrix, vmatrix)


