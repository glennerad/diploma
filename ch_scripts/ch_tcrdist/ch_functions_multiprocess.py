import pandas as pd
from multiprocessing import Pool

def get_groups(sample_tcrs, epitope_col, tcr_col, chains, nbrdist_percentiles, distance_params, rep_dists, multiproc, db_tcrs=None, ch_distance_type='default'):
    print('Start computing distances')


    if 'default' in ch_distance_type:
        sign = 1
        import ch_tcr_distances_default as ch_tcr_distances
    else:
        sign = -1  # just little hack for distance computation. Biopython is looking for alignments with higheest score. So I made score < 0 and dist = -score
        import ch_tcr_distances_custom as ch_tcr_distances

    iti = 0

    if db_tcrs is None:
        db_tcrs = sample_tcrs

    numepitopes = len(pd.unique(sample_tcrs[epitope_col]))
    if numepitopes < 1:
        return pd.DataFrame(data=None, columns=tcr_col), tcr_col

    def compute_dists(index):
        ntcr = sample_tcrs.loc[index,] #tcrr
        edists = []
        for tindex, tcr in epigroup.iterrows():
            if tcr['tcr_info'] == ntcr['tcr_info']:
                continue
            edists.append(ch_tcr_distances.compute_distance(tcr['tcr_info'], ntcr['tcr_info'],
                                                            chains, distance_params, rep_dists=rep_dists))
        edists.sort()
        for nbrdist_percentile in nbrdist_percentiles:
            nbrdist = ch_tcr_distances.sort_and_compute_nbrdist_from_distances(edists, nbrdist_percentile,
                                                                               dont_sort=True)
            wtd_nbrdist = ch_tcr_distances.sort_and_compute_weighted_nbrdist_from_distances(edists,
                                                                                            nbrdist_percentile,
                                                                                            dont_sort=True)
            ntcr['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = '{:.3f}'.format(
                nbrdist)
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
            sample_tcrs['{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
            sample_tcrs['{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile)] = ''
            tcr_col.append('{}_{}_nbrdist{}'.format(epitope, chains, nbrdist_percentile))
            tcr_col.append('{}_{}_wtd_nbrdist{}'.format(epitope, chains, nbrdist_percentile))

        pool = Pool(multiproc)
        results = pool.map(compute_dists, [i for i in range(len(sample_tcrs.index))])
        pool.close()
        pool.join()
        all_tcrs = pd.concat(results, axis=1).T

    return all_tcrs, tcr_col