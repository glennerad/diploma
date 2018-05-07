import os
import sys
import re
import copy
import argparse
import pandas as pd
from Bio import SeqIO

from ch_base import *
import ch_db_functions as db_funcs
import ch_imgt_subscripts as ch_imgt
import csv

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from vdjdb_full and imgt.fa and then "
                                             "it will get sequences of CDR1, CDR2 and CDR3 by imgt anchor positions.")

parser.add_argument("-i", nargs=1, type=str, default="../vdjdb-2018-04-09/vdjdb.txt",
                    help="vdjdb path.")
parser.add_argument("-i2", nargs=1, type=str, default="../imgt_work/imgt_all.fasta", #"../vdjdb-2017-06-13/imgt.fa",
                    help="imgt.fasta path.")
parser.add_argument("-o", nargs=1, type=str, default="../imgt_work/",
                    help="Output path.")
parser.add_argument("-o2", nargs=1, type=str, default="vdjdb_imgted.txt",
                    help="Output file name.")
parser.add_argument("-addfix", nargs=1, type=bool, default=False,
                    help="Add fix columns?.")
parser.add_argument("-full", nargs=1, type=bool, default=False,
                    help="use vdjdb_full.txt or vdjdb.txt?")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
input_file = os.path.abspath(args.i)

if type(args.i2) is list:
    args.i2 = args.i2[0]
input_file2 = os.path.abspath(args.i2)

if type(args.o) is list:
    args.o = args.o[0]
output_path = os.path.abspath(args.o)

if type(args.addfix) is list:
    args.addfix = args.addfix[0]

if type(args.full) is list:
    args.full = args.full[0]

if args.full is True:
    use_full = True
else:
    use_full = False

if args.addfix is True:
    addcol = True
else:
    addcol = False

crdir(output_path)

if type(args.o2) is list:
    args.o2 = args.o2[0]
output_file = args.o2


#http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
#anchors and mouse_anchors are in ch_base

colnames_default = ['cdr1.alpha', 'cdr2.alpha', 'cdr2.5.alpha', 'cdr3.alpha', 'v.alpha', 'j.alpha',
                    'cdr1.beta', 'cdr2.beta', 'cdr2.5.beta', 'cdr3.beta', 'v.beta', 'j.beta',
                    'species', 'mhc.a', 'mhc.b', 'mhc.class',
                    'epitope.seq', 'epitope',
                    'antigen.species', 'reference.id',  # no meta.epitope.id/subject.id/clone.id
                    'meta.structure.id']
dictimgt = ch_imgt.get_sliced_dictimgt(input_file2)


def get_partial_melted_vdjdb(inpfile, to_tsv=None, colnames=colnames_default):
    colnames=copy.deepcopy(colnames)

    inpcols = ['complex.id', 'gene',
               'cdr3', 'v.segm', 'j.segm',  # -> .{alpha/beta}
               'species', 'mhc.a', 'mhc.b', 'mhc.class',
               'antigen.epitope', 'antigen.gene',
               'antigen.species', 'reference.id', 'meta']

    if addcol is True:
        inpcols += ['cdr3fix']

    annot = pd.read_table(inpfile, usecols=inpcols)[inpcols]

    # get sequences with 3D structures
    pdbs = copy.deepcopy(annot[(~annot['meta'].str.contains('"structure.id": ""'))
                               & (~annot['meta'].str.contains('"structure.id": "NA"'))])
    pdbs['meta.structure.id'] = pdbs['meta'].str.replace(r'.*"structure.id": "(.{,4})".*', r'\1').str.lower()

    # get sequences without 3D structures
    annot = copy.deepcopy(annot[(annot['meta'].str.contains('"structure.id": ""'))
                                | (annot['meta'].str.contains('"structure.id": "NA"'))])
    annot['meta.structure.id'] = ''

    # merge df together
    annot = pdbs.append(annot, ignore_index=True).drop(columns=['meta'])

    annot['cdr1'], annot['cdr2'], annot['cdr2.5'] = '', '', ''

    annotcols = ['complex.id', 'gene',
                 'cdr1', 'cdr2', 'cdr2.5', 'cdr3', 'v.segm', 'j.segm',  # -> .{alpha/beta}
                 'species', 'mhc.a', 'mhc.b', 'mhc.class',
                 'antigen.epitope', 'antigen.gene',
                 'antigen.species', 'reference.id', 'meta.structure.id']

    change_cols = ['cdr1', 'cdr2', 'cdr2.5', 'cdr3', 'v', 'j']

    if addcol is True:
        colnames += ['cdr3fix.alpha', 'cdr3fix.beta']
        annot['cdr3fix'] = annot['cdr3fix'].str.replace(r'.*"cdr3": "(.*)", "cdr3_old":.*', r'\1')
        annotcols += ['cdr3fix']
        change_cols += ['cdr3fix']

    annot = annot[annotcols]
    annot = annot.rename(index=str,
                         columns={'v.segm': 'v', 'j.segm': 'j', 'antigen.epitope': 'epitope.seq',
                                  'antigen.gene': 'epitope'})

    unpaired = annot[annot['complex.id'] == 0]
    unpaired = db_funcs.ch_partial_melt(unpaired, melt_column='gene',
                                     change_columns=change_cols,
                                     suffixdict={'TRA':'alpha', 'TRB':'beta'}, method='append')
    annot = annot[annot['complex.id'] != 0]

    annot = db_funcs.ch_partial_melt(annot, melt_column='gene',
                                     change_columns=change_cols,
                                     suffixdict={'TRA':'alpha', 'TRB':'beta'})

    annot = annot.append(unpaired, ignore_index=True)[colnames]

    if to_tsv is not None:
        annot.to_csv(to_tsv, sep='\t', index=None)
    return annot


def make_vdjdb_imgted(vdjdbpath, outfile, colnames=colnames_default, type='full'):
    colnames = copy.deepcopy(colnames)
    if addcol is True:
        colnames += ['cdr3fix.alpha', 'cdr3fix.beta']

    result_table = []

    if type == 'full':
        species, mhc_a, mhc_b, mhc_class, epitope_seq, epitope, antigen_spec, ref_id, struct_id = \
            'species', 'mhc.a', 'mhc.b', 'mhc.class', 'antigen.epitope', 'antigen.gene', 'antigen.species', 'reference.id', 'meta.structure.id'
    elif type == 'not_full':
        species, mhc_a, mhc_b, mhc_class, epitope_seq, epitope, antigen_spec, ref_id, struct_id = \
            'species', 'mhc.a', 'mhc.b', 'mhc.class', 'epitope.seq', 'epitope', 'antigen.species', 'reference.id', 'meta.structure.id'

    with open(vdjdbpath) as inp:
        reader = csv.DictReader(inp, delimiter='\t')
        # header = next(reader)
        for row in reader:
            info = row
            cdralpha = ['', '', '', info['cdr3.alpha'], info['v.alpha'], info['j.alpha']]
            cdrbeta = ['', '', '', info['cdr3.beta'], info['v.beta'], info['j.beta']]
            otherinfo = [info[species], info[mhc_a], info[mhc_b], info[mhc_class], info[epitope_seq],
                         info[epitope], info[antigen_spec], info[ref_id], info[struct_id]]
            if addcol is True:
                otherinfo += [info['cdr3fix.alpha'], info['cdr3fix.beta']]
            info[species] = info[species].lower()
            if info[species] in dictimgt:
                if info['v.alpha'] in dictimgt[info[species]]:
                    if info[species] == 'musmusculus':  # use another anchors
                        outdict, outdict2 = ch_imgt.changeseq(dictimgt[info[species]][info['v.alpha']], mouse_anchors)
                    else:
                        outdict, outdict2 = ch_imgt.changeseq(dictimgt[info[species]][info['v.alpha']], anchors)

                    cdralpha[0], cdralpha[1], cdralpha[2] = outdict['CDR1'], outdict['CDR2'], outdict['CDR2.5']

                if info['v.beta'] in dictimgt[info[species]]:
                    outdict, outdict2 = ch_imgt.changeseq(dictimgt[info[species]][info['v.beta']], anchors)

                    cdrbeta[0], cdrbeta[1], cdrbeta[2] = outdict['CDR1'], outdict['CDR2'], outdict['CDR2.5']

            result_table.extend('\t'.join(cdralpha + cdrbeta + otherinfo) + '\n')
    with open(outfile, 'w') as out:
        out.write('\t'.join(colnames)+'\n')
        out.writelines(result_table)


if use_full is True: #MUCH FASTER...
    make_vdjdb_imgted(input_file, os.path.join(output_path, output_file), colnames=colnames_default, type='full')
else:
    get_partial_melted_vdjdb(input_file, to_tsv=os.path.join(BASE_HOMEDIR, 'imgt_work/vdjdb.partial.melted.txt'))
    make_vdjdb_imgted(os.path.join(BASE_HOMEDIR, 'imgt_work/vdjdb.partial.melted.txt'),
                      os.path.join(output_path, output_file),
                      colnames=colnames_default, type='not_full')
