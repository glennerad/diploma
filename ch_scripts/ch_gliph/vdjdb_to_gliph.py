import os
import pandas as pd
import argparse

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from vdjdb and convert it to gliph table.")

parser.add_argument("-i", nargs=1, type=str, default="/Users/AlekseyYeliseev/Desktop/AdaptiveImm/imgt_work/vdjdb_imgted_filtered_M.D_and_chains_B_organism_HomoSapiens.txt",
                    help="vdjdb path.")
parser.add_argument("-o", nargs=1, type=str, default="vdjdb_gliph_table_M.D_B.txt",
                    help="output path.")
parser.add_argument("-c", nargs=1, type=str, default="B",
                    help="chains.")

args = parser.parse_args()

if type(args.i) is list:
    args.i = args.i[0]
if type(args.o) is list:
    args.o = args.o[0]
if type(args.c) is list:
    args.c = args.c[0]

vdjdbpath = args.i
gliphpath = args.o
chains = args.c

#CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	PatientCounts
#CAADTSSGANVLTF	TRBV30	TRBJ2-6	CALSDEDTGRRALTF	TRAV19		TRAJ5	09/02171
#CAATGGDRAYEQYF	TRBV2	TRBJ2-7	CAASSGANSKLTF	TRAV13-1	TRAJ56	03/04922
#CAATQQGETQYF	TRBV2	TRBJ2-5	CAASYGGSARQLTF	TRAV13-1	TRAJ22	02/02591
#CACVSNTEAFF	TRBV28	TRBJ1-1	CAGDLNGAGSYQLTF	TRAV25		TRAJ28	PBMC8631
#CAGGKGNSPLHF	TRBV2	TRBJ1-6	CVVLRGGSQGNLIF	TRAV12-1	TRAJ42	02/02071
#CAGQILAGSDTQYF	TRBV6-4	TRBJ2-3	CATASGNTPLVF	TRAV17		TRAJ29	09/00181
#CAGRTGVSTDTQYF	TRBV5-1	TRBJ2-3	CAVTPGGGADGLTF	TRAV41		TRAJ45	02/02591
#CAGYTGRANYGYTF	TRBV2	TRBJ1-2	CVVNGGFGNVLHC	TRAV12-1	TRAJ35	01/08733

#needed = CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ, PatientCounts

#vdjdb = cdr3.alpha, j.alpha, v.alpha, cdr3.beta, j.beta, v.beta, patientcounts

def vdjdb_to_gliph(vdjdbpath, gliphpath, chains = "BA"):
    gliphcols = []
    vdjdbcols = []
    if len(chains) == 2 and chains.startswith('A'):
        chains = 'BA'
    for chain in chains:
        if chain == "A":
            gliphcols += ['CDR3a', 'TRAV', 'TRAJ']
            vdjdbcols += ['cdr3.alpha', 'v.alpha', 'j.alpha']
        elif chain == 'B':
            gliphcols += ['CDR3b', 'TRBV', 'TRBJ']
            vdjdbcols += ['cdr3.beta', 'v.beta', 'j.beta']

    vdjdb = pd.read_table(vdjdbpath, usecols=vdjdbcols).reindex_axis(vdjdbcols, axis=1)
    vdjdb.columns = gliphcols
    vdjdb = vdjdb.dropna()
    for col in gliphcols:
        vdjdb[col] = vdjdb[col].str.replace(r'\*\d*', '')
    vdjdb.to_csv(gliphpath, sep='\t', index=None)

vdjdb_to_gliph(vdjdbpath, gliphpath, chains)