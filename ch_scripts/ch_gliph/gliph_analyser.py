import argparse

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will get information from gliph-convergence groups and calculate TP, TN, FP, FN")

parser.add_argument("-v", nargs=1, type=str, default="/Users/AlekseyYeliseev/Desktop/AdaptiveImm/imgt_work/vdjdb_imgted_filtered_M.D_and_chains_B_organism_HomoSapiens_epitope_info.txt",
                    help="vdjdb path.")
parser.add_argument("-i", nargs=1, type=str, default="vdjdb_gliph_table_M.D_B-convergence-groups.txt",
                    help="gliph input path.")

args = parser.parse_args()

if type(args.v) is list:
    args.v = args.v[0]
if type(args.i) is list:
    args.i = args.i[0]


vdjdbpath = args.v
gliphpath = args.i


with open(vdjdbpath, 'r') as inp:
    info = inp.read().strip().split('\n')
    vdjdb_epitopes = {i.split('\t')[1]: set(i.split('\t')[2].split(' ')) for i in info}


with open(gliphpath, 'r') as inp:
    info = inp.read().strip().split('\n')
    gliph_epitopes = {i.split('\t')[1]: set(i.split('\t')[2].split(' ')) for i in info}

gliph_epitopes_info = {} #{group_name:{epitope:[here, total], 'total':numcdrs}}

totalnum = 0
for gepitope in gliph_epitopes:
    numcdrs = len(gliph_epitopes[gepitope])
    totalnum += numcdrs
    gliph_epitopes_info[gepitope] = {'total': numcdrs}
    for epitope in vdjdb_epitopes:
        numintersect = len(gliph_epitopes[gepitope] & vdjdb_epitopes[epitope])
        if numintersect > 0:
            gliph_epitopes_info[gepitope][epitope] = [numintersect, len(vdjdb_epitopes[epitope])]


tp = 0
fp = 0
tn = 0
fn = 0

for gepitope in gliph_epitopes_info:
    total = gliph_epitopes_info[gepitope]['total']
    for epitope in gliph_epitopes_info[gepitope]:
        if epitope == 'total': continue
        locepit = gliph_epitopes_info[gepitope][epitope]
        tp += locepit[0]*(locepit[0]-1)
        tn += locepit[0]*(totalnum-locepit[1])
        fp += locepit[0]*(total-locepit[0])
        fn += locepit[0]*(locepit[1]-locepit[0])


ttp = 0
ttn = 0
for epitope in vdjdb_epitopes:
    total = len(vdjdb_epitopes[epitope])
    ttp += total*(total-1)
    ttn += total*(totalnum-total)


with open('{}_analysis.txt'.format(gliphpath[:gliphpath.rfind('.')]), 'w') as out:
    out.write('TP: {}; TN: {}; FP: {}; FN: {}\n'.format(tp, tn, fp, fn))
    out.write('Precision: {:.4f}; Recall: {:.4f}\n'.format(tp/(tp+fp), tp/(tp+fn)))
    out.write('Specificity: {:.4f}\n'.format(tn/(tn+fp)))