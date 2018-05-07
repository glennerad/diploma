import argparse
import os

#==============================
curr_version = 1.0
parser = argparse.ArgumentParser(description="This script will run all needed scripts to get TP/TN/FP/FN")

parser.add_argument("-i", nargs=1, type=str, default="/Users/AlekseyYeliseev/Desktop/AdaptiveImm/imgt_work/vdjdb_imgted_filtered_M.D.txt",
                    help="input file path.")
parser.add_argument("-c", nargs=1, type=str, default="B",
                    help="chains.")
parser.add_argument("-o", nargs=1, type=str, default="/Users/AlekseyYeliseev/Desktop/AdaptiveImm/imgt_work/vdjdb_imgted_filtered_M.D_and_chains_B_organism_HomoSapiens.txt",
                    help="output file path.")
parser.add_argument("-org", nargs=1, type=str, default="HomoSapiens",
                    help="Organism.")
parser.add_argument("-path_to_scripts", nargs=1, type=str, default="/Users/AlekseyYeliseev/Desktop/AdaptiveImm/ch_scripts/",
                    help="Path to scripts.")

args = parser.parse_args()


if type(args.i) is list:
    args.i = args.i[0]
input_file = args.i

if type(args.o) is list:
    args.o = args.o[0]
output_file = args.o

if type(args.c) is list:
    args.c = args.c[0]
chains = args.c

if type(args.org) is list:
    args.org = args.org[0]
organism = args.org

if type(args.path_to_scripts) is list:
    args.path_to_scripts = args.path_to_scripts[0]
path_to_scripts = args.path_to_scripts

#PART1 ch_vdjdb_filter_argscr
print("PART1")
cmd = 'python3 {}/ch_vdjdb_filter_argscr.py -i {} -c {} -o {} -org {}'.format(path_to_scripts, input_file, chains, output_file, organism)
os.system(cmd)

#PART2
print("PART2")
input_file = output_file
vdjdb_file = "{}_epitope_info.txt".format(input_file[:input_file.rfind('.')])
cmd = 'python3 {}/ch_analyse_epitopes.py -i {} -c {}'.format(path_to_scripts, input_file, chains)
os.system(cmd)

#PART3
print("PART3")
output_file = '{}_gliph_table.txt'.format(input_file[:output_file.rfind('.')])
cmd = 'python3 vdjdb_to_gliph.py -i {} -c {} -o {}'.format(input_file, chains, output_file)
os.system(cmd)

#PART4
print("PART4")
cmd = 'perl bin/gliph-group-discovery.pl --tcr={}'.format(output_file)
os.system(cmd)

#PART5
print("PART5")
input_file = '{}-convergence-groups.txt'.format(output_file[:output_file.rfind('.')])
cmd = 'python3 gliph_analyser.py -v {} -i {}'.format(vdjdb_file, input_file)
os.system(cmd)