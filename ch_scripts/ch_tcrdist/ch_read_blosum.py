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


import os
import pandas as pd
import numpy as np
import copy
class Found(Exception): pass


occurence = {'A': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'C': 0,
             'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}


def read_blosum(inppath):
    scores = {'matrix': pd.DataFrame(float(0), index=sorted(occurence), columns=sorted(occurence))}
    with open(inppath, 'r') as inp:
        for line in inp:
            if line.startswith('#'):
                continue
            else:
                if line.strip().startswith('A'):
                    upline = line.strip().split()
                    for line in inp:
                        info = line.strip().split()
                        if info[0] in occurence:
                            for i in range(len(upline)):
                                if upline[i] in occurence:
                                    scores['matrix'][info[0]][upline[i]] = int(info[i+1])
    return(scores)


def read_blosum_default(matrix, path_to_dir='../../blosum_matrices'):
    blosum_path = os.path.join(path_to_dir, 'BLOSUM{}'.format(matrix))
    return(read_blosum(blosum_path))


def getmax(matrix):
    maxmatrix = copy.deepcopy(matrix)
    np.fill_diagonal(maxmatrix.values, 0)
    return(maxmatrix.values.max())


def compute_cdrmatrix(importmatrix, lim=None, sign=-1):
    if lim is None:
        lim = getmax(importmatrix['matrix'])
    cdrdistmatrix = {}
    for a in amino_acids:
        for b in amino_acids:
            bab = importmatrix['matrix'][a][b]
            if a == b:
                cdrdistmatrix[(a, b)] = 0
            else:  ## different
                if bab < 0:
                    cdrdistmatrix[(a, b)] = sign*lim
                else:
                    cdrdistmatrix[(a, b)] = sign*(lim - bab)
    return (cdrdistmatrix)


def compute_matrix_shugai(inpfile, params=None):
    inp = pd.read_table(inpfile)
    inp.columns = ['aa_1', 'aa_2', 'cdr3_score']
    if params is None:
        outdict = dict(zip(zip(inp.aa_1, inp.aa_2), -1*inp.cdr3_score))
    else:
        print('params {} yet to be implemented'.format(params))
        return 0
    for aa in occurence:
        outdict[('C', aa)] = 0
        outdict[(aa, 'C')] = 0
    return outdict


blosum = {('S', 'W'): -3, ('G', 'G'): 6, ('E', 'M'): -2, ('A', 'N'): -2, ('A', 'Y'): -2, ('W', 'Q'): -2, ('V', 'N'): -3, ('F', 'K'): -3, ('G', 'E'): -2, ('E', 'D'): 2, ('W', 'P'): -4, ('I', 'T'): -1, ('F', 'D'): -3, ('K', 'V'): -2, ('C', 'Y'): -2, ('G', 'D'): -1, ('T', 'N'): 0, ('W', 'W'): 11, ('S', 'S'): 4, ('K', 'C'): -3, ('E', 'F'): -3, ('N', 'L'): -3, ('A', 'K'): -1, ('Q', 'P'): -1, ('F', 'G'): -3, ('D', 'S'): 0, ('C', 'V'): -1, ('V', 'T'): 0, ('H', 'P'): -2, ('P', 'V'): -2, ('I', 'Q'): -3, ('F', 'V'): -1, ('W', 'T'): -2, ('H', 'F'): -1, ('P', 'D'): -1, ('Q', 'R'): 1, ('D', 'Q'): 0, ('K', 'Q'): 1, ('D', 'F'): -3, ('V', 'W'): -3, ('T', 'C'): -1, ('A', 'F'): -2, ('T', 'H'): -2, ('A', 'Q'): -1, ('Q', 'T'): -1, ('V', 'F'): -1, ('F', 'C'): -2, ('C', 'R'): -3, ('V', 'P'): -2, ('H', 'T'): -2, ('E', 'L'): -3, ('F', 'R'): -3, ('I', 'G'): -4, ('C', 'Q'): -3, ('Y', 'V'): -1, ('T', 'A'): 0, ('T', 'V'): 0, ('Q', 'V'): -2, ('S', 'K'): 0, ('K', 'K'): 5, ('E', 'N'): 0, ('N', 'T'): 0, ('A', 'H'): -2, ('A', 'C'): 0, ('V', 'S'): -2, ('Q', 'H'): 0, ('H', 'S'): -1, ('Q', 'Y'): -1, ('P', 'N'): -2, ('I', 'Y'): -1, ('P', 'G'): -2, ('F', 'N'): -3, ('H', 'N'): 1, ('K', 'H'): -1, ('N', 'W'): -4, ('S', 'Y'): -2, ('W', 'N'): -4, ('D', 'Y'): -3, ('E', 'Q'): 2, ('K', 'Y'): -2, ('S', 'G'): 0, ('Y', 'S'): -2, ('G', 'R'): -2, ('A', 'L'): -1, ('A', 'G'): 0, ('T', 'K'): -1, ('T', 'P'): -1, ('M', 'V'): 1, ('Q', 'L'): -2, ('E', 'S'): 0, ('H', 'W'): -2, ('I', 'D'): -3, ('K', 'F'): -3, ('N', 'A'): -2, ('T', 'I'): -1, ('Q', 'N'): 0, ('K', 'W'): -3, ('S', 'C'): -1, ('Y', 'Y'): 7, ('G', 'V'): -3, ('L', 'V'): 1, ('A', 'R'): -1, ('M', 'R'): -1, ('Y', 'L'): -1, ('D', 'C'): -3, ('P', 'P'): 7, ('D', 'H'): -1, ('Q', 'Q'): 5, ('I', 'V'): 3, ('P', 'F'): -4, ('I', 'A'): -1, ('F', 'F'): 6, ('K', 'T'): -1, ('L', 'T'): -1, ('S', 'Q'): 0, ('W', 'F'): 1, ('D', 'A'): -2, ('E', 'Y'): -2, ('K', 'A'): -1, ('Q', 'S'): 0, ('A', 'D'): -2, ('L', 'R'): -2, ('T', 'S'): 1, ('A', 'V'): 0, ('M', 'N'): -2, ('Q', 'D'): 0, ('E', 'P'): -1, ('V', 'V'): 4, ('D', 'N'): 1, ('I', 'S'): -2, ('P', 'M'): -2, ('H', 'D'): -1, ('I', 'L'): 2, ('K', 'N'): 0, ('L', 'P'): -3, ('Y', 'I'): -1, ('N', 'I'): -3, ('T', 'Q'): -1, ('Q', 'F'): -3, ('S', 'M'): -1, ('E', 'R'): 0, ('Q', 'W'): -2, ('G', 'N'): 0, ('L', 'Y'): -1, ('L', 'N'): -3, ('A', 'S'): 1, ('D', 'T'): -1, ('S', 'T'): 1, ('P', 'S'): -1, ('V', 'R'): -3, ('D', 'K'): -1, ('P', 'H'): -2, ('H', 'C'): -3, ('Q', 'I'): -3, ('H', 'H'): 8, ('I', 'I'): 4, ('L', 'W'): -2, ('L', 'L'): 4, ('D', 'R'): -2, ('S', 'I'): -2, ('D', 'I'): -3, ('E', 'A'): -1, ('K', 'I'): -3, ('Q', 'K'): 1, ('T', 'D'): -1, ('A', 'W'): -3, ('Y', 'R'): -2, ('M', 'F'): 0, ('S', 'P'): -1, ('H', 'Q'): 0, ('Y', 'N'): -2, ('I', 'P'): -3, ('E', 'C'): -4, ('H', 'G'): -2, ('P', 'E'): -1, ('Q', 'M'): 0, ('H', 'L'): -3, ('L', 'S'): -2, ('L', 'H'): -3, ('N', 'Q'): 0, ('T', 'Y'): -2, ('K', 'G'): -2, ('S', 'E'): 0, ('Y', 'E'): -2, ('W', 'R'): -3, ('V', 'M'): 1, ('N', 'R'): 0, ('G', 'F'): -3, ('F', 'Y'): 3, ('L', 'Q'): -2, ('M', 'Y'): -1, ('A', 'P'): -1, ('S', 'N'): 1, ('C', 'L'): -1, ('L', 'F'): 0, ('D', 'W'): -4, ('S', 'L'): -2, ('P', 'R'): -2, ('P', 'K'): -1, ('Y', 'G'): -3, ('C', 'K'): -3, ('H', 'K'): -1, ('Q', 'A'): -1, ('I', 'F'): 0, ('K', 'D'): -1, ('N', 'C'): -3, ('L', 'D'): -4, ('Y', 'K'): -2, ('S', 'A'): 1, ('W', 'V'): -3, ('E', 'I'): -3, ('V', 'I'): 3, ('Q', 'C'): -3, ('T', 'G'): -2, ('T', 'L'): -1, ('L', 'M'): 2, ('A', 'T'): 0, ('C', 'H'): -3, ('P', 'Y'): -3, ('S', 'H'): -1, ('H', 'Y'): 2, ('E', 'K'): 1, ('C', 'G'): -3, ('I', 'C'): -1, ('Q', 'E'): 2, ('K', 'R'): 2, ('T', 'E'): -1, ('L', 'K'): -2, ('M', 'W'): -1, ('N', 'Y'): -2, ('N', 'H'): 1, ('V', 'E'): -2, ('Q', 'G'): -2, ('Y', 'D'): -3, ('F', 'Q'): -3, ('G', 'Y'): -3, ('L', 'I'): 2, ('M', 'Q'): 0, ('R', 'A'): -1, ('C', 'D'): -3, ('S', 'V'): -2, ('D', 'D'): 6, ('S', 'D'): 0, ('P', 'C'): -3, ('C', 'C'): 9, ('W', 'K'): -3, ('I', 'N'): -3, ('K', 'L'): -2, ('N', 'K'): 0, ('L', 'G'): -4, ('M', 'S'): -1, ('R', 'C'): -3, ('R', 'D'): -2, ('V', 'A'): 0, ('W', 'I'): -3, ('T', 'T'): 5, ('F', 'M'): 0, ('L', 'E'): -3, ('M', 'M'): 5, ('R', 'E'): 0, ('W', 'H'): -2, ('S', 'R'): -1, ('E', 'W'): -3, ('P', 'Q'): -1, ('H', 'A'): -2, ('Y', 'A'): -2, ('E', 'H'): 0, ('R', 'F'): -3, ('I', 'K'): -3, ('N', 'E'): 0, ('T', 'M'): -1, ('T', 'R'): -1, ('M', 'T'): -1, ('G', 'S'): 0, ('L', 'C'): -1, ('R', 'G'): -2, ('Y', 'M'): -1, ('N', 'F'): -3, ('Y', 'Q'): -1, ('N', 'P'): -2, ('R', 'H'): 0, ('W', 'M'): -1, ('C', 'N'): -3, ('V', 'L'): 1, ('F', 'I'): 0, ('G', 'Q'): -2, ('L', 'A'): -1, ('M', 'I'): 1, ('R', 'I'): -3, ('W', 'L'): -2, ('D', 'G'): -1, ('D', 'L'): -4, ('I', 'R'): -3, ('C', 'M'): -1, ('H', 'E'): 0, ('Y', 'W'): 2, ('G', 'P'): -2, ('W', 'C'): -2, ('M', 'P'): -2, ('N', 'S'): 1, ('G', 'W'): -2, ('M', 'K'): -1, ('R', 'K'): 2, ('D', 'E'): 2, ('K', 'E'): 1, ('R', 'L'): -2, ('A', 'I'): -1, ('V', 'Y'): -1, ('W', 'A'): -3, ('Y', 'F'): 3, ('T', 'W'): -2, ('V', 'H'): -3, ('F', 'E'): -3, ('M', 'E'): -2, ('R', 'M'): -1, ('E', 'T'): -1, ('H', 'R'): 0, ('P', 'I'): -3, ('F', 'T'): -2, ('C', 'I'): -1, ('H', 'I'): -3, ('G', 'T'): -2, ('I', 'H'): -3, ('R', 'N'): 0, ('C', 'W'): -2, ('W', 'G'): -2, ('N', 'M'): -2, ('M', 'L'): 2, ('G', 'K'): -2, ('M', 'G'): -3, ('K', 'S'): 0, ('E', 'V'): -2, ('N', 'N'): 6, ('V', 'K'): -2, ('R', 'P'): -2, ('A', 'M'): -1, ('W', 'E'): -3, ('F', 'W'): 1, ('C', 'F'): -2, ('V', 'D'): -3, ('F', 'A'): -2, ('G', 'I'): -4, ('M', 'A'): -1, ('R', 'Q'): 1, ('C', 'T'): -1, ('W', 'D'): -4, ('H', 'V'): -3, ('S', 'F'): -2, ('P', 'T'): -1, ('F', 'P'): -4, ('C', 'E'): -4, ('H', 'M'): -2, ('I', 'E'): -3, ('G', 'H'): -2, ('R', 'R'): 5, ('K', 'P'): -1, ('C', 'S'): -1, ('D', 'V'): -3, ('M', 'H'): -2, ('M', 'C'): -1, ('R', 'S'): -1, ('D', 'M'): -3, ('E', 'E'): 5, ('K', 'M'): -1, ('V', 'G'): -3, ('R', 'T'): -1, ('A', 'A'): 4, ('V', 'Q'): -2, ('W', 'Y'): 2, ('F', 'S'): -2, ('G', 'M'): -3, ('C', 'P'): -3, ('E', 'G'): -2, ('I', 'W'): -3, ('P', 'A'): -1, ('F', 'L'): 0, ('C', 'A'): 0, ('G', 'L'): -4, ('R', 'V'): -3, ('T', 'F'): -2, ('Y', 'P'): -3, ('M', 'D'): -3, ('G', 'C'): -3, ('R', 'W'): -3, ('N', 'D'): 1, ('N', 'V'): -3, ('V', 'C'): -1, ('A', 'E'): -1, ('Y', 'H'): 2, ('D', 'P'): -1, ('G', 'A'): 0, ('R', 'Y'): -2, ('P', 'W'): -4, ('Y', 'C'): -2, ('P', 'L'): -3, ('F', 'H'): -1, ('I', 'M'): 1, ('Y', 'T'): -2, ('N', 'G'): 0, ('W', 'S'): -3}

## convert from a similarity to a "distance", not necessarily a true metric as reported by the DEV lines in __main__ output
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

bsd4 = {}

for a in amino_acids:
    for b in amino_acids:
        bab = blosum[(a,b)]
        if a==b:
            bsd4[(a,b)] = 0
        else: ## different
            assert bab<4
            if bab<0:
                bsd4[(a,b)] = 4
            else:
                bsd4[(a,b)] = 4-bab
bsd4neg = {}
for i in bsd4:
    bsd4neg[i] = -bsd4[i]

custommatrix = compute_cdrmatrix(read_blosum('../../blosum_matrices/BLOSUM62'))