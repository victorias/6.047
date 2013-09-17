#!/usr/bin/env python

import sys

base_idx = { 'A' : 0, 'G' : 1, 'C' : 2, 'T' : 3 }
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3

def seqalignDP(seq1,seq2,subst_matrix,gap_pen):
    """return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
       Note: gap_pen should be positive (it is subtracted)
    """
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
    for i in range(1,len(seq1)+1):
        F[i][0] = 0 - i*gap_pen
        TB[i][0] = PTR_GAP2 # indicates a gap in seq2
    for j in range(1,len(seq2)+1):
        F[0][j] = 0 - j*gap_pen
        TB[0][j] = PTR_GAP1 # indicates a gap in seq1

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            seq1Letter = seq1[i-1]
            seq2Letter = seq2[j-1]
            match = F[i-1][j-1] + subst_matrix[base_idx[seq1Letter]][base_idx[seq2Letter]]
            delete = F[i-1][j] - gap_pen
            insert = F[i][j-1] - gap_pen
            
            m = max(match, delete, insert)
            if m == match:
                TB[i][j] = PTR_BASE
            elif m == delete:
                TB[i][j] = PTR_GAP2
            elif m == insert:
                TB[i][j] = PTR_GAP1
            F[i][j] = m

    return F[len(seq1)][len(seq2)], F, TB

def traceback(seq1,seq2,TB):
    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)

    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i=i-1
            j=j-1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j=j-1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i=i-1
        else: assert False

    return s1,s2

def readSeq(filename):
    """reads in a FASTA sequence"""

    stream = open(filename)
    seq = []

    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())

    return "".join(seq)

S = [
    # A  G   C   T
    [3, -1, -2, -2], # A
    [-1, 3, -2, -2], # G
    [-2, -2, 3, -1], # C
    [-2, -2, -1, 3]  # T
    ]
gap_pen = 4

edit_S = [
    # A  G   C   T
    [0, -1, -1, -1], # A
    [-1, 0, -1, -1], # G
    [-1, -1, 0, -1], # C
    [-1, -1, -1, 0]  # T
    ]
edit_gap_pen = 1

def main():
    # parse commandline
    if len(sys.argv) < 3:
        print "you must call program as: python ps1-seqalign.py <FASTA 1> <FASTA 2>"
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = seqalignDP(seq1,seq2,S,gap_pen)

    print >> sys.stderr, score

    s1, s2 = traceback(seq1,seq2,TB)
    print s1
    print s2


def main2():


    seq1 = "AGGTGAT"
    seq2 = "AGTAA"

    score, F, TB = seqalignDP(seq1,seq2,edit_S,edit_gap_pen)

    print >> sys.stderr, score

    s1, s2 = traceback(seq1,seq2,TB)
    print s1
    print s2
    
    print F
    print TB


def mainEdit():
    # parse commandline
    if len(sys.argv) < 3:
        print "you must call program as: python ps1-seqalign.py <FASTA 1> <FASTA 2>"
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = seqalignDP(seq1,seq2,edit_S,edit_gap_pen)

    print >> sys.stderr, score

    s1, s2 = traceback(seq1,seq2,TB)
    print s1
    print s2


if __name__ == "__main__":
    mainEdit()
