#!/usr/bin/python
# Will Tian, Jul 20, 2015
# Refactored by Etienne Richan April, 2018


from os.path import join

import numpy as np

from defs import INPUT_DIR


# This function takes information from the "strands" file in input and places it in one line starting with the name of the pdb
# It basically lists all the "strands" files to be studied together

def run_pretest(pdbs, max_strand_len=None):
    if max_strand_len == None:
        print 'No max_strand_len provided, using default value 12'
        max_strand_len = 12

    with open('bmp51_len' + str(max_strand_len).zfill(2) + '.pretest', 'w') as f:
        for pdb in pdbs:
            f.write(pdb)
            strands = np.loadtxt('../inputs/' + pdb + '/' + pdb + '.strands').astype(int)
            for i in range(len(strands)):
                if i % 2 == 0:
                    start = strands[i][0]
                    end = start + max_strand_len - 1
                    if end > strands[i][1]:
                        end = strands[i][1]
                else:
                    end = strands[i][1]
                    start = end - max_strand_len + 1
                    if start < strands[i][0]:
                        start = strands[i][0]
                f.write(' ' + str(start) + ' ' + str(end))
            f.write('\n')


#Return a dict of pdbs where each entry is a list of tuples of each start/end pair
def get_pretest_dict(pdbs, max_strand_len=None):
    if max_strand_len == None:
        print 'No max_strand_len provided, using default value 12'
        max_strand_len = 12
    pdb_dict = {}
    for pdb in pdbs:
        pdb_dict[pdb]= []
        strands = np.loadtxt(join(INPUT_DIR, *[pdb, pdb + '.strands'])).astype(int)
        for i in range(len(strands)):
            if i % 2 == 0:
                start = strands[i][0]
                end = start + max_strand_len - 1
                if end > strands[i][1]:
                    end = strands[i][1]
            else:
                end = strands[i][1]
                start = end - max_strand_len + 1
                if start < strands[i][0]:
                    start = strands[i][0]
            pdb_dict[pdb].append((start,end))
    return pdb_dict