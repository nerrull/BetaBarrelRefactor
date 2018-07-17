#!/usr/bin/python
import glob
import numpy as np
from os.path import join
from defs import BARREL_OUTPUT_DIR, REINDEX_DIR


def correct_and_trim(pdbs):
    for resol in ['ca', 'bb', 'sc']:
        for pdb in pdbs:
            fn = join(BARREL_OUTPUT_DIR, *[resol, pdb + '_reidx_ext.pdb'])

            foutn = join(BARREL_OUTPUT_DIR, *[resol, pdb + '.pdb'])
            reindexmap = np.loadtxt(join(REINDEX_DIR,  pdb + '.map')).astype(int)
            index = [str(i).rjust(4) for i in reindexmap[:, 0].tolist()]
            seqindex = [str(i).rjust(4) for i in reindexmap[:, 1].tolist()]

            idxdict = dict(zip(index, seqindex))

            # print pdb
            # print idxdict
            with open(fn) as f:
                lines = f.readlines()

            with open(foutn, 'w') as fout:
                for line in lines:
                    if len(line) > 26:
                        if line[22:26] in idxdict:
                            fout.write(line[:22] + idxdict[line[22:26]] + line[26:])
                    else:
                        fout.write(line)
