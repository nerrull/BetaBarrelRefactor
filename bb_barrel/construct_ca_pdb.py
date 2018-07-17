#!/usr/bin/python

from scripts.barrel_2zz import Barrel, write_pdb
from os.path import join
from defs import CONSTRUCTION_DIR, CA_DIR, REINDEX_DIR

def construct_ca(pdbs):
    for pdb in pdbs:
        # best params for 2zz model
        bb = Barrel(pdb, 3.345, 4.83, 0.85, 0.22, join(CONSTRUCTION_DIR,'results'), True)
        write_pdb(bb.balls, join(CA_DIR, pdb + '_reidx_ext.pdb'), bb.strandlens, True, bb.reindexmap,
                  join(REINDEX_DIR, pdb + '.map'))
