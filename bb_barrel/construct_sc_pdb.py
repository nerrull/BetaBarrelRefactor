#!/usr/bin/python
import glob
import os
from os.path import join
from defs import SCWRL_PATH, BARREL_OUTPUT_DIR, BB_DIR, SC_DIR

def construct_sc(pdbs):
    print 'Sidechain construction depends on Scwrl4.'
    print 'You can find the software through this link http://dunbrack.fccc.edu/scwrl4/'

    print(SCWRL_PATH)
    if os.path.exists(SCWRL_PATH):
        path = SCWRL_PATH
    else:
        print 'please enter the path of your Scwrl4 executable (eg. /path/to/scwrl4/ ):'
        path=raw_input()

    for pdb in pdbs:
        fn = join(BB_DIR, pdb + "_reidx_ext.pdb")
        foutn = join (SC_DIR, pdb+'_reidx_ext.pdb')
        os.system(path+'/Scwrl4 -h -i '+fn+' -o '+ foutn )
