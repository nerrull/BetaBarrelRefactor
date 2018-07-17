#!/usr/bin/python
import glob
import os
from os.path import join
from defs import CA_DIR, TMPPDB_DIR, BB_DIR, BBQ_DIR


def construct_bb(pdbs):
    cwd = os.getcwd()
    for pdb in pdbs:
        ca_file = join(CA_DIR, pdb + "_reidx_ext.pdb")
        print pdb
        with open(ca_file) as fin:
            lines = fin.readlines()
        for i in range(len(lines)):
            lines[i] = lines[i].replace(' CA   ', '  CA  ')
        prebbqfn = join(TMPPDB_DIR, pdb + '_prebbq.pdb')
        with open(prebbqfn, 'w') as fout:
            fout.writelines(lines)

        # run bbq
        postbbqfn = join(BB_DIR, pdb + '_reidx_ext.pdb')
        # print( '/usr/bin/java -classpath /home/wtian7/projects/bb-barrel/scripts/BBQ:/home/wtian7/projects/bb-barrel/scripts/BBQ/jbcl.jar BBQ -d=/home/wtian7/projects/bb-barrel/scripts/BBQ/q_50_xyz.dat -r=/home/wtian7/projects/bb-barrel/' +prebbqfn+' > /home/wtian7/projects/bb-barrel/'+postbbqfn )
        # os.system( '/usr/bin/java -classpath /home/wtian7/projects/bb-barrel/scripts/BBQ:/home/wtian7/projects/bb-barrel/scripts/BBQ/jbcl.jar BBQ -d=/home/wtian7/projects/bb-barrel/scripts/BBQ/q_50_xyz.dat -r=/home/wtian7/projects/bb-barrel/' +prebbqfn+' > /home/wtian7/projects/bb-barrel/'+postbbqfn )
        bbqcall ='/usr/bin/java -classpath ' +BBQ_DIR +":" + join(BBQ_DIR, 'jbcl.jar')+ ' BBQ -d='+ join(BBQ_DIR, 'q_50_xyz.dat') + ' -r=' + prebbqfn + ' > ' +  postbbqfn
        #print(bbqcall)
        os.system( bbqcall)
