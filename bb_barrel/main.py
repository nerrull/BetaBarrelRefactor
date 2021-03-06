#!/usr/bin/python
import os
from os.path import join, exists
from os import mkdir
import time
from shutil import rmtree, copy
from construct_ca_pdb import construct_ca
from construct_bb_pdb import construct_bb
from construct_sc_pdb import construct_sc
from correct_index_and_trim import correct_and_trim
from defs import CONSTRUCTION_DIR, BARREL_OUTPUT_DIR, BARREL_OUTPUT_PDB_DIR, BASE_DIR, CA_DIR, SC_DIR, BB_DIR, REINDEX_DIR,TMPPDB_DIR, TMP_DIR, RESULTS_DIR

def clean():
    results_dir = join(CONSTRUCTION_DIR, 'results')
    rmtree(results_dir)
    os.mkdir(results_dir)
    for dir in [ CA_DIR, SC_DIR, BB_DIR, TMPPDB_DIR, REINDEX_DIR]:
        rmtree(dir)
        os.mkdir(dir)

def copy_outputs(pdb, level):

    pdb_path = join(SC_DIR, pdb + ".pdb")
    pdb_path_2 = join(SC_DIR, pdb+ "_reidx_ext.pdb")

    out_path = join(BARREL_OUTPUT_PDB_DIR, pdb + "_l0{}.pdb".format(level))
    out_path_2 = join(BARREL_OUTPUT_PDB_DIR, pdb + "_ext_l0{}.pdb".format(level))

    copy(pdb_path, out_path)
    copy(pdb_path_2, out_path_2)

    return (out_path, out_path_2)

def makedirs():
    dirs = [CONSTRUCTION_DIR, RESULTS_DIR, BARREL_OUTPUT_PDB_DIR, BARREL_OUTPUT_DIR, BASE_DIR, CA_DIR, SC_DIR, BB_DIR, REINDEX_DIR, TMPPDB_DIR, TMP_DIR]
    for dir in dirs:
        if not exists(dir):
            mkdir(dir)

def generate_barrel_structure(pdbs,level):

    makedirs()

    print ("Cleaning up")
    clean()

    print ('preparing input from register prediciton results...')
    os.chdir(CONSTRUCTION_DIR)
    os.system('bash process_all.sh')
    os.chdir(BASE_DIR)

    print 'constructing Ca atoms...'
    construct_ca(pdbs)
    time.sleep(0.1)

    print 'constructing backbone atoms from Ca atoms using BBQ algorithm...'
    construct_bb(pdbs)
    time.sleep(0.1)

    print 'constructing sidechains using Scwrl4...'
    construct_sc(pdbs)
    time.sleep(0.1)

    print 'cleaning tmp info to generate final structures...'
    correct_and_trim(pdbs)
    time.sleep(0.1)

    output_pdbs = copy_outputs(pdbs[0], level)
    return output_pdbs


if __name__ =="__main__":
    pdbs = ["unkn"]
    generate_barrel_structure(pdbs, 0)

#TODO FIGURE OUT WHY IT DOESNT WORK in the old version
