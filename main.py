from os.path import exists, splitext,basename
import argparse
import time

from utils import loadSearchPaths, check_input_files, autodetectLevel
from bb_register import run_register_prediction
from bb_barrel import generate_barrel_structure
from process_all import prepare_input
from LoopModeller import prepModellerFiles

lvls ='''
1 : Small BMPs (N < 16) without inplugs or outclamps\n
2 : Small BMPs (N < 16) with inplugs or outclamps\n
3 : Medium oligomeric BMPs (16 <= N < 20)\n
4 : Medium monomeric BMPs (16 <= N < 20)\n
5 : Large BMPs (N >= 20)
'''


#Add subdirectories to python path
loadSearchPaths()

parser =argparse.ArgumentParser()
parser.add_argument("-f","--fastafile", metavar= "F", type = str, help="Path to the FASTA file containing the sequences to predict")
parser.add_argument("-l","--level", metavar= "L", type = int, default=None, help="The level to run the prediction at (default auto-detects the appropriate level): -1 to run all levels. \n" + lvls)


if __name__ ==  "__main__":
    args =parser.parse_args()
    fasta_file = args.fastafile
    if not exists(fasta_file):
        print("Couldn't find the fasta file : \n{}".format(fasta_file))

    pdbs = prepare_input(fasta_file, scramble=True)

    if not check_input_files(pdbs):
        exit(0)

    print ("Running ")
    level = args.level

    modellerCommands = []
    for pdb in pdbs:
        levels = [level]
        if level == None:
            levels =autodetectLevel(pdb)
        elif level ==-1:
            levels =[1,2,3,4,5]

        for l in levels:
            print ("Predicting structure for level {}".format(l))
            run_register_prediction([pdb], l)
            time.sleep(0.5)

            output_files = generate_barrel_structure([pdb], l)
            time.sleep(1.)
            extended_pdb = output_files[1]

            prepModellerFiles(pdb,extended_pdb)
            pdb_name = splitext(basename(extended_pdb))[0]
            modellerCommands.append("python ./LoopModeller.py -b {0} -l {1} -n 100".format(pdb_name, l) )

    print("To add loops with MODELLER to all generated PDBS, run these commands")
    for command in modellerCommands:
        print(command)
    print("To inpect PDBs look in /output")


