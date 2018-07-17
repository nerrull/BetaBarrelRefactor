from utils import loadSearchPaths, check_input_files
from bb_register import run_register_prediction
from bb_barrel import generate_barrel_structure
import argparse



''' Level descriptions: '''
'''
1 : Small BMPs (N < 16) without inplugs or outclamps
2 : Small BMPs (N < 16) with inplugs or outclamps
3 : Medium oligomeric BMPs (16 <= N < 20)
4 : Medium monomeric BMPs (16 <= N < 20)
5 : Large BMPs (N >= 20)
'''

parser =argparse.ArgumentParser()
parser.add_argument("--fastafile", metavar= "F", type = str, help="Path to the FASTA file containing ")


#Add subdirectories to python path
loadSearchPaths()

# Arguments
pdbs = ["A9WGN5_scrambled"]

if not check_input_files(pdbs):
    exit(0)

level =4
run_register_prediction(pdbs, level)
generate_barrel_structure(pdbs, level)
