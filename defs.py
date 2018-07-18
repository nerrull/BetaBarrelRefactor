import os
from os.path import join, basename, dirname

defs_file_path = os.path.realpath(__file__)

BASE_DIR =dirname(defs_file_path)
INPUT_DIR = join(BASE_DIR, "inputs")
OUTPUT_DIR = join(BASE_DIR, "output")

# Paths for regiter predition
REGISTER_DIR  = join(BASE_DIR, "bb_register")
ODDS_DIR = join(join(REGISTER_DIR, "pred_combine"), "odds")
HAMMAD_DIR = join(ODDS_DIR, "hammad25")
TMP_DIR =  join(REGISTER_DIR, "tmp_files")


EXEC_FILE  = join(REGISTER_DIR, *["pred_combine","exec", "new_pred_reg_bestscore_ecparam110030080.out"])

# Paths for barrel construction
BARREL_DIR = join(BASE_DIR, "bb_barrel")
BARREL_OUTPUT_DIR = join(BARREL_DIR, "output")
CONSTRUCTION_DIR = join(BARREL_DIR, "construction_inputs")

CA_DIR = join(BARREL_OUTPUT_DIR, "ca")
BB_DIR = join(BARREL_OUTPUT_DIR, "bb")
SC_DIR = join(BARREL_OUTPUT_DIR, "sc")
TMPPDB_DIR = join(BARREL_OUTPUT_DIR, "tmppdbs")
REINDEX_DIR = join(BARREL_OUTPUT_DIR, "reindexmap")
BBQ_DIR  = join(BARREL_DIR, *["scripts", "BBQ"])



SCWRL_PATH =  join(BASE_DIR, "scwrl")

l01 = ['1bxw', '1qj8', '1p4t', '2f1t', '1thq', '2erv', '2lhf', '2mlh', '3dzm', '1qd6', '2f1c', '1k24', '1i78', '2wjr', '4pr7', 'unkn']
l02 = ['1t16', '1uyn', '1tly', '3aeh', '3bs0', '3dwo', '3fid', '3kvn', '4e1s']
l03 = ['2mpr', '1a0s', '2omf', '2por', '1prn', '1e54', '2o4v', '3vzt', '4k3c', '4k3b', '4c4v', '4n75', '3emn']
l04 = ['2qdz', '2ynk', '3rbh', '3syb', '3szv', '4c00', '4gey', '3emn']
l05 = ['1fep', '2fcp', '1kmo', '1nqe', '1xkw', '2vqi', '3csl', '3rfz', '3v8x', '4q35']
