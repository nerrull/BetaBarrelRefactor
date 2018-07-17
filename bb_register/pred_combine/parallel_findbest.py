#!/usr/bin/python


import sys

sys.dont_write_bytecode = True
import os
from os.path import join
from defs import ODDS_DIR, EXEC_FILE, INPUT_DIR, TMP_DIR


execprog = './exec/new_pred_reg_bestscore_ecparam110030080.out'


def find_best_combination(test_file, score_file, level):

    output_file = join(TMP_DIR, "results.l0{}".format(level))

    command = "{} {} {} {} {} {} > {} &".format(EXEC_FILE, test_file, score_file,ODDS_DIR,INPUT_DIR, level, output_file)
    print("Running :{}".format(command))
    os.system(command)

