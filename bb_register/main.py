import sys
from os import walk, getcwd, mkdir
from os.path import join, exists
from shutil import rmtree

from pretest import get_pretest_dict, prepare_test_dict
from registration import get_scores
from pred_combine import find_best_combination
from shear_adjustment import get_shear_adjustments
from defs import TMP_DIR,EXEC_FILE, BASE_DIR


def cleartmp():
    if exists(TMP_DIR):
        rmtree(TMP_DIR)
    mkdir(TMP_DIR)


def run_register_prediction(pdbs, level):

    print ("Clearing temp dir")
    cleartmp()

    print ('Beta barrel membrane protein register prediction')
    # prepare the information for the test
    print ('Preparing input for register prediction...')

    pdb_dict = get_pretest_dict(pdbs)
    test_dict, test_file_path = prepare_test_dict(pdb_dict)

    # Get registration scores?
    ec_raw_dir = join(BASE_DIR,'sc.psicov')
    weights = [110, 30, 80]
    log = 0
    theta =0

    score_output_file= join(TMP_DIR, "scores.ecs")
    scores = get_scores(test_dict, weights, log, theta, score_output_file, ec_raw_dir)

    if not exists(EXEC_FILE):
        print("Need to build the register prediction code, go to bb-register/pred-combine in a console and then type 'make'")

    print ('Predicting registers...')
    #Todo add options to test all levels, or select one
    find_best_combination(test_file_path,score_output_file, level )

    print 'Optimizing local register for better global shear...'
    outfile = "shear_adjustments"
    get_shear_adjustments(test_file_path, TMP_DIR, outfile)


# Do this so we can run our code from the console
def LoadSearchPaths():
    for i, j, y in walk(join(getcwd(), "..")):
        if (str(i).find('__pycache__') == -1 and str(i).find('.vscode') == -1):
            sys.path.append(i)

if __name__ =="__main__":
    LoadSearchPaths()
    pdbs = ["unkn"]
    run_register_prediction(pdbs)