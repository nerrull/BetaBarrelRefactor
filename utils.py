

import sys
from os import walk, getcwd
from os.path import join, exists
from defs import INPUT_DIR

def check_input_files(pdbs):
    required_ext = [".res", ".strands", ".tmregs", ".tmstrands"]
    for pdb in pdbs:
        input_path = join(INPUT_DIR, pdb)
        if not exists(input_path):
            print "Missing input directory for {}".format(pdb)
            return False
        else:
            for ext in required_ext:
                f = pdb + ext
                required_file = join(input_path, f)
                if not exists(required_file):
                    print "Missing file for  {}: {}".format(pdb, f)
                    return False
    return True


# Do this so we can run our code from the console
def loadSearchPaths():
    for i, j, y in walk(getcwd()):
        if (str(i).find('__pycache__') == -1 and str(i).find('.vscode') == -1):
            sys.path.append(i)