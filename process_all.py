import sys
import re
import random
from os import mkdir
from os.path import exists, join
from defs import INPUT_DIR, OUTPUT_DIR

aad = dict()
aad["A"] = 0
aad["R"] = 1
aad["N"] = 2
aad["D"] = 3
aad["C"] = 4
aad["Q"] = 5
aad["E"] = 6
aad["G"] = 7
aad["H"] = 8
aad["I"] = 9
aad["L"] = 10
aad["K"] = 11
aad["M"] = 12
aad["F"] = 13
aad["P"] = 14
aad["S"] = 15
aad["T"] = 16
aad["W"] = 17
aad["Y"] = 18
aad["V"] = 19


def convert_string(seq, out_fn):
    with open(out_fn, "w") as f:
        for c in seq:
            c = c.upper()
            if c in aad:
                f.write("{0}\n".format(aad[c]))
            else:
                raise ValueError("{0} is not a valid aminoacid letter!".format(c))


def get_seqs(fn):
    seq = ""
    tm = ""

    with open(fn) as f:
        l = f.readline()
        if re.search('^>', l):

            l = f.readline().rstrip()
            seq += l
            l = f.readline().rstrip()
            while (len(l) > 0):
                if re.search('^>', l):
                    break
                seq += l
                l = f.readline().rstrip()
            l = f.readline().rstrip()
            tm += l
            l = f.readline().rstrip()
            while (len(l) > 0):
                tm += l
                l = f.readline().rstrip()
        else:
            raise ValueError("Wrong file format for the sequences\n")

    return (seq, tm)


def read_sequence(f):
    seq = ""
    l = f.readline().rstrip()
    while (len(l) > 0):
        if re.search('^>', l):
            break
        seq += l
        l = f.readline().rstrip()
    return  f,seq, l

def get_sequences(fn):
    returnDict = {}

    with open(fn) as f:
        l = f.readline()
        if re.search('^>', l):

            while len(l) >0:
                name = re.search('\|([A-Za-z_0-9]+)\s*$', l).group(1)
                f, seq , _= read_sequence(f)
                f, tm, l = read_sequence(f)
                returnDict[name] = (seq, tm)
        else:
            raise ValueError("Wrong file format for the sequences\n")

    return returnDict



def write_regs(seq, tm, out_fn):
    str_fn = "{0}.strands".format(out_fn)
    tmstr_fn = "{0}.tmstrands".format(out_fn)
    regs_fn = "{0}.tmregs".format(out_fn)
    with open(str_fn, "w") as f:
        with open(tmstr_fn, "w") as ftm:
            with open(regs_fn, "w") as regs:
                if len(seq) != len(tm):
                    raise ValueError("Different lengths for sequence and TM regions")
                tm_flag = False
                start = None
                for i, (a, b) in enumerate(zip(seq, tm)):
                    if tm_flag and b != "M":
                        tm_flag = False
                        f.write("{0} {1}\n".format(start, i))
                        ftm.write("{0} {1}\n".format(start, i))
                        regs.write("0\n")
                    if b == "M" and not tm_flag:
                        start = i
                        tm_flag = True
                if tm_flag:
                    f.write("{0} {1}\n".format(start, len(seq)))
                    ftm.write("{0} {1}\n".format(start, len(seq)))
                    regs.write("0\n")


def fisher_yates(seq):
    l = list(seq)
    for i in range(len(l) - 1):
        j = random.randint(0, i)  # random number 0 <= j <= i
        tmp = l[j]
        l[j] = l[i]
        l[i] = tmp
    seq = "".join(l)
    return seq


def prepare_input(fasta_file, scramble = False):
    sequenceDict= get_sequences(fasta_file)

    for sequence_name, (seq, tm) in sequenceDict.iteritems():
        out_dir = join(INPUT_DIR, sequence_name)
        if not exists(out_dir):
            mkdir(out_dir)
        write_regs(seq, tm, "{0}/{1}".format(out_dir, sequence_name))
        convert_string(seq, "{0}/{1}.res".format(out_dir, sequence_name))

        with open("{0}/{1}.seq".format(out_dir, sequence_name), "w") as f:
            f.write(seq)
            f.write("\n")

        #Generate scrambled control sequence
        if scramble:
            out_dir = join(INPUT_DIR, sequence_name+"_scrambled")
            if not exists(out_dir):
                mkdir(out_dir)
            write_regs(seq, tm, "{0}/{1}_scrambled".format(out_dir, sequence_name))
            scramble = fisher_yates(seq)
            convert_string(scramble, "{0}/{1}_scrambled.res".format(out_dir, sequence_name))

            with open("{0}/{1}_scrambled.seq".format(out_dir, sequence_name), "w") as f:
                f.write(scramble)
                f.write("\n")

    returnKeys =sequenceDict.keys()
    if scramble:
        scrambled_pdbs = [pdb +"_scrambled" for pdb in sequenceDict.keys()]
        returnKeys = sequenceDict.keys() + scrambled_pdbs

    return returnKeys


if __name__=="__main__":
    prepare_input("fasta_test.txt")


