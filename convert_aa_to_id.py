import sys

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
                raise ValueError("{0} is not a valide aminoacid letter!".format(c))

if __name__ == "__main__":
    if (len(sys.argv) != 3):
        raise ValueError("I need 2 args:  the sequence file and the output name\n")
    in_fn = sys.argv[1]
    out_fn = sys.argv[2]
    with open(in_fn) as f:
        line = f.readline().rstrip()
        convert_string(line, out_fn)