import glob
import os
import sys
import numpy as np
from defs import INPUT_DIR
from os.path import join

sys.dont_write_bytecode = True


# construct a barrel. for an N strands barrel, it will be like
# [ 
#	strand N-1 range
#	strand 0   range
#	strand 1   range
#	......
#	......
#	strand N-2 range
#	strand N-1 range
#	strand 0   range
# ]
# the range follows seq id order but not peri to extra
# Notice: it is not the barrel, just the ranges
def get_strand_ranges(pdb):
    ranges = []
    path = join(join(INPUT_DIR, pdb), "{}.strands".format(pdb))
    f = open(path)
    lines = f.readlines()
    f.close()
    for line in lines:
        start, end = line.split()
        ranges.append((int(start), int(end) + 1))
    ranges.append(ranges[0])
    ranges = [ranges[-2]] + ranges
    return ranges


def get_tmstrand_ranges(pdb):
    ranges = []
    f = open('../bb-pipe/inputs/' + pdb + '/' + pdb + '.tmstrands')
    lines = f.readlines()
    f.close()
    for line in lines:
        start, end = line.split()
        ranges.append((int(start), int(end) + 1))
    ranges.append(ranges[0])
    ranges = [ranges[-2]] + ranges
    return ranges


# construct a circular barrel
# this function is similar with  get_strand_ranges(pdb)
def get_circular_barrel(pdb):
    strand_ranges = get_strand_ranges(pdb)
    barrel = []
    for i in range(len(strand_ranges)):
        if i % 2 == 0:
            strand = list(reversed(range(strand_ranges[i][0], strand_ranges[i][1])))
        else:
            strand = range(strand_ranges[i][0], strand_ranges[i][1])
        barrel.append(strand)
    return barrel


# get a dict of interstrand ecs
# the dict is like
# {
#	(strand0, strand1) : { (resi, resj):prob }
#	(strand1, strand2) : { (resi, resj):prob }
#	......
#	(strandN-1, strand0) : { (resi, resj):prob }
# }
def get_neighbor_ec_dict(pdb, ecfn):
    strand_ranges = get_strand_ranges(pdb)
    strandn = len(strand_ranges) - 2

    # ec_dict = { (i,i+1):{} for i in range(strandn-1) }
    ec_dict = {}
    for i in range(strandn - 1):
        ec_dict[i, i + 1] = {}
    ec_dict[strandn - 1, 0] = {}

    if os.path.exists(ecfn):
        ecf = open(ecfn)
        lines = ecf.readlines()
        ecf.close()

        # parse ec data
        for line in lines:
            splits = line.split()
            res1 = int(splits[0])
            res2 = int(splits[1])
            if res1 > res2:
                res1, res2 = res2, res1
            prob = float(splits[-1])

            # check if the contact is a neighbour strand contact
            for i in range(1, len(strand_ranges) - 1):
                if res1 >= strand_ranges[i][0] and res1 < strand_ranges[i][1]:
                    if res2 >= strand_ranges[i + 1][0] and res2 < strand_ranges[i + 1][1]:
                        strand1 = i - 1
                        strand2 = i % (len(strand_ranges) - 2)
                        ec_dict[strand1, strand2][res1, res2] = prob
                    elif res2 >= strand_ranges[i - 1][0] and res2 < strand_ranges[i - 1][1]:
                        strand1 = i - 1
                        strand2 = (i - 2) % (len(strand_ranges) - 2)
                        if strand1 == 0 and strand2 != 1:
                            strand1, strand2 = strand2, strand1
                            res1, res2 = res2, res1
                        ec_dict[strand1, strand2][res1, res2] = prob

    return ec_dict


if __name__ == "__main__":
    # get_neighbor_ec_dict('1bxw', ['ec_results/ec_coinfold/1bxwA.pred','ec_results/ec_metapsicov/1bxwA.pred','ec_results/ec_psicov/1bxw.psicov'], 'avg')
    a = get_neighbor_ec_dict('1bxw', 'ec_results/ec_coinfold.hammad/1bxwA.pred', 'avg')
    for k in sorted(a.keys()):
        print k
        for i in a[k]:
            print '\t', i, a[k][i]
