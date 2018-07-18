#!/usr/bin/python

''' this script is to calculate the reg-ec scores
'''

import argparse
import glob
import sys
import math
from parse_ec import get_circular_barrel, get_strand_ranges, get_neighbor_ec_dict
import regpy
from os.path import join

def get_ecs(pdb, peris, ecfn, weights, theta, log):
    ## get ec dict
    ec_dict = get_neighbor_ec_dict(pdb, ecfn)
    ## construct circular barrel
    barrel = get_circular_barrel(pdb)

    rtn = {}
    ## scan all the strands
    for i in range(1, len(barrel) - 1):
        rtn[i] = []
        offset1 = barrel[i][0] - peris[i - 1]
        try:
            offset2 = barrel[i + 1][0] - peris[i]
        except IndexError:
            offset2 = barrel[i + 1][0] - peris[0]
        ## enumerate registration
        for reg in range(-len(barrel[i]) + 1, len(barrel[i + 1])):
            strand1, strand2 = regpy.construct_hairpin(barrel[i], barrel[i + 1], reg)
            ## adjustment of reg to interface
            if i % 2 == 0:
                reg = reg + offset1 + offset2
            else:
                reg = reg - offset1 - offset2

            ## get distant contact list
            dist_cont_list = []
            dist_cont_list.append(regpy.get_dist_contacts(strand1, strand2, 0))
            dist_cont_list.append(regpy.get_dist_contacts(strand1, strand2, 1))
            dist_cont_list.append(regpy.get_dist_contacts(strand1, strand2, 2))

            score = 0
            contnum = 0
            for j in range(len(dist_cont_list)):
                conts = dist_cont_list[j]
                for cont in conts:
                    ## strand indexes are i-1, i instead of i, i+1
                    try:
                        prob = ec_dict[i - 1, i % (len(barrel) - 2)][cont]
                        if prob < theta:
                            continue
                        if log == 1:
                            score -= weights[j] * math.log(prob)
                        else:
                            score += weights[j] * prob
                        contnum += 1
                    except KeyError:
                        pass

            rtn[i].append([reg, score, contnum])
    if log == 2:
        minp = 1e-10
        for i in range(1, len(barrel) - 1):
            tot = 0.0
            for dummy, score, dummy in rtn[i]:
                tot += score
            if tot == 0.0:
                continue
            for j in range(len(rtn[i])):
                if rtn[i][j][1] == 0.0:
                    rtn[i][j][1] = math.log(minp / tot)
                else:
                    rtn[i][j][1] = math.log(rtn[i][j][1] / tot)
    return rtn


##################################################################
############################ MAIN ################################
##################################################################
# Write the scores to the tmp file and return them as a dict
def get_scores(test_dict, weights, log, theta, output_file, ec_raw_dir):
    dummy = sum(weights)
    wtag = '_' + str(weights[0]).zfill(3) + str(weights[1]).zfill(3) + str(weights[2]).zfill(3)

    if log == 1:
        ltag = '_lp'
    elif log == 2:
        ltag = '_lnp'
    else:
        ltag = '_p'

    ttag = '_%1.2E' % theta

    score_dict = {}
    fout = open(output_file, 'w')
    for pdb in test_dict.keys():
        print ("------------------------PDB : {}".format(pdb))
        print ("------------------------ec_raw_dir : {}".format(ec_raw_dir))

        ec_fn = join(ec_raw_dir, "empty.psicov")
        psicov_files = glob.glob(join(ec_raw_dir, pdb) + '*')
        if len(psicov_files) !=0:
            ec_fn = psicov_files[0]
        else:
            print("------------------------ No psicov file found for {}".format(pdb))

        ecs = get_ecs(pdb, test_dict[pdb]["peris"], ec_fn, weights, theta, log)
        score_dict[pdb] = ecs

        fout.write(pdb + '\n')
        for i in sorted(ecs.keys()):
            fout.write(str(len(ecs[i])))
            for itm in ecs[i]:
                fout.write('\t' + str(itm[0]) + '\t' + str(itm[1]))
            fout.write('\n')
    fout.close()


