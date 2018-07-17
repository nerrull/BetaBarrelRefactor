#!/usr/bin/python

import os
import random
from os.path import join

import numpy as np
import pandas

from defs import HAMMAD_DIR, INPUT_DIR, TMP_DIR


# TODO compare changes with original code
def get_true_facings(pdb, periids):
    cens = np.loadtxt(join(INPUT_DIR, *[pdb, pdb + '.cen'])).astype(int)
    # assert len(cens)==len(periids)
    with open(join(INPUT_DIR, *[pdb, pdb + '.sidechain'])) as f:
        lines = f.readlines()

    true_facings = {}
    for line in lines:
        split = line.split()
        periid = int(split[0])
        fac = split[1]
        true_facings[periid] = fac

    true_cen_facings = []
    for cen in cens:
        if cen in true_facings:
            true_cen_facings.append(true_facings[cen])
        elif cen + 2 in true_facings:
            true_cen_facings.append(true_facings[cen + 2])
        else:
            true_cen_facings.append(true_facings[cen - 2])

    # determine true peri facings based on true cen facings
    true_peri_facings = []
    for i in range(len(periids)):
        if (periids[i] - cens[i]) % 2 != 0:
            if true_cen_facings[i] == 'IN':
                true_peri_facings.append('OUT')
            else:
                true_peri_facings.append('IN')
        else:
            true_peri_facings.append(true_cen_facings[i])

    return true_peri_facings


# dirn = '../pred_combine/odds/meishan60/'
dirn = HAMMAD_DIR
oddsei = np.loadtxt(join(dirn, 'ExtraIn.odds'))
oddseo = np.loadtxt(join(dirn, 'ExtraOut.odds'))
oddspi = np.loadtxt(join(dirn, 'PeriIn.odds'))
oddspo = np.loadtxt(join(dirn, 'PeriOut.odds'))
oddsci = np.loadtxt(join(dirn, 'CoreIn.odds'))
oddsco = np.loadtxt(join(dirn, 'CoreOut.odds'))


def determine_facings_odds(pdb, periids, extraids):
    res = np.loadtxt(join(INPUT_DIR, *[pdb, pdb + '.res'])).astype(int)
    predicted_peri_facings = []
    for i in range(len(periids)):
        if periids[i] < extraids[i]:
            strand = range(periids[i], extraids[i] + 1, 1)
        else:
            strand = range(periids[i], extraids[i] - 1, -1)

        for j in range(len(strand)):
            strand[j] = res[strand[j] - 1]

        bestscore = -10000000000000
        for peri_num, core_num, extra_num in [(2, 5, 2), (2, 5, 3), (2, 6, 2), (3, 5, 2), (2, 6, 3), (3, 6, 2),
                                              (3, 5, 3), (3, 6, 3), (2, 5, 1), (2, 4, 2), (1, 5, 2), (2, 4, 1),
                                              (1, 5, 1), (1, 4, 2)]:
            for peri_facing_out in [0, 1]:
                score = 0.0
                for j in range(peri_num):
                    try:
                        strand[j]
                    except IndexError:
                        continue
                    if j % 2 == peri_facing_out:
                        score += oddspi[strand[j]]
                    else:
                        score += oddspo[strand[j]]
                for j in range(peri_num, peri_num + core_num):
                    try:
                        strand[j]
                    except IndexError:
                        continue
                    if j % 2 == peri_facing_out:
                        score += oddsci[strand[j]]
                    else:
                        score += oddsco[strand[j]]
                for j in range(peri_num + core_num, peri_num + core_num + extra_num):
                    try:
                        strand[j]
                    except IndexError:
                        continue
                    if j % 2 == peri_facing_out:
                        score += oddsei[strand[j]]
                    else:
                        score += oddseo[strand[j]]
                if score > bestscore:
                    bestscore = score
                    best_peri_facing_out = peri_facing_out
        if best_peri_facing_out:
            predicted_peri_facings.append('OUT')
        else:
            predicted_peri_facings.append('IN')
    return predicted_peri_facings


getfep_score_dict = {0: 0, 1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 7: -1, 8: -1, 9: 0, 10: 0, 11: -1, 12: 0, 13: 0,
                     14: -1, 15: -1, 16: -1, 17: 0, 18: 0, 19: 0}


def determine_facings_GeTFEP(pdb, periids, extraids):
    res = np.loadtxt('../inputs/' + pdb + '/' + pdb + '.res').astype(int)
    predicted_peri_facings = []

    for i in range(len(periids)):
        if periids[i] < extraids[i]:
            strand = range(periids[i], extraids[i] + 1, 1)
        else:
            strand = range(periids[i], extraids[i] - 1, -1)

        for j in range(len(strand)):
            strand[j] = res[strand[j] - 1]

        periout_score = 0
        for j in [0, 2, 4]:
            try:
                periout_score += getfep_score_dict[strand[j]]
            except:
                pass
        for j in [1, 3, 5]:
            try:
                periout_score -= getfep_score_dict[strand[j]]
            except:
                pass

        if periout_score >= 0:
            predicted_peri_facings.append('OUT')
        else:
            predicted_peri_facings.append('IN')

    return predicted_peri_facings


def get_possible_targetshears(pdbdata, pdb):
    N = len(pdbdata)
    if N <= 8:
        targetshears = [N, N + 2]
    elif N <= 10:
        targetshears = [N + 2]
    elif N <= 12:
        targetshears = [N + 2, N + 4]
    elif N <= 14:
        targetshears = [N, N + 2]
    elif N <= 16:
        targetshears = [N + 4]
    elif N <= 18:
        targetshears = [N + 2, N + 4]
    elif N <= 24:
        targetshears = [N + 2]
    elif N <= 26:
        targetshears = [N + 4]


def get_common_targetshears(pdbdata, pdb):
    N = len(pdbdata)
    if N <= 12:
        targetshears = [N + 2]
    elif N <= 14:
        targetshears = [N]
    elif N <= 18:
        targetshears = [N + 4]
    elif N <= 24:
        targetshears = [N + 2]
    elif N <= 26:
        targetshears = [N + 2]
        # targetshears = [N+4]
    return targetshears


def get_group_targetshears(pdbdata, pdb):
    print("SHEAR PDB PLEASE : {}".format(pdb))
    N = len(pdbdata)
    ## grpshear
    l01 = ['1bxw', '1qj8', '1p4t', '2f1t', '1thq', '2erv', '2lhf', '2mlh', '3dzm', '1qd6', '2f1c', '1k24', '1i78',
           '2wjr', '4pr7', 'unkn']
    l02 = ['1t16', '1uyn', '1tly', '3aeh', '3bs0', '3dwo', '3fid', '3kvn', '4e1s']
    l03 = ['2mpr', '1a0s', '2omf', '2por', '1prn', '1e54', '2o4v', '3vzt', '4k3c', '4k3b', '4c4v', '4n75']
    l04 = ['2qdz', '2ynk', '3rbh', '3syb', '3szv', '4c00', '4gey']
    l05 = ['1fep', '2fcp', '1kmo', '1nqe', '1xkw', '2vqi', '3csl', '3rfz', '3v8x', '4q35']
    if pdb in l01:
        targetshears = [N + 2]
    elif pdb in l02:
        # targetshears = [N]
        targetshears = [N + 2]
        # targetshears = [N,N+2]
    elif pdb in l03:
        targetshears = [N + 4]
    elif pdb in l04:
        targetshears = [N + 4]
    elif pdb in l05:
        targetshears = [N + 2]
    return targetshears


def show_tagged_data(pdbdata):
    newdata = pdbdata.tolist()
    xcount = 0
    scount = 0
    for i in range(len(newdata)):
        # newdata[i] = newdata[i].tolist()
        tag = ''
        if newdata[i][4] == newdata[i][6] and newdata[i][3] != newdata[i][6]:
            tag = '**'
            scount += 1
        elif newdata[i][3] != newdata[i][6]:
            tag = 'xx'
            xcount += 1
        newdata[i].append(tag)
    newdata = np.array(newdata)
    labs = ['diff', 'score0', 'score1', 'preg0', 'preg1', 'pfac', 'treg', 'tfac', 'strdid', 'seqid', 'tag']
    newdata = pandas.DataFrame(newdata, columns=labs)
    print newdata
    return scount, xcount


def random_improvement(pdbdata, imprv):
    N = len(pdbdata)
    newdata = []
    for i in range(N):
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            pdbdata[i]
        if pred0_reg != truereg:
            if random.random() <= imprv:
                score_diff = 0
                pred0_score = 1000
                pred1_score = 1000
                pred0_reg = truereg
                pred1_reg = truereg
        newdata.append(
            (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri))
    newdata = np.array(sorted(newdata))
    return newdata


def parity_first_shear_adjustment(pdbdata, targetshears):
    N = len(pdbdata)
    currshear = sum(pdbdata[:, 3].astype(int))

    newdata = []
    for i in range(N):
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            pdbdata[i]
        score_diff = float(score_diff)
        pred0_reg = int(pred0_reg)
        pred1_reg = int(pred1_reg)
        regmod01 = (pred0_reg - pred1_reg) % 2

        if currshear < targetshears[0]:
            if pred0_reg < pred1_reg:
                ## adjust shear parity first unless top 2 regs can make the shear larger
                if i < 2 or (currshear % 2 == regmod01 and targetshears[-1] - currshear >= pred1_reg - pred0_reg):
                    # #adjust shear parity first
                    # if (currshear%2==regmod01 and targetshears[-1]-currshear>=pred1_reg-pred0_reg):
                    currshear += pred1_reg - pred0_reg
                    pred0_reg = pred1_reg
                    pred0_score = pred1_score
                    score_diff = 0
        elif currshear > targetshears[-1]:
            if pred0_reg > pred1_reg:
                # #adjust shear parity first unless top 2 regs can make the shear smaller
                if i < 2 or (currshear % 2 == regmod01 and targetshears[0] - currshear <= pred1_reg - pred0_reg):
                    # #adjust shear parity first
                    # if (currshear%2==regmod01 and targetshears[0]-currshear<=pred1_reg-pred0_reg):
                    currshear += pred1_reg - pred0_reg
                    pred0_reg = pred1_reg
                    pred0_score = pred1_score
                    score_diff = 0

        newdata.append(
            (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri))

    newdata = np.array(sorted(newdata))
    return newdata


def brute_shear_adjustment(alldata, targetshears):
    N = len(alldata)
    currshear = sum(alldata[:, 3].astype(int))

    newdata = []
    for i in range(N):
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            alldata[i]
        score_diff = float(score_diff)
        pred0_reg = int(pred0_reg)
        pred1_reg = int(pred1_reg)
        peri = int(peri)

        if currshear < targetshears[0]:
            if pred0_reg < pred1_reg:
                if i < 2 or targetshears[-1] - currshear >= pred1_reg - pred0_reg:
                    currshear += pred1_reg - pred0_reg
                    pred0_reg = pred1_reg
                    pred0_score = pred1_score
                    score_diff = 0
        elif currshear > targetshears[-1]:
            if pred0_reg > pred1_reg:
                if i < 2 or targetshears[0] - currshear <= pred1_reg - pred0_reg:
                    currshear += pred1_reg - pred0_reg
                    pred0_reg = pred1_reg
                    pred0_score = pred1_score
                    score_diff = 0

        newdata.append(
            (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri))

    newdata = np.array(sorted(newdata))
    return newdata


def get_final(alldata):
    N = len(alldata)
    newdata = []
    for i in range(N):
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            alldata[i]
        strandid = int(strandid)
        newdata.append((strandid, pred0_reg, predface, truereg, trueface, peri))
    newdata = sorted(newdata)
    return np.array(newdata)


## correct facing according to reg
##    what it actually does:
##    if 2 consuctive pred_reg0 needs to be corrected by pred_reg1 according to the pred_fac, it is probably that the pred_fac is wrong
def correct_fac(pdbdata):
    newdata = []
    N = len(pdbdata)
    facingcorrecting_cond1_strandids = []
    facingcorrecting_cond2_strandids = []

    for i in range(N):
        ip1 = (i + 1) % N
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            pdbdata[i]
        dummy, dummy, dummy, dummy, dummy, predface_ip1, dummy, dummy, dummy, dummy = pdbdata[ip1]
        pred0_reg = int(pred0_reg)
        pred1_reg = int(pred1_reg)
        truereg = int(truereg)
        if predface == predface_ip1:
            if pred0_reg % 2 != 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                ### facing correction condition 1
                # print '*wew', i, predface, trueface
                # print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
                facingcorrecting_cond1_strandids.append(i)
            ###------------------------------------------
            elif pred0_reg % 2 == 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                pass
            ###------------------------------------------
            elif pred0_reg % 2 != 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                ### facing correction condition 2
                # print '*hew', i, predface, trueface
                # print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
                facingcorrecting_cond2_strandids.append(i)
                ###------------------------------------------
        else:  ## predface != predface_ip1
            if pred0_reg % 2 == 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                ### facing correction condition 2
                # print '*how', i, predface, trueface
                # print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
                facingcorrecting_cond2_strandids.append(i)
            ###------------------------------------------
            elif pred0_reg % 2 != 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                pass
            ###------------------------------------------
            elif pred0_reg % 2 == 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                ### facing correction condition 1
                # print '*wow', i, predface, trueface
                # print ' ', (i+1)%N, pdbdata[(i+1)%N][5], pdbdata[(i+1)%N][7]
                facingcorrecting_cond1_strandids.append(i)
                ###------------------------------------------
        newdata.append(
            (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri))

    ## correct facings here
    for i in range(len(newdata)):
        im1 = (i - 1) % N
        if i in facingcorrecting_cond1_strandids and im1 in facingcorrecting_cond1_strandids:
            if newdata[i][5] == 'OUT':
                newdata[i][5] == 'IN'
            else:
                newdata[i][5] == 'OUT'
                # if i in facingcorrecting_cond2_strandids and im1 in facingcorrecting_cond2_strandids:
    # if newdata[i][5] == 'OUT':
    #			newdata[i][5] == 'IN'
    #		else:
    #			newdata[i][5] == 'OUT'

    return np.array(newdata)


def filter_by_fac(pdbdata):
    newdata = []
    N = len(pdbdata)

    for i in range(N):
        ip1 = (i + 1) % N
        score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri = \
            pdbdata[i]
        dummy, dummy, dummy, dummy, dummy, predface_ip1, dummy, dummy, dummy, dummy = pdbdata[ip1]
        pred0_reg = int(pred0_reg)
        pred1_reg = int(pred1_reg)
        truereg = int(truereg)
        if predface == predface_ip1:
            if pred0_reg % 2 != 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                ### filter regs using facing
                score_diff = 0
                pred0_score = pred1_score
                pred0_reg = pred1_reg
                # print pdb,i,'wew'
                pass
            ###------------------------------------------
            elif pred0_reg % 2 == 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                ### filter regs using facing
                score_diff = 0
                pred1_score = pred0_score
                pred1_reg = pred0_reg
                pass
            ###------------------------------------------
            elif pred0_reg % 2 != 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                pass
                ###------------------------------------------
        else:  ## predface != predface_ip1
            if pred0_reg % 2 == 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                pass
            ###------------------------------------------
            elif pred0_reg % 2 != 0 and pred1_reg % 2 == 0:
                ###------------------------------------------
                ### filter regs using facing
                score_diff = 0
                pred1_score = pred0_score
                pred1_reg = pred0_reg
                # print pdb,i,'how'
                pass
            ###------------------------------------------
            elif pred0_reg % 2 == 0 and pred1_reg % 2 != 0:
                ###------------------------------------------
                ### in this case, correcting pred0 does more bad than good, according to manual counted results
                # score_diff = 0
                # pred0_score = pred1_score
                # pred0_reg = pred1_reg
                # print pdb,i,'wow'
                pass
                ###------------------------------------------
        newdata.append(
            (score_diff, pred0_score, pred1_score, pred0_reg, pred1_reg, predface, truereg, trueface, strandid, peri))

    return np.array(newdata)


'''

		print 'Usage:',sys.argv[0],'<test_file> <score_file_dir> [optimiaztion#] [output_file]'
		#print 'Usage:',sys.argv[0],'<test_file> <score_file> [optimiaztion#] [output_file]'
		print 'Optimization #:'
		print '0             : fac_crct'
		print '1             : fac_crct + fac_fltr'
		print '2             : fac_crct + fac_fltr + parity'
		print '3             : fac_crct + parity'
		print '4             : fac_crct + fac_fltr + brute'
		print '5             : fac_crct + brute'
		print '6             : fac_crct + fac_fltr + parity + brute'
		print '7             : fac_crct + parity + brute'
		print '101~114       : rand_impr + fac_crct + fac_fltr + parity'
		print '                for 1xx, randomly correct xx% wrong predictions first, then do the others'
		print 'o/w (default) : original'
		sys.exit(0)
'''


# Note by Etienne
# We had to scalpel out all the code that compared the predicted results to the ground truth
# because we don't have the ground truth information for new protein structures

def shear_adjust(testfn, scoredir, optnum, out_file):
    os.system('cat ' + join(scoredir, '*.l0?') + ' > ' + join(scoredir, 'all.scores'))
    scorefn = join(scoredir, 'all.scores')

    out_file = join(TMP_DIR, out_file + '.rlt')
    fout = open(out_file, 'w')

    ### read test file, and predict facing
    ### testdict restores peris used to predict regs and the corresponding true regs
    ### sizedict restores number of strands
    ### perifacing_dict restores facings predicted
    with open(testfn) as f:
        lines = f.readlines()
    sizedict = {}
    testdict = {}
    perifacing_dict = {}

    correct_facing_odds = 0
    correct_facing_getfep = 0
    for line in lines:
        split = line.split()
        pdb = split[0]
        strandnum = int(split[1])
        sizedict[pdb] = strandnum
        testdict[pdb] = []
        peris = []
        extras = []
        for i in range(2, strandnum * 3 + 2, 3):
            peri = int(split[i])
            extra = int(split[i + 1])
            peris.append(peri)
            extras.append(extra)
            trueperireg = int(split[i + 2])
            testdict[pdb].append((peri, trueperireg))

        pred_facing = determine_facings_odds(pdb, peris, extras)  # accuracy 723/736
        # pred_facing = determine_facings_GeTFEP(pdb, peris, extras) # accuracy 689/736
        true_facing = pred_facing  # get_true_facings(pdb, peris)
        perifacing_dict[pdb] = zip(pred_facing, true_facing)
    sizedict['3emn'] = 18  # special treatment for vdac as its strand number is odd

    ###################################
    ### adjustment starts from here ###
    ###################################

    correct_facing0 = 0  # num before adjustment
    correct_facing1 = 0  # num after adjustment
    correct_reg0 = 0  # num before adjustment
    correct_reg1 = 0  # num after adjustment
    total_strand = 0
    badshear = []

    print scoredir, optnum
    #### read score file
    with open(scorefn) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith('## '):
                pdb = line.split()[1]

                ## The following loop read info of the pdb from score files and compile them with the info get above
                pdbdata_ori = []
                for strandid in range(sizedict[pdb]):
                    line = f.readline()
                    split = line.split()
                    strandpairs = []
                    for pair in split:
                        perireg, score = pair.split(':')
                        perireg = int(perireg)
                        score = float(score)
                        strandpairs.append((score, perireg))
                    cand1, cand0 = sorted(strandpairs)[-2:]
                    trueperireg = testdict[pdb][strandid][1]
                    peri = testdict[pdb][strandid][0]
                    perifacing = perifacing_dict[pdb][strandid][0]
                    trueperifacing = perifacing_dict[pdb][strandid][1]
                    ## original data without adjustment
                    ## 0,			1,		2,		3,			4,			5,			6,			7,,			8,			9
                    ## score_diff,	score0,	score1,	pred_reg0,	pred_reg1,	pred_fac,	true_reg,	true_face,	strand_id,	peri_seqid
                    pdbdata_ori.append((float('%.3f' % ((cand0[0] - cand1[0]) / abs(cand0[0] + cand1[0]))),
                                        float('%.3f' % (cand0[0])), float('%.3f' % (cand1[0])), cand0[1], cand1[1],
                                        perifacing, trueperireg, trueperifacing, strandid, peri))
                # all these tried c0+c1, abs(c0+c1), c0, 1
                # best is abs(c0+c1)
                # pdbdata_ori.append( (cand0[0]-cand1[0],cand0[0],cand1[0],cand0[1],cand1[1],perifacing,trueperireg,trueperifacing,strandid,peri) )

                ################################################
                ### process info of this pdb start from here ###
                ################################################

                ## get target shear
                targetshears = get_common_targetshears(pdbdata_ori, pdb)
                # targetshears = get_group_targetshears(pdbdata_ori,pdb)

                ## get statistics of original prediction
                pdbdata_ori = np.array(pdbdata_ori)
                correct_facing0 += sum(pdbdata_ori[:, 5] == pdbdata_ori[:, 7])
                correct_reg0 += sum(pdbdata_ori[:, 3] == pdbdata_ori[:, 6])
                total_strand += len(pdbdata_ori)

                ## test
                # labs = [ 'diff', 'scr0', 'scr1', 'preg0', 'preg1', 'pfac', 'treg', 'tfac', 'strdi', 'peri' ]
                # print pandas.DataFrame(pdbdata_ori, columns=labs)
                # print '-'*30

                shear0 = sum(pdbdata_ori[:, 3].astype(int))

                if optnum == 0:
                    ## optnum 0 : no adjustment, just a little bit facing correction
                    pdbdata0 = np.array(pdbdata_ori)
                    tofinal = correct_fac(pdbdata0)

                elif optnum == 1:
                    ## optnum 1 : facing filtering
                    pdbdata1 = np.array(pdbdata_ori)
                    pdbdata1 = correct_fac(pdbdata1)
                    tofinal = filter_by_fac(pdbdata1)

                elif optnum == 2:
                    ## optnum 2 : facing filtering + parity first adjustment
                    pdbdata2 = np.array(pdbdata_ori)
                    pdbdata2 = correct_fac(pdbdata2)
                    pdbdata2 = filter_by_fac(pdbdata2)
                    pdbdata2 = np.array(sorted(pdbdata2.tolist()))
                    tofinal = parity_first_shear_adjustment(pdbdata2, targetshears)

                elif optnum == 3:
                    ## optnum 3 : parity first adjustment
                    pdbdata3 = np.array(pdbdata_ori)
                    pdbdata3 = correct_fac(pdbdata3)
                    pdbdata3 = np.array(sorted(pdbdata3.tolist()))
                    tofinal = parity_first_shear_adjustment(pdbdata3, targetshears)

                elif optnum == 4:
                    ## optnum 4 : facing filtering + brute adjustment
                    pdbdata4 = np.array(pdbdata_ori)
                    pdbdata4 = correct_fac(pdbdata4)
                    pdbdata4 = filter_by_fac(pdbdata4)
                    pdbdata4 = np.array(sorted(pdbdata4.tolist()))
                    tofinal = brute_shear_adjustment(pdbdata4, targetshears)

                elif optnum == 5:
                    ## optnum 5 : brute adjustment
                    pdbdata5 = np.array(pdbdata_ori)
                    pdbdata5 = correct_fac(pdbdata5)
                    pdbdata5 = np.array(sorted(pdbdata5.tolist()))
                    tofinal = brute_shear_adjustment(pdbdata5, targetshears)

                elif optnum == 6:
                    ## optnum 6 : facing filtering + parity first adjustment + brute adjustment
                    pdbdata6 = np.array(pdbdata_ori)
                    pdbdata6 = correct_fac(pdbdata6)
                    pdbdata6 = filter_by_fac(pdbdata6)
                    pdbdata6 = np.array(sorted(pdbdata6.tolist()))
                    pdbdata6 = parity_first_shear_adjustment(pdbdata6, targetshears)
                    tofinal = brute_shear_adjustment(pdbdata6, targetshears)

                elif optnum == 7:
                    ## optnum 7 : parity first adjustment + brute adjustment
                    pdbdata7 = np.array(pdbdata_ori)
                    pdbdata7 = correct_fac(pdbdata7)
                    pdbdata7 = np.array(sorted(pdbdata7.tolist()))
                    pdbdata7 = parity_first_shear_adjustment(pdbdata7, targetshears)
                    tofinal = brute_shear_adjustment(pdbdata7, targetshears)

                elif optnum.isdigit() and int(optnum) >= 99:
                    ## optnum 1xx : random improvement + facing filtering + parity first adjustment + brute adjustment
                    imprv = (float(optnum) - 100) / 14.0
                    imprv = min(1.0, imprv)
                    pdbdata100 = np.array(pdbdata_ori)
                    pdbdata100 = random_improvement(pdbdata100, imprv)
                    pdbdata100 = correct_fac(pdbdata100)
                    pdbdata100 = filter_by_fac(pdbdata100)
                    pdbdata100 = np.array(sorted(pdbdata100.tolist()))
                    pdbdata100 = parity_first_shear_adjustment(pdbdata100, targetshears)
                    tofinal = brute_shear_adjustment(pdbdata100, targetshears)

                else:
                    ## optnum m1 : no adjustment
                    pdbdatam1 = np.array(pdbdata_ori)
                    tofinal = pdbdatam1

                ## test
                # labs = [ 'diff', 'scr0', 'scr1', 'preg0', 'preg1', 'pfac', 'treg', 'tfac', 'strdi', 'peri' ]
                # print pandas.DataFrame(tofinal[tofinal[:,8].argsort()], columns=labs)
                # print '-'*30
                # shear1 = sum(tofinal[:,3].astype(int))
                # print pdb, shear0, shear1, targetshears[0]

                ## get statistics of adjusted prediction
                correct_facing1 += sum(tofinal[:, 5] == tofinal[:, 7])
                correct_reg1 += sum(tofinal[:, 3] == tofinal[:, 6])
                shear1 = sum(tofinal[:, 3].astype(int))
                if shear1 != targetshears[0]:
                    badshear.append((pdb, shear1, targetshears[0]))

                finaldata = get_final(tofinal)

                if fout is not None:
                    for i in range(len(finaldata)):
                        fout.write(pdb)
                        for itm in finaldata[i]:
                            fout.write(' ' + str(itm))
                        fout.write('\n')
                else:
                    # print pdb, len(finaldata), sum(finaldata[:,1].astype(int))
                    # labs = [ 'strdid', 'predreg', 'predfac', 'truereg', 'truefac', 'peri' ]
                    # print pandas.DataFrame(finaldata, columns=labs)
                    # print '='*40
                    pass
                    ## test
                    # print '='*40

    print total_strand
    print correct_facing0, correct_facing1
    print correct_reg0, correct_reg1
    print np.array(badshear)

    fout.close()
