import sys

sys.dont_write_bytecode = True
import Bio.PDB
import math
import numpy as np
import copy
from dsspdata import DSSPData
import AminoAcid as AA
import exceptions
from os.path import join
from defs import INPUT_DIR, CONSTRUCTION_DIR
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


class Ball(list):
    p_ballid, p_resseqid, p_coordx, p_coordy, p_coordz, p_aaid, p_facing = range(7)
    Facing = enum('IN', 'OUT', 'SP')

    def getballid(self):
        return self[Ball.p_ballid]

    def getcoord(self):
        return np.array([self[Ball.p_coordx], self[Ball.p_coordy], self[Ball.p_coordz]])

    def setcoord(self, coord):
        self[Ball.p_coordx] = coord[0]
        self[Ball.p_coordy] = coord[1]
        self[Ball.p_coordz] = coord[2]

    def getseqid(self):
        return self[Ball.p_resseqid]

    def gettype(self):
        return self[Ball.p_aaid]

    def getfacing(self):
        return self[Ball.p_facing]

    def __str__(self):
        return str(self[Ball.p_ballid]) + " " + str(self[Ball.p_resseqid]) + " " + \
               str(self.getcoord()) + " " + \
               AA.index_to_one(self[Ball.p_aaid]) + " " + str(self[Ball.p_facing])


class Barrel:
    extra_ball_num = 4

    def __init__(self, pdb, A=3.30, B=4.72, dr=0.51, dw=0.2, inputdirn=None, extraballs=False,
                 np1_right=False):  # zigzag
        self.pdb = pdb

        # params to build barrel
        self.A = A  # intrastrand Ca distance
        self.B = B  # interstrand Ca distance
        self.dr = dr  # anterior-posterior zigzag radius
        self.dw = dw  # left-right zigzag radius

        # A0=3.79 # avg intra length of dataset. A0 is the actual intrastrand length, A is just a part of A0
        # self.A = 2*math.sqrt(A0*A0/4-dw*dw-dr*dr)

        self.inputdirn = inputdirn

        self.strandlens = []
        self.reindexmap = []

        self.balls = []
        self.extraballs = extraballs
        self.np1_right = np1_right

        self.radius = 0  # radius of the barrel

        self._construct_ideal_balls_()

    def _construct_ideal_balls_(self):

        # load residues
        strands_file = join(INPUT_DIR,*[self.pdb, self.pdb+'.strands'])
        residues_file = join(INPUT_DIR,*[self.pdb, self.pdb+'.res'])
        regs_file = join(CONSTRUCTION_DIR,*['results','regs',self.pdb+'.regs'])
        peris_file = join(CONSTRUCTION_DIR,*['results','peris',self.pdb+'.peris'])
        facings_file = join(CONSTRUCTION_DIR,*['results','facings',self.pdb+'.facings'])

        residues = np.loadtxt(residues_file).astype(int).tolist()
        residues = [200] + residues

        # load regs and peris got from reg adjustment
        periregs = np.loadtxt(regs_file).astype(int)
        peris = np.loadtxt(peris_file).astype(int)
        # load facings  got from reg adjustment
        with open(facings_file) as f:
            lines = f.readlines()
        firstfacings = [line.strip() for line in lines]

        # load strandends
        strandends = np.loadtxt(strands_file).astype(int)

        # the strand ranges here may not be the ranges used for reg prediction
        # correct reg according to the difference between strandends to construct barrel and peris used in reg pred
        # after this loop, array periregs will be the predicted regs of strandends
        for strdi in range(len(strandends)):
            strdim1 = (strdi - 1) % len(strandends)
            if strdi % 2 == 0:
                periregs[strdi] -= strandends[strdi][0] - peris[strdi]
                periregs[strdim1] += strandends[strdi][0] - peris[strdi]
            else:
                periregs[strdi] += strandends[strdi][1] - peris[strdi]
                periregs[strdim1] -= strandends[strdi][1] - peris[strdi]

        # correct facing according to the difference between strandends to construct barrel and peris used in reg pred
        for strdi in range(len(strandends)):
            if strdi % 2 == 0:
                if (strandends[strdi][0] - peris[strdi]) % 2 != 0:
                    if firstfacings[strdi] == 'OUT':
                        firstfacings[strdi] = 'IN'
                    else:
                        firstfacings[strdi] = 'OUT'
            else:
                if (strandends[strdi][1] - peris[strdi]) % 2 != 0:
                    if firstfacings[strdi] == 'OUT':
                        firstfacings[strdi] = 'IN'
                    else:
                        firstfacings[strdi] = 'OUT'

        # add extra residues for bbq
        for strdi in range(len(strandends)):
            strandends[strdi][0] -= Barrel.extra_ball_num
            strandends[strdi][1] += Barrel.extra_ball_num

        # construct facing arrays for all residues (including extra residues)
        facings = []
        for fac in firstfacings:
            if fac == 'OUT':
                facings.append([Ball.Facing.OUT])
            else:
                facings.append([Ball.Facing.IN])
        for strdi in range(len(strandends)):
            for resi in range(strandends[strdi][1] - strandends[strdi][0]):
                if facings[strdi][resi] == Ball.Facing.OUT:
                    facings[strdi].append(Ball.Facing.IN)
                else:
                    facings[strdi].append(Ball.Facing.OUT)

        peripositions = np.cumsum(np.hstack(([0], -periregs)))

        N = len(strandends)  # strand num
        A = self.A  # intrastrand Ca distance
        B = self.B  # interstrand Ca distance
        S = sum(periregs)  # shear number

        ## circle formula
        # a = math.sqrt( (N*B)**2+(S*A)**2 ) / 2.0 / math.pi # tilt angle
        # theta = math.asin(S*A/2.0/math.pi/a) # radius

        ## polygan formula
        theta = math.atan(S * A / (N * B))  # tilt angle
        a = B / (2 * math.sin(math.pi / N) * math.cos(theta))  # radius

        self.radius = a
        b = a / math.tan(theta)  # vertical speed
        c = math.sqrt(a * a + b * b)
        delta = 2 * math.pi * a * a / (
        c * N)  # offset on the neigbouring strand to ensure inter H-bond is perpendicular to the strand

        currid = 0
        ids = []
        seqids = []
        restypes = []
        cacoords = []

        # construct the barrel
        for strdi in range(N):
            # seq ids
            if strdi % 2 == 0:
                seqids.append(range(strandends[strdi][0], strandends[strdi][1] + 1))
            else:
                seqids.append(range(strandends[strdi][0], strandends[strdi][1] + 1)[::-1])
            ids.append(range(currid, currid + len(seqids[strdi])))
            currid += len(seqids[strdi])
            # residue types
            restypes.append([])
            for seqid in seqids[strdi]:
                try:
                    restypes[strdi].append(AA.index_to_one(residues[seqid]))
                except:
                    restypes[strdi].append('C')

            cacoords.append([])
            for resi in range(strandends[strdi][1] - strandends[strdi][0] + 1):
                # zigzag deviation
                if facings[strdi][resi] == Ball.Facing.OUT:
                    dr = self.dr
                    # righthand side of out facing residue is always SH
                    # lefthand side NH
                    if strdi % 2 == 0:
                        dw = self.dw
                    else:
                        dw = -self.dw
                    if self.np1_right:  # test TODO
                        if strdi % 2 == 0:
                            dw = -self.dw
                        else:
                            dw = self.dw
                else:
                    dr = -self.dr
                    if strdi % 2 == 0:
                        dw = -self.dw
                    else:
                        dw = self.dw
                    if self.np1_right:  # test TODO
                        if strdi % 2 == 0:
                            dw = self.dw
                        else:
                            dw = -self.dw

                s = (peripositions[strdi] + resi) * A + strdi * delta
                x = (a + dr) * math.cos(s / c - 2 * math.pi * strdi / N);
                y = (a + dr) * math.sin(s / c - 2 * math.pi * strdi / N);
                if self.np1_right:  # test TODO
                    s = (peripositions[strdi] + resi) * A + (N - strdi) * delta
                    x = (a + dr) * math.sin(s / c + 2 * math.pi * strdi / N);
                    y = (a + dr) * math.cos(s / c + 2 * math.pi * strdi / N);

                z = b * s / c;

                xn1 = (a + dr) * (- math.cos(s / c - 2 * math.pi * strdi / N) + math.cos(
                    (s + delta) / c - 2 * math.pi * (strdi + 1) / N));
                yn1 = (a + dr) * (- math.sin(s / c - 2 * math.pi * strdi / N) + math.sin(
                    (s + delta) / c - 2 * math.pi * (strdi + 1) / N));
                zn1 = b * delta / c
                if (strdi % 2 == 1 and facings[strdi][resi] == Ball.Facing.OUT) or (
                            strdi % 2 == 0 and facings[strdi][resi] == Ball.Facing.IN):
                    xn1 = (a + dr) * (- math.cos(s / c - 2 * math.pi * (strdi - 1) / N) + math.cos(
                        (s + delta) / c - 2 * math.pi * strdi / N));
                    yn1 = (a + dr) * (- math.sin(s / c - 2 * math.pi * (strdi - 1) / N) + math.sin(
                        (s + delta) / c - 2 * math.pi * strdi / N));
                    zn1 = b * delta / c

                if self.np1_right:  # test TODO
                    xn1 = (a + dr) * (- math.sin(s / c + 2 * math.pi * strdi / N) + math.sin(
                        (s + delta) / c + 2 * math.pi * (strdi + 1) / N));
                    yn1 = (a + dr) * (- math.cos(s / c + 2 * math.pi * strdi / N) + math.cos(
                        (s + delta) / c + 2 * math.pi * (strdi + 1) / N));
                    if (strdi % 2 == 1 and facings[strdi][resi] == Ball.Facing.OUT) or (
                                strdi % 2 == 0 and facings[strdi][resi] == Ball.Facing.IN):
                        xn1 = (a + dr) * (- math.sin(s / c + 2 * math.pi * (strdi - 1) / N) + math.sin(
                            (s + delta) / c + 2 * math.pi * strdi / N));
                        yn1 = (a + dr) * (- math.cos(s / c + 2 * math.pi * (strdi - 1) / N) + math.cos(
                            (s + delta) / c + 2 * math.pi * strdi / N));

                n1norm = math.sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1)
                xn1 = xn1 / n1norm
                yn1 = yn1 / n1norm
                zn1 = zn1 / n1norm
                x += dw * xn1
                y += dw * yn1
                z += dw * zn1

                cacoords[strdi].append(np.array([x, y, z]))

            self.strandlens.append(len(ids[strdi]))

        for i in range(len(ids)):
            for j in range(len(ids[i])):
                ## following line is for model/param selections
                # ball = Ball([ ids[i][j], seqids[i][j], cacoords[i][j][0], cacoords[i][j][1], cacoords[i][j][2], AA.one_to_index(restypes[i][j]), facings[i][j] ])
                ## store ballids instead of seqids. needs to be correted after bbq
                ## if using seqids, bbq will have problems
                if i % 2 != 0:
                    ball = Ball([ids[i][j], ids[i][len(ids[i]) - j - 1], cacoords[i][j][0], cacoords[i][j][1],
                                 cacoords[i][j][2], AA.one_to_index(restypes[i][j]), facings[i][j]])
                    if j >= Barrel.extra_ball_num and j < len(ids[i]) - Barrel.extra_ball_num:
                        self.reindexmap.append((ids[i][len(ids[i]) - j - 1], seqids[i][j]))
                else:
                    ball = Ball([ids[i][j], ids[i][j], cacoords[i][j][0], cacoords[i][j][1], cacoords[i][j][2],
                                 AA.one_to_index(restypes[i][j]), facings[i][j]])
                    if j >= Barrel.extra_ball_num and j < len(ids[i]) - Barrel.extra_ball_num:
                        self.reindexmap.append((ids[i][j], seqids[i][j]))
                self.balls.append(ball)


def write_pdb(balls, fn, strandlens, extraballs=False, reindexmap=None, mapfn=None):
    tmpballs = []
    pos = 0
    for i in range(len(strandlens)):
        strandlen = strandlens[i]
        tmpstrandballs = []
        for j in range(strandlen):
            tmpstrandballs.append(balls[pos])
            pos += 1
        # remove the first and the last balls from the strands, which are extra balls
        if not extraballs:
            tmpstrandballs.pop(0)
            tmpstrandballs.pop(0)
            tmpstrandballs.pop(0)
            tmpstrandballs.pop(0)
            tmpstrandballs.pop(-1)
            tmpstrandballs.pop(-1)
            tmpstrandballs.pop(-1)
            tmpstrandballs.pop(-1)

        if i % 2 == 1:
            tmpstrandballs = tmpstrandballs[::-1]
        tmpballs += tmpstrandballs

    chain = Bio.PDB.Chain.Chain('A')
    for i in range(len(tmpballs)):
        try:
            res_id = (' ', tmpballs[i][Ball.p_resseqid], ' ')
            restype = AA.index_to_three(tmpballs[i][Ball.p_aaid])
            residue = Bio.PDB.Residue.Residue(res_id, restype, ' ')
            cacoord = tmpballs[i].getcoord()
            atom = Bio.PDB.Atom.Atom('CA', cacoord, 0, 0, ' ', 'CA', tmpballs[i][Ball.p_resseqid], 'C')
            residue.add(atom)
            chain.add(residue)
        except:
            res_id = ('A', tmpballs[i][Ball.p_resseqid], ' ')
            restype = AA.index_to_three(tmpballs[i][Ball.p_aaid])
            residue = Bio.PDB.Residue.Residue(res_id, restype, ' ')
            cacoord = tmpballs[i].getcoord()
            atom = Bio.PDB.Atom.Atom('CA', cacoord, 0, 0, ' ', 'CA', tmpballs[i][Ball.p_resseqid], 'C')
            residue.add(atom)
            chain.add(residue)
    model = Bio.PDB.Model.Model(1)
    model.add(chain)
    structure = Bio.PDB.Structure.Structure("ref")
    structure.add(model)
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(fn, write_end=False)

    if reindexmap is not None and mapfn is not None:
        np.savetxt(mapfn, reindexmap, fmt='%d')
