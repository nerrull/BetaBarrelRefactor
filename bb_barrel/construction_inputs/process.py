#!/usr/bin/python
import glob
import os
import sys

if len(sys.argv) != 2:
    print 'Usage: ' + sys.argv[0] + ' <pred_results folder>'
    sys.exit(0)

dirn = sys.argv[1]
regdirn = dirn + '/regs'
if not os.path.isdir(regdirn):
    os.makedirs(regdirn)
peridirn = dirn + '/peris'
if not os.path.isdir(peridirn):
    os.makedirs(peridirn)
facdirn = dirn + '/facings'
if not os.path.isdir(facdirn):
    os.makedirs(facdirn)

regdict = {}
peridict = {}
facdict = {}

with open(glob.glob(dirn + '/*.rlt')[0]) as f:
    lines = f.readlines()

for line in lines:
    pdb, strandid, peripredreg, peripredfacing, peritruereg, peritruefacing, peri = line.split()
    peripredreg = int(peripredreg)
    peritruereg = int(peritruereg)
    peri = int(peri)

    # peripredreg = peritruereg
    # peripredfacing = peritruefacing

    if pdb in regdict:
        regdict[pdb].append(peripredreg)
        peridict[pdb].append(peri)
        facdict[pdb].append(peripredfacing)
    else:
        regdict[pdb] = [peripredreg]
        peridict[pdb] = [peri]
        facdict[pdb] = [peripredfacing]

for pdb in regdict.keys():
    with open(dirn + '/facings/' + pdb + '.facings', 'w') as f:
        for peripredfacing in facdict[pdb]:
            f.write(peripredfacing + '\n')

    with open(dirn + '/regs/' + pdb + '.regs', 'w') as f:
        for peripredreg in regdict[pdb]:
            f.write(str(peripredreg) + '\n')

    with open(dirn + '/peris/' + pdb + '.peris', 'w') as f:
        for peri in peridict[pdb]:
            f.write(str(peri) + '\n')
