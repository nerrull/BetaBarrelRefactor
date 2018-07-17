#!/usr/bin/python

#######################################################
#
# The script reads a set of PDB files and generates
# a database for backbone reconstruction
#
######################################################

import commands
import sys
import re
import os
import glob

#---------------------- SET UP SOME OPTIONS -----------------------------
# Where PDB files are?
PDB_DIR='../prot_set'
# File for quadrilaterals 
Q_OUT='q_15_xyz.q'
# File for database
Q_DAT='q_15_xyz.dat'
# JAVA command (include CLASSPATH here if necessary)
JAVA='java'
# Maximum R-factor allowed
MAX_RF='15'

#---------------------- THE SCRIPT ITSELF  -----------------------------

f2=open(Q_OUT,"wt")

s = glob.glob(PDB_DIR+'/*.pdb')
for i in s:
  cmd=JAVA+" BBQ -q="+i+" -max_r_factor="+MAX_RF+" 2>> ERROR_LOG"
  print cmd
  status,output = commands.getstatusoutput(cmd);
  f2.write(output+"\n")



cnt = {}

cx = {}
cy = {}
cz = {}

ox = {}
oy = {}
oz = {}

nx = {}
ny = {}
nz = {}

f=open(Q_OUT,"rt")
for line in f:
  if len(line) <92 : 
    continue
  key = line[0:6].strip()
  if cnt.has_key(key) :
    cnt[key] = cnt[key] + 1
    cx[key] = cx[key] + float(line[26:32].strip())
    cy[key] = cy[key] + float(line[33:39].strip())
    cz[key] = cz[key] + float(line[40:46].strip())

    ox[key] = ox[key] + float(line[49:55].strip())
    oy[key] = oy[key] + float(line[56:62].strip())
    oz[key] = oz[key] + float(line[63:69].strip())

    nx[key] = nx[key] + float(line[72:78].strip())
    ny[key] = ny[key] + float(line[79:85].strip())
    nz[key] = nz[key] + float(line[86:92].strip())
  else :
    cnt[key] = 1
    cx[key] = float(line[26:32].strip())
    cy[key] = float(line[33:39].strip())
    cz[key] = float(line[40:46].strip())

    ox[key] = float(line[49:55].strip())
    oy[key] = float(line[56:62].strip())
    oz[key] = float(line[63:69].strip())

    nx[key] = float(line[72:78].strip())
    ny[key] = float(line[79:85].strip())
    nz[key] = float(line[86:92].strip())

l = cnt.keys()
f4=open(Q_DAT,"wt")
for i in l :
  print >> f4, "%(key)6d 0.00 0.00   0.00 C %(cx)6.3f %(cy)6.3f %(cz)6.3f O %(ox)6.3f %(oy)6.3f %(oz)6.3f N %(nx)6.3f %(ny)6.3f %(nz)6.3f %(cnt)4d" % \
	{'key': int(i),'cnt':cnt[i],'cx':cx[i]/cnt[i],'cy':cy[i]/cnt[i],'cz':cz[i]/cnt[i]
	,'ox':ox[i]/cnt[i],'oy':oy[i]/cnt[i],'oz':oz[i]/cnt[i]
	,'nx':nx[i]/cnt[i],'ny':ny[i]/cnt[i],'nz':nz[i]/cnt[i]}
#  print i,cnt[i],cx[i]/cnt[i],cy[i]/cnt[i],cz[i]/cnt[i]
f2.close();
f4.close();
f.close();
