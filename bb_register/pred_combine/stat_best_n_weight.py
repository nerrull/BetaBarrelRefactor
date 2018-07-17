#!/usr/bin/python

import numpy as np
def best_of_file(fn):
	f = open(fn)
	lines = f.readlines()
	f.close()
	best = 0
	bestweights = []
	for line in lines:
		split = line.split()
		try:
			correct = int(split[6])
			weights = [float(w) for w in split[:6]]
		except:
			continue
		if correct > best:
			best = correct
			bestweights = [weights]
		elif correct == best:
			bestweights.append(weights)
	return best, np.array(bestweights)


def correct_of_ec_weight_dir(dirn):
	correct = 0
	dirweights = []
	for fn in sorted(glob.glob(dirn+'/*.l??')):
		filebest, filebestweights = best_of_file(fn)
		correct += filebest
		dirweights.append(filebestweights)
	return correct, dirweights

def correct_of_type(typen):
	stat = []
	for dirn in sorted(glob.glob(typen+'/*')):
		correct, dirweights = correct_of_ec_weight_dir(dirn)
		stat.append((correct, dirn[dirn.find('/')+1 : ], dirweights))
	stat = sorted(stat)
	return stat


import sys
sys.dont_write_bytecode = True
import glob

if len(sys.argv)!=2:
	print 'Usage: '+sys.argv[0]+' <result folder>'
	sys.exit(0)

ll = [01,02,03,04,05]

correct = correct_of_type(sys.argv[1])
count = 0
for count_n_file_n_ws in reversed(correct):
	print 'tot correct # :', count_n_file_n_ws[0], '||', count_n_file_n_ws[1][:count_n_file_n_ws[1].find('/')], '||', count_n_file_n_ws[1][count_n_file_n_ws[1].rfind('/')+6:count_n_file_n_ws[1].rfind('/')+15]
	for i in range(5):
		print '\t# lt '+str(ll[i])+' strds [ec, ub, neg, wstrong, wvdw, wweak] : '
		outstat=[]
		for ws in count_n_file_n_ws[2][i]:
			ws = np.array(ws)
			outstat.append(ws)
			#print '\t\t', ws
		outstat = np.array(outstat)
		#for i in range(outstat.shape[-1]):
		#	print np.min(outstat[:,i])-0.01, np.max(outstat[:,i])+0.01

	count+=1
	if count==400:
		break
	print '~'*40
	print '~'*40

