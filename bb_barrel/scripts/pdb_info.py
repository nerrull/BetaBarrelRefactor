def read_pdb_info(fn):
	f = open(fn)
	lines = f.readlines()
	f.close()
	info = {}
	keys = lines[0].strip('#').split()
	for line in lines:
		if line[0]=='#':
			continue
		vals = line.split()
		pdb = vals[0]
		info[pdb] = {}
		for i in range(1,len(keys)):
			try:
				info[pdb][keys[i]]=vals[i]
			except IndexError:
				info[pdb][keys[i]]=''
	return info

