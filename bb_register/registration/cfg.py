def read_cfg(fn):
	f = open(fn)
	lines = f.readlines()
	f.close()
	cfgdict = {}
	for line in lines:
		strip = line.strip()
		if strip[0]=='#':
			continue
		try:
			key, val = line.split(':',1)
			key, val = key.strip(), val.strip()
			cfgdict[key]=val
		except ValueError: # no ':' in line
			pass
	return cfgdict
