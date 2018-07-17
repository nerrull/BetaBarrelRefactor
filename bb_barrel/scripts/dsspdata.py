import re

class DSSPData:
	def __init__(self, fn):
		self.res_dict = {}
		self.records = {}
		self._parseDSSP(fn)

	def _parseDSSP(self, fn):
		with open(fn, 'r') as dsspf:
			start=False
			for line in dsspf:
				if not start and re.search('#', line):
					start=True
					continue
	
				if start :
					#if not line[5:10].strip().isdigit():
					try:
						int(line[5:10].strip())
					except:
						continue
					num     = int(line[0:5].strip())
					resnum  = int(line[5:10].strip())
					inscode = line[10:11].strip() 
					chain   = line[11:12].strip()
					aa      = line[12:14].strip()
					struct  = line[14:25]
					bp1     = int(line[25:29].strip())
					bp2     = int(line[29:33].strip())
					shlabel = line[33:34].strip()
					acc     = line[34:38].strip()
					h_nho1  = line[38:50].strip()
					h_ohn1  = line[50:61].strip()
					h_nho2  = line[61:72].strip()
					h_ohn2  = line[72:83].strip()
					tco     = line[83:91].strip()
					kappa   = line[91:97].strip()
					alpha   = line[97:103].strip() 
					phi     = line[103:109].strip() 
					psi     = line[109:115].strip() 
					xca     = float(line[115:122].strip())
					yca     = float(line[122:129].strip())
					zca     = float(line[129:136].strip())

					self.res_dict[ (chain, resnum) ] = num
					self.records[ num ] = {	'num'     : num,
											'resnum'  : resnum  ,
											'inscode' : inscode ,
											'chain'   : chain   ,
											'aa'      : aa      ,
											'struct'  : struct  ,
											'bp1'     : bp1     ,
											'bp2'     : bp2     ,
											'shlabel' : shlabel ,
											'acc'     : acc     ,
											'h_nho1'  : h_nho1  ,
											'h_ohn1'  : h_ohn1  ,
											'h_nho2'  : h_nho2  ,
											'h_ohn2'  : h_ohn2  ,
											'tco'     : tco     ,
											'kappa'   : kappa   ,
											'alpha'   : alpha   ,
											'phi'     : phi     ,
											'psi'     : psi     ,
											'xca'     : xca     ,
											'yca'     : yca     ,
											'zca'     : zca     }

	def get(self, chain, resnum, key):
		return self.records[ self.res_dict[(chain, resnum)] ][key]

	def get_chain_res(self, chain):
		resids = []
		for num in self.records.keys():
			if self.records[num]['chain'] == chain:
				resids.append(self.records[num]['resnum'])
		return sorted(resids)

