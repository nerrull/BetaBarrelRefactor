import sys
sys.dont_write_bytecode = True

# for a given hairpin, get the list of interstrand contacting residues in the given distance.
def get_dist_contacts(strand1, strand2, dist):
	contacts = set()
	for i in range(len(strand1)):
		try:
			contact = strand1[i], strand2[i+dist]
			if not '_' in contact:
				contacts.add(contact)
		except IndexError:
			pass
		try:
			if i-dist<0: # negative index has meaning in python, so deal with separately here
				continue
			contact = strand1[i], strand2[i-dist]
			if not '_' in contact:
				contacts.add(contact)
		except IndexError:
			pass
	return sorted(list(contacts))


# construct a hairpin according to the given strands and the registration.
# position of index 0 of the strands are the starting points of the strands on the periplasmic side
def construct_hairpin(strand1, strand2, reg):
	if reg<0:
		strd1 = strand1
		strd2 = ['_']*(-reg)+strand2
	else:
		strd1 = ['_']*reg+strand1
		strd2 = strand2
	return strd1, strd2


# read registrations from the input test file
def read_true_regs(testfn):
	f = open(testfn)
	lines = f.readlines()
	f.close()
	regdict = {}
	for line in lines:
		split = line.split()
		pdb = split[0]
		strandn = int(split[1])
		regdict[pdb] = []
		for i in range(4,strandn*3+2,3):
			regdict[pdb].append(int(split[i]))
	return regdict
