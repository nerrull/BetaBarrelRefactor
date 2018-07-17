one_list = [
	'A', 'R', 'N', 'D', 'C',
	'Q', 'E', 'G', 'H', 'I',
	'L', 'K', 'M', 'F', 'P',
	'S', 'T', 'W', 'Y', 'V'
	]

three_list = [ 
	'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
	'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
	'LEU', 'LYS', 'MET', 'PHE', 'PRO',
	'SER', 'THR', 'TRP', 'TYR', 'VAL'
	]

index_list = range(20)

# strandard
one_to_three_dict = dict( zip(one_list, three_list) )

# strandard
three_to_one_dict = dict( zip(three_list, one_list) )
# add non strandard
three_to_one_dict[ 'MSE' ] = 'M'

one_to_index_dict = dict( zip(one_list, index_list) )
index_to_one_dict = dict( zip(index_list, one_list) )

def three_to_one(key):
	return three_to_one_dict[key]
def three_to_index(key):
	return one_to_index_dict[ three_to_one_dict[key] ]
def one_to_three(key):
	return one_to_three_dict[key]
def one_to_index(key):
	return one_to_index_dict[key]
def index_to_one(key):
	return index_to_one_dict[key]
def index_to_three(key):
	return one_to_three_dict[ index_to_one_dict[key] ]
