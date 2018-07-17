#!/usr/bin/python

import os
from new_reg_adjust import shear_adjust
#opts = ['0','1','2','3','4','5','6', '7', 'm1']
# opts = ['2']
#
# dirn = '../pred_combine/results/bmp51_len12/pred_110030080_p_0.00E+00/'
# for opt in opts:
# 	print '-'*20, opt, '-'*20
# 	os.system('./new_reg_adjust.py ../testfiles/bmp51_len12/bmp51_len12.test '+dirn+' '+opt+' w110030080_o'+opt)



#in the reference code only opt level 2 is used soooooooo...
def get_shear_adjustments( testfile_path, prediction_results_path, outfile, opt=2):
	shear_adjust(testfile_path, prediction_results_path, opt, outfile )
