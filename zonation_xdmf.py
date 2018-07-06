# SUMMARY:      zonation_xdmf.py
# USAGE:        Generate XDMF file for eSTOMP zonation file
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    July-2018
# DESCRIPTION:  
# DESCRIP-END.
# COMMENTS:     only deal cartesian sturtured grids
#
# Last Change: 2018-07-02


import h5py as h5
import numpy as np
import re

simu_dir = '/pic/scratch/song884/bcomplex/model/'

input_file = open(simu_dir + "input", "r")

estomp_input = input_file.readlines()
input_file.close()

# remove comments and blank lines in input deck
estomp_input = [re.split('[# ! \n]', x)[0] for x in estomp_input]
estomp_input = [x for x in estomp_input if x]

grid_line = [i for i, s in enumerate(estomp_input) if "~Grid" in s][0]
