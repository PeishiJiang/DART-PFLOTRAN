"""Generate the prior ensemble in NetCDF format by reading (1) the model output in HDF; (2) PFLOTRAN parameter.h5; (3) the list of variable names to be updated/analyzed/assimilated."""

# Author: Peishi Jiang

import os
import re
import sys
import numpy as np
from netCDF4 import num2date, date2num, Dataset

###############################
# Parameters
###############################
pflotran_out_file  = sys.argv[1]
pflotran_para_file = sys.argv[2]
dart_prior_file    = sys.argv[3]
nens               = int(sys.argv[4])
pflotran_var_set   = sys.argv[5:]

# Parse the PFLOTRAN variables to be updated/analyzed/assimilated
p = re.compile('[A-Z_]+')
pflotran_var_set = [p.search(v).group() for v in pflotran_var_set]

ens_set = np.arange(1,nens+1)
# Get the file names of all ensembles for PFLOTRAN output
# pflotran_para_file_set = [re.sub(r"[\*]+",str(ens),pflotran_para_file) for ens in ens_set]
pflotran_out_file_set = [re.sub(r"\[ENS\]",str(ens),pflotran_out_file) for ens in ens_set]

#TODO check the existences of these files

# Get the file names of all ensembles for DART restart file
dart_prior_file_set = [re.sub(r"\[ENS\]",str(ens),dart_prior_file) for ens in ens_set]

###############################
# Read the model output in HDF
###############################


###############################
# Read the parameter.h5
###############################


###############################
# Write the prior to NetCDF file
###############################
