"""Update the PFLOTRAN input files based on the DART posterior output data"""

# Author: Peishi Jiang

import os
import re
import sys
import h5py
import shutil
import f90nml
import numpy as np
from netCDF4 import num2date, date2num, Dataset

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)


###############################
# TODO
# Update time in PFLOTRAN.in
###############################


###############################
# TODO
# Update the parameter values in
# parameter_prior.h5 based on DART
# posterior output
###############################
