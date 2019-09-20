"""Prepare the namelist for convert_nc.f90 file"""

# Author: Peishi Jiang

import os
import sys

########################
# Parameters
########################
nc_file       = sys.argv[1]
dart_obs_file = sys.argv[2]
nml_file      = "../work/input.nml.convert_nc"

########################
# Create the namelist
########################
# Remove the current one
if os.path.exists(nml_file):
    os.remove(nml_file)

with open(nml_file, "w") as f:
    f.write("&convert_nc_nml\n")
    f.write("    netcdf_file = '"+nc_file+"'\n")
    f.write("    out_file  = '"+dart_obs_file+"'\n")
    f.write(" /\n")
