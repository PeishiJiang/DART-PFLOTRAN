"""Convert the NetCDF file to DART observation file."""

# Author: Peishi Jiang

import re
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

###############################
# Parameters
###############################
obs_type_qty_ind = sys.arg(1)
obs_nc   = sys.arg(2)
obs_dart = sys.arg(3)

###############################
# Read the NetCDF file
###############################
root_nc = Dataset(obs_nc, 'r')

# Get dimensions
dimensions = root_nc.dimensions.keys()

# Get observation variable name list
all_var = root_nc.variables.keys()
obs_var_list = list(set(all_var)-set(dimensions))
nvar    = len(obs_var)

# Get dimension values
time = root_nc.variables['time']
x_loc = root_nc.variables['x_location']
y_loc = root_nc.variables['y_location']
z_loc = root_nc.variables['z_location']
nt, nx, ny, nz = len(time), len(x_loc), len(y_loc), len(z_loc)

# Get variable values
for i in range(len(nvar)):
    var_val = root_nc.variables[obs_var[i]]

# Compute the total number of observations
num_obs = nt*nx*ny*nz*nvar

root_nc.close()


###############################
# Read the map with DART quantities
###############################
p = re.compile('[A-Z_]+')
obs_dart_map = {}
with open(obs_type_qty_ind, 'r') as f:
    # variable names
    line = f.readline()[:-2]
    obs_var_all = line.split(',')
    # DART variable quantity names
    line = f.readline()
    obs_qty_all = line.split(',')
    # DART variable quantity indices
    line = f.readline()
    obs_qty_ind_all = line.split(',')
    # Save them in a dictionary
    obs_qty_info_all = list(zip(obs_qty_all, obs_qty_ind_all))
    obs_dart_map = dict(zip(obs_var_all, obs_qty_info_all))


###############################
# Write the DART observation file
###############################
# TODO
with open(obs_dart, 'w') as f:
    # Write the head
    f.write('obs_sequence\n')
    f.write('obs_kind_definitions\n')
    f.write(str(nvar)+'\n')
    for i in range(nvar):
        var_name = obs_var_list[i]
        qty_ind  = obs_dart_map[var_name][1]
        f.write(str(qty_ind)+" "+var_name+"\n")
    f.write("num_copies:   1  num_qc:   1\n")
    f.write("num_obs:   "+str(num_obs)+"  max_num_obs:   "+str(num_obs)+"\n")


    # Write each observation
    pass
