"""Convert the HDF5 file generated by PFLOTRAN to NetCDF format."""

# Author: Peishi Jiang

import numpy as np
import h5py
from netCDF4 import Dataset

# Let's define some utility functions

# Parameters
hdf_file = 'pflotran.h5'
nc_file  = 'pflotran_test.nc'

########################
# Read the HDF5 file
########################
f = h5py.File(hdf_file, "r")

# Get the grids/coordinates of the domain
coordinates = f["Coordinates"]
x_set = coordinates['X [m]'][:][:-1]
y_set = coordinates['Y [m]'][:][:-1]
z_set = coordinates['Z [m]'][:][:-1]

# Get the state values at each time step
time_set = list(f.keys())
time_set = [t for t in time_set if t.startswith("Time")]
state_set = dict.fromkeys(time_set)
for t in time_set:
    dataset = f[t]
    state_set[t] = {}
    for varname in list(dataset.keys()):
        state_set[t][varname] = dataset[varname][:]
        print(dataset[varname][:].shape)

# ntime    = len(time_set)
ntime    = 1
nx,ny,nz = len(x_set), len(y_set), len(z_set)

var_set  = list(state_set[t].keys())

state_set2 = dict.fromkeys(var_set)
for var in var_set:
    state_set2[var] = np.zeros([nx,ny,nz])
    state_set2[var][:,:,:] = state_set[time_set[ntime-1]][var]
    # print(state_set[time_set[ntime-1]][var])
    # state_set2[var] = np.zeros([ntime,nx,ny,nz])
    # for i in range(ntime):
        # t = time_set[i]
        # state_set2[var][i,:,:,:] = state_set[t][var]

time_set  = [t.split()[1:] for t in state_set.keys()]
time_vset = [float(t[0]) for t in time_set]
time_unit = time_set[0][1]

###############################
# Write it into NetCDF format
###############################
root_nc = Dataset(nc_file, 'w')

# Create the dimension
# location = root_nc.createDimension()
xloc     = root_nc.createDimension('x_location', nx)
yloc     = root_nc.createDimension('y_location', ny)
zloc     = root_nc.createDimension('z_location', nz)
time     = root_nc.createDimension('time', ntime)
# member   = root_nc.createDimension('member', 1)

# Write the values
times = root_nc.createVariable('time', 'f8', ('time',))
xloc = root_nc.createVariable('x_location', 'f8', ('x_location',))
yloc = root_nc.createVariable('y_location', 'f8', ('y_location',))
zloc = root_nc.createVariable('z_location', 'f8', ('z_location',))

# Create the attributes
times.units = time_unit
times.calendar = 'none'
times[:] = time_vset[:ntime]

# Write coordinates values
xloc[:], yloc[:], zloc[:] = x_set, y_set, z_set

var_dict_nc = dict.fromkeys(var_set)
for var in var_set:
    # Separate variable name from unit
    varinfo = var.split()
    if len(varinfo) == 1:
        varn, varunit = varinfo[0], ''
    elif len(varinfo) == 2:
        varn, varunit = varinfo[0], varinfo[1]
    else:
        print(varinfo)
        raise Exception('Wrong!')
    # Make the variable name uppercase
    varn = varn.upper()
    # Save the variable into netcdf
    # vargrp               = root_nc.createVariable(varn, 'f8', ('time', 'member', 'x_location','y_location','z_location'))
    # vargrp               = root_nc.createVariable(varn, 'f8', ('x_location','y_location','z_location'))
    # Save the variables in the opposite ordering so that the Fortran-NetCDF can read them in a correct way
    # See more explanations at: https://github.com/Unidata/netcdf4-python/issues/337
    vargrp               = root_nc.createVariable(varn, 'f8', ('z_location','y_location','x_location'))
    vargrp.unit          = varunit
    var_dict_nc[varn]    = vargrp
    # var_dict_nc[varn][:] = state_set2[var]
    var_dict_nc[varn][:] = np.asfortranarray(state_set2[var].T)
    print(var_dict_nc[varn][:][1,:,0])

root_nc.close()
