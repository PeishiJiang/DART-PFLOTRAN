"""Clip the observation NetCDF file based on the required spatial and temporal domains."""

# Author: Peishi Jiang

import os
import sys
import f90nml
import warnings
import numpy as np
from netCDF4 import Dataset

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

# The original observation NetCDF file
obs_nc_original = configs["file_cfg"]["obs_nc_original_file"]        # The converted observation file
obs_nc          = configs["file_cfg"]["obs_nc_file"]        # The converted observation file

# The spatial domain
xmin, xmax = configs["obs_space_cfg"]["obs_space_xlimit"]
ymin, ymax = configs["obs_space_cfg"]["obs_space_ylimit"]
zmin, zmax = configs["obs_space_cfg"]["obs_space_zlimit"]

# The temporal domain
first_obs_time = configs["time_cfg"]["first_obs_time_size"]
last_obs_time  = configs["time_cfg"]["last_obs_time_size"]


def oned_clip(data, dmin, dmax):
    """A function for clipping 1D-array data based on dmin and dmax."""
    if (dmin == -99999) and (dmax == -99999):
        clipped_ind = np.ones(len(data), dtype=bool)
    elif (dmin == -99999) and (dmax != -99999):
        clipped_ind = data <= dmax
    elif (dmin != -99999) and (dmax == -99999):
        clipped_ind = data >= dmin
    else:
        clipped_ind = (data >= dmin) & (data <= dmax)

    return data[clipped_ind], clipped_ind


###############################
# Read the original data
###############################
if os.path.isfile(obs_nc):
    warnings.warn("The following file exists and thus is deleted: %s" % obs_nc)

root_nc_original = Dataset(obs_nc_original, 'r')
root_nc          = Dataset(obs_nc, 'w')


###############################
# Clip the temporal domain
###############################
time_original = root_nc_original.variables['time']
time_clipped,  time_clipped_ind = oned_clip(time_original[:], first_obs_time, last_obs_time)
ntime_clipped = np.sum(time_clipped_ind)


###############################
# Clip the spatial domain
###############################
xloc_original = root_nc_original.variables['x_location']
yloc_original = root_nc_original.variables['y_location']
zloc_original = root_nc_original.variables['z_location']

_,  xloc_clipped_ind = oned_clip(xloc_original[:], xmin, xmax)
_,  yloc_clipped_ind = oned_clip(yloc_original[:], ymin, ymax)
_,  zloc_clipped_ind = oned_clip(zloc_original[:], zmin, zmax)
loc_clipped_ind = xloc_clipped_ind & yloc_clipped_ind & zloc_clipped_ind
xloc_clipped = xloc_original[:][loc_clipped_ind]
yloc_clipped = yloc_original[:][loc_clipped_ind]
zloc_clipped = zloc_original[:][loc_clipped_ind]
nloc_clipped = np.sum(loc_clipped_ind)


###############################
# Create the spatial and temporal dimensions for the clipped NetCDF file
###############################
time = root_nc.createDimension('time', ntime_clipped)
loc  = root_nc.createDimension('location', nloc_clipped)

# Write the values
times = root_nc.createVariable('time', 'f8', ('time',))
xloc = root_nc.createVariable('x_location', 'f8', ('location',))
yloc = root_nc.createVariable('y_location', 'f8', ('location',))
zloc = root_nc.createVariable('z_location', 'f8', ('location',))
# locs  = root_nc.createVariable('location', 'f8', ('location',3))

# Create the attributes
zloc.units, zloc.type = zloc_original.units, zloc_original.type
yloc.units, yloc.type = yloc_original.units, yloc_original.type
xloc.units, xloc.type = xloc_original.units, xloc_original.type

times.calendar = time_original.calendar
times.units = time_original.units
times.type = time_original.type

# Write values
times[:] = time_clipped
xloc[:] = xloc_clipped
yloc[:] = yloc_clipped
zloc[:] = zloc_clipped


###############################
# Clip the variables based on the spatial and temporal domains
###############################
for varn in root_nc_original.variables.keys():
    var_original = root_nc_original.variables[varn]
    # Skip the dimension variables
    if var_original.type == "dimension_value":
        continue
    # If not, read the variable and clip it
    vargrp = root_nc.createVariable(
            varname=varn,
            datatype=var_original.datatype,
            dimensions=('location', 'time'),
            fill_value=var_original._FillValue
        )
    vargrp.unit = var_original.unit
    vargrp.type = var_original.type
    vargrp[:]   = var_original[loc_clipped_ind, time_clipped_ind]

root_nc_original.close()
root_nc.close()
