"""Convert the temperature csv data to NetCDF file."""

# Author: Peishi Jiang

# Note that the original observation data format is case by case.
# Here, we consider csv file with and NetCDF format.

import sys
import numpy as np
import pandas as pd
from netCDF4 import num2date, date2num, Dataset
from datetime import datetime, timedelta

###############################
# Parameters
###############################
obs_original = sys.argv[1]
obs_nc       = sys.argv[2]

###############################
# Read the csv file
###############################
obs_pd = pd.read_csv(obs_original)

# Get time and vertical depths
time_set  = obs_pd['# Datetime'].values
z_set     = [float(z[:-2])/100 for z in obs_pd.keys()[1:]]
ntime, nz = obs_pd.shape
nz        = nz-1
dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in time_set]

# Get the temperature values
temperature = obs_pd[obs_pd.keys()[1:]].values

###############################
# Write it into NetCDF format
###############################
root_nc = Dataset(obs_nc, 'w')

# Create the dimension
# TODO
# The 3D cartesian locations need to be revised later on
xloc     = root_nc.createDimension('x_location', 1)
yloc     = root_nc.createDimension('y_location', 1)
zloc     = root_nc.createDimension('z_location', nz)
time     = root_nc.createDimension('time', ntime)

# Write the values
times = root_nc.createVariable('time', 'f8', ('time',))
xloc = root_nc.createVariable('x_location', 'f8', ('x_location',))
yloc = root_nc.createVariable('y_location', 'f8', ('y_location',))
zloc = root_nc.createVariable('z_location', 'f8', ('z_location',))

# TODO
# Add the units following the 'CF' conventions later on.
# For example,
# double time(time) ;
#        time:units = "days since 1986-01-01 00:00:00" ;
#        time:calendar = "gregorian" ;
# Create the attributes
zloc.units  = 'm'
yloc.units  = 'm'
xloc.units  = 'm'
times.units = 'minutes since 2017-04-01 00:00:00'
times.calendar = 'gregorian'
times[:] = date2num(dates,units=times.units,calendar=times.calendar)

# Write coordinates values
xloc[:], yloc[:], zloc[:] = [0], [0], z_set

# Write the variable information
vargrp      = root_nc.createVariable('TEMPERATURE', 'f8', ('time','z_location','y_location','x_location'))
vargrp.unit = 'C'
# vargrp[:][:,:,0,0] = temperature
vargrp[:] = temperature

root_nc.close()
