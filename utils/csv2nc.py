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
nloc      = nz*1*1
dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in time_set]

# Get the temperature values
temperature = obs_pd[obs_pd.keys()[1:]].values
print(temperature.shape)

###############################
# Write it into NetCDF format
###############################
root_nc = Dataset(obs_nc, 'w')

# Create the dimension
# Here, the data points are saved in two dimensions:
#    -- time    -- location
time = root_nc.createDimension('time', ntime)
loc  = root_nc.createDimension('location', nloc)

# Write the values
times = root_nc.createVariable('time', 'f8', ('time',))
xloc = root_nc.createVariable('x_location', 'f8', ('location',))
yloc = root_nc.createVariable('y_location', 'f8', ('location',))
zloc = root_nc.createVariable('z_location', 'f8', ('location',))
# locs  = root_nc.createVariable('location', 'f8', ('location',3))

# TODO
# Add the units following the 'CF' conventions later on.
# For example,
# double time(time) ;
#        time:units = "days since 1986-01-01 00:00:00" ;
#        time:calendar = "gregorian" ;
# Create the attributes
zloc.units, zloc.type  = 'm', 'dimension_value'
yloc.units, yloc.type  = 'm', 'dimension_value'
xloc.units, xloc.type  = 'm', 'dimension_value'
times.units    = 'seconds since 2017-04-01 00:00:00'
# times.units    = 'None'
times.calendar = 'gregorian'
times.type     = 'dimension_value'
times[:] = date2num(dates,units=times.units,calendar=times.calendar)

# Write coordinates values
xloc[:] = np.zeros(nloc)
yloc[:] = np.zeros(nloc)
zloc[:] = z_set

# TODO add missing_value
# Write the variable information
vargrp      = root_nc.createVariable(
               varname='TEMPERATURE',
               datatype='f8',
               # dimensions=('time','location'),
               dimensions=('location','time'),
               fill_value=-99999)
vargrp.unit = 'C'
vargrp[:]   = temperature.T
vargrp.type = 'observation_value'

root_nc.close()

print("Finished converting raw observation in NetCDF format...")
