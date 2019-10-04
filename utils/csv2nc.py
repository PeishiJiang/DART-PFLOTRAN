"""Convert the temperature csv data to NetCDF file."""

# Author: Peishi Jiang

# Note that the original observation data format is case by case.
# Here, we consider csv file with and NetCDF format.

import sys
import f90nml
import numpy as np
import pandas as pd
from netCDF4 import num2date, date2num, Dataset
from datetime import datetime, timedelta

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

obs_original = configs["file_cfg"]["obs_original_file"]        # The original observation file
obs_nc       = configs["file_cfg"]["obs_nc_file"]        # The converted observation file
# spinup_time  = float(configs["time_cfg"]["spinup_length"]) # The spinup time (day)
assim_start_str = configs["time_cfg"]["assim_start"]     # The map between the start of observation and spinup time

# Get the reference time
# ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")
ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")

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

dates_ref = [t-ref_time for t in dates]
# dates_ref_values = [t.days+float(t.seconds)/86400.+spinup_time for t in dates_ref]
dates_ref_values = [t.days+float(t.seconds)/86400. for t in dates_ref]

# Get the temperature values
temperature = obs_pd[obs_pd.keys()[1:]].values
# print(temperature.shape)

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

# Create the attributes
zloc.units, zloc.type  = 'm', 'dimension_value'
yloc.units, yloc.type  = 'm', 'dimension_value'
xloc.units, xloc.type  = 'm', 'dimension_value'
# times.units    = 'seconds since 2017-04-01 00:00:00'
# times.units    = 'days since 2017-04-01 00:00:00'
# times.calendar = 'gregorian'
# times.type     = 'dimension_value'
# times[:] = date2num(dates,units=times.units,calendar=times.calendar)

times.calendar = 'None'
times.units    = 'days'
times.type     = 'dimension_value'
times[:]       = dates_ref_values

# Write coordinates values
xloc[:] = np.zeros(nloc)
yloc[:] = np.zeros(nloc)
zloc[:] = z_set

# Write the variable information
# Missing value is assigned as -99999 as the fill_value
vargrp      = root_nc.createVariable(
               varname='TEMPERATURE',
               datatype='f8',
               # dimensions=('time','location'),
               dimensions=('location', 'time'),
               fill_value=-99999)
vargrp.unit = 'C'
vargrp[:]   = temperature.T
vargrp.type = 'observation_value'

root_nc.close()

print("Finished converting raw observation in NetCDF format...")
