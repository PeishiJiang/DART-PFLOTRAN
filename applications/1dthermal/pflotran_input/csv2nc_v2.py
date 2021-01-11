"""Convert the temperature csv data to NetCDF file."""

# Author: Peishi Jiang

# Note that the original observation data format is case by case.
# Here, we consider csv file and NetCDF format.

import os
import sys
import f90nml
import warnings
import numpy as np
import pandas as pd
from netCDF4 import num2date, date2num, Dataset
from datetime import datetime, timedelta

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
obs_original    = 'temperature_1d_true_3.30-5.1_fmt_3col.csv'
obs_nc          = 'temperature_1d_true_3.30-5.1_fmt_3col.nc'
model_start_str = "2017-03-30 00:00:00"

# Get the reference time
# ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")
ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")

missing_value = -99999

###############################
# Read the csv file
###############################
obs_pd = pd.read_csv(obs_original)

# Get time and vertical depths
print(obs_pd.keys())
time_set  = obs_pd['# Datetime'].values
z_set     = [float(z[:-2])/100 for z in obs_pd.keys()[1:]]
ntime, nz = obs_pd.shape
nz        = nz-1
nloc      = nz*1*1
dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in time_set]

dates_ref = [t-ref_time for t in dates]
# dates_ref_values = [t.days+float(t.seconds)/86400.+spinup_time for t in dates_ref]
dates_ref_values = [t.days+float(t.seconds)/86400. for t in dates_ref]

# TODO: Get the true values based on the required spatial information
# TODO: Here, we assume the observation errors are uncorrelated.
# Get the true and perturbed temperature values
temperature = obs_pd[obs_pd.keys()[1:]].values
ntime,  nz  = temperature.shape

# print(temperature.shape)

###############################
# Write it into NetCDF format
###############################
if os.path.isfile(obs_nc):
    warnings.warn("The following file already exists: %s" % obs_nc)
    os.remove(obs_nc)

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

# print(times[:])

# Write coordinates values
xloc[:] = np.zeros(nloc)
yloc[:] = np.zeros(nloc)
zloc[:] = z_set

# Write the perturbed observations of variable information
# Missing value is assigned as -99999 as the fill_value
vargrp  = root_nc.createVariable(
    varname='TEMPERATURE',
    datatype='f8',
    # dimensions=('time','location'),
    dimensions=('location', 'time'),
    fill_value=missing_value)
vargrp.unit = 'C'
vargrp[:]   = temperature.T

root_nc.close()