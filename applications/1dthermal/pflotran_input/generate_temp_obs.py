
"""Generate the observation data from the true data"""

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
true_error_type = "obsolute"
true_error = 0.05

missing_value = -99999

###############################
# Read the csv file
###############################
true_pd = pd.read_csv("temperature_true.csv")

# Get time and vertical depths
time_set  = true_pd['# Datetime'].values
z_set     = [float(z[:-2])/100 for z in true_pd.keys()[1:]]
ntime, nz = true_pd.shape
nz        = nz-1
nloc      = nz*1*1
dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in time_set]


###############################
# Get the true temperature values
###############################
temperature = true_pd[true_pd.keys()[1:]].values
ntime,  nz  = temperature.shape
if true_error_type == "obsolute":
    errors_var = true_error / 3. * np.ones([ntime, nz])
elif true_error_type == "relative":
    errors_var = true_error / 3. * temperature
else:
    raise Exception("Unknown observation error type %s" % true_error_type)

###############################
# Get the observed temperature values
###############################
obs_pd = true_pd.copy()
temperature_obs = temperature + np.random.normal(0,1,[ntime,nz])*errors_var
obs_pd[obs_pd.keys()[1:]] = temperature_obs

obs_pd.to_csv("temperature_obs_error005.csv", index=False)
