"""Plot the data assimilation results."""

# Author: Peishi Jiang

import os
import re
import sys
import f90nml
import warnings
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import gridspec


###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

# DART prior and posterior files
dart_posterior_file = configs["file_cfg"]["dart_posterior_nc_file"]
dart_prior_file     = configs["file_cfg"]["dart_prior_nc_file"]
dart_prior_template = configs["file_cfg"]["dart_prior_template_file"]

# Variables
obs_set             = configs["obspara_set_cfg"]["obs_set"]
para_set            = configs["obspara_set_cfg"]["para_set"]
if isinstance(obs_set, str):
    obs_set = [obs_set]
if isinstance(para_set, str):
    para_set = [para_set]
pflotran_var_set  = obs_set + para_set
nvar              = len(pflotran_var_set)

# Model time, assimilation start time, number of ensemble
assim_start_time    = configs["time_cfg"]["assim_start"]
model_time          = float(configs["time_cfg"]["current_model_time"])   # days
model_time_list     = configs["time_cfg"]["model_time_list"]
exceeds_obs_time    = configs["time_cfg"]["exceeds_obs_time"]
assim_window        = float(configs["da_cfg"]["assim_window_size"])
nens                = configs["da_cfg"]["nens"]
ntime               = len(model_time_list)
model_start_time    = model_time_list[0]
model_end_time      = model_time_list[-1] + assim_window

# Observation file in NetCDF
obs_nc = configs["file_cfg"]["obs_nc_file"]

if not exceeds_obs_time:
    warnings.warn("The data assimilation is not completed!")


###############################
# Get the model spatial domain
###############################
root_template = Dataset(dart_prior_template, 'r')
x_loc = root_template.variables['x_location'][:]
y_loc = root_template.variables['y_location'][:]
z_loc = root_template.variables['z_location'][:]
nx, ny, nz = len(x_loc), len(y_loc), len(z_loc)
root_template.close()

###############################
# Read in the prior and posterior data
###############################
prior     = np.zeros([nvar, nens, ntime, nz, ny, nx])
posterior = np.zeros([nvar, nens, ntime, nz, ny, nx])
for i in range(nens):
    ens = i + 1
    for j in range(ntime):
        model_time_ind = j + 1
        # Prior data
        prior_nc_file = re.sub(r"\[ENS\]", str(ens)+"_time"+str(model_time_ind),
                               dart_prior_file)
        root_prior = Dataset(prior_nc_file, 'r')
        for k in range(nvar):
            varn = pflotran_var_set[k]
            prior_var = root_prior.variables[varn][:]
            # print(prior_var.shape, nz, ny, nx)
            prior[k, i, j, :, :, :] = prior_var
        root_prior.close()
        # Posterior data
        posterior_nc_file = re.sub(r"\[ENS\]", str(ens)+"_time"+str(model_time_ind),
                                   dart_posterior_file)
        root_posterior = Dataset(posterior_nc_file, 'r')
        for k in range(nvar):
            varn = pflotran_var_set[k]
            posterior_var = root_posterior.variables[varn][:]
            posterior[k, i, j, :, :, :] = posterior_var
        root_posterior.close()


###############################
# Read in the observation data
###############################
obs_value_set      = dict.fromkeys(obs_set)
obs_value_set_used = dict.fromkeys(obs_set)
root_obs = Dataset(obs_nc, 'r')
obs_time_set = root_obs.variables['time'][:]
obs_used_ind = (obs_time_set >= model_start_time) & (obs_time_set <= model_end_time)
obs_time_set_used = obs_time_set[obs_used_ind]
for i in range(len(obs_set)):
    obs_var = obs_set[i]
    obs_value_set[obs_var]      = root_obs.variables[obs_var][:]
    obs_value_set_used[obs_var] = obs_value_set[obs_var][:, obs_used_ind]
root_obs.close()


###############################
# Plot -- mean throughout the whole domain
###############################
fig = plt.figure(num=1, dpi=150, figsize=(12, 6*nvar))
gs  = gridspec.GridSpec(nvar, 2, width_ratios=[1,1])

for i in range(nvar):
    varn = pflotran_var_set[i]
    # Plot the prior
    ax = plt.subplot(gs[i, 0])
    for j in range(nens):
        prior_ens_mean = np.mean(prior[i, j, :, :, :, :], axis=(1, 2, 3))
        ax.plot(model_time_list, prior_ens_mean, color='grey',
                linewidth=0.5, linestyle=':')
    prior_mean = np.mean(prior[i, :, :, :, :, :], axis=(0, 2, 3, 4))
    ax.plot(model_time_list, prior_mean, color='red', linewidth=2)
    if varn in obs_set:
        obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
        ax.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=2)
    ax.set_xlabel('Time (day)')
    ax.set_ylabel(varn)
    ax.set_title('Prior')

    # Plot the posterior
    ax = plt.subplot(gs[i, 1])
    for j in range(nens):
        posterior_ens_mean = np.mean(posterior[i, j, :, :, :, :], axis=(1, 2, 3))
        ax.plot(model_time_list, posterior_ens_mean, color='grey',
                linewidth=0.5, linestyle=':')
    posterior_mean = np.mean(posterior[i, :, :, :, :, :], axis=(0, 2, 3, 4))
    ax.plot(model_time_list, posterior_mean, color='red', linewidth=2)
    if varn in obs_set:
        obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
        ax.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=2)
    ax.set_xlabel('Time (day)')
    ax.set_ylabel(varn)
    ax.set_title('Posterior')

plt.savefig('test_da_plot2.png')
# plt.show()

###############################
# TODO
# Plot -- along with a given axis
###############################
