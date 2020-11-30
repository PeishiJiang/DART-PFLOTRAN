"""This file is used for updating the observation ensemble posterior from rerunning the model with parameter posterior as input."""

# Author: Peishi Jiang

import os
# Put the following before importing numpy is to avoid the following error
# OpenBLAS blas_thread_init: pthread_create failed for thread 19 of 64: Resource temporarily unavailable
# See more at: https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import re
import sys
import h5py
import shutil
import f90nml
import subprocess
import numpy as np
from math import ceil, log10
from netCDF4 import num2date, date2num, Dataset

from parse_pflotran_files_utils import pflotran_files

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

dart_data_dir            = configs["other_dir_cfg"]["dart_data_dir"]
update_pflotran_input_py = configs["file_cfg"]["update_pflotran_input_file"]
run_pflotran_py          = configs["file_cfg"]["run_pflotran_file"]

pflotran_out_file  = configs["file_cfg"]["pflotran_out_file"]
pflotran_para_file = configs["file_cfg"]["pflotran_para_file"]
dart_input_list    = configs["file_cfg"]["dart_input_list_file"]
dart_output_list   = configs["file_cfg"]["dart_output_list_file"]
# obs_set            = configs["obspara_set_cfg"]["obs_set"]
obs_pflotran_set   = configs["obspara_set_cfg"]["obs_pflotran_set"]
para_set           = configs["obspara_set_cfg"]["para_set"]
# para_take_log      = configs["obspara_set_cfg"]["para_take_log"]
model_time         = float(configs["time_cfg"]["current_model_time"])  # days
model_time_list    = configs["time_cfg"]["model_time_list"]
nens               = configs["da_cfg"]["nens"]
required           = configs["da_cfg"]["obs_ens_posterior_from_model"]
spinup_time        = configs["time_cfg"]["spinup_length"] # days

# Get the start and end time of the current assimilation window
assim_end_days    = configs["da_cfg"]["assim_end_days"]
assim_end_seconds = configs["da_cfg"]["assim_end_seconds"]
assim_start_days    = configs["da_cfg"]["assim_start_days"]
assim_start_seconds = configs["da_cfg"]["assim_start_seconds"]
start_obs_sec = assim_start_days * 86400 + assim_start_seconds
end_obs_sec   = assim_end_days * 86400 + assim_end_seconds

# Get the list of all required PFLOTRAN variables
if not isinstance(obs_pflotran_set, list):
    obs_pflotran_set = [obs_pflotran_set]
    
if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

# Get some constants
missing_value  = 99999

# Get the start and end time of the current assimilation window
assim_window_fixed = configs["da_cfg"]["assim_window_fixed"]
if assim_window_fixed:
    assim_window = float(configs["da_cfg"]["assim_window_size"])  # days
else:
    assim_window_list = configs["da_cfg"]["assim_window_list"]  # days
    if not isinstance(assim_window_list, list):
        assim_window_list = [assim_window_list]
    assim_window = assim_window_list[len(model_time_list)-1]
assim_end_days    = configs["da_cfg"]["assim_end_days"]
assim_end_seconds = configs["da_cfg"]["assim_end_seconds"]
assim_start_days    = configs["da_cfg"]["assim_start_days"]
assim_start_seconds = configs["da_cfg"]["assim_start_seconds"]
start_obs_sec = assim_start_days * 86400 + assim_start_seconds
end_obs_sec   = assim_end_days * 86400 + assim_end_seconds
assim_start, assim_end = start_obs_sec / 86400, end_obs_sec / 86400


###############################
# Check whether the observation ens posterior should be obtained from rerunning the model
###############################
if not required:
    ###############################
    # Remove the temporary prior files
    ###############################
    last_obs_time = configs['time_cfg']['last_obs_time_size']
    # DART prior
    if last_obs_time <= assim_end:
        dart_prior_file_set = []
        if not os.path.isfile(dart_input_list):
            raise Exception("The DART input list file does not exist in the path: %s" % dart_input_list)
        with open(dart_input_list, "r") as f:
            dart_prior_file_list = f.readlines()
        for i in range(len(dart_prior_file_list)):
            file_name_old = dart_prior_file_list[i]
            file_name     = file_name_old.replace("\n", "")
            subprocess.run("cd {0}; rm {1}".format(dart_data_dir, file_name), shell=True, check=True)

    exit()

print("Now update the observation ensemble posterior from rerunning the model.")


# raise Exception('Stop')
###############################
# Update the model parameters from DART posterior
###############################
configs["da_cfg"]["update_obs_ens_posterior_now"] = True
configs.write(config_nml, force=True)
subprocess.run("python {} {}".format(update_pflotran_input_py, config_nml), shell=True, check=True)


###############################
# Rerun the model
###############################
subprocess.run("python {} {}".format(run_pflotran_py, config_nml), shell=True, check=True)


###############################
# Check whether the DART posterior files exit.
###############################
# DART posterior
dart_posterior_file_set = []
if not os.path.isfile(dart_output_list):
    raise Exception(
        "The DART output list file does not exist in the path: %s" %
        dart_output_list)
with open(dart_output_list, "r") as f:
    dart_posterior_file_list = f.readlines()
for i in range(len(dart_posterior_file_list)):
    file_name_old = dart_posterior_file_list[i]
    file_name     = file_name_old.replace("\n", "")
    dart_posterior_file_set.append(file_name)


###############################
# Read the model output
###############################
pflotran_file_reader = pflotran_files(config_nml)
dart_state_dict, ntime_state, time_state, nloc_state, x_loc_state, y_loc_state, z_loc_state, cell_ids = \
    pflotran_file_reader.read_pflotran_output(assim_start, assim_end, missing_value)

###############################
# Assign the model output to observation ensemble posterior
###############################
for i in range(nens):
    dart_posterior_file = dart_posterior_file_set[i]
    # pf_fname            = pflotran_out_file_set[i]

    # dart_var_dict = dict.fromkeys(obs_set)

    # # Read the PFLOTRAN output
    # f_out = h5py.File(pf_fname, "r")

    # # Get the grids/coordinates of the domain
    # coordinates = f_out["Coordinates"]
    # # The following options are choosing the center of the grid
    # x_loc_state = coordinates['X [m]'][:]
    # y_loc_state = coordinates['Y [m]'][:]
    # z_loc_state = coordinates['Z [m]'][:]
    # x_loc_state = [(a + b) / 2 for a, b in zip(x_loc_state[:-1], x_loc_state[1:])]
    # y_loc_state = [(a + b) / 2 for a, b in zip(y_loc_state[:-1], y_loc_state[1:])]
    # z_loc_state = [(a + b) / 2 for a, b in zip(z_loc_state[:-1], z_loc_state[1:])]
    # nx_state, ny_state, nz_state = len(x_loc_state), len(y_loc_state), len(z_loc_state)

    # # Get all the time steps
    # time_set_o = np.array([t for t in list(f_out.keys()) if t.startswith("Time")])
    # time_set   = np.array([t.split()[1:] for t in time_set_o])
    # time_vset  = np.array([float(t[0]) for t in time_set])
    # time_unit  = time_set[0][1]

    # # Get the time steps within the assimilation window
    # time_set_assim_ind = (time_vset > start_obs_sec) & (time_vset <= end_obs_sec)
    # time_vset_assim    = time_vset[time_set_assim_ind]
    # time_set_o_assim   = time_set_o[time_set_assim_ind]
    # time_set_assim     = time_set[time_set_assim_ind]

    # ntime_state = len(time_vset_assim)
    # nloc_state = nx_state*ny_state*nz_state

    # # Initialize the dart_var_dict
    # for varn in obs_set:
    #     dart_var_dict[varn] = {"value": np.zeros([ntime_state, nloc_state]), "unit": ""}

    # # Get the state/parameter/variable values required in pflotran_var_set
    # for j in range(ntime_state):
    #     time_o           = time_set_o_assim[j]
    #     dataset          = f_out[time_o]
    #     pl_out_var_set_o = list(dataset.keys())
    #     # pl_out_var_dict  = dict()
    #     for v in pl_out_var_set_o:
    #         # Get the variable name and unit from the original variable name
    #         varinfo = v.split()
    #         if len(varinfo) == 1:
    #             varn, varunit = varinfo[0].upper(), ''
    #         elif len(varinfo) == 2:
    #             varn, varunit = varinfo[0].upper(), varinfo[1]
    #         else:
    #             raise Exception('Invalid variable name %s!' % v)
    #         # Check if the variable is required by pflotran_var_set
    #         if varn in obs_set:
    #             # dart_var_dict[varn]["value"].append(dataset[v][:]) 
    #             # TODO: make sure the shapes between dataset and dart_var_dict[varn]["value"] are compatible.
    #             dart_var_dict[varn]["value"][j,:] = dataset[v][:].flatten()
    #             dart_var_dict[varn]["unit"] = varunit 

    # f_out.close()


    # Open the posterior data
    nc_posterior = Dataset(dart_posterior_file, 'r+')

    # Get the ensemble index
    # This is because the order of the file names in dart_posterior_file_list
    # may not follow 1,2,3,...,ens
    ens_num = int(nc_posterior.variables["member"][:])
    ens_ind = ens_num - 1

    if ens_ind != i:
        raise Exception("The ensemble member does not match! Check your DART posterior file list.")

    for j in range(len(obs_pflotran_set)):  # For each observation variable ...
        varn = obs_pflotran_set[j]
        nc_posterior.variables[varn][:] = dart_state_dict[varn]["value"][i,:,:]

    # Close the posterior NetCDF file
    nc_posterior.close()

###############################
# Remove the temporary prior files
###############################
last_obs_time = configs['time_cfg']['last_obs_time_size']
# DART prior
if last_obs_time <= assim_end:
    dart_prior_file_set = []
    if not os.path.isfile(dart_input_list):
        raise Exception("The DART input list file does not exist in the path: %s" % dart_input_list)
    with open(dart_input_list, "r") as f:
        dart_prior_file_list = f.readlines()
    for i in range(len(dart_prior_file_list)):
        file_name_old = dart_prior_file_list[i]
        file_name     = file_name_old.replace("\n", "")
        subprocess.run("cd {0}; rm {1}".format(dart_data_dir, file_name), shell=True, check=True)


# Once it is done, switch the update back to False
configs["da_cfg"]["update_obs_ens_posterior_now"] = False
configs.write(config_nml, force=True)