"""Update the PFLOTRAN input files based on the DART posterior output data"""

# Author: Peishi Jiang

import os
import sys
import h5py
import f90nml
import numpy as np
from scipy.stats import truncnorm
from netCDF4 import num2date, date2num, Dataset

# TODO: For now, a single value for each variable is assumed, it should be more generic in the future by considering the 3D case.

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_in_file         = configs["file_cfg"]["pflotran_in_file"]
pflotran_para_file       = configs["file_cfg"]["pflotran_para_file"]
pflotran_checkpoint_file = configs["file_cfg"]["pflotran_checkpoint_file"]
dart_posterior_file      = configs["file_cfg"]["dart_posterior_nc_file"]
dart_output_list         = configs["file_cfg"]["dart_output_list_file"]
spinup_length            = configs["time_cfg"]["spinup_length"]
model_time_list          = configs["time_cfg"]["model_time_list"]
current_model_time       = configs["time_cfg"]["current_model_time"]
enks_mda_iteration_step  = configs["da_cfg"]["enks_mda_iteration_step"]
assim_window_fixed       = configs["da_cfg"]["assim_window_fixed"]

try:
    update_obs_ens_posterior = configs["da_cfg"]["update_obs_ens_posterior"]
except:
    update_obs_ens_posterior = False

if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

# Get the current model end time
if assim_window_fixed:
    assim_window_days    = configs["da_cfg"]["assim_window_days"]
    assim_window_seconds = configs["da_cfg"]["assim_window_seconds"]
    assim_window_size    = configs["da_cfg"]["assim_window_size"]
    one_sec                    = 1./86400.  # one second (fractional days)
    current_model_time_end     = current_model_time + (assim_window_size - one_sec) / 2
else:
    assim_window_list = configs["da_cfg"]["assim_window_list"]
    assim_time_sofar = np.sum(assim_window_list[:len(model_time_list)])
    current_model_time_end = spinup_length + assim_time_sofar

current_model_time_end_sec = current_model_time_end * 86400
spinup_length_sec          = spinup_length * 86400

para_set      = configs["obspara_set_cfg"]["para_set"]
para_take_log = configs["obspara_set_cfg"]["para_take_log"]
para_min_set  = configs["obspara_set_cfg"]["para_min_set"]
para_max_set  = configs["obspara_set_cfg"]["para_max_set"]
para_mean_set = configs["obspara_set_cfg"]["para_mean_set"]
para_std_set  = configs["obspara_set_cfg"]["para_std_set"]
para_dist_set = configs["obspara_set_cfg"]["para_dist_set"]
rescaled      = configs["obspara_set_cfg"]["para_prior_rescaled"]
para_resampled_set = configs["obspara_set_cfg"]["para_resampled_set"]
nens               = configs["da_cfg"]["nens"]

if not isinstance(para_resampled_set, list):
    para_resampled_set = [para_resampled_set]

if not isinstance(para_set, list):
    para_set      = [para_set]
    para_min_set  = [para_min_set]
    para_max_set  = [para_max_set]
    para_mean_set = [para_mean_set]
    para_std_set  = [para_std_set]
    para_dist_set = [para_dist_set]

# Check whether this is the first time to update pflotran input files
first_time_update = True if len(model_time_list) == 1 else False

# Check whether this is the first time to update pflotran input files
second_time_update = True if len(model_time_list) == 2 else False

# Check whether the spinup was conducted before
is_spinup_length_zero = True if spinup_length == 0 else False


###############################
# Update the following in PFLOTRAN.in
# (1) The model running time
# (2) The restart file config
###############################
# If it is the first iteration, revise the pflotran.in file
if enks_mda_iteration_step == 1 and not update_obs_ens_posterior:
    # Read the current PFLOTRAN.in information
    if not os.path.isfile(pflotran_in_file):
        raise Exception("The pflotran in file does not exist in the path: %s" %
                        pflotran_in_file)
    with open(pflotran_in_file, 'r') as f:
        pflotranin = f.readlines()

    # Delete the current file
    os.remove(pflotran_in_file)

    # Update the model run time based on the assimilation window
    with open(pflotran_in_file, 'w') as f:
        for i, s in enumerate(pflotranin):
            if "FINAL_TIME" in s:
                # pflotranin[i] = "  FINAL_TIME {} sec".format(spinup_length_sec + current_model_time_end_sec) + "\n"
                pflotranin[i] = "  FINAL_TIME {} sec".format(current_model_time_end_sec) + "\n"

            if "SUBSURFACE_FLOW" in s and "MODE" in pflotranin[i + 1] and first_time_update and not is_spinup_length_zero:
                pflotranin.insert(i + 2, "        OPTIONS \n")
                pflotranin.insert(i + 3, "            REVERT_PARAMETERS_ON_RESTART \n")
                pflotranin.insert(i + 4, "        / \n")
            elif "SUBSURFACE_FLOW" in s and "MODE" in pflotranin[i + 1] and second_time_update and is_spinup_length_zero:
                pflotranin.insert(i + 2, "        OPTIONS \n")
                pflotranin.insert(i + 3, "            REVERT_PARAMETERS_ON_RESTART \n")
                pflotranin.insert(i + 4, "        / \n")

            if "CHECKPOINT" in s and "/" in pflotranin[i + 1] and first_time_update and not is_spinup_length_zero:
                pflotranin.insert(i + 2, "  RESTART \n")
                pflotranin.insert(i + 3, "    FILENAME " + pflotran_checkpoint_file + " \n")
                pflotranin.insert(i + 4, "    REALIZATION_DEPENDENT \n")
                pflotranin.insert(i + 5, "  / \n")
            elif "CHECKPOINT" in s and "/" in pflotranin[i + 1] and second_time_update and is_spinup_length_zero:
                pflotranin.insert(i + 2, "  RESTART \n")
                pflotranin.insert(i + 3, "    FILENAME " + pflotran_checkpoint_file + " \n")
                pflotranin.insert(i + 4, "    REALIZATION_DEPENDENT \n")
                pflotranin.insert(i + 5, "  / \n")

            if "SNAPSHOT_FILE" in s and first_time_update:
                pflotranin.insert(i + 1, "   PERIODIC TIME 300.0d0 sec \n")

        f.writelines(pflotranin)


###############################
# If it is the first iteration at the initial model time. No further posterior data is required.
###############################
if len(model_time_list) == 1 and enks_mda_iteration_step == 1 and not update_obs_ens_posterior:
    print("It is the first iteration at the initial model time. No further conversion from DART posterior is needed.")
    exit()

###############################
# Get the list of DART posterior files
###############################
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
# TODO Modify the parameter from a single value to 3D later on
# Update the parameter values in
# parameter_prior.h5 based on DART
# posterior output
###############################
if not os.path.isfile(pflotran_para_file):
    raise Exception("The file does not exist in the path: %s" % pflotran_para_file)
f_para = h5py.File(pflotran_para_file, "r+")

posterior = np.zeros([len(para_set), nens])
for i in range(nens):  # For each ensemble...
    dart_posterior_file = dart_posterior_file_set[i]
    # Open the posterior data
    nc_posterior = Dataset(dart_posterior_file, 'r')

    # Get the ensemble index
    # This is because the order of the file names in dart_posterior_file_list
    # may not follow 1,2,3,...,ens
    ens_num = int(nc_posterior.variables["member"][:])
    ens_ind = ens_num - 1

    if ens_ind != i:
        raise Exception("The ensemble member does not match! Check your DART posterior file list.")

    for j in range(len(para_set)):  # For each model parameter...
        varn = para_set[j]
        # Read the posterior data of para_set
        # posterior_data = nc_posterior.variables[varn][:][0, 0, 0]
        posterior_data = np.mean(nc_posterior.variables[varn][:])

        posterior[j, i] = posterior_data
        # # Replace it with the values in f_para
        # f_para[varn][:][ens_ind] = posterior_data

    # Close the posterior NetCDF file
    nc_posterior.close()

# Replace it with the values in f_para
for j in range(len(para_set)):
    varn     = para_set[j]
    var_dist = para_dist_set[j]

    var_min , var_max = para_min_set[j], para_max_set[j]
    var_mean, var_std = para_mean_set[j], para_std_set[j]

    # resample the prior if it is required and it is not the time for updating observation ensemble posterior
    if varn in para_resampled_set and not update_obs_ens_posterior:
        var_mean_nc = np.mean(posterior[j, :])
        var_std_nc  = np.std(posterior[j, :])
        # Generate the ensemble
        if rescaled:  # if rescaling the posterior at the previous time step is required
            posterior[j, :] = (posterior[j, :] - var_mean_nc) / var_std_nc * var_std + var_mean_nc

        else:  # if resampling is required.
            if var_dist.lower() == 'normal':
                posterior[j, :] = np.random.normal(var_mean_nc, var_std, nens)
                # Exclude those values outside of [minv, maxv]
                if var_min != -99999:
                    posterior[j, :][posterior[j, :] < var_min] = var_min
                if var_max != 99999:
                    posterior[j, :][posterior[j, :] > var_max] = var_max

            elif var_dist.lower() == 'lognormal':
                # logmean = np.exp(var_mean_nc + var_std**2 / 2.)
                # logstd  = np.exp(2 * var_mean_nc + var_std**2) * (np.exp(var_std**2) - 1)
                # posterior[j, :]  = np.random.lognormal(logmean, logstd)
                logvalues  = np.random.normal(var_mean_nc, var_std, nens)
                if var_min != -99999:
                    logvalues[logvalues < var_min] = var_min
                if var_max != 99999:
                    logvalues[logvalues > var_max] = var_max
                posterior[j, :] = np.power(10, logvalues)

            # elif var_dist.lower() == 'truncated_normal':
            #     posterior[j, :] = truncnorm.rvs(var_min, var_max, loc=var_mean_nc, scale=var_std, size=nens)

            elif var_dist.lower() == 'uniform':
                posterior[j, :] = np.random.uniform(var_min, var_max, nens)
                # Exclude those values outside of [minv, maxv]
                if var_min != -99999:
                    posterior[j, :][posterior[j, :] < var_min] = var_min
                if var_max != 99999:
                    posterior[j, :][posterior[j, :] > var_max] = var_max

            elif var_dist.lower() == 'test':
                posterior[j, :] = posterior[j, :]
                # Exclude those values outside of [minv, maxv]
                if var_min != -99999:
                    posterior[j, :][posterior[j, :] < var_min] = var_min
                if var_max != 99999:
                    posterior[j, :][posterior[j, :] > var_max] = var_max

            else:
                raise Exception("unknown distribution %s" % var_dist)


    else:
        # Exclude those values outside of [minv, maxv]
        # if var_dist.lower() == 'uniform' or var_dist.lower() == 'truncated_normal':
        posterior[j, :][posterior[j, :] < var_min] = var_min
        posterior[j, :][posterior[j, :] > var_max] = var_max

    f_para[varn][:] = posterior[j, :]
    # f_para[varn][:] = posterior[j, :]

# Close the parameter_prior.h5 file
f_para.close()
