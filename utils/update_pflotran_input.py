"""Update the PFLOTRAN input files based on the DART posterior output data"""

# Author: Peishi Jiang

import os
# Put the following before importing numpy is to avoid the following error
# OpenBLAS blas_thread_init: pthread_create failed for thread 19 of 64: Resource temporarily unavailable
# See more at: https://stackoverflow.com/questions/52026652/openblas-blas-thread-init-pthread-create-resource-temporarily-unavailable
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import re
import sys
import h5py
import f90nml
import subprocess
import numpy as np
from scipy.stats import truncnorm
from math import ceil, log10
from netCDF4 import num2date, date2num, Dataset

from parse_pflotran_files_utils import pflotran_files

# TODO: For now, a single value for each variable is assumed, it should be more generic in the future by considering the 3D case.

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

pflotran_in_file               = configs["file_cfg"]["pflotran_in_file"]
pflotran_para_file             = configs["file_cfg"]["pflotran_para_file"]
pflotran_para_backup_file      = configs["file_cfg"]["pflotran_para_backup_file"]
use_para_initial_at_nth_window = configs["da_cfg"]["use_para_initial_at_nth_window"]
pflotran_restart_file          = configs["file_cfg"]["pflotran_restart_file"]
dart_prior_file                = configs["file_cfg"]["dart_prior_nc_file"]
dart_posterior_file            = configs["file_cfg"]["dart_posterior_nc_file"]
dart_input_list                = configs["file_cfg"]["dart_input_list_file"]
dart_output_list               = configs["file_cfg"]["dart_output_list_file"]
use_obs_tecfile_for_prior      = configs["file_cfg"]["use_obs_tecfile_for_prior"]
save_immediate_mda_result      = configs["file_cfg"]["save_immediate_mda_result"]

spinup_length           = configs["time_cfg"]["spinup_length"]
model_time_list         = configs["time_cfg"]["model_time_list"]
current_model_time      = configs["time_cfg"]["current_model_time"]
enks_mda_iteration_step = configs["da_cfg"]["enks_mda_iteration_step"]
assim_window_fixed      = configs["da_cfg"]["assim_window_fixed"]
total_iterations        = configs["da_cfg"]["enks_mda_total_iterations"]
iteration_step          = configs["da_cfg"]["enks_mda_iteration_step"]

try:
    update_obs_ens_posterior_now = configs["da_cfg"]["update_obs_ens_posterior_now"]
except:
    update_obs_ens_posterior_now = False

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
    if not isinstance(assim_window_list, list):
        assim_window_list = [assim_window_list]
    assim_time_sofar = np.sum(assim_window_list[:len(model_time_list)])
    current_model_time_end = spinup_length + assim_time_sofar

current_model_time_end_sec = current_model_time_end * 86400
spinup_length_sec          = spinup_length * 86400
current_model_time_end_hr  = current_model_time_end * 24
spinup_length_hr           = spinup_length * 24

para_set           = configs["obspara_set_cfg"]["para_set"]
nens               = configs["da_cfg"]["nens"]

if not isinstance(para_set, list):
    para_set = [para_set]

# Check whether this is the first time to update pflotran input files
first_time_update = True if len(model_time_list) == 1 else False

# Check whether this is the first time to update pflotran input files
second_time_update = True if len(model_time_list) == 2 else False

# Check whether the spinup was conducted before
is_spinup_length_zero = True if spinup_length == 0 else False

# Get some constants
# TODO: ntimestep will change when more restarts are involved.
ntimestep      = int(configs["da_cfg"]["ntimestep"])
nens           = configs["da_cfg"]["nens"]
ndigit_time    = int(ceil(log10(ntimestep))) + 1
ndigit_ens     = int(ceil(log10(nens))) + 1
ens_set        = np.arange(1, nens + 1)  # the set of ensembles
# Get the right iteration step for posterior file
if iteration_step == 1:
    model_time_ind = len(model_time_list)-1  # the ith model step
    post_file_iteration_step = total_iterations
else:
    model_time_ind = len(model_time_list)  # the ith model step
    post_file_iteration_step = iteration_step-1

###############################
# Update the following in PFLOTRAN.in
# (1) The model running time
# (2) The restart file config
# (3) TODO: print out time
###############################
# If it is the first iteration, revise the pflotran.in file
if enks_mda_iteration_step == 1 and not update_obs_ens_posterior_now:
    # Read the current PFLOTRAN.in information
    if not os.path.isfile(pflotran_in_file):
        raise Exception("The pflotran in file does not exist in the path: %s" %
                        pflotran_in_file)

    # Read the pflotran input deck first to get what is inside, regarding
    # REVERT_PARAMETERS_ON_RESTART in OPTIONS
    # CHECKPOINT
    # RESTART FILENAME
    # RESTART REALIZATION_DEPENDENT
    is_revert_parameters_on_restart_on                 = False
    is_checkpoint_on, checkpoint_section_nrow          = False, 1
    is_restart_on, is_restart_realization_dependent_on = False, False
    restart_filename_nrow                              = 0
    is_snap_on, snap_periodtime_nrow, snap_time_nrow   = False, 0, 0
    is_obs_on, obs_periodtime_nrow, obs_time_nrow      = False, 0, 0
    with open(pflotran_in_file, 'r') as f:
        pflotranin = f.readlines()
        for i, s in enumerate(pflotranin):
            if "REVERT_PARAMETERS_ON_RESTART" in s and "#" not in s:
                is_revert_parameters_on_restart_on = True
            
            if "CHECKPOINT" in s and "#" not in s:
                is_checkpoint_on, k = True, 1
                # while "/" not in pflotranin[i+k]: k += 1
                while not pflotranin[i+k].strip().startswith("/"): k += 1
            
            if "RESTART" in s.split() and "#" not in s:
                is_restart_on, k = True, 1
                # while "/" not in pflotranin[i+k]:
                while not pflotranin[i+k].strip().startswith("/"):
                    if "REALIZATION_DEPENDENT" in pflotranin[i+k] and "#" not in pflotranin[i+k]: is_restart_realization_dependent_on = True
                    if "FILENAME" in pflotranin[i+k] and "#" not in pflotranin[i+k]: restart_filename_nrow = k
                    k += 1
            
            if "SNAPSHOT_FILE" in s and "#" not in s:
                is_snap_on, k = True, 1
                # while "/" not in pflotranin[i+k]:
                while not pflotranin[i+k].strip().startswith("/"):
                    if "TIMES" in pflotranin[i+k].split() and "#" not in pflotranin[i+k]: snap_time_nrow = k
                    if set(["PERIODIC","TIME"]).issubset(set(pflotranin[i+k].split())) and "#" not in pflotranin[i+k]: snap_periodtime_nrow = k
                    k += 1

            if "OBSERVATION_FILE" in s and "#" not in s:
                is_snap_on, k = True, 1
                # while "/" not in pflotranin[i+k]:
                while not pflotranin[i+k].strip().startswith("/"):
                    if "TIMES" in pflotranin[i+k].split() and "#" not in pflotranin[i+k]: obs_time_nrow = k
                    if set(["PERIODIC","TIME"]).issubset(set(pflotranin[i+k].split())) and "#" not in pflotranin[i+k]: obs_periodtime_nrow = k
                    k += 1

    # Read the file again and delete the current file
    with open(pflotran_in_file, 'r') as f:
        pflotranin = f.readlines()
    os.remove(pflotran_in_file)

    # Update the model run time based on the assimilation window
    with open(pflotran_in_file, 'w') as f:
        for i, s in enumerate(pflotranin):
            #  --------- Fixed the FINAL_TIME according to the current assimilation window ----------
            if "FINAL_TIME" in s:
                # pflotranin[i] = "  FINAL_TIME {} sec".format(spinup_length_sec + current_model_time_end_sec) + "\n"
                # pflotranin[i] = "  FINAL_TIME {} sec".format(current_model_time_end_sec) + "\n"
                pflotranin[i] = "  FINAL_TIME {} h".format(current_model_time_end_hr) + "\n"

            #  --------- Turn on REVERT_PARAMETERS_ON_RESTART if not ----------
            if not is_revert_parameters_on_restart_on:
                if "SUBSURFACE_FLOW" in s and "MODE" in pflotranin[i + 1] and first_time_update and not is_spinup_length_zero:
                    pflotranin.insert(i + 2, "        OPTIONS \n")
                    pflotranin.insert(i + 3, "            REVERT_PARAMETERS_ON_RESTART \n")
                    pflotranin.insert(i + 4, "        / \n")
                elif "SUBSURFACE_FLOW" in s and "MODE" in pflotranin[i + 1] and second_time_update and is_spinup_length_zero:
                    pflotranin.insert(i + 2, "        OPTIONS \n")
                    pflotranin.insert(i + 3, "            REVERT_PARAMETERS_ON_RESTART \n")
                    pflotranin.insert(i + 4, "        / \n")

            #  --------- Add RESTART card ----------
            # if RESTART is not there, add it after CHECKPOINT card
            if not is_restart_on and "CHECKPOINT" in s and "#" not in s:
                if (first_time_update and not is_spinup_length_zero) or (second_time_update and is_spinup_length_zero):
                    pflotranin.insert(i + checkpoint_section_nrow + 1, "  RESTART \n")
                    if "[ENS]" in pflotran_restart_file:
                        pflotran_restart_file = re.sub(r"R\[ENS\]", "", pflotran_restart_file)
                        pflotranin.insert(i + checkpoint_section_nrow + 2, "    FILENAME " + pflotran_restart_file + " \n")
                        pflotranin.insert(i + checkpoint_section_nrow + 3, "    REALIZATION_DEPENDENT \n")
                        pflotranin.insert(i + checkpoint_section_nrow + 4, "  / \n")
                    else:
                        pflotranin.insert(i + checkpoint_section_nrow + 2, "    FILENAME " + pflotran_restart_file + " \n")
                        pflotranin.insert(i + checkpoint_section_nrow + 3, "  / \n")
            # if RESTART is on, check whether realization dependent is on and add it if not
            if is_restart_on and "RESTART" in s.split() and len(model_time_list) != 1 and "#" not in s:
                if restart_filename_nrow == 0:
                    raise Exception("There is no FILENAME in the original RESTART card in pflotran!")
                if "[ENS]" in pflotran_restart_file and not is_restart_realization_dependent_on:
                    pflotran_restart_file = re.sub(r"R\[ENS\]", "", pflotran_restart_file)
                    # pflotranin.insert(i + restart_filename_nrow, "    FILENAME " + pflotran_restart_file + " \n")
                    pflotranin[i + restart_filename_nrow] = "    FILENAME " + pflotran_restart_file + " \n"
                    pflotranin.insert(i + restart_filename_nrow + 1, "    REALIZATION_DEPENDENT \n")
                elif "[ENS]" in pflotran_restart_file:
                    pflotran_restart_file = re.sub(r"R\[ENS\]", "", pflotran_restart_file)
                    # pflotranin.insert(i + restart_filename_nrow, "    FILENAME " + pflotran_restart_file + " \n")
                    pflotranin[i + restart_filename_nrow] = "    FILENAME " + pflotran_restart_file + " \n"
                else:
                    pflotranin[i + restart_filename_nrow] = "    FILENAME " + pflotran_restart_file + " \n"
                    # pflotranin.insert(i + restart_filename_nrow, "    FILENAME " + pflotran_restart_file + " \n")

            #  --------- Revised the SNAPSHOT_FILE ----------
            # TODO: irregular time point should be added later on
            # TODO: flexible edits on time periods should be done here.
            if not use_obs_tecfile_for_prior and "SNAPSHOT_FILE" in s:
                # pflotranin.insert(i + snap_periodtime_nrow, "   PERIODIC TIME 300.0d0 sec \n")
                if snap_periodtime_nrow != 0:
                    pflotranin[i + snap_periodtime_nrow] = "   PERIODIC TIME 300.0d0 sec \n"
                else:
                    pflotranin.insert(i + 1, "   PERIODIC TIME 300.0d0 sec \n")
            elif use_obs_tecfile_for_prior and "OBSERVATION_FILE" in s:
                # pflotranin.insert(i + 1, "   PERIODIC TIME 300.0d0 sec \n")
                # pflotranin.insert(i + obs_periodtime_nrow, "   PERIODIC TIME 1.0d0 h \n")
                if obs_periodtime_nrow != 0:
                    pflotranin[i + obs_periodtime_nrow] = "   PERIODIC TIME 1.0d0 h \n"
                else:
                    pflotranin.insert(i + 1, "   PERIODIC TIME 1.0d0 h \n")

        f.writelines(pflotranin)


###############################
# If it is the first iteration at the initial model time. No further posterior data is required.
###############################
if len(model_time_list) == 1 and enks_mda_iteration_step == 1 and not update_obs_ens_posterior_now:
    print("It is the first iteration at the initial model time. No further conversion from DART posterior is needed.")
    exit()


###############################
# Get the lists of netcdf file names for DART inputs and outputs
###############################
# Get the file names of all ensembles for DART restart file
dart_prior_file_set = [
    # re.sub(r"\[ENS\]", str(ens) + "_time" + str(model_time_ind), dart_prior_file)
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file)
    for ens in ens_set
]
dart_prior_file_set = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_prior_file_each)
    for dart_prior_file_each in dart_prior_file_set
]
dart_posterior_file_set = [
    # re.sub(r"\[ENS\]", str(ens) + "_time" + str(model_time_ind), dart_posterior_file)
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_file)
    for ens in ens_set
]
dart_posterior_file_set = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_posterior_file_each)
    for dart_posterior_file_each in dart_posterior_file_set
]
if save_immediate_mda_result:
    dart_prior_file_set = [
        re.sub("[.]nc", "_iter"+str(1)+".nc", dart_prior_file_each)
        for dart_prior_file_each in dart_prior_file_set
    ]
    dart_posterior_file_set = [
        # TODO: use re
        # re.sub("[.]nc", "_iter"+str(total_iterations)+".nc", dart_posterior_file_each)
        re.sub("[.]nc", "_iter"+str(post_file_iteration_step)+".nc", dart_posterior_file_each)
        for dart_posterior_file_each in dart_posterior_file_set
    ]

# print(dart_posterior_file_set)
# print(dart_prior_file_set)
# exit()
# ###############################
# # Get the list of DART prior and posterior files
# ###############################
# # prior
# dart_prior_file_set = []
# if not os.path.isfile(dart_input_list):
#     raise Exception("The DART output list file does not exist in the path: %s" % dart_input_list)
# with open(dart_input_list, "r") as f:
#     dart_prior_file_list = f.readlines()

# for i in range(len(dart_prior_file_list)):
#     file_name_old = dart_prior_file_list[i]
#     file_name     = file_name_old.replace("\n", "")
#     dart_prior_file_set.append(file_name)

# # posterior
# dart_posterior_file_set = []
# if not os.path.isfile(dart_output_list):
#     raise Exception("The DART output list file does not exist in the path: %s" % dart_output_list)
# with open(dart_output_list, "r") as f:
#     dart_posterior_file_list = f.readlines()

# for i in range(len(dart_posterior_file_list)):
#     file_name_old = dart_posterior_file_list[i]
#     file_name     = file_name_old.replace("\n", "")
#     dart_posterior_file_set.append(file_name)


###############################
# Check whether the parameter file should be adopted from
# the initial.
###############################
if use_para_initial_at_nth_window == len(model_time_list) and not update_obs_ens_posterior_now and enks_mda_iteration_step == 1:
    subprocess.run("mv {} {}".format(pflotran_para_backup_file, pflotran_para_file), shell=True, check=True)
    sys.exit()


###############################
# TODO Modify the parameter from a single value to 3D later on
# Get the prior and posterior of the parameters
###############################
if not os.path.isfile(pflotran_para_file):
    raise Exception("The file does not exist in the path: %s" % pflotran_para_file)
# f_para = h5py.File(pflotran_para_file, "r+")

posterior = dict.fromkeys(para_set)
prior     = dict.fromkeys(para_set)
for i in range(nens):  # For each ensemble...
    dart_prior_file     = dart_prior_file_set[i]
    dart_posterior_file = dart_posterior_file_set[i]

    # Prior
    with Dataset(dart_prior_file, 'r') as nc_prior:
        # Get the ensemble index
        # This is because the order of the file names in dart_prior_file_list may not follow 1,2,3,...,ens
        ens_num = int(nc_prior.variables["member"][:])
        ens_ind = ens_num - 1

        if ens_ind != i:
            raise Exception("The ensemble member does not match! Check your DART prior file list.")

        # Get the spatial and temporal dimension and further initialize prior dictionary
        # Get the cell ids and physical location of these parameters
        if i == 0:
            # Spatio-temporal dimensions
            ntime_para = nc_prior.dimensions["para_time"].size
            if ntime_para != 1: raise Exception("The parameters in each assimilation window should be temporally independent.")
            # nloc_para  = nc_prior.dimensions["para_location"].size
            # for varn in para_set:
            #     # prior[varn]["value"] = np.zeros([nens, ntime_para, nloc_para])
            #     prior[varn] = {"value": np.zeros([nens, nloc_para])}
            # # Get the cell ids
            # cell_ids = nc_prior.variables['para_cell_ids'][:]
            # # Get the physical locations of the parameters
            # xloc_para = nc_prior.variables['para_x_location'][:]
            # yloc_para = nc_prior.variables['para_y_location'][:]
            # zloc_para = nc_prior.variables['para_z_location'][:]
            # para_loc_set = np.array([xloc_para, yloc_para, zloc_para])

        # Get the value for each model parameter
        for varn in para_set:
            prior_data = nc_prior.variables[varn][:]
            # prior_data = prior_data.filled(np.nan)

            # Cell ids and physical locations
            if i == 0:
                cell_ids     = nc_prior.variables['para_cell_ids'][:]
                xloc_para    = nc_prior.variables['para_x_location'][:]
                yloc_para    = nc_prior.variables['para_y_location'][:]
                zloc_para    = nc_prior.variables['para_z_location'][:]
                para_loc_set = np.array([xloc_para, yloc_para, zloc_para])

                if prior_data.mask != False:
                    cell_ids_varn     = cell_ids[~prior_data.mask[0,:]]
                    para_loc_set_varn = para_loc_set[:,~prior_data.mask[0,:]]
                else:
                    cell_ids_varn = cell_ids
                    para_loc_set_varn = para_loc_set
                nloc_varn         = len(cell_ids_varn)

                prior[varn] = {"value": np.zeros([nens, nloc_varn]), 
                               "cell_ids": cell_ids_varn, "loc_set": para_loc_set_varn}
            prior_data        = prior_data.data[~prior_data.mask]
            prior[varn]["value"][i,:] = prior_data


    # Posterior
    with Dataset(dart_posterior_file, 'r') as nc_posterior:
        # Get the ensemble index
        # This is because the order of the file names in dart_posterior_file_list may not follow 1,2,3,...,ens
        ens_num = int(nc_posterior.variables["member"][:])
        ens_ind = ens_num - 1

        if ens_ind != i:
            raise Exception("The ensemble member does not match! Check your DART posterior file list.")

        # Get the spatial and temporal dimension and further initialize posterior dictionary
        if i == 0:
            ntime_para = nc_posterior.dimensions["para_time"].size
            if ntime_para != 1: raise Exception("The parameters in each assimilation window should be temporally independent.")
            # nloc_para  = nc_posterior.dimensions["para_location"].size
            # for varn in para_set:
            #     # posterior[varn]["value"] = np.zeros([nens, ntime_para, nloc_para])
            #     # posterior[varn]["value"] = np.zeros([nens, nloc_para])
            #     posterior[varn] = {"value": np.zeros([nens, nloc_para])}

        # Get the value for each model parameter
        for varn in para_set:
            posterior_data = nc_posterior.variables[varn][:]
            # posterior_data = posterior_data.filled(np.nan)

            # Cell ids and physical locations
            if i == 0:
                cell_ids     = nc_posterior.variables['para_cell_ids'][:]
                xloc_para    = nc_posterior.variables['para_x_location'][:]
                yloc_para    = nc_posterior.variables['para_y_location'][:]
                zloc_para    = nc_posterior.variables['para_z_location'][:]
                para_loc_set = np.array([xloc_para, yloc_para, zloc_para])

                if posterior_data.mask != False:
                    cell_ids_varn     = cell_ids[~posterior_data.mask[0,:]]
                    para_loc_set_varn = para_loc_set[:,~posterior_data.mask[0,:]]
                else:
                    cell_ids_varn = cell_ids
                    para_loc_set_varn = para_loc_set
                nloc_varn         = len(cell_ids_varn)

                posterior[varn] = {"value": np.zeros([nens, nloc_varn]), 
                                   "cell_ids": cell_ids_varn, "loc_set": para_loc_set_varn}
            posterior_data                = posterior_data.data[~posterior_data.mask]
            posterior[varn]["value"][i,:] = posterior_data


###############################
# Update the posterior parameters in PFLOTRAN parameter file
###############################
pflotran_file_parser = pflotran_files(config_nml)
pflotran_file_parser.update_parameters_in_pflotran_parameter_file(prior, posterior)
# pflotran_file_parser.update_parameters_in_pflotran_parameter_file(prior, posterior, cell_ids, para_loc_set)

# # Replace it with the values in f_para
# for j in range(len(para_set)):
#     varn     = para_set[j]
#     var_dist = para_dist_set[j]

#     var_min , var_max = para_min_set[j], para_max_set[j]
#     var_mean, var_std = para_mean_set[j], para_std_set[j]

#     # resample the prior if 
#     # (1) it is required and 
#     # (2) it is not the time for updating observation ensemble posterior
#     # (3) it is not during ES-MDA iteration 
#     if varn in para_resampled_set and not update_obs_ens_posterior_now and enks_mda_iteration_step == 1:
#         var_mean_nc = np.mean(posterior[j, :])
#         var_std_nc  = np.std(posterior[j, :])
#         # Generate the ensemble
#         if rescaled:  # if rescaling the posterior at the previous time step is required
#             posterior[j, :] = (posterior[j, :] - var_mean_nc) / var_std_nc * var_std + var_mean_nc
#             if var_dist.lower() == 'lognormal':
#                 posterior[j, :] = np.power(10, logvalues)

#         else:  # if resampling is required.
#             if var_dist.lower() == 'normal':
#                 posterior[j, :] = np.random.normal(var_mean_nc, var_std, nens)
#                 # Exclude those values outside of [minv, maxv]
#                 if var_min != -99999:
#                     posterior[j, :][posterior[j, :] < var_min] = var_min
#                 if var_max != 99999:
#                     posterior[j, :][posterior[j, :] > var_max] = var_max

#             elif var_dist.lower() == 'lognormal':
#                 # logmean = np.exp(var_mean_nc + var_std**2 / 2.)
#                 # logstd  = np.exp(2 * var_mean_nc + var_std**2) * (np.exp(var_std**2) - 1)
#                 # posterior[j, :]  = np.random.lognormal(logmean, logstd)
#                 logvalues  = np.random.normal(var_mean_nc, var_std, nens)
#                 if var_min != -99999:
#                     logvalues[logvalues < var_min] = var_min
#                 if var_max != 99999:
#                     logvalues[logvalues > var_max] = var_max
#                 posterior[j, :] = np.power(10, logvalues)

#             # elif var_dist.lower() == 'truncated_normal':
#             #     posterior[j, :] = truncnorm.rvs(var_min, var_max, loc=var_mean_nc, scale=var_std, size=nens)

#             elif var_dist.lower() == 'uniform':
#                 posterior[j, :] = np.random.uniform(var_min, var_max, nens)
#                 # Exclude those values outside of [minv, maxv]
#                 if var_min != -99999:
#                     posterior[j, :][posterior[j, :] < var_min] = var_min
#                 if var_max != 99999:
#                     posterior[j, :][posterior[j, :] > var_max] = var_max

#             elif var_dist.lower() == 'test':
#                 posterior[j, :] = posterior[j, :]
#                 # Exclude those values outside of [minv, maxv]
#                 if var_min != -99999:
#                     posterior[j, :][posterior[j, :] < var_min] = var_min
#                 if var_max != 99999:
#                     posterior[j, :][posterior[j, :] > var_max] = var_max

#             else:
#                 raise Exception("unknown distribution %s" % var_dist)


#     else:
#         if var_dist.lower() == 'lognormal':
#             posterior[j, :] = np.power(10, posterior[j, :])
#         # # Exclude those values outside of [minv, maxv]
#         # # if var_dist.lower() == 'uniform' or var_dist.lower() == 'truncated_normal':
#         # if var_dist.lower() == 'lognormal':
#         #     posterior[j, :][posterior[j, :] < var_min] = var_min
#         #     posterior[j, :][posterior[j, :] > var_max] = var_max
#         #     posterior[j, :] = np.power(10, posterior[j, :])
#         # else:
#         #     posterior[j, :][posterior[j, :] < var_min] = var_min
#         #     posterior[j, :][posterior[j, :] > var_max] = var_max

#     f_para[varn][:] = posterior[j, :]
#     # f_para[varn][:] = posterior[j, :]

# # Close the parameter_prior.h5 file
# f_para.close()