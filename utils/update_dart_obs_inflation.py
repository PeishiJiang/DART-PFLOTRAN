"""This file is used for updating the variance of observation in DART based on the inflation coefficient defined in EnKS-MDA."""

import os
import sys
import f90nml
import subprocess


###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

# Read in the required configurations
input_nml_file            = configs["file_cfg"]["input_nml_file"]
convert_nc_exe            = configs["file_cfg"]["convert_nc_exe"]
current_model_time        = configs["time_cfg"]["current_model_time"]
model_time_list           = configs["time_cfg"]["model_time_list"]
last_obs_time_days        = configs["time_cfg"]["last_obs_time_days"]
last_obs_time_seconds     = configs["time_cfg"]["last_obs_time_seconds"]
exceeds_obs_time          = configs["time_cfg"]["exceeds_obs_time"]
assim_start_days          = configs["da_cfg"]["assim_start_days"]
assim_start_seconds       = configs["da_cfg"]["assim_start_seconds"]
assim_end_days            = configs["da_cfg"]["assim_end_days"]
assim_end_seconds         = configs["da_cfg"]["assim_end_seconds"]
assim_window_days         = configs["da_cfg"]["assim_window_days"]
assim_window_seconds      = configs["da_cfg"]["assim_window_seconds"]
enks_mda_iteration_step   = configs["da_cfg"]["enks_mda_iteration_step"]
enks_mda_alpha            = configs["da_cfg"]["enks_mda_alpha"]
enks_mda_total_iterations = configs["da_cfg"]["enks_mda_total_iterations"]

if not isinstance(enks_mda_alpha, list):
    enks_mda_alpha = [enks_mda_alpha]


###############################
# Update the convert_nc_nml in input.nml
###############################
# Read in the template file
nml = f90nml.read(input_nml_file)

# Update the convert_nc_nml
nml["convert_nc_nml"]["inflation_alpha"] = enks_mda_alpha[enks_mda_iteration_step-1]

# Save input.nml
nml.write(input_nml_file, force=True)


###############################
# Generate the corresponding DART obs data
###############################
subprocess.run(convert_nc_exe, shell=True, check=True)


###############################
# Update the iteration step if it is not the last iteration
###############################
if enks_mda_iteration_step <= enks_mda_total_iterations:
    configs["da_cfg"]["enks_mda_iteration_step"] = enks_mda_iteration_step + 1
    configs.write(config_nml_file, force=True)