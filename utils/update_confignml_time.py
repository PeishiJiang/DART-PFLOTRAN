"""Update the model/obs window times one based on the
assimilation window in the config.nml and input.nml, and
change the MDA iteration to 1 again."""

# Author: Peishi Jiang

import sys
import f90nml

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

# Read in the required configurations
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
enks_mda_total_iterations = configs["da_cfg"]["enks_mda_total_iterations"]

# Get the last_obs_time in days
last_obs_time = last_obs_time_days + float(last_obs_time_seconds) / 86400.


###############################
# Change the EnKS-MDA iteration step back to 1
###############################
if enks_mda_iteration_step - 1 != enks_mda_total_iterations:
    raise Exception("The current iteration step {} is not the same as the total iterations {}".format(enks_mda_iteration_step, enks_mda_total_iterations))
configs["da_cfg"]["enks_mda_iteration_step"] = 1


###############################
# Check if the updated current_model_time exceeds the last_obs_time
###############################
current_model_time += assim_window_days + float(assim_window_seconds) / 86400.
if current_model_time >= (last_obs_time - 1e-8):
    configs["time_cfg"]["exceeds_obs_time"] = True

    ###############################
    # Save them
    ###############################
    configs.write(config_nml_file, force=True)

else:
    ###############################
    # Update the model/obs window times in config.nml
    ###############################
    # Update the model time
    configs["time_cfg"]["current_model_time"] = current_model_time

    # Update the list model time
    if not isinstance(model_time_list, list):
        model_time_list = [model_time_list]
    model_time_list.append(current_model_time)
    configs["time_cfg"]["model_time_list"] = model_time_list

    # print(configs["time_cfg"]["current_model_time"])
    # Update the observation window time
    assim_start_days    += assim_window_days
    assim_start_seconds += assim_window_seconds
    assim_end_days      += assim_window_days
    assim_end_seconds   += assim_window_seconds

    assim_start_time = assim_start_days + float(assim_start_seconds) / 86400.
    assim_end_time   = assim_end_days + float(assim_end_seconds) / 86400.

    # If the assimilation start or end times are larger than the last observation time,
    # set it/them as the last observation time
    if assim_start_time > last_obs_time:
        assim_start_days    = last_obs_time_days
        assim_start_seconds = last_obs_time_seconds
    if assim_end_time > last_obs_time:
        assim_end_days    = last_obs_time_days
        assim_end_seconds = last_obs_time_seconds

    configs["da_cfg"]["assim_start_days"]    = assim_start_days
    configs["da_cfg"]["assim_start_seconds"] = assim_start_seconds
    configs["da_cfg"]["assim_end_days"]      = assim_end_days
    configs["da_cfg"]["assim_end_seconds"]   = assim_end_seconds

    ###############################
    # Update the obs window times in input.nml for both smoother_nml and convert_nc_nml
    ###############################
    input_nml_file = configs["file_cfg"]["input_nml_file"]
    input_nml      = f90nml.read(input_nml_file)

    input_nml["smoother_nml"]["first_obs_days"]    = assim_start_days
    input_nml["smoother_nml"]["first_obs_seconds"] = assim_start_seconds
    input_nml["smoother_nml"]["last_obs_days"]     = assim_end_days
    input_nml["smoother_nml"]["last_obs_seconds"]  = assim_end_seconds

    input_nml["convert_nc_nml"]["obs_start_day"]    = assim_start_days
    input_nml["convert_nc_nml"]["obs_start_second"] = assim_start_seconds
    input_nml["convert_nc_nml"]["obs_end_day"]      = assim_end_days
    input_nml["convert_nc_nml"]["obs_end_second"]   = assim_end_seconds

    ###############################
    # Save them
    ###############################
    configs.write(config_nml_file, force=True)
    input_nml.write(input_nml_file, force=True)
