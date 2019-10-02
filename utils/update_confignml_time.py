"""Update the model/obs window times one based on the
assimilation window in the config.nml and input.nml"""

# Author: Peishi Jiang

import sys
import f90nml

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)


###############################
# Update the model/obs window times in config.nml
###############################
current_model_time   = configs["time_cfg"]["current_model_time"]
assim_start_days     = configs["da_cfg"]["assim_start_days"]
assim_start_seconds  = configs["da_cfg"]["assim_start_seconds"]
assim_end_days       = configs["da_cfg"]["assim_end_days"]
assim_end_seconds    = configs["da_cfg"]["assim_end_seconds"]
assim_window_days    = configs["da_cfg"]["assim_window_days"]
assim_window_seconds = configs["da_cfg"]["assim_window_seconds"]

# Update the model time
configs["time_cfg"]["current_model_time"] += assim_window_days + float(assim_window_seconds)/86400.

print(configs["time_cfg"]["current_model_time"])
# Update the observation window time
assim_start_days    += assim_window_days
assim_start_seconds += assim_window_seconds
assim_end_days      += assim_window_days
assim_end_seconds   += assim_window_seconds

configs["da_cfg"]["assim_start_days"]    = assim_window_days
configs["da_cfg"]["assim_start_seconds"] = assim_window_seconds
configs["da_cfg"]["assim_end_days"]      = assim_window_days
configs["da_cfg"]["assim_end_seconds"]   = assim_window_seconds

###############################
# Update the obs window times in input.nml
###############################
input_nml_file = configs["file_cfg"]["input_nml_file"]
input_nml      = f90nml.read(input_nml_file)

input_nml["convert_nc_nml"]["obs_start_day"]    = assim_start_days
input_nml["convert_nc_nml"]["obs_start_second"] = assim_start_seconds
input_nml["convert_nc_nml"]["obs_end_day"]      = assim_end_days
input_nml["convert_nc_nml"]["obs_end_second"]   = assim_end_seconds

###############################
# Save them
###############################
configs.write(config_nml_file, force=True)
input_nml.write(input_nml_file, force=True)
