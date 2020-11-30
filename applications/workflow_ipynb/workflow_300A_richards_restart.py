# This script is used for restarting DART-PFLOTRAN.
# By restarting DART-PFLOTRAN, we are actually referring to continue the data assimilation starting from the end of the last run
#
# Therefore, the following should be updated:
# (1) the last observation time
# (2) the numbers of cores for running DART and PFLOTRAN.

# %%
import os
import time
import f90nml
import subprocess


# %%
# Main directory names
app_par_dir  = "/global/cscratch1/sd/peishi89/300A-analysis"
app_dir_name = "Richards_av_std1.5_restart"
app_dir      = os.path.join(app_par_dir, app_dir_name)          # The application folder name
# dart_dir     = "/Users/jian449/Codes/DART/manhattan"
# dart_pf_dir  = "/Users/jian449/Codes/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name


# %%
# The numbers of cores for DART and PFLOTRAN
ncore_da = 10  # The number of MPI cores used by DART
ncore_pf = 380  # The number of MPI cores used by PFLOTRAN


# %%
# Config the additional simulation/assimilation period time
additional_simulation_time_days = 2
additional_simulation_time_seconds = 0
additional_simulation_time_size = additional_simulation_time_days+float(additional_simulation_time_seconds)/86400. # day
# last_obs_time_days = 360
# last_obs_time_seconds = 129600
# last_obs_time_size = 361.5


# %%
#--------------- The following is not required to be changed by users --------------------

# %%
# Check the existence of the config.nml and input.nml files
config_nml = os.path.join(app_dir, "work/config.nml")
input_nml  = os.path.join(app_dir, "work/input.nml")
if not os.path.isfile(config_nml):
    raise Exception("The configuration file does not exist under: {}".format(config_nml))
else:
    configs = f90nml.read(config_nml)
if not os.path.isfile(input_nml):
    raise Exception("The input namelist file does not exist under: {}".format(input_nml))


# %%
# Change the locations of the all the files saved in the original config_file according to the new defined locations if necessary 
# (this is used for the case that the application folder is moved/copied from another location)
dart_pf_dir = configs['main_dir_cfg']['dart_pf_dir']
change_path_name_files = os.path.join(dart_pf_dir, "utils/change_file_paths_in_confignml.py")
subprocess.run("python {} {} {} {}".format(change_path_name_files, config_nml, dart_pf_dir, app_dir), shell=True, check=True)
# Reload the configuration
configs = f90nml.read(config_nml)


# %%
# Check whether the previous run is done
if configs['time_cfg']['exceeds_obs_time']:
    configs['time_cfg']['exceeds_obs_time'] = False
else:
    raise Exception("The previous assimilation does noe end normally. Please check!")


# %%
# Conduct the changed configuration
# Numbers of cores
configs['exe_cfg']['ncore_da'] = ncore_da
configs['exe_cfg']['ncore_pf'] = ncore_pf
# Time
configs['time_cfg']['last_obs_time_days'] += additional_simulation_time_days
configs['time_cfg']['last_obs_time_seconds'] += additional_simulation_time_seconds
configs['time_cfg']['last_obs_time_size'] += additional_simulation_time_size
# configs['time_cfg']['enks_mda_'] += additional_simulation_time_size
configs.write(config_nml, force=True)
# Update other time configurations accordingly using update_confignml_time.py
update_confignml_time_file = configs['file_cfg']['update_confignml_time_file']
subprocess.run("python {} {} {}".format(update_confignml_time_file, config_nml, input_nml), shell=True, check=True)


# %%
# Reclip the observations
print("\n")
print("------------------------------------------------------------")
print("Clip the NetCDF file based on the defined spatial and temporal domains...")
clip_obs_nc = configs["file_cfg"]["clip_obs_nc_file"]
subprocess.run("python {} {}".format(clip_obs_nc, config_nml), shell=True, check=True)


# %%
# In the following, run the shell script to couple DART and PFLOTRAN.
dart_work_dir      = configs['other_dir_cfg']["dart_work_dir"]
run_DART_PFLOTRAN  = configs['file_cfg']["run_da_csh"]
concatenate_output = configs['file_cfg']["concatenate_dart_output_file"]
start_time         = time.time()
# exit()

# %%
# Run the script
print("\n")
print("------------------------------------------------------------")
print("Assimilation starts here...")
subprocess.run("cd {}; csh {} {} {}".format(dart_work_dir, run_DART_PFLOTRAN, input_nml, config_nml), shell=True, check=True)


# %%
print("\n")
print("------------------------------------------------------------")
print("Concatenate the prior and posterior at all times ...")
subprocess.run("python {} {}".format(concatenate_output, config_nml), shell=True, check=True)


# %%
end_time = time.time()
print("The total time usage of running DART and PFLOTRAN is %.3f (second): " % (end_time-start_time))