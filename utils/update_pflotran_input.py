"""Update the PFLOTRAN input files based on the DART posterior output data"""

# Author: Peishi Jiang

import os
import sys
import h5py
import f90nml
import numpy as np
from netCDF4 import num2date, date2num, Dataset

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_in_file    = configs["file_cfg"]["pflotran_in_file"]
pflotran_para_file  = configs["file_cfg"]["pflotran_para_file"]
pflotran_checkpoint_file = configs["file_cfg"]["pflotran_checkpoint_file"]
dart_posterior_file = configs["file_cfg"]["dart_posterior_nc_file"]
dart_output_list    = configs["file_cfg"]["dart_output_list_file"]
model_time_list     = configs["time_cfg"]["model_time_list"]

if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

para_set            = configs["obspara_set_cfg"]["para_set"]
para_min_set        = configs["obspara_set_cfg"]["para_min_set"]
para_max_set        = configs["obspara_set_cfg"]["para_max_set"]
para_mean_set       = configs["obspara_set_cfg"]["para_mean_set"]
para_std_set        = configs["obspara_set_cfg"]["para_std_set"]
para_dist_set       = configs["obspara_set_cfg"]["para_dist_set"]
nens                = configs["da_cfg"]["nens"]

if not isinstance(para_set, list):
    para_set      = [para_set]
    para_min_set  = [para_min_set]
    para_max_set  = [para_max_set]
    para_mean_set = [para_mean_set]
    para_std_set  = [para_std_set]
    para_dist_set = [para_dist_set]

# Check whether this is the first time to update pflotran input files
first_time_update = True if len(model_time_list) <= 2 else False

###############################
# Update the following in PFLOTRAN.in
# (1) The model running time
# (2) The restart file config
###############################
# Get the current model time and assimilation window from config.nml
assim_window_days    = configs["da_cfg"]["assim_window_days"]
assim_window_seconds = configs["da_cfg"]["assim_window_seconds"]
assim_window = assim_window_seconds + assim_window_days*86400  # seconds

# Read the current PFLOTRAN.in information
if not os.path.isfile(pflotran_in_file):
    raise Exception("The pflotran in file does not exist in the path: %s" % pflotran_in_file)
with open(pflotran_in_file, 'r') as f:
    pflotranin = f.readlines()

# Delete the current file
os.remove(pflotran_in_file)

# Update the model run time based on the assimilation window
with open(pflotran_in_file, 'w') as f:
    for i, s in enumerate(pflotranin):
        if "FINAL_TIME" in s:
            pflotranin[i] = "  FINAL_TIME {} sec".format(assim_window)+"\n"
        if "SUBSURFACE_FLOW" in s and "MODE" in pflotranin[i+1] and first_time_update:
            pflotranin.insert(i+2, "        OPTIONS \n")
            pflotranin.insert(i+3, "            REVERT_PARAMETERS_ON_RESTART \n")
            pflotranin.insert(i+4, "        / \n")
        if "CHECKPOINT" in s and "/" in pflotranin[i+1] and first_time_update:
            pflotranin.insert(i+2, "  RESTART \n")
            pflotranin.insert(i+3, "    FILENAME " + pflotran_checkpoint_file + " \n")
            pflotranin.insert(i+4, "    REALIZATION_DEPENDENT \n")
            pflotranin.insert(i+5, "    RESET_TO_TIME_ZERO \n")
            pflotranin.insert(i+6, "  / \n")

    f.writelines(pflotranin)


###############################
# Get the list of DART posterior files
###############################
dart_posterior_file_set = []
if not os.path.isfile(dart_output_list):
    raise Exception("The DART output list file does not exist in the path: %s" % dart_output_list)
with open(dart_output_list, "r") as f:
    dart_posterior_file_list = f.readlines()

for i in range(len(dart_posterior_file_list)):
    file_name_old = dart_posterior_file_list[i]
    file_name = file_name_old.replace("\n", "")
    dart_posterior_file_set.append(file_name)


###############################
# TODO Modify the parameter to 3D later on
# Update the parameter values in
# parameter_prior.h5 based on DART
# posterior output
###############################
if not os.path.isfile(pflotran_para_file):
    raise Exception("The file does not exist in the path: %s" % pflotran_para_file)
f_para = h5py.File(pflotran_para_file, "r+")

posterior = np.zeros([len(para_set), nens])
for i in range(nens):   # For each ensemble...
    dart_posterior_file = dart_posterior_file_set[i]
    # Open the posterior data
    nc_posterior = Dataset(dart_posterior_file, 'r')

    # Get the ensemble index
    # This is because the order of the file names in dart_posterior_file_list
    # may not follow 1,2,3,...,ens
    ens_num = int(nc_posterior.variables["member"][:])
    ens_ind = ens_num - 1

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
    varn = para_set[j]
    f_para[varn][:] = posterior[j, :]
    # f_para[varn][:] = posterior[j, :]

# Close the parameter_prior.h5 file
f_para.close()
