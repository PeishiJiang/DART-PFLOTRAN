"""This file is used for updating the observation ensemble posterior from rerunning the model with parameter posterior as input."""

# Author: Peishi Jiang

import os
import re
import sys
import h5py
import shutil
import f90nml
import subprocess
import numpy as np
from math import ceil, log10
from netCDF4 import num2date, date2num, Dataset

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

update_pflotran_input_py = configs["file_cfg"]["update_pflotran_input_file"]
run_pflotran_py          = configs["file_cfg"]["run_pflotran_file"]

pflotran_out_file        = configs["file_cfg"]["pflotran_out_file"]
pflotran_para_file       = configs["file_cfg"]["pflotran_para_file"]
dart_output_list         = configs["file_cfg"]["dart_output_list_file"]
obs_set                  = configs["obspara_set_cfg"]["obs_set"]
para_set                 = configs["obspara_set_cfg"]["para_set"]
model_time               = float(configs["time_cfg"]["current_model_time"])  # days
model_time_list          = configs["time_cfg"]["model_time_list"]
nens                     = configs["da_cfg"]["nens"]
required                 = configs["da_cfg"]["obs_ens_posterior_from_model"]

# Get the list of all required PFLOTRAN variables
if not isinstance(obs_set, list):
    obs_set = [obs_set]

###############################
# Check whether the observation ens posterior should be obtained from rerunning the model
###############################
if not required:
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
# Check whether the PFLOTRAN output HDF and DART posterior files exit.
###############################
# PFLOTRAN output
ens_set        = np.arange(1, nens + 1)  # the set of ensembles
pflotran_out_file_set = [
    re.sub(r"\[ENS\]", str(ens), pflotran_out_file) for ens in ens_set
]
for f in pflotran_out_file_set:
    if not os.path.isfile(f):
        raise Exception("The PFLOTRAN output file %s does not exits!" % f)

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
# Assign the model output to observation ensemble posterior
###############################
for i in range(nens):
    dart_posterior_file = dart_posterior_file_set[i]
    pf_fname            = pflotran_out_file_set[i]

    dart_var_dict = dict.fromkeys(obs_set)

    # Read the PFLOTRAN output
    f_out = h5py.File(pf_fname, "r")
    # Get the last time step
    time_set_o = [t for t in list(f_out.keys()) if t.startswith("Time")]
    time_set   = [t.split()[1:] for t in time_set_o]
    time_vset  = [float(t[0]) for t in time_set]

    last_time, last_time_ind = np.max(time_vset), np.argmax(time_vset)
    time_unit, last_time_o   = time_set[last_time_ind][1], time_set_o[last_time_ind]

    # Get the state/parameter/variable values required in pflotran_var_set
    dataset          = f_out[last_time_o]
    pl_out_var_set_o = list(dataset.keys())
    pl_out_var_dict  = dict()
    for v in pl_out_var_set_o:
        # Get the variable name and unit from the original variable name
        varinfo = v.split()
        if len(varinfo) == 1:
            varn, varunit = varinfo[0].upper(), ''
        elif len(varinfo) == 2:
            varn, varunit = varinfo[0].upper(), varinfo[1]
        else:
            raise Exception('Invalid variable name %s!' % v)
        pl_out_var_dict[varn] = {"unit": varunit, "original_name": v}
        # Check if the variable is required by pflotran_var_set
        if varn in obs_set:
            dart_var_dict[varn] = dataset[v][:]
    f_out.close()


    # Open the posterior data
    nc_posterior = Dataset(dart_posterior_file, 'r+')

    # Get the ensemble index
    # This is because the order of the file names in dart_posterior_file_list
    # may not follow 1,2,3,...,ens
    ens_num = int(nc_posterior.variables["member"][:])
    ens_ind = ens_num - 1

    if ens_ind != i:
        raise Exception("The ensemble member does not match! Check your DART posterior file list.")

    for j in range(len(obs_set)):  # For each observation variable ...
        varn = obs_set[j]
        nc_posterior.variables[varn][:] = dart_var_dict[varn]

    # Close the posterior NetCDF file
    nc_posterior.close()


# Once it is done, switch the update back to False
configs["da_cfg"]["update_obs_ens_posterior_now"] = False
configs.write(config_nml, force=True)