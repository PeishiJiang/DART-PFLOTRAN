"""This file is used for running PFLOTRAN."""

# Author: Peishi Jiang

import re
import os
import sys
import glob
import f90nml
import subprocess
import numpy as np


###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_exe          = configs["exe_cfg"]["pflotran_exe"]
mpirun                = configs["exe_cfg"]["mpi_exe_pf"]
ncore_pf              = configs["exe_cfg"]["ncore_pf"]
pflotran_in_file      = configs["file_cfg"]["pflotran_in_file"]
pflotran_in_dir       = configs["other_dir_cfg"]["pflotran_in_dir"]
pflotran_out_dir      = configs["other_dir_cfg"]["pflotran_out_dir"]
pflotran_out_prefix   = configs["file_cfg"]["pflotran_out_prefix"]
pflotran_obs_file     = configs["file_cfg"]["pflotran_obs_file"]
pflotran_out_file     = configs["file_cfg"]["pflotran_out_file"]
pflotran_log_file     = configs["file_cfg"]["pflotran_log_file"]
pflotran_restart_file = configs["file_cfg"]["pflotran_restart_file"]
nreaz                 = configs["da_cfg"]["nens"]
ngroup                = configs["exe_cfg"]["ngroup_pf"]
iteration_step        = configs["da_cfg"]["enks_mda_iteration_step"]
total_iterations      = configs["da_cfg"]["enks_mda_total_iterations"]
is_spinup_done        = configs["time_cfg"]["is_spinup_done"]
model_time            = float(configs["time_cfg"]["current_model_time"])  # days

update_obs_ens_posterior_required = configs["da_cfg"]["obs_ens_posterior_from_model"]
try:
    update_obs_ens_posterior_now = configs["da_cfg"]["update_obs_ens_posterior_now"]
except:
    update_obs_ens_posterior_now = False

# model_time_list = configs["time_cfg"]["model_time_list"]
# print(model_time_list)
# if isinstance(model_time_list, list):
#     print(ncore_pf, model_time_list)
#     raise Exception("Stop!")


###############################
# Move to PFLOTRAN input folder
###############################
os.chdir(pflotran_in_dir)


###############################
# Run forward simulation
##############################
# print("{} -n {} {} -pflotranin {} -output_prefix {} -stochastic -num_realizations {} -num_groups {}".format(
#         mpirun, ncore_pf, pflotran_exe, pflotran_in_file, pflotran_out_prefix, nreaz, ngroup
#     ))
print("PFLOTRAN simulation starts ...")
if ncore_pf > 1:
    # subprocess.run("{} -n {} {} -pflotranin {} -output_prefix {} -stochastic -num_realizations {} -num_groups {} -screen_output off".format(
    subprocess.run("{} -n {} {} -pflotranin {} -output_prefix {} -stochastic -num_realizations {} -num_groups {}".format(
        mpirun, ncore_pf, pflotran_exe, pflotran_in_file, pflotran_out_prefix, nreaz, ngroup
    ), shell=True, check=True)
else:
    subprocess.run("{} -pflotranin {} -stochastic -output_prefix {} -num_realizations {} -num_groups {} -screen_output off".format(
        pflotran_exe, pflotran_in_file, pflotran_out_prefix, nreaz, ngroup
    ), shell=True, check=True)

print("PFLOTRAN simulation finishes ...")


###############################
# Rename different PFLOTRAN files
###############################
# PFLOTRAN output file
pflotran_out_file = re.sub(r"\[ENS\]", "*", pflotran_out_file)
pflotran_out_file = os.path.basename(pflotran_out_file)

# PFLOTRAN observation file
pflotran_obs_file = re.sub(r"\[ENS\]", "*", pflotran_obs_file)
pflotran_obs_file = re.sub(r"\[ANY\]", "*", pflotran_obs_file)
pflotran_obs_file = os.path.basename(pflotran_obs_file)

# PFLOTRAN restart file
pflotran_restart_file = re.sub(r"\[ENS\]", "*", pflotran_restart_file)
pflotran_restart_file = os.path.basename(pflotran_restart_file)

# PFLOTRAN log file
pflotran_log_file = re.sub(r"\[ENS\]", "*", pflotran_log_file)
pflotran_log_file = os.path.basename(pflotran_log_file)


###############################
# Move the PFLOTRAN outputs files to output folder
# - Move all of them if it is the first iteration given a time step 
# - Move only .h5 and .out and remove the new .chk if else 
###############################
# If the observation posterior needs to be updated from rerunning the model, increment the total_iteration by one.
if update_obs_ens_posterior_required:
    total_iterations = total_iterations + 1

# print(pflotran_out_file)
# print(pflotran_restart_file)
# print(pflotran_obs_file)
# print(is_spinup_done, update_obs_ens_posterior_required, update_obs_ens_posterior_now, iteration_step, total_iterations)
if is_spinup_done:
    # If iteration is done and model state is not required to be updated
    if not update_obs_ens_posterior_required and iteration_step == total_iterations:
        # copy the latest output from either snapshot file from pflotran_in to pflotran_out
        subprocess.run("cd {0}; cp {1} {2}".format(pflotran_in_dir, pflotran_out_file, pflotran_out_dir), shell=True, check=True)
        # move the observation file from pflotran_in to pflotran_out if exits
        if len(glob.glob(pflotran_obs_file)) != 0:
            subprocess.run("cd {0}; mv {1} {2}".format(pflotran_in_dir, pflotran_obs_file, pflotran_out_dir), shell=True, check=True)
        # move the latest checkout/restart file from pflotran_in to pflotran_out
        subprocess.run("cd {0}; mv {1} {3}; mv {2} {3}".format(pflotran_in_dir, pflotran_restart_file, pflotran_log_file, pflotran_out_dir), shell=True, check=True)
        # subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(
        #     pflotran_in_dir, pflotran_out_dir), shell=true, check=true)

    # If iteration is done and mode state needs update (note that iteration_step would be total_iteration + 1)
    elif update_obs_ens_posterior_required and iteration_step == total_iterations:
        # when model state has been updated:
        if update_obs_ens_posterior_now:
            # print("Cleaning after updating the model states from rerunning the model using posterior parameters ...")
            # copy the latest output from either snapshot file from pflotran_in to pflotran_out
            subprocess.run("cd {0}; cp {1} {2}".format(pflotran_in_dir, pflotran_out_file, pflotran_out_dir), shell=True, check=True)
            # move the observation file from pflotran_in to pflotran_out if exits
            if len(glob.glob(pflotran_obs_file)) != 0:
                subprocess.run("cd {0}; mv {1} {2}".format(pflotran_in_dir, pflotran_obs_file, pflotran_out_dir), shell=True, check=True)
            # move the latest checkout/restart file from pflotran_in to pflotran_out
            subprocess.run("cd {0}; mv {1} {3}; mv {2} {3}".format(pflotran_in_dir, pflotran_restart_file, pflotran_log_file, pflotran_out_dir), shell=True, check=True)
            # subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(
            #     pflotran_in_dir, pflotran_out_dir), shell=True, check=True)
        else:
            raise Exception("Model can't be run happen when iteration is done and update_obs_ens_posterior_now is false!")

    # If iteration is not done
    elif (not update_obs_ens_posterior_required and iteration_step < total_iterations) or \
         (update_obs_ens_posterior_required and iteration_step < total_iterations):
        # It's important to remove the unused restart file first!
        # remove the latest checkout/restart file
        subprocess.run("cd {0}; rm {1}; mv {2} {3}".format(pflotran_in_dir, pflotran_restart_file, pflotran_log_file, pflotran_out_dir), shell=True, check=False)
        # copy the latest output from either snapshot file from pflotran_in to pflotran_out
        subprocess.run("cd {0}; cp {1} {2}".format(pflotran_in_dir, pflotran_out_file, pflotran_out_dir), shell=True, check=False)
        # move the observation file from pflotran_in to pflotran_out if exits
        if len(glob.glob(pflotran_obs_file)) != 0:
            subprocess.run("cd {0}; mv {1} {2}".format(pflotran_in_dir, pflotran_obs_file, pflotran_out_dir), shell=True, check=False)
        # subprocess.run("cd {0}; rm pflotran*.chk; cp pflotran*.h5 {1}; mv pflotran*.out {1}".format(
        #     pflotran_in_dir, pflotran_out_dir), shell=True, check=True)

    else:
        print(iteration_step, update_obs_ens_posterior_required)
        raise Exception("The following situation shouldn't happen -- iteration step: {}; update_obs_ens_posterior_required: {}".format(iteration_step, update_obs_ens_posterior_required))

else:
    # copy the latest output from either snapshot file from pflotran_in to pflotran_out
    subprocess.run("cd {0}; cp {1} {2}".format(pflotran_in_dir, pflotran_out_file, pflotran_out_dir), shell=True, check=True)
    # move the observation file from pflotran_in to pflotran_out if exits
    if len(glob.glob(pflotran_obs_file)) != 0:
        subprocess.run("cd {0}; mv {1} {2}".format(pflotran_in_dir, pflotran_obs_file, pflotran_out_dir), shell=True, check=True)
    # move the latest checkout/restart file from pflotran_in to pflotran_out
    subprocess.run("cd {0}; mv {1} {3}; mv {2} {3}".format(pflotran_in_dir, pflotran_restart_file, pflotran_log_file, pflotran_out_dir), shell=True, check=True)
    # subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(
    #     pflotran_in_dir, pflotran_out_dir), shell=True, check=True)