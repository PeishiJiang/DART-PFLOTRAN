"""This file is used for running PFLOTRAN."""

# Author: Peishi Jiang

import os
import sys
import f90nml
import subprocess
import numpy as np


###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_exe     = configs["exe_cfg"]["pflotran_exe"]
mpirun           = configs["exe_cfg"]["mpi_exe_pf"]
ncore_pf         = configs["exe_cfg"]["ncore_pf"]
pflotran_in_file = configs["file_cfg"]["pflotran_in_file"]
pflotran_in_dir  = configs["other_dir_cfg"]["pflotran_in_dir"]
pflotran_out_dir = configs["other_dir_cfg"]["pflotran_out_dir"]
nreaz            = configs["da_cfg"]["nens"]
ngroup           = configs["exe_cfg"]["ngroup_pf"]
iteration_step   = configs["da_cfg"]["enks_mda_iteration_step"]
total_iterations = configs["da_cfg"]["enks_mda_total_iterations"]
is_spinup_done   = configs["time_cfg"]["is_spinup_done"]
model_time       = float(configs["time_cfg"]["current_model_time"])  # days

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
# Run forward simulation
###############################
if ncore_pf > 1:
    subprocess.run("{} -n {} {} -pflotranin {} -stochastic -num_realizations {} -num_groups {} -screen_output off".format(
        mpirun, ncore_pf, pflotran_exe, pflotran_in_file, nreaz, ngroup
    ), shell=True, check=True)
else:
    subprocess.run("{} -pflotranin {} -stochastic -num_realizations {} -num_groups {} -screen_output off".format(
        pflotran_exe, pflotran_in_file, nreaz, ngroup
    ), shell=True, check=True)


###############################
# Remove the original PFLOTRAN outputs files in output folder
# - Remove all of them if it is the first iteration given a time step 
# - Remove only .h5 and .out if else 
###############################
# if is_spinup_done:
#     # if iteration_step == total_iterations + 1:  ## the final iteration
#     if iteration_step == total_iterations:  ## the final iteration
#         # subprocess.run("cd {}; rm pflotran*.h5; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#         # subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#         if update_obs_ens_posterior_required and not update_obs_ens_posterior_now:
#             subprocess.run("cd {}; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#             # subprocess.run("cd {}; rm pflotran*.out; rm pflotran*.h5".format(pflotran_out_dir), shell=True, check=True)
#         elif update_obs_ens_posterior_required and update_obs_ens_posterior_now:
#             # print("Come here")
#             subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#             # subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out; rm pflotran*.h5".format(pflotran_out_dir), shell=True, check=True)
#         else:
#             # subprocess.run("cd {}; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#             subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#             # subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out; rm pflotran*.h5".format(pflotran_out_dir), shell=True, check=True)

#     else:  ## the other iterations
#         subprocess.run("cd {}; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=True)
#         # subprocess.run("cd {}; rm pflotran*.out; rm pflotran*.h5".format(pflotran_out_dir), shell=True, check=True)

# else:  
#     subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=False)


###############################
# Move the PFLOTRAN outputs files to output folder
# - Move all of them if it is the first iteration given a time step 
# - Move only .h5 and .out and remove the new .chk if else 
###############################
# If the observation posterior needs to be updated from rerunning the model, increment the total_iteration by one.
if update_obs_ens_posterior_required:
    total_iterations = total_iterations + 1

if is_spinup_done:
    # If iteration is done and model state is not required to be updated
    if not update_obs_ens_posterior_required and iteration_step == total_iterations:
        # copy the latest .h5 from pflotran_in to pflotran_out
        # move the latest checkout point file from pflotran_in to pflotran_out
        subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(pflotran_in_dir, pflotran_out_dir), shell=True, check=True)

    # If iteration is done and mode state needs update (note that iteration_step would be total_iteration + 1)
    elif update_obs_ens_posterior_required and iteration_step == total_iterations:
        if update_obs_ens_posterior_now:
        # when model state has been updated:
        # move the latest checkout point file in pflotran_in to pflotran_out
        # copy the latest .h5 from pflotran_in to pflotran_out
            subprocess.run("cd {0}; touch hi_{1}.txt".format(pflotran_out_dir, model_time), shell=True, check=True)
            subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(pflotran_in_dir, pflotran_out_dir), shell=True, check=True)
        else:
            raise Exception("Model can't be run happen when iteration is done and update_obs_ens_posterior_now is false!")

    # If iteration is not done
    # remove the latest checkout point file in pflotran_in to pflotran_out
    # copy the latest .h5 from pflotran_in to pflotran_out
    elif (not update_obs_ens_posterior_required and iteration_step < total_iterations) or \
         (update_obs_ens_posterior_required and iteration_step < total_iterations):
        subprocess.run("cd {0}; rm pflotran*.chk; cp pflotran*.h5 {1}; mv pflotran*.out {1}".format(pflotran_in_dir, pflotran_out_dir), shell=True, check=True)

    else:
        print(iteration_step, update_obs_ens_posterior_required)
        raise Exception("The following situation shouldn't happen -- iteration step: {}; update_obs_ens_posterior_required: {}".format(iteration_step, update_obs_ens_posterior_required))

else:
    # subprocess.run("cd {}; rm pflotran*.chk; rm pflotran*.out".format(pflotran_out_dir), shell=True, check=False)
    subprocess.run("cd {0}; cp pflotran*.h5 {1}; mv pflotran*.chk {1}; mv pflotran*.out {1}".format(pflotran_in_dir, pflotran_out_dir), shell=True, check=True)