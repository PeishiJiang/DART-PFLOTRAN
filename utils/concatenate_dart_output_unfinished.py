"""This python file is used for concatenate dart output for the unfinished applications."""

import os
import re
import sys
import f90nml
import subprocess
from math import ceil, log10
import numpy as np

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs = f90nml.read(config_nml)

app_dir                     = configs["main_dir_cfg"]["app_dir"]
dart_prior_file             = configs["file_cfg"]["dart_prior_nc_file"]
dart_posterior_file         = configs["file_cfg"]["dart_posterior_nc_file"]
dart_prior_all_file         = configs["file_cfg"]["dart_prior_nc_all_file"]
dart_posterior_all_file     = configs["file_cfg"]["dart_posterior_nc_all_file"]
dart_prior_all_ens_file     = configs["file_cfg"]["dart_prior_nc_all_ens_file"]
dart_posterior_all_ens_file = configs["file_cfg"]["dart_posterior_nc_all_ens_file"]
keep_each_ens_file          = configs["file_cfg"]["keep_each_ens_file"]
ntimestep                   = int(configs["da_cfg"]["ntimestep"])
model_time_list             = configs["time_cfg"]["model_time_list"]
nens                        = configs["da_cfg"]["nens"]

# iteration_step            = configs["da_cfg"]["enks_mda_iteration_step"]
# total_iterations          = configs["da_cfg"]["enks_mda_total_iterations"]
total_iterations          = configs["da_cfg"]["enks_mda_iteration_step"] - 1
save_immediate_mda_result = configs["file_cfg"]["save_immediate_mda_result"]

dart_inout_dir = os.path.join(app_dir, "dart_inout")

# Convert the model_time_list to a list (model_time_list = 0 in the first model tim)
if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

# ndigit = np.ceil(np.log10(ntimestep), dtype=int)
ndigit_time = int(ceil(log10(ntimestep))) + 1
ndigit_ens = int(ceil(log10(nens))) + 1
# Use the real number of time steps here
ntimestep = len(model_time_list)

# Go to the directory
os.chdir(dart_inout_dir)

###############################
# Concatenate all the time steps per ensemble
###############################
for i in range(nens):
    ens = i + 1

    # When the immediate mda results are not saved...
    if not save_immediate_mda_result:
        # Get the first prior and posterior files
        dart_prior_temp1    = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file)
        dart_prior_temp     = re.sub(r"\[TIME\]", "'" + "." * ndigit_time + "'", dart_prior_temp1)
        dart_prior_temp_rm  = re.sub(r"\[TIME\]", "*", dart_prior_temp1)
        dart_posterior_temp1    = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_file)
        dart_posterior_temp     = re.sub(r"\[TIME\]", "'" + "." * ndigit_time + "'", dart_posterior_temp1)
        dart_posterior_temp_rm  = re.sub(r"\[TIME\]", "*", dart_posterior_temp1)

        dart_prior_temp     = os.path.basename(dart_prior_temp)
        dart_posterior_temp = os.path.basename(dart_posterior_temp)

        dart_prior_all     = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_all_file)
        dart_posterior_all = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_all_file)

        # Remove them first if exists
        subprocess.run("rm {}".format(dart_prior_all), shell=True) 
        subprocess.run("rm {}".format(dart_posterior_all), shell=True)

        # Concatanate the prior data
        # subprocess.run("ls | grep {}".format(dart_prior_temp), shell=True, check=True)
        subprocess.run("ls | grep {} | ncrcat -o {}".format(dart_prior_temp, dart_prior_all), shell=True, check=True)

        # Concatanate the posterior data
        subprocess.run("ls | grep {} | ncrcat -o {}".format(dart_posterior_temp, dart_posterior_all), shell=True, check=True)

        # Remove each individual file if needed
        if not keep_each_ens_file:    
            subprocess.run("rm {}".format(dart_prior_temp_rm), shell=True)
            subprocess.run("rm {}".format(dart_posterior_temp_rm), shell=True)
    
    # When the immediate mda results are saved...
    else:
        for i in range(total_iterations):
            iteration_step = i + 1

            # Get the first prior and posterior files
            dart_prior_temp1    = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file)
            dart_prior_temp1    = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_prior_temp1)
            dart_prior_temp     = re.sub(r"\[TIME\]", "'" + "." * ndigit_time + "'", dart_prior_temp1)
            dart_prior_temp_rm  = re.sub(r"\[TIME\]", "*", dart_prior_temp1)
            dart_posterior_temp1    = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_file)
            dart_posterior_temp1    = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_posterior_temp1)
            dart_posterior_temp     = re.sub(r"\[TIME\]", "'" + "." * ndigit_time + "'", dart_posterior_temp1)
            dart_posterior_temp_rm  = re.sub(r"\[TIME\]", "*", dart_posterior_temp1)

            dart_prior_temp     = os.path.basename(dart_prior_temp)
            dart_posterior_temp = os.path.basename(dart_posterior_temp)

            dart_prior_all     = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_all_file)
            dart_prior_all     = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_prior_all)
            dart_posterior_all = re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_all_file)
            dart_posterior_all = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_posterior_all)

            # Remove them first if exists
            subprocess.run("rm {}".format(dart_prior_all), shell=True) 
            subprocess.run("rm {}".format(dart_posterior_all), shell=True)

            # Concatanate the prior data
            # subprocess.run("ls | grep {}".format(dart_prior_temp), shell=True, check=True)
            subprocess.run("ls | grep {} | ncrcat -o {}".format(dart_prior_temp, dart_prior_all), shell=True, check=True)

            # Concatanate the posterior data
            subprocess.run("ls | grep {} | ncrcat -o {}".format(dart_posterior_temp, dart_posterior_all), shell=True, check=True)

            # Remove each individual file if needed
            if not keep_each_ens_file:    
                subprocess.run("rm {}".format(dart_prior_temp_rm), shell=True)
                subprocess.run("rm {}".format(dart_posterior_temp_rm), shell=True)


###############################
# Concatenate all the ensemble
###############################
# When the immediate mda results are not saved...
if not save_immediate_mda_result:
    # Remove them first if exists
    subprocess.run("rm {}".format(dart_prior_all_ens_file), shell=True) 
    subprocess.run("rm {}".format(dart_posterior_all_ens_file), shell=True)

    dart_prior_all_temp     = re.sub(r"\[ENS\]", "'" + "." * ndigit_ens + "'", dart_prior_all_file)
    dart_posterior_all_temp = re.sub(r"\[ENS\]", "'" + "." * ndigit_ens + "'", dart_posterior_all_file)

    dart_prior_all_temp     = os.path.basename(dart_prior_all_temp)
    dart_posterior_all_temp = os.path.basename(dart_posterior_all_temp)

    subprocess.run("ls | grep {} | ncecat -o {}".format(dart_prior_all_temp, dart_prior_all_ens_file), shell=True, check=True)
    subprocess.run("ls | grep {} | ncecat -o {}".format(dart_posterior_all_temp, dart_posterior_all_ens_file), shell=True, check=True)
    subprocess.run("ncrename -d record,ensemble {}".format(dart_prior_all_ens_file), shell=True, check=False)
    subprocess.run("ncrename -d record,ensemble {}".format(dart_posterior_all_ens_file), shell=True, check=False)

    # Remove each individual file if needed
    dart_prior_all_temp_rm     = re.sub(r"\[ENS\]", "*", dart_prior_all_file)
    dart_posterior_all_temp_rm = re.sub(r"\[ENS\]", "*", dart_posterior_all_file)
    if not keep_each_ens_file:    
        subprocess.run("rm {}".format(dart_prior_all_temp_rm), shell=True)
        subprocess.run("rm {}".format(dart_posterior_all_temp_rm), shell=True)

# When the immediate mda results are saved...
else:
    for i in range(total_iterations):
        iteration_step = i + 1

        dart_prior_all_ens_iter     = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_prior_all_ens_file)
        dart_posterior_all_ens_iter = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_posterior_all_ens_file)

        # Remove them first if exists
        subprocess.run("rm {}".format(dart_prior_all_ens_iter), shell=True) 
        subprocess.run("rm {}".format(dart_posterior_all_ens_iter), shell=True)

        dart_prior_all_temp1     = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_prior_all_file)
        dart_posterior_all_temp1 = re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_posterior_all_file)

        dart_prior_all_temp     = re.sub(r"\[ENS\]", "'" + "." * ndigit_ens + "'", dart_prior_all_temp1)
        dart_posterior_all_temp = re.sub(r"\[ENS\]", "'" + "." * ndigit_ens + "'", dart_posterior_all_temp1)

        dart_prior_all_temp     = os.path.basename(dart_prior_all_temp)
        dart_posterior_all_temp = os.path.basename(dart_posterior_all_temp)

        subprocess.run("ls | grep {} | ncecat -o {}".format(dart_prior_all_temp, dart_prior_all_ens_iter), shell=True, check=True)
        subprocess.run("ls | grep {} | ncecat -o {}".format(dart_posterior_all_temp, dart_posterior_all_ens_iter), shell=True, check=True)
        subprocess.run("ncrename -d record,ensemble {}".format(dart_prior_all_ens_iter), shell=True, check=False)
        subprocess.run("ncrename -d record,ensemble {}".format(dart_posterior_all_ens_iter), shell=True, check=False)

        # Remove each individual file if needed
        dart_prior_all_temp_rm     = re.sub(r"\[ENS\]", "*", dart_prior_all_temp1)
        dart_posterior_all_temp_rm = re.sub(r"\[ENS\]", "*", dart_posterior_all_temp1)
        if not keep_each_ens_file:    
            subprocess.run("rm {}".format(dart_prior_all_temp_rm), shell=True)
            subprocess.run("rm {}".format(dart_posterior_all_temp_rm), shell=True)


###############################
# Concatenate all the iterations if needed
###############################
if save_immediate_mda_result:
    dart_prior_iter_temp     = re.sub("[.]nc", "_iter'.'.nc", dart_prior_all_ens_file)
    dart_posterior_iter_temp = re.sub("[.]nc", "_iter'.'.nc", dart_posterior_all_ens_file)

    dart_prior_iter_temp     = os.path.basename(dart_prior_iter_temp)
    dart_posterior_iter_temp = os.path.basename(dart_posterior_iter_temp)

    # Concatanate the prior data
    if os.path.isfile(dart_prior_all_ens_file):
        subprocess.run("rm {}".format(dart_prior_all_ens_file), shell=True)
    subprocess.run("ls | grep {} | ncecat -o {}".format(dart_prior_iter_temp, dart_prior_all_ens_file), shell=True, check=True)
    subprocess.run("ncrename -d record,iteration {}".format(dart_prior_all_ens_file), shell=True, check=False)

    # Concatanate the posterior data
    if os.path.isfile(dart_posterior_all_ens_file):
        subprocess.run("rm {}".format(dart_posterior_all_ens_file), shell=True)
    subprocess.run("ls | grep {} | ncecat -o {}".format(dart_posterior_iter_temp, dart_posterior_all_ens_file), shell=True, check=True)
    subprocess.run("ncrename -d record,iteration {}".format(dart_posterior_all_ens_file), shell=True, check=False)

    # Remove each individual file if needed
    dart_prior_iter_temp_rm     = re.sub("[.]nc", "_iter*.nc", dart_prior_all_ens_file)
    dart_posterior_iter_temp_rm = re.sub("[.]nc", "_iter*.nc", dart_posterior_all_ens_file)
    if not keep_each_ens_file:    
        subprocess.run("rm {}".format(dart_prior_iter_temp_rm), shell=True)
        subprocess.run("rm {}".format(dart_posterior_iter_temp_rm), shell=True)