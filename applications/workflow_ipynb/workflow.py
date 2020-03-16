# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
#   # Overview
#   The **objective** of this notebook is to present the workflow of conducting data assimilation on [PFLOTRAN](https://www.pflotran.org/) by using [DART](https://www.image.ucar.edu/DAReS/DART/). Briefly, the procedures are as follows:
#   - [x] [Configuration](#parameter): define directories, file locations, and other parameters
#   - [x] [PFLOTRAN preparation](#pflotran_prepare): generate PFLOTRAN input files
#   - [x] [PFLOTRAN model spin-up](#pflotran_spinup): conduct model spin-up
#   - [x] [DART files preparation](#dart_prepare): add new DART quantities, prepare DART input namelists, prepare DART prior data, prepare observations in DART format, and check ```model_mod``` interface
#   - [x] [Generate all the executable files](#dart_executables): generate all the executables, convert observations in DART format, check ```model_mod``` interface, and test the data assimilation engine
#   - [x] [Run DART and PFLOTRAN](#run_dart_pflotran): run the shell script for integrating DART and PFLOTRAN model
# 
#   Here, we perform inverse modeling on a 1D thermal model for illustration. The model assimilates temperature observation to update its parameters (i.e., flow flux).
# %% [markdown]
#   <a id='parameter'></a>
#   # Step 1: Configuration

# %%
import os
import re
import sys
import time
import shutil
import pickle
import f90nml
import subprocess
import numpy as np
from math import floor, ceil
from datetime import datetime, timedelta
# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')

# %% [markdown]
#   ****************
#   **Define the locations of MPI, PFLOTRAN, application folder and DART-PFLOTRAN interface folder**
# 
#   It is suggested that <span style="background-color:lightgreen">mpi_exe</span> is defined based on the mpi utility (e.g., mpirun) installed by PFLOTRAN.
# 
#   **Important:** You must make sure the <span style="background-color:lightgreen">mpi_exe</span> has the same MPI system as the settings ```MPIFC``` and ```MPILD``` in ```mkmf.template```.

# %%
# MPI settings
# mpi_exe_da  = '/usr/local/bin/mpirun'  # The location of mpirun
mpi_exe_da  = '/Users/jian449/Codes/petsc/arch-darwin-c-opt/bin/mpirun'
mpi_exe_pf  = '/Users/jian449/Codes/petsc/arch-darwin-c-opt/bin/mpirun'
# mpi_exe_da  = '/software/petsc_v3.11.3/arch-linux2-c-opt/bin/mpirun'  # The location of mpirun
# mpi_exe_pf  = '/software/petsc_v3.11.3/arch-linux2-c-opt/bin/mpirun'
ncore_da = 4  # The number of MPI cores used by DART
ncore_pf = 4  # The number of MPI cores used by PFLOTRAN
ngroup_pf= 4  # The number of group used by stochastic running in PFLOTRAN

# PFLOTRAN executable
# pflotran_exe  = '/global/project/projectdirs/m1800/pin/pflotran-haswell/src/pflotran/pflotran'
pflotran_exe  = '/Users/jian449/Codes/pflotran/src/pflotran/pflotran'
# pflotran_exe  = '/software/pflotran/src/pflotran/pflotran'

# Main directory names
temp_app_dir = "/Users/jian449/Codes/DART/manhattan/models/pflotran/applications/1dthermal"          # The template for application folder
app_dir      = "/Users/jian449/Codes/DART/manhattan/models/pflotran/applications/1dthermal_parastate"          # The application folder name
dart_dir     = "/Users/jian449/Codes/DART/manhattan"
dart_pf_dir  = "/Users/jian449/Codes/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
# temp_app_dir = "/home/jian449/DART/manhattan/models/pflotran/applications/template"          # The template for application folder
# app_dir      = "/home/jian449/DART/manhattan/models/pflotran/applications/1dthermal_test_1month_4mda_v2"          # The application folder name
# dart_dir     = "/home/jian449/DART/manhattan"
# dart_pf_dir  = "/home/jian449/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
# temp_app_dir = "/global/cscratch1/sd/peishi89/DART_PFLOTRAN_APP/applications/template"          # The template for application folder
# app_dir      = "/global/cscratch1/sd/peishi89/DART_PFLOTRAN_APP/applications/1dthermal"          # The application folder name
# dart_dir     = "/global/homes/p/peishi89/DART/manhattan"
# dart_pf_dir  = "/global/homes/p/peishi89/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
#temp_app_dir = os.path.abspath("../template" )          # The template for application folder
#app_dir      = os.path.abspath("../1dthermal/")          # The application folder name
#dart_dir     = os.path.abspath("../../../../")
#dart_pf_dir  = os.path.join(dart_dir, "models/pflotran")     # The dart pflotran utitlity folder name

# configs = {}
configs = f90nml.namelist.Namelist()
configs["main_dir_cfg"] = {"app_dir": app_dir, "dart_dir": dart_dir, "dart_pf_dir": dart_pf_dir}
configs["exe_cfg"]      = {"pflotran_exe": pflotran_exe, "mpi_exe_da": mpi_exe_da, "mpi_exe_pf": mpi_exe_pf,
                           "ncore_pf": ncore_pf, "ncore_da": ncore_da, "ngroup_pf": ngroup_pf}

# %% [markdown]
#   ****************
#   **Generate the application directory if it does not exit**

# %%
# Create the application directory if it does not exists
if not os.path.isdir(app_dir):
    shutil.copytree(temp_app_dir, app_dir)
    

# %% [markdown]
#   ****************
#   **Load all the required file paths/names from ```file_paths.nml```**
# 
#   ```file_paths.nml``` defines the relative locations of all the files (e.g., utility files, shell scripts, data files, etc) used by this application

# %%
sys.path.append(dart_pf_dir)
from utils.read_filepaths_nml import read_filepaths_nml


# %%
dirs_cfg, files_cfg      = read_filepaths_nml(app_dir=app_dir, dart_pf_dir=dart_pf_dir)
files_cfg["keep_each_ens_file"] = False   # Whether each ensemble nc file is kept when the DA ends
configs["other_dir_cfg"] = dirs_cfg
configs["file_cfg" ]     = files_cfg
config_file              = files_cfg["config_file"]

# %% [markdown]
#   ****************
#   **Clean up the folders under the app_dir**

# %%
pflotran_in_dir   = dirs_cfg["pflotran_in_dir"]
pflotran_out_dir  = dirs_cfg["pflotran_out_dir"]
dart_inout_dir    = dirs_cfg["dart_data_dir"]
# subprocess.run("cd {}; rm *".format(pflotran_in_dir), shell=True, check=False)
subprocess.run("cd {}; rm -f *".format(pflotran_out_dir), shell=True, check=True)
subprocess.run("cd {}; ls | xargs rm".format(dart_inout_dir), shell=True, check=False)
# subprocess.run("cd {}; ls | xargs -r rm".format(dart_inout_dir), shell=True, check=True)

# %% [markdown]
#   ****************
#   **Specify the following types of variables used in DART**
#   - the observation data to be assimilated
#   - the PFLOTRAN parameters to be analyzed
#   - the statistics of the parameters: (1) the value range; (2) the mean and std for randomly sampling; (3) the distribution to be sampled (i.e., normal, log normal, and uniform distributions)
#   - Whether the prior at one time window is rescaled from the posterior from the previous window
#   - The list of parameters whose prior would be resampled based on the mean of the corresponding posterior at the previous time step.

# %%
# Observation data to be assimilated
obs_set  = ['TEMPERATURE']

# The PFLOTRAN parameters to be analyzed
# para_set = ['FLOW_FLUX','POROSITY','THERMAL_CONDUCTIVITY']
para_set = ['FLOW_FLUX']
para_homogeneous = True

# The statistics of the parameters
# The index order follows para_set
# para_min_set  = [-10.0, 0.01, 0.5]  # The minimum values (-99999 means no lower bound limit)
# para_max_set  = [10.0, 0.7, 2.5]  # The maximum values (99999 means no upper bound limit)
# para_mean_set = [0.0, 0.3, 1.5]  # The mean values
# para_std_set  = [0.5, 0.1, 0.5]  # The standard deviation values
# para_dist_set = ["normal", "normal", "normal"]  # The assumed distribution to be sampled
para_min_set  = [-5.0]  # The minimum values (-99999 means no lower bound limit)
para_max_set  = [5.0]  # The maximum values (99999 means no upper bound limit)
para_mean_set = [0.0]  # The mean values
para_std_set  = [0.5]  # The standard deviation values
para_dist_set = ["normal"]  # The assumed distribution to be sampled
para_prior_rescaled = False  # Whether the prior is generated by rescaling the posterior

para_resampled_set = ['FLOW_FLUX']   # The parameters to be resampled at each time step
# para_resampled_set = ['']   # The parameters to be resampled at each time step

configs["obspara_set_cfg"] = {"obs_set": obs_set, "para_set": para_set, 
                              "para_prior_rescaled": para_prior_rescaled, "para_homogeneous": para_homogeneous,
                              "para_min_set": para_min_set, "para_max_set": para_max_set,
                              "para_mean_set": para_mean_set, "para_std_set": para_std_set,
                              "para_dist_set": para_dist_set, "para_resampled_set": para_resampled_set}

# %% [markdown]
#   ****************
#   **Specify the spatial domains of the observation data to be assimilated**
# 
#   The limits of x, y, and z are bounded by the minimum and maximum boundaries through (min, max). If the limit is not specified, -99999 and 99999 are assigned for the lower and upper bounds, respectively.

# %%
obs_space_xlimit = [-99999, 99999]
obs_space_ylimit = [-99999, 99999]
obs_space_zlimit = [-0.50, -0.03]
configs["obs_space_cfg"] = {"obs_space_xlimit": obs_space_xlimit, "obs_space_ylimit": obs_space_ylimit,
                            "obs_space_zlimit": obs_space_zlimit}

# %% [markdown]
#   ****************
#   **Specify the data assimilation time window**

# %%
da_cfg = {}

# Assimilation time window time_step_days+time_step_seconds
# Assimilation window
da_cfg["assim_window_days"]    = 0     # assimilation time window/step (day)
da_cfg["assim_window_seconds"] = 3600  # assimilation time window/step  (second)
da_cfg["assim_window_size"] = da_cfg["assim_window_days"]+float(da_cfg["assim_window_seconds"])/86400. # day

# %% [markdown]
#   ****************
#   **Specify the temporal information**
#   - model spinup time/start time
#   - the current model time (considered as the middle time of the current assimilation time window)
#   - the list of model time
#   - the list of starting assimilation time (the first starting from the end of spinup)
#   - the map between the observation start time and model start time (considered as the end of spinup)
#   - the first and the last observation times
# 
#   **note that** model start time is considered after the spinup

# %%
time_cfg = {}
one_sec = 1./86400. # fraction of day
# Model spinup length
time_cfg["spinup_length_days"] = 0    # number of days for spinup
time_cfg["spinup_length_seconds"] = 7200    # number of seconds for the remaining spinup
time_cfg["spinup_length"]  = time_cfg["spinup_length_days"]+float(time_cfg["spinup_length_seconds"])/86400. # spinup time (day)
time_cfg["is_spinup_done"] = False  # whether spinup is conducted

# Model start time
time_cfg["current_model_time"] = time_cfg["spinup_length"] + (da_cfg["assim_window_size"] + one_sec)/2  # model start time should be the middle of the first assimilation window (after spinup)
# time_cfg["model_time_list"]    = [0.]   # the list of model time
time_cfg["model_time_list"]    = [time_cfg["current_model_time"]]   # the list of model time

# Map between observation assimilation start time and model start time
obs_start   = datetime(2017,4,1,0,0,0)
assim_start = obs_start + timedelta(days=time_cfg["spinup_length"]) # assimilation time should be after the model spinup
time_cfg["model_start"] = obs_start.strftime("%Y-%m-%d %H:%M:%S")
time_cfg["assim_start"] = assim_start.strftime("%Y-%m-%d %H:%M:%S")

# The maximum time for observation
time_cfg["first_obs_time_days"]    = time_cfg["spinup_length_days"]
time_cfg["first_obs_time_seconds"] = time_cfg["spinup_length_seconds"]
time_cfg["first_obs_time_size"] = time_cfg["first_obs_time_days"]+float(time_cfg["first_obs_time_seconds"])/86400. # day
time_cfg["last_obs_time_days"]    = time_cfg["first_obs_time_days"]
# time_cfg["last_obs_time_seconds"] = time_cfg["first_obs_time_seconds"] + 3600*24*29
time_cfg["last_obs_time_seconds"] = time_cfg["first_obs_time_seconds"] + 3600*3
time_cfg["last_obs_time_size"] = time_cfg["last_obs_time_days"]+float(time_cfg["last_obs_time_seconds"])/86400. # day

# Whether the model time exceeds the last observation
time_cfg["exceeds_obs_time"] = time_cfg["current_model_time"] >= time_cfg["last_obs_time_size"]

# Save them to configs
configs["time_cfg"] = time_cfg

# %% [markdown]
#   ****************
#   **Config the data assimilation**
#   - observation error and error type
#   - number of realizations
#   - the start and end time of the current assimilation time window (obtained based on the previous temporal information settings)
#   - decision on whether the observation posteriors are recomputed based on the updated parameters

# %%
# Observation error, number of ensembles
da_cfg["obs_reso"]  = 300.0      # observation resolution (second)
da_cfg["nens"]      = 100        # number of ensembles
da_cfg["obs_error"] = 0.05       # the observation error
da_cfg["obs_error_type"] = "absolute" # the type of observation error (i.e., relative and absolute)

# The start and end time of the current assimilation time window
da_cfg["assim_start_days"]    = int(floor(time_cfg["current_model_time"] - da_cfg["assim_window_size"]/2))
da_cfg["assim_start_seconds"] = int((time_cfg["current_model_time"] - da_cfg["assim_window_size"]/2 - da_cfg["assim_start_days"])*86400+1)
da_cfg["assim_end_days"]      = int(floor(time_cfg["current_model_time"] + da_cfg["assim_window_size"]/2))
da_cfg["assim_end_seconds"]   = int((time_cfg["current_model_time"] + da_cfg["assim_window_size"]/2 - da_cfg["assim_end_days"])*86400)

# Compute the number of time steps
# print(np.ceil((time_cfg["last_obs_time_size"] - time_cfg["first_obs_time_size"]) / da_cfg["assim_window_size"]))
da_cfg["ntimestep"] = ceil((time_cfg["last_obs_time_size"] - time_cfg["first_obs_time_size"]) / da_cfg["assim_window_size"])

# Decide whether the observation posterior should be recomputed from the model output based on the parameter posterior
da_cfg["obs_ens_posterior_from_model"] = True

# %% [markdown]
#   ****************
#   **Config the MDA settings**

# %%
# The inflation settings used in EnKS-MDA (the alpha value)
# da_cfg["enks_mda_alpha"] = [4., 4., 4., 4.]  # Note that the summation of the inverse of alpha should be one
# da_cfg["enks_mda_alpha"] = [3., 3., 3.]  # Note that the summation of the inverse of alpha should be one
da_cfg["enks_mda_alpha"] = [2., 2.]  # Note that the summation of the inverse of alpha should be one
# da_cfg["enks_mda_alpha"] = [1.]  # Note that the summation of the inverse of alpha should be one
da_cfg["enks_mda_iteration_step"] = 1  # the ith iteration (1 for the first iteration)
da_cfg["enks_mda_total_iterations"] = len(da_cfg["enks_mda_alpha"])  # Note that the summation of the inverse of alpha should be one

# Check whether the sum of the inverse of enks_mda_alpha is one
alpha_inv_sum = sum([1./alpha for alpha in da_cfg["enks_mda_alpha"]])
if alpha_inv_sum - 1 > 1e-8:
    raise Exception("The sum of the inverse of enks_mda_alpha should be one!")

# Save them to configs
configs["da_cfg"] = da_cfg

# %% [markdown]
#   ****************
#   **Save all the configurations in pickle**

# %%
configs.write(config_file, force=True)

# %% [markdown]
#   <a id='pflotran_prepare'></a>
#   # Step 2: PFLOTRAN preparation
#   *Here, we use Kewei's 1D thermal model as an example for generating PFLOTRAN input card and parameter.h5.*
# 
#   In this section, the following procedures are performed:
#   - generate PFLOTRAN input deck file ```PFLOTRAN.in```
#   - generate the parameter files in HDF 5, ```parameter_prior.h5```, used by PFLOTRAN input deck file
# 
#   **Note that**
#   - ```PFLOTRAN.in``` for each DA scenario should be prepared by users.

# %%
prep_pflotran_inputdeck, pflotran_in = files_cfg["prep_pflotran_inputdeck_file"], files_cfg["pflotran_in_file"]
prep_pflotran_parameterprior = files_cfg["prep_pflotran_para_file"]

# %% [markdown]
#   ****************
#   **Generate the ensembles of PFLOTRAN prior**
# 
#   **Run code**
#   - Run: ```prepare_pflotran_parameterprior.py```

# %%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_pflotran_parameterprior" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare PFLOTRAN ensemble parameter values...")
subprocess.run("python {} {}".format(prep_pflotran_parameterprior, config_file), shell=True, check=True)

# %% [markdown]
#   ****************
#   **Generate PFLOTRAN input deck file**
# 
#   **Run code**
#   - Run: ```prepare_pflotran_inputdeck.py```

# %%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_pflotran_inputdeck" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare PFLOTRAN input deck...")
subprocess.run("python {} {}".format(prep_pflotran_inputdeck, config_file), shell=True, check=True)

# %% [markdown]
#   <a id='pflotran_spinup'></a>
#   # Step 3: PFLOTRAN model spin-up
#   Take in the ```pflotran.in``` and ```parameter.h5``` files and conduct the model spin-up by running ```pflotran.sh``` file. The ```pflotran.sh``` is a simple shell script executing ensemble simulation of PFLOTRAN by using MPI.
# 
#   **Run the code**
#   - Run: ```run_pflotran.py```
# 

# %%
# pflotran_sh, pflotran_out_dir = files_cfg["pflotran_sh_file"], dirs_cfg["pflotran_out_dir"]
run_pflotran, pflotran_out_dir = files_cfg["run_pflotran_file"], dirs_cfg["pflotran_out_dir"]


# %%
print("\n")
print("------------------------------------------------------------")
print("Model spinup and run the forward simulation in the first assimilation time window...")
# subprocess.run("{} {}".format(pflotran_sh, config_file), shell=True, check=True)
subprocess.run("python {} {}".format(run_pflotran, config_file), shell=True, check=True)

# %% [markdown]
#   ****************
#   **Once the model spinup finishes, modify the corresponding configuration entry**

# %%
configs["time_cfg"]["is_spinup_done"] = True
configs.write(config_file, force=True)

# %% [markdown]
#   <a id='dart_prepare'></a>
#   # Step 4: DART files preparation
#   In this section, the following procedures are performed:
#   - [TO BE REVISED]: generate the template for DART generic variable quantity files (i.e., ```DEFAULT_obs_kind_mod.F90``` and ```obs_def_pflotran_mod.f90```);
#   - generate the DART input namelists;
#   - generate DART prior NetCDF data ```prior_ensemble_[ENS].nc``` from PFLOTRAN's parameter and outputs;
#   - generate DART posterior NetCDF files (*sharing the same variable names and dimensions as the prior NetCDF files but without the data values*);
#   - convert the observation file to DART observation format;
#   - check ```model_mod.F90``` based on current setting by using the ```check_model_mod``` provided by DART.
# %% [markdown]
#   <a id='dart_generic_prepare'></a>
#   ## Generate the templates for DART generic variable quantity files
#   - Run: ```list2dartqty.py``` to sequentially generate/modify
#       - a mapping between PFLOTRAN variales and DART generic quantities in ```obs_def_pflotran_mod.F90```
#       - the default DART generic quantity definition file ```DEFAULT_obs_kind_mod.F90```
# 

# %%
to_dartqty, obs_type_file = files_cfg["to_dartqty_file"], files_cfg["obs_type_file"]


# %%
# get_ipython().run_cell_magic('script', 'bash -s "$to_dartqty" "$config_file"', 'python $1 $2 $3 $4')
print("\n")
print("------------------------------------------------------------")
print("Add PFLOTRAN variables to DEFAULT_obs_kind_mod.F90 if necessary...")
subprocess.run("python {} {}".format(to_dartqty, config_file), shell=True, check=True)

# %% [markdown]
#   ## Generate  DART input namelists in ```input.nml```
# 
#   The ```input.nml``` file is generated based on a template ```input.nml.template``` by modifying the following namelist entries:
# 
#   ```input.nml.template``` $\rightarrow$ ```input.nml```
# 
#   **Namelists from DART**
#   - [obs_kind_nml](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/modules/observations/obs_kind_mod.html#Namelist): namelist for controling what observation types are to be assimilated
#   - [preprocess_nml](https://www.image.ucar.edu/DAReS/Codes/DART/manhattan/assimilation_code/programs/preprocess/preprocess): namelist of the DART-supplied preprocessor program which creates observation kind and observation definition modules from a set of other specially formatted Fortran 90 files
#   - and others (see DART documentation for details.)
# 
#   **Self-defined namelists**
#   - smoother_nml: a self-defined namelist of the main module for driving ensemble data assimilations modified from DART [filter_nml](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/modules/assimilation/filter_mod.html#Namelist)
#   - model_nml: a self-defined namelist for providing the basic information in the model
#       - template_file: the template prior NetCDF file for ```model_mod.F90``` to digest the spatial information of the model
#       - nvar_para: the number of parameter variables
#       - para_var_names: the original variable names for model parameters
#       - para_var_qtynames: the corresponding DART variable quantities for model parameters
#       - nvar_state: the number of state variables
#       - state_var_names: the original variable names for model states
#       - state_var_qtynames: the corresponding DART variable quantities for model states
#       - debug: whether the debug mode is turned on
#       - max_time_diff_seconds: the maximum time difference allowed for finding the closest ensemble for observations (second)
#       
#   - convertnc_nml: a self-defined namelist for providing the NetCDF observation file name and the DART observation file name used in ```convert_nc.f90```
#       - netcdf_file: the location of the NetCDF file containing the observation data
#       - out_file: the location of the DART observation file
# 
#   **Note that**
#   - There are more namelists or other items in the above namelist in input.nml.template. Users can edit the below python dictionary ```inputnml``` to include their modifications.
#   - Users can also include more namelists provided by DART by modifying ```inputnml```.
# %% [markdown]
#   ***************
#   **Assemble all the namelists in input.nml**

# %%
# Parameters for different namelists in input.nml
smoother_nml = {"input_state_file_list":files_cfg["dart_input_list_file"],
              "output_state_file_list":files_cfg["dart_output_list_file"],
              "ens_size":da_cfg["nens"],
              "num_output_state_members":da_cfg["nens"],
              "obs_sequence_in_name":files_cfg["obs_dart_file"],
              "first_obs_days": da_cfg["assim_start_days"],
              "first_obs_seconds": da_cfg["assim_start_seconds"],
              "last_obs_days": da_cfg["assim_end_days"],
              "last_obs_seconds": da_cfg["assim_end_seconds"],
              "output_mean": False,
              "output_sd": False}
obs_kind_nml = {"assimilate_these_obs_types":obs_set}
assim_tools_nml = {"filter_kind":2}
model_nml = {"nvar_state":len(obs_set),
             "state_var_names":obs_set,
             "state_var_qtynames":['QTY_PFLOTRAN_'+v for v in obs_set],
             "nvar_para":len(para_set),
             "para_var_names":para_set,
             "para_var_qtynames":['QTY_PFLOTRAN_'+v for v in para_set],
             "max_time_diff_seconds": 10,
             "debug": False,
             "template_file":files_cfg["dart_prior_template_file"]}
preprocess_nml = {"input_files":files_cfg["obs_type_file"],
                  "input_obs_kind_mod_file":files_cfg["def_obs_kind_file"]}
convertnc_nml = {"netcdf_file": files_cfg["obs_nc_file"],
                 "out_file": files_cfg["obs_dart_file"],
                 "obs_start_day": da_cfg["assim_start_days"],
                 "obs_start_second": da_cfg["assim_start_seconds"],
                 "obs_end_day": da_cfg["assim_end_days"],
                 "obs_end_second":da_cfg["assim_end_seconds"],
                 "inflation_alpha":da_cfg["enks_mda_alpha"][da_cfg["enks_mda_iteration_step"]-1]}
modelmodcheck_nml = {"input_state_files": files_cfg["dart_prior_template_file"]}
inputnml = {"smoother_nml":smoother_nml,
            "obs_kind_nml":obs_kind_nml,
            "assim_tools_nml":assim_tools_nml,
            "model_nml":model_nml,
            "preprocess_nml":preprocess_nml,
            "convert_nc_nml":convertnc_nml,
            "model_mod_check_nml":modelmodcheck_nml}


configs["inputnml_cfg"] = inputnml

# Save the configurations
configs.write(config_file, force=True)
# with open(config_pickle, 'wb') as f:
#     pickle.dump(configs, f)

# %% [markdown]
#   ***************
#   **Run the code**
#   - Run: ```prepare_inputnml.py```

# %%
prep_inputnml = files_cfg["prep_inputnml_file"]


# %%
# get_ipython().run_cell_magic('script', 'bash -s  "$prep_inputnml" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Generate input.nml file...")
subprocess.run("python {} {}".format(prep_inputnml, config_file), shell=True, check=True)

# %% [markdown]
#   <a id='observationconvertion'></a>
#   ## Prepare the observation conversion to DART observation format
#   In this section, we prepare the process of converting the observation data to DART format. We first convert observation data in raw format into NetCDF format. Then, a fortran script is prepared for the conversion from the NetCDF to to DART format. The structure of NetCDF file for recording observation file.
# 
#   | NetCDF dimensions |           NetCDF variables          |
#   |:-----------------:|:-----------------------------------:|
#   | time: ntime       | time: shape(time)                   |
#   | location: nloc    | location: shape(location)           |
#   |                   | physical variable: shape(time,nloc) |
# 
#   **Note that**
#   - if the time calendar follows *gregorian*, the time unit should be entered as ```seconds since YYYY-MM-DD HH:MM:SS```. Otherwise, put the time calender as *None* and time unit as ```second``` (make sure convert your measurement times to seconds).
# %% [markdown]
#   ***************
#   **Convert the raw csv temperature observations to NetCDF file**
#   - Run: ```csv2nc.py```
# 
#   **Note that**
#   - Here, we illustrate the conversion from csv as an example. However, the conversion to netCDF format can be arbitrary based on the raw data format. 
#  

# %%
csv_to_nc, obs_nc_original = files_cfg["csv_to_nc_file"], files_cfg["obs_nc_original_file"]


# %%
# get_ipython().run_cell_magic('script', 'bash -s "$csv_to_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Convert raw observation data to NetCDF file...")
subprocess.run("python {} {}".format(csv_to_nc, config_file), shell=True, check=True)

# %% [markdown]
#   ***************
#   **Clip the NetCDF file based on the defined spatial and temporal domains**
# 
#   The NetCDF file generated in the previous step is further processed by selecting data observed in the required spatial and temporal domains
#   - Run: ```clip_obs_nc.py```

# %%
clip_obs_nc, obs_nc = files_cfg["clip_obs_nc_file"], files_cfg["obs_nc_file"]


# %%
# get_ipython().run_cell_magic('script', 'bash -s "$clip_obs_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Clip the NetCDF file based on the defined spatial and temporal domains...")
subprocess.run("python {} {}".format(clip_obs_nc, config_file), shell=True, check=True)

# %% [markdown]
#   ***************
#   **Prepare the ```convert_nc.f90``` based on the list of observation variables**
#   - Run: ```prepare_convert_nc.py```
#   - Code input arguments:
#       - <span style="background-color:lightgreen">obs_nc</span>: filename for the observation NetCDF file

# %%
prep_convert_nc, convert_nc_file = files_cfg["prep_convert_nc_file"], files_cfg["convert_nc_file"]


# %%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_convert_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare the convert_nc.f90 based on the list of observation variables to be assimilated...")
subprocess.run("python {} {}".format(prep_convert_nc, config_file), shell=True, check=True)

# %% [markdown]
#   <a id='dart_executables'></a>
#   # Step 5: Generate all the executable files
#   Now, we compile all the executables from ```mkmf_*```. The following executables are generated here:
#   - ```preprocess```: for preprocessing the [prepared DART generic variable quantity files prepared](#dart_generic_prepare)
#   - ```convert_nc```: for [converting the observations from NetCDF to DART format](#observationconvertion)
#   - ```model_mod_check```: for checking ```model_mod.F90``` interface file
#   - ```smoother```: for conducting the DART data assimilation
# %% [markdown]
#   ## Generate the executables
#   - Run: ```quickbuild.csh```
#   - Code input arguments:
#       - <span style="background-color:lightgreen">app_work_dir</span>: location of the application work folder

# %%
dart_work_dir, app_work_dir = dirs_cfg["dart_work_dir"], dirs_cfg["app_work_dir"]
quickbuild = files_cfg["quickbuild_csh"]


# %%
print("\n")
print("------------------------------------------------------------")
print("Generate all the executables...")
subprocess.run("cd {}; csh {} {} -mpi".format(dart_work_dir, quickbuild, app_work_dir), shell=True, check=True)
# subprocess.run("cd {}; csh {} {}".format(dart_work_dir, quickbuild, app_work_dir), shell=True, check=True)

# %% [markdown]
#   ## Check ```model_mod.F90``` interface file
#   - Run: ```model_mod_check```

# %%
model_mod_check = files_cfg["model_mod_check_exe"]

# %% [markdown]
#   <a id='run_dart_pflotran'></a>
#   # Step 6: Run DART and PFLOTRAN
#   In this section, run the shell script to couple DART and PFLOTRAN.

# %%
dart_work_dir     = dirs_cfg["dart_work_dir"]
run_DART_PFLOTRAN = files_cfg["run_da_csh"]
concatenate_output = files_cfg["concatenate_dart_output_file"]
inputnml_file     = files_cfg["input_nml_file"]
start_time        = time.time()


# %%
print("\n")
print("------------------------------------------------------------")
print("Assimilation starts here...")
subprocess.run("cd {}; csh {} {} {}".format(dart_work_dir, run_DART_PFLOTRAN, inputnml_file, config_file), shell=True, check=True)


# %%
print("\n")
print("------------------------------------------------------------")
print("Concatenate the prior and posterior at all times ...")
subprocess.run("python {} {}".format(concatenate_output, config_file), shell=True, check=True)


# %%
end_time = time.time()
print("The total time usage of running DART and PFLOTRAN is %.3f (second): " % (end_time-start_time))



