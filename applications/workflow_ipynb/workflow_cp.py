# To add a new cell, type '#%%'
# To add a new markdown cell, type '#%% [markdown]'
#%% [markdown]
# # Overview
# The **objective** of this notebook is to present the workflow of conducting data assimilation on [PFLOTRAN](https://www.pflotran.org/) by using [DART](https://www.image.ucar.edu/DAReS/DART/). Briefly, the procedures are as follows:
# - [x] [Configuration](#parameter): define directories, file locations, and other parameters
# - [x] [PFLOTRAN preparation](#pflotran_prepare): generate PFLOTRAN input files
# - [x] [PFLOTRAN model spin-up](#pflotran_spinup): conduct model spin-up
# - [x] [DART files preparation](#dart_prepare): add new DART quantities, prepare DART input namelists, prepare DART prior data, prepare observations in DART format, and check ```model_mod``` interface
# - [x] [Generate all the executable files](#dart_executables): generate all the executables, convert observations in DART format, check ```model_mod``` interface, and test the filter
# - [x] [Run DART and PFLOTRAN](#run_dart_pflotran): run the shell script for integrating DART filter and PFLOTRAN model
# 
# Here, we perform inverse modeling on a 1D thermal model for illustration. The model assimilates temperature observation to update its parameters (i.e., flow flux, porosity, and thermal conductivity). For now, the ensemble Kalman filter (EnKF) is used for assimilation.
#%% [markdown]
# <a id='parameter'></a>
# # Step 1: Configuration

#%%
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

#%% [markdown]
# ****************
# **Define the locations of MPI, PFLOTRAN, application folder and DART-PFLOTRAN interface folder**
# 
# It is suggested that <span style="background-color:lightgreen">mpi_exe</span> is defined based on the mpi utility (e.g., mpirun) installed by PFLOTRAN.
# 
# **Important:** You must make sure the <span style="background-color:lightgreen">mpi_exe</span> has the same MPI system as the settings ```MPIFC``` and ```MPILD``` in ```mkmf.template```.

#%%
# MPI settings
mpi_exe_da  = '/usr/local/bin/mpirun'  # The location of mpirun
mpi_exe_pf  = '/Users/jian449/Codes/petsc/arch-darwin-c-opt/bin/mpirun'
ncore_da = 2                       # The number of MPI cores used by DART
ncore_pf = 10                      # The number of MPI cores used by PFLOTRAN
ngroup_pf= 10                       # The number of group used by stochastic running in PFLOTRAN

# PFLOTRAN executable
# pflotran_exe  = '/global/project/projectdirs/m1800/pin/pflotran-haswell/src/pflotran/pflotran'
pflotran_exe  = '/Users/jian449/Codes/pflotran/src/pflotran/pflotran'

# Main directory names
temp_app_dir = "/Users/jian449/Codes/DART/manhattan/models/pflotran/applications/template"          # The template for application folder
app_dir      = "/Users/jian449/Codes/DART/manhattan/models/pflotran/applications/1dthermal"          # The application folder name
dart_dir     = "/Users/jian449/Codes/DART/manhattan"
dart_pf_dir  = "/Users/jian449/Codes/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
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

#%% [markdown]
# ****************
# **Generate the application directory if it does not exit**

#%%
# Create the application directory if it does not exists
if not os.path.isdir(app_dir):
    shutil.copytree(temp_app_dir, app_dir)

#%% [markdown]
# ****************
# **Load all the required file paths/names from ```file_paths.nml```**
# 
# ```file_paths.nml``` defines the relative locations of all the files (e.g., utility files, shell scripts, data files, etc) used by this application

#%%
sys.path.append(dart_pf_dir)
from utils.read_filepaths_nml import read_filepaths_nml


#%%
dirs_cfg, files_cfg      = read_filepaths_nml(app_dir=app_dir, dart_pf_dir=dart_pf_dir)
configs["other_dir_cfg"] = dirs_cfg
configs["file_cfg" ]     = files_cfg
config_file              = files_cfg["config_file"]

#%% [markdown]
# ****************
# **Specify the following types of variables used in DART**
# - the observation data to be assimilated
# - the PFLOTRAN parameters to be analyzed
# - the statistics of the parameters: (1) the value range; (2) the mean and std for randomly sampling; (3) the distribution to be sampled (i.e., normal, truncated normal, and uniform distributions).
# - The list of parameters whose prior would be resampled based on the mean of the corresponding posterior at the previous time step.

#%%
# Observation data to be assimilated
obs_set  = ['TEMPERATURE']

# The PFLOTRAN parameters to be analyzed
# para_set = ['FLOW_FLUX','POROSITY','THERMAL_CONDUCTIVITY']
para_set = ['FLOW_FLUX']

# The statistics of the parameters
# The index order follows para_set
# para_min_set  = [-10.0, 0.01, 0.5]  # The minimum values (-99999 means no lower bound limit)
# para_max_set  = [10.0, 0.7, 2.5]  # The maximum values (99999 means no upper bound limit)
# para_mean_set = [0.0, 0.3, 1.5]  # The mean values
# para_std_set  = [0.5, 0.1, 0.5]  # The standard deviation values
# para_dist_set = ["normal", "normal", "normal"]  # The assumed distribution to be sampled
para_min_set  = [-10.0]  # The minimum values (-99999 means no lower bound limit)
para_max_set  = [10.0]  # The maximum values (99999 means no upper bound limit)
para_mean_set = [0.0]  # The mean values
para_std_set  = [0.5]  # The standard deviation values
para_dist_set = ["normal"]  # The assumed distribution to be sampled

para_resampled_set = ['FLOW_FLUX']   # The parameters to be resampled at each time step
# para_resampled_set = ['']   # The parameters to be resampled at each time step

configs["obspara_set_cfg"] = {"obs_set": obs_set, "para_set": para_set,
                              "para_min_set": para_min_set, "para_max_set": para_max_set,
                              "para_mean_set": para_mean_set, "para_std_set": para_std_set,
                              "para_dist_set": para_dist_set, "para_resampled_set": para_resampled_set}

#%% [markdown]
# ****************
# **Specify the spatial domains of the observation data to be assimilated**
# 
# The limits of x, y, and z are bounded by the minimum and maximum boundaries through (min, max). If the limit is not specified, -99999 and 99999 are assigned for the lower and upper bounds, respectively.

#%%
obs_space_xlimit = [-99999, 99999]
obs_space_ylimit = [-99999, 99999]
obs_space_zlimit = [-0.5, -0.04]
configs["obs_space_cfg"] = {"obs_space_xlimit": obs_space_xlimit, "obs_space_ylimit": obs_space_ylimit,
                            "obs_space_zlimit": obs_space_zlimit}

#%% [markdown]
# ****************
# **Specify the temporal information**
# - model spinup time/start time
# - the map between the begin of observation assimilation and model start time
# - the list of model time or the list of starting assimilation time (starting from zero)
# 
# **note that** model start time is considered after the spinup

#%%
time_cfg = {}
# Model spinup length
time_cfg["spinup_length"]  = 0.1    # spinup time (day)
time_cfg["is_spinup_done"] = False  # whether spinup is conducted

# Model start time
time_cfg["current_model_time"] = 0.     # model start time zero (after spinup)
time_cfg["model_time_list"]    = [0.]   # the list of model time

# Map between assimilation start time and model start time
obs_start   = datetime(2017,4,1,0,0,0) 
assim_start = obs_start + timedelta(days=time_cfg["spinup_length"]) # assimilation time should be after the model spinup
time_cfg["assim_start"] = assim_start.strftime("%Y-%m-%d %H:%M:%S")

# The maximum time for observation
time_cfg["first_obs_time_days"]    = 0
time_cfg["first_obs_time_seconds"] = 0
time_cfg["first_obs_time_size"] = time_cfg["first_obs_time_days"]+float(time_cfg["first_obs_time_seconds"])/86400. # day
time_cfg["last_obs_time_days"]    = 0
time_cfg["last_obs_time_seconds"] = 3600*2
time_cfg["last_obs_time_size"] = time_cfg["last_obs_time_days"]+float(time_cfg["last_obs_time_seconds"])/86400. # day

# Whether the model time exceeds the last observation
time_cfg["exceeds_obs_time"] = time_cfg["current_model_time"] >= time_cfg["last_obs_time_size"]

# Save them to configs
configs["time_cfg"] = time_cfg

#%% [markdown]
# ****************
# **Define the data assimilation configurations**

#%%
da_cfg = {}
# More need to be added...
# And later on, these DA setting can be saved in a txt or pickel file for further loading
da_cfg["obs_reso"]  = 300.0      # observation resolution (second)
da_cfg["nens"]      = 100         # number of ensembles
da_cfg["obs_error"] = 0.02       # the observation error
da_cfg["obs_error_type"] = "relative" # the type of observation error (i.e., relative and absolute)

# Assimilation time window time_step_days+time_step_seconds
# Assimilation window
da_cfg["assim_window_days"]    = 0     # assimilation time window/step (day)
da_cfg["assim_window_seconds"] = 3600  # assimilation time window/step  (second)
da_cfg["assim_window_size"] = da_cfg["assim_window_days"]+float(da_cfg["assim_window_seconds"])/86400. # day

# Assimilation start and end time
da_cfg["assim_start_days"]    = int(floor(time_cfg["current_model_time"]))
da_cfg["assim_start_seconds"] = int((time_cfg["current_model_time"] - da_cfg["assim_start_days"])*86400)
da_cfg["assim_end_days"]      = int(floor(time_cfg["current_model_time"]+da_cfg["assim_window_size"]))
da_cfg["assim_end_seconds"]   = int((time_cfg["current_model_time"]+da_cfg["assim_window_size"] - 
                                     da_cfg["assim_end_days"])*86400-1)

# Compute the number of time steps
# print(np.ceil((time_cfg["last_obs_time_size"] - time_cfg["first_obs_time_size"]) / da_cfg["assim_window_size"]))
da_cfg["ntimestep"] = ceil((time_cfg["last_obs_time_size"] - time_cfg["first_obs_time_size"]) / da_cfg["assim_window_size"])

# The inflation settings used in EnKS-MDA (the alpha value)
da_cfg["enks_mda_alpha"] = [4, 4, 4, 4]  # Note that the summation of the inverse of alpha should be one
da_cfg["enks_mda_iteration_step"] = 1  # the ith iteration (1 for the first iteration)
da_cfg["enks_mda_total_iterations"] = len(da_cfg["enks_mda_alpha"])  # Note that the summation of the inverse of alpha should be one

# Save them to configs
configs["da_cfg"] = da_cfg

#%% [markdown]
# ****************
# **Save all the configurations in pickle**

#%%
configs.write(config_file, force=True)

#%% [markdown]
# ****************
# **Clean up the work, PFLOTRAN output, and DART in/out directories**
shutil.rmtree(dirs_cfg["dart_data_dir"], ignore_errors=True)
shutil.rmtree(dirs_cfg["pflotran_out_dir"], ignore_errors=True)
shutil.rmtree(dirs_cfg["app_work_dir"], ignore_errors=True)
os.mkdir(dirs_cfg["dart_data_dir"])
os.mkdir(dirs_cfg["pflotran_out_dir"])
os.mkdir(dirs_cfg["app_work_dir"])
# subprocess.run("rm -f {}/*".format(dirs_cfg["dart_data_dir"]), shell=True, check=True)
# subprocess.run("rm -f {}/*".format(dirs_cfg["pflotran_out_dir"]), shell=True, check=True)
# subprocess.run("rm -f {}/*".format(dirs_cfg["app_work_dir"]), shell=True, check=True)

#%%
# Save it in a temperory pickle file or namelist???
configs.write(config_file, force=True)

#%% [markdown]
# <a id='pflotran_prepare'></a>
# # Step 2: PFLOTRAN preparation
# *Here, we use Kewei's 1D thermal model as an example for generating PFLOTRAN input card and parameter.h5.*
# 
# In this section, the following procedures are performed:
# - generate PFLOTRAN input deck file ```PFLOTRAN.in```
# - generate the parameter files in HDF 5, ```parameter_prior.h5```, used by PFLOTRAN input deck file
# 
# **Note that**
# - ```PFLOTRAN.in``` for each DA scenario should be prepared by users.

#%%
prep_pflotran_inputdeck, pflotran_in = files_cfg["prep_pflotran_inputdeck_file"], files_cfg["pflotran_in_file"]
prep_pflotran_parameterprior = files_cfg["prep_pflotran_para_file"]

#%% [markdown]
# ****************
# **Generate the ensembles of PFLOTRAN prior**
# 
# **Run code**
# - Run: ```prepare_pflotran_parameterprior.py```
# - Code input arguments (loaded from the configuration file):
#     - <span style="background-color:lightgreen">pflotran_para</span>: filename for ```parameter_prior.h5```
#     - <span style="background-color:lightgreen">obs_resolution, obs_error, nens, spinup_length, spinup</span>: data assimilation settings (i.e., observation timestep, observation error, number of ensemble, whether it is spinup, **to be revised**)

#%%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_pflotran_parameterprior" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare PFLOTRAN ensemble parameter values...")
subprocess.run("python {} {}".format(prep_pflotran_parameterprior, config_file), shell=True, check=True)

#%% [markdown]
# ****************
# **Generate PFLOTRAN input deck file**
# 
# **Run code**
# - Run: ```prepare_pflotran_inputdeck.py```
# - Code input arguments (loaded from the configuration file):
#     - <span style="background-color:lightgreen">pflotran_in</span>: filename for ```pflotran.in```
#     - <span style="background-color:lightgreen">pflotran_para</span>: filename for ```parameter_prior.h5```
#     - <span style="background-color:lightgreen">obs_resolution, obs_error, nens, spinup_length, spinup</span>: data assimilation settings (i.e., observation timestep, observation error, number of ensemble, whether it is spinup, **to be revised**)

#%%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_pflotran_inputdeck" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare PFLOTRAN input deck...")
subprocess.run("python {} {}".format(prep_pflotran_inputdeck, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$pflotran_in"', 'cat $1')

#%% [markdown]
# <a id='pflotran_spinup'></a>
# # Step 3: PFLOTRAN model spin-up
# Take in the ```pflotran.in``` and ```parameter.h5``` files and conduct the model spin-up by running ```pflotran.sh``` file. The ```pflotran.sh``` is a simple shell script executing ensemble simulation of PFLOTRAN by using MPI.
# 
# **Run the code**
# - Run: ```prepare_pflotran_inpara.py```
# - Code input arguments (loaded from the configuration file):
#     - <span style="background-color:lightgreen">pflotran_exe</span>: location of the executable PFLOTRAN
#     - <span style="background-color:lightgreen">pflotran_in</span>: filename for ```pflotran.in```
#     - <span style="background-color:lightgreen">pflotran_out_dir</span>: directory of PFLOTRAN output
#     - <span style="background-color:lightgreen">nens</span>: number of ensemble
#     - <span style="background-color:lightgreen">mpi_exe, ncore</span>: location of mpirun and number of cpu cores

#%%
pflotran_sh, pflotran_out_dir = files_cfg["pflotran_sh_file"], dirs_cfg["pflotran_out_dir"]


#%%
print("\n")
print("------------------------------------------------------------")
print("Model spinup and run the forward simulation in the first assimilation time window...")
subprocess.run("{} {}".format(pflotran_sh, config_file), shell=True, check=True)

#%% [markdown]
# ****************
# **Once the model spinup finishes, modify the corresponding configuration entry**

#%%
configs["time_cfg"]["is_spinup_done"] = True
configs.write(config_file, force=True)


#%%
# cd $1
# ls *.h5

#%% [markdown]
# <a id='dart_prepare'></a>
# # Step 4: DART files preparation
# In this section, the following procedures are performed:
# - generate the template for DART generic variable quantity files (i.e., ```DEFAULT_obs_kind_mod.F90``` and ```obs_def_pflotran_mod.f90```);
# - generate the DART input namelists;
# - generate DART prior NetCDF data ```prior_ensemble_[ENS].nc``` from PFLOTRAN's parameter and outputs;
# - generate DART posterior NetCDF files (*sharing the same variable names and dimensions as the prior NetCDF files but without the data values*);
# - convert the observation file to DART observation format;
# - check ```model_mod.F90``` based on current setting by using the ```check_model_mod``` provided by DART.
#%% [markdown]
# <a id='dart_generic_prepare'></a>
# ## Generate the templates for DART generic variable quantity files
# - Run: ```list2dartqty.py``` to sequentially generate
#     - a mapping between PFLOTRAN variales and DART generic quantities in ```obs_def_pflotran_mod.F90```
#     - the default DART generic quantity definition file ```DEFAULT_obs_kind_mod.F90```
# - Code input arguments:
#     - <span style="background-color:lightgreen">obs_type</span>: filename for ```DEFAULT_obs_kind_mod.F90```
#     - <span style="background-color:lightgreen">def_obs_kind</span>: filename for ```obs_def_pflotran_mod.F90```
#     - <span style="background-color:lightgreen">pflotran_parastate_set</span>: a list of variables required to be assimilated

#%%
to_dartqty, obs_type_file = files_cfg["to_dartqty_file"], files_cfg["obs_type_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$to_dartqty" "$config_file"', 'python $1 $2 $3 $4')
print("\n")
print("------------------------------------------------------------")
print("Add PFLOTRAN variables to DEFAULT_obs_kind_mod.F90 if necessary...")
subprocess.run("python {} {}".format(to_dartqty, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$obs_type_file"', 'cat $1')

#%% [markdown]
# ## Generate  DART input namelists in ```input.nml```
# 
# The ```input.nml``` file is generated based on a template ```input.nml.template``` by modifying the following namelist entries:
# 
# ```input.nml.template``` $\rightarrow$ ```input.nml```
# 
# |filter_nml|obs_kind_nml|preprocess_nml|model_nml|convertnc_nml|
# |:--:|:--:|:--:|:--:|:--:|
# | input_state_file_list, output_state_file_list, ens_size, async, adv_ens_command, obs_sequence_in_name | assimilate_these_obs_types | input_files, input_obs_kind_mod_file | time_step_days, time_step_seconds, nvar, var_names, template_file, var_qtynames | netcdf_file, out_file |
# 
# **Namelists from DART**
# - [filter_nml](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/modules/assimilation/filter_mod.html): namelist of the main module for driving ensemble filter assimilations
# - [obs_kind_nml](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/modules/observations/obs_kind_mod.html#Namelist): namelist for controling what observation types are to be assimilated
# - [preprocess_nml](https://www.image.ucar.edu/DAReS/Codes/DART/manhattan/assimilation_code/programs/preprocess/preprocess): namelist of the DART-supplied preprocessor program which creates observation kind and observation definition modules from a set of other specially formatted Fortran 90 files
# 
# **Self-defined namelists**
# - model_nml: a self-defined namelist for providing the basic information in the model
#     - time_step_days, time_step_seconds: the assimilation time window
#     - template_file: the template prior NetCDF file for ```model_mod.F90``` to digest the spatial information of the model
#     - var_names: the original variable names
#     - var_qtynames: the corresponding DART variable quantities
#     - nvar: the number of variables
# - convertnc_nml: a self-defined namelist for providing the NetCDF observation file name and the DART observation file name used in ```convert_nc.f90```
#     - netcdf_file: the location of the NetCDF file containing the observation data
#     - out_file: the location of the DART observation file
# 
# **Note that**
# - There are more namelists or other items in the above namelist in input.nml.template. Users can edit the below python dictionary ```inputnml``` to include their modifications.
# - Users can also include more namelists provided by DART by modifying ```inputnml```.
#%% [markdown]
# ***************
# **Assemble all the namelists in input.nml**

#%%
# Parameters for different namelists in input.nml
filter_nml = {"input_state_file_list":files_cfg["dart_input_list_file"],
              "output_state_file_list":files_cfg["dart_output_list_file"],
              "ens_size":da_cfg["nens"],
              "num_output_state_members":da_cfg["nens"],
              "obs_sequence_in_name":files_cfg["obs_dart_file"]}
#               "obs_window_days":obs_window_days,
#               "obs_window_seconds":obs_window_seconds}
obs_kind_nml = {"assimilate_these_obs_types":obs_set}
model_nml = {"time_step_days":da_cfg["assim_window_days"],
             "time_step_seconds":da_cfg["assim_window_seconds"],
             "nvar":len(obs_set)+len(para_set),
             "var_names":obs_set+para_set,
             "template_file":files_cfg["dart_prior_template_file"],
             "var_qtynames":['QTY_PFLOTRAN_'+v for v in obs_set]+['QTY_PFLOTRAN_'+v for v in para_set]}
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
inputnml = {"filter_nml":filter_nml,
            "obs_kind_nml":obs_kind_nml,
            "model_nml":model_nml,
            "preprocess_nml":preprocess_nml,
            "convert_nc_nml":convertnc_nml,
            "model_mod_check_nml":modelmodcheck_nml}


configs["inputnml_cfg"] = inputnml

# Save the configurations
configs.write(config_file, force=True)
# with open(config_pickle, 'wb') as f:
#     pickle.dump(configs, f)

#%% [markdown]
# ***************
# **Run the code**
# - Run: ```prepare_inputnml.py```
# - Code input arguments:
#     - <span style="background-color:lightgreen">input_nml</span>: the ```input.nml``` namelist file
#     - <span style="background-color:lightgreen">input_nml_dict</span>: the ```inputnml.p``` pickle file

#%%
prep_inputnml = files_cfg["prep_inputnml_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s  "$prep_inputnml" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Generate input.nml file...")
subprocess.run("python {} {}".format(prep_inputnml, config_file), shell=True, check=True)

#%% [markdown]
# ## Convert the model output to DART prior NetCDF and  generate the preliminary DART posterior NetCDF file
# - The structure of ```prior_ensemble_[ENS].nc``` and ```posterior_ensemble_[ENS].nc``` files (```[ENS]``` refers to the ensemble number):
# 
# | NetCDF dimensions |                      NetCDF variables                      |
# |:-----------------:|:----------------------------------------------------------:|
# | time: 1           | time: shape(time)                                          |
# | x_location: nx    | x_location: shape(x_location)                              |
# | y_location: ny    | y_location: shape(y_location)                              |
# | z_location: nz    | z_location: shape(z_location)                              |
# | member: 1         | member: shape(member)                                      |
# |                   | physical variable: shape(x_location,y_location,z_location) |
# 
# **Note that** 
# - required by DART, each ```prior_R[ENS].nc``` file only includes the state/parameter values of one ensemble member at one given time. 
# - For the time, we set the initial time as 0, with time units converted *day* (requied by DART's ```read_model_time``` subroutine). 
# - Also, it is different from the definition for the [observation NetCDF](#observationconvertion), because ```prior_R[ENS].nc``` aims for the structured cartesian grids while the observation NetCDF aims for a general case.
# 
# **Run the code**
# - Run: ```prepare_prior_nc.py``` to generate 
#     - the DART prior input file ```prior_ensemble_[ENS].nc```
#     - the DART posterior output file ```prior_ensemble_[ENS].nc``` (*sharing the same variable names and dimensions as the prior files but without the variable values*)
#     - the prior template file (copied from ```prior_ensemble_1.nc```) used by ```input.nml```
#     - the dart_input_list and dart_output_list used by DART
# - Code input arguments:
#     - <span style="background-color:lightgreen">pflotran_out</span>: filename ```R[ENS].h5``` from PFLOTRAN model output
#     - <span style="background-color:lightgreen">pflotran_para</span>: pflotran parameter HDF file ```parameter.h5```
#     - <span style="background-color:lightgreen">dart_prior_nc</span>: filename ```prior_R[ENS].nc``` for the prior input file for DART
#     - <span style="background-color:lightgreen">dart_input_list</span>: filename for recording the list of dart_prior_nc
#     - <span style="background-color:lightgreen">nens</span>: number of ensemble
#     - <span style="background-color:lightgreen">spinup</span>: whether it is spinup (if yes, the time is set to zero; otherwise, the time is read from ```R[ENS].h5```)
#     - <span style="background-color:lightgreen">pflotran_parastate_set</span>: a list of variables to be assimilated

#%%
prep_prior_nc, dart_prior_template = files_cfg["prep_prior_nc_file"], files_cfg["dart_prior_template_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_prior_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Convert PFLOTRAN ensemble outputs and parameters to DART prior NetCDF file...")
subprocess.run("python {} {}".format(prep_prior_nc, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$dart_prior_template"', 'ncdump -h $1')


#%% [markdown]
# <a id='observationconvertion'></a>
# ## Prepare the observation conversion to DART observation format
# In this section, we prepare the process of converting the observation data to DART format. We first convert observation data in raw format into NetCDF format. Then, a fortran script is prepared for the conversion from the NetCDF to to DART format. The structure of NetCDF file for recording observation file.
# 
# | NetCDF dimensions |           NetCDF variables          |
# |:-----------------:|:-----------------------------------:|
# | time: 1           | time: shape(time)                   |
# | location: nloc    | location: shape(location)           |
# |                   | physical variable: shape(time,nloc) |
# 
# **Note that** 
# - if the time calendar follows *gregorian*, the time unit should be entered as ```seconds since YYYY-MM-DD HH:MM:SS```. Otherwise, put the time calender as *None* and time unit as ```second``` (make sure convert your measurement times to seconds).
#%% [markdown]
# ***************
# **Convert the raw csv temperature observations to NetCDF file**
# - Run: ```csv2nc.py```
# - Code input arguments:
#     - <span style="background-color:lightgreen">obs_original</span>: filename for the original observed temperature file
#     - <span style="background-color:lightgreen">obs_nc</span>: filename for the observation NetCDF file
#     - <span style="background-color:lightgreen">assim_start_str</span>: the reference time to set zero

#%%
csv_to_nc, obs_nc_original = files_cfg["csv_to_nc_file"], files_cfg["obs_nc_original_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$csv_to_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Convert raw observation data to NetCDF file...")
subprocess.run("python {} {}".format(csv_to_nc, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$obs_nc_original"', 'ncdump -h $1')


#%% [markdown]
# ***************
# **Clip the NetCDF file based on the defined spatial and temporal domains**
# 
# The NetCDF file generated in the previous step is further processed by selecting data observed in the required spatial and temporal domains
# - Run: ```clip_obs_nc.py```

#%%
clip_obs_nc, obs_nc = files_cfg["clip_obs_nc_file"], files_cfg["obs_nc_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$clip_obs_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Clip the NetCDF file based on the defined spatial and temporal domains...")
subprocess.run("python {} {}".format(clip_obs_nc, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$obs_nc"', 'ncdump -h $1')


#%% [markdown]
# ***************
# **Prepare the ```convert_nc.f90``` based on the list of observation variables**
# - Run: ```prepare_convert_nc.py```
# - Code input arguments:
#     - <span style="background-color:lightgreen">obs_nc</span>: filename for the observation NetCDF file

#%%
prep_convert_nc, convert_nc_file = files_cfg["prep_convert_nc_file"], files_cfg["convert_nc_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$prep_convert_nc" "$config_file"', 'python $1 $2')
print("\n")
print("------------------------------------------------------------")
print("Prepare the convert_nc.f90 based on the list of observation variables to be assimilated...")
subprocess.run("python {} {}".format(prep_convert_nc, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$convert_nc_file"', 'head $1')


#%% [markdown]
# <a id='dart_executables'></a>
# # Step 5: Generate all the executable files
# Now, we compile all the executables from ```mkmf_*```. The following executables are generated here:
# - ```preprocess```: for preprocessing the [prepared DART generic variable quantity files prepared](#dart_generic_prepare)
# - ```convert_nc```: for [converting the observations from NetCDF to DART format](#observationconvertion)
# - ```model_mod_check```: for checking ```model_mod.F90``` interface file
# - ```filter```: for conducting the [DART data assimilation](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/filter/filter.html)
#%% [markdown]
# ## Generate the executables
# - Run: ```quickbuild.csh```
# - Code input arguments:
#     - <span style="background-color:lightgreen">app_work_dir</span>: location of the application work folder

#%%
dart_work_dir, app_work_dir = dirs_cfg["dart_work_dir"], dirs_cfg["app_work_dir"]
quickbuild = files_cfg["quickbuild_csh"]


#%%
# print("\n")
# print("------------------------------------------------------------")
# print("Generate all the executables...")
# subprocess.run("cd {}; csh {} {} -mpi".format(dart_work_dir, quickbuild, app_work_dir), shell=True, check=True)

#%% [markdown]
# ## Convert the observation file in NetCDF to DART format
# - Run: ```convert_nc```

#%%
convert_nc, obs_dart = files_cfg["convert_nc_exe"], files_cfg["obs_dart_file"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$app_work_dir" "$convert_nc"', 'cd $1\n$2')
print("\n")
print("------------------------------------------------------------")
print("Conver the NetCDF observation files to DART format...")
subprocess.run("cd {}; {}".format(app_work_dir, convert_nc), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$obs_dart"', 'head -n100 $1')

#%% [markdown]
# ## Check ```model_mod.F90``` interface file
# - Run: ```model_mod_check```

#%%
model_mod_check = files_cfg["model_mod_check_exe"]


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$app_work_dir" "$model_mod_check"', 'cd $1\n$2')


#%% [markdown]
# <a id='run_dart_pflotran'></a>
# # Step 6: Run DART and PFLOTRAN
# In this section, run the shell script to couple DART and PFLOTRAN.

#%%
dart_work_dir     = dirs_cfg["dart_work_dir"]
run_DART_PFLOTRAN = files_cfg["run_filter_csh"]
concatenate_output = files_cfg["concatenate_dart_output_file"]
inputnml_file     = files_cfg["input_nml_file"]
start_time        = time.time()


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$dart_work_dir" "$inputnml_file" "$config_file"', 'cd $1\ncsh run_DART_PFLOTRAN.csh $2 $3')
print("\n")
print("------------------------------------------------------------")
print("Assimilation starts here...")
subprocess.run("cd {}; csh {} {} {}".format(dart_work_dir, run_DART_PFLOTRAN, inputnml_file, config_file), shell=True, check=True)


#%%
# get_ipython().run_cell_magic('script', 'bash -s "$dart_work_dir" "$inputnml_file" "$config_file"', 'cd $1\ncsh run_DART_PFLOTRAN.csh $2 $3')
print("\n")
print("------------------------------------------------------------")
print("Concatenate the prior and posterior at all times ...")
subprocess.run("python {} {}".format(concatenate_output, config_file), shell=True, check=True)

#%%
end_time = time.time()
print("The total time usage of running DART and PFLOTRAN is %.3f (second): " % (end_time-start_time))

#%% [markdown]
# # Some tests
#%% [markdown]
# ## Test ```update_confignml_time.py```

#%%
# update_confignml_time, inputnml_file = files_cfg["update_confignml_time_file"], files_cfg["input_nml_file"]


#%%
# python $1 $2

#%% [markdown]
# ## Test ```filter```
# - Run: ```filter```

#%%
# filter_exe = files_cfg["filter_exe"]


#%%
# cd $1
# $2


#%%


