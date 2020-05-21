"""Generate the prior ensemble in NetCDF format by reading (1) the model output; (2) PFLOTRAN parameter.h5; (3) the list of variable names to be updated/analyzed/assimilated."""

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

from parse_pflotran_files_utils import pflotran_files

###############################
# Configurations
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

dart_data_dir       = configs["other_dir_cfg"]["dart_data_dir"]
pflotran_out_file   = configs["file_cfg"]["pflotran_out_file"]
pflotran_para_file  = configs["file_cfg"]["pflotran_para_file"]
material_id_file    = configs["file_cfg"]["material_id_file"]
dart_prior_file     = configs["file_cfg"]["dart_prior_nc_file"]
dart_prior_file_temp = configs["file_cfg"]["dart_prior_nc_file_temp"]
dart_input_list     = configs["file_cfg"]["dart_input_list_file"]
dart_posterior_file = configs["file_cfg"]["dart_posterior_nc_file"]
dart_output_list    = configs["file_cfg"]["dart_output_list_file"]
dart_prior_template = configs["file_cfg"]["dart_prior_template_file"]
model_time          = float(configs["time_cfg"]["current_model_time"])  # days
model_time_list     = configs["time_cfg"]["model_time_list"]
ntimestep           = int(configs["da_cfg"]["ntimestep"])
nens                = configs["da_cfg"]["nens"]
spinup_time         = configs["time_cfg"]["spinup_length"]

# Parameters and observations configuration
# obs_set               = configs["obspara_set_cfg"]["obs_set"]
obs_pflotran_set      = configs["obspara_set_cfg"]["obs_pflotran_set"]
para_set              = configs["obspara_set_cfg"]["para_set"]
para_take_log_set     = configs["obspara_set_cfg"]["para_take_log_set"]
para_homogeneous      = configs["obspara_set_cfg"]["para_homogeneous"]
# para_material_id_set  = configs["obspara_set_cfg"]["para_material_id_set"]
# para_hdf_dataset_name_set = configs["obspara_set_cfg"]["para_hdf_dataset_name_set"]
# para_isotropic_set        = configs["obspara_set_cfg"]["para_isotropic_set"]

iteration_step            = configs["da_cfg"]["enks_mda_iteration_step"]
total_iterations          = configs["da_cfg"]["enks_mda_total_iterations"]
save_immediate_mda_result = configs["file_cfg"]["save_immediate_mda_result"]

# Convert the model_time_list to a list (model_time_list = 0 in the first model tim)
if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]
# if not isinstance(para_material_id_set, list):
#     para_material_id_set = [para_material_id_set]

# Get the start and end time of the current assimilation window
assim_window_fixed = configs["da_cfg"]["assim_window_fixed"]
if assim_window_fixed:
    assim_window = float(configs["da_cfg"]["assim_window_size"])  # days
else:
    assim_window_list = configs["da_cfg"]["assim_window_list"]  # days
    if not isinstance(assim_window_list, list):
        assim_window_list = [assim_window_list]
    assim_window = assim_window_list[len(model_time_list)-1]
assim_end_days    = configs["da_cfg"]["assim_end_days"]
assim_end_seconds = configs["da_cfg"]["assim_end_seconds"]
assim_start_days    = configs["da_cfg"]["assim_start_days"]
assim_start_seconds = configs["da_cfg"]["assim_start_seconds"]
start_obs_sec = assim_start_days * 86400 + assim_start_seconds
end_obs_sec   = assim_end_days * 86400 + assim_end_seconds
assim_start, assim_end = start_obs_sec / 86400, end_obs_sec / 86400
# start_obs    , end_obs     = model_time - assim_window / 2. + spinup_time, model_time + assim_window / 2. + spinup_time
# start_obs    , end_obs     = model_time - assim_window / 2., model_time + assim_window / 2.
# start_obs_sec, end_obs_sec = start_obs * 86400,              end_obs * 86400

# ndigit = np.ceil(np.log10(ntimestep), dtype=int)
ndigit_time = int(ceil(log10(ntimestep))) + 1
ndigit_ens = int(ceil(log10(nens))) + 1

# Get the list of all required PFLOTRAN variables
# if isinstance(obs_set, str):
#     obs_set = [obs_set]
if isinstance(obs_pflotran_set, str):
    obs_pflotran_set = [obs_pflotran_set]
if isinstance(para_set, str):
    para_set = [para_set]
    para_take_log_set = [para_take_log_set]
# pflotran_var_set = obs_set + para_set

# Get some constants
one_sec        = 1. / 86400.  # one second in fractional days
ens_set        = np.arange(1, nens + 1)  # the set of ensembles
model_time_ind = len(model_time_list)  # the ith model step
missing_value  = 99999

###############################
# Parse the PFLOTRAN variables to be updated/analyzed/assimilated
###############################
p                = re.compile('[A-Z_]+')
obs_pflotran_set = [p.search(v).group() for v in obs_pflotran_set]
para_set         = [p.search(v).group() for v in para_set]
# obs_set          = [p.search(v).group() for v in obs_set]
# pflotran_var_set = [p.search(v).group() for v in pflotran_var_set]


###############################
# Get the lists of netcdf file names for DART inputs and outputs
###############################
# Get the file names of all ensembles for DART restart file
dart_prior_file_set_temp = [
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file_temp)
    for ens in ens_set
]
dart_prior_file_set_temp = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_prior_file_each)
    for dart_prior_file_each in dart_prior_file_set_temp
]
dart_posterior_file_set = [
    # re.sub(r"\[ENS\]", str(ens) + "_time" + str(model_time_ind), dart_posterior_file)
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_posterior_file)
    for ens in ens_set
]
dart_posterior_file_set = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_posterior_file_each)
    for dart_posterior_file_each in dart_posterior_file_set
]
if save_immediate_mda_result:
    dart_posterior_file_set = [
        # TODO: use re
        re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_posterior_file_each)
        for dart_posterior_file_each in dart_posterior_file_set
    ]

# Save the list of prior.nc to dart_input_list
if os.path.isfile(dart_input_list):
    os.remove(dart_input_list)
with open(dart_input_list, "w") as f:
    for fname in dart_prior_file_set_temp:
        f.write(fname + "\n")

# Save the list of posterior.nc to dart_output_list
if os.path.isfile(dart_output_list):
    os.remove(dart_output_list)
with open(dart_output_list, "w") as f:
    for fname in dart_posterior_file_set:
        f.write(fname + "\n")


###############################
# Initialize a dictionary for saving the state/parameter values
###############################
# dart_var_dict = dict.fromkeys(pflotran_var_set)
dart_para_dict = dict.fromkeys(para_set)
dart_state_dict = dict.fromkeys(obs_pflotran_set)


###############################
# Initialize PFLOTRAN file reader
###############################
pflotran_file_reader = pflotran_files(config_nml)


###############################
# Obtain the prior ensemble of model states
###############################
dart_state_dict, ntime_state, time_state, nloc_state, x_loc_state, y_loc_state, z_loc_state, cell_ids_state = \
    pflotran_file_reader.read_pflotran_output(assim_start, assim_end, missing_value)


###############################
# Obtain the prior ensemble of parameters to be updated
###############################
if not para_homogeneous:
    dart_para_dict, ntime_para, time_para, nloc_para, x_loc_para, y_loc_para, z_loc_para, cell_ids_para = \
        pflotran_file_reader.read_pflotran_parameter(missing_value=missing_value)
else:
    dart_para_dict, ntime_para, time_para, nloc_para, cell_ids_para = \
        pflotran_file_reader.read_pflotran_parameter(missing_value=missing_value)
    x_loc_para, y_loc_para, z_loc_para = x_loc_state.mean(), y_loc_state.mean(), z_loc_state.mean()


###############################
# Convert the prior state/parameter information to NetCDF file and
# Create an empty NetCDF file for each posterior ensemble
###############################
for i in range(nens):

    ens = ens_set[i]

    # print("Converting state/parameter into NetCDF file for ensemble %d..." % ens)

    nc_fname_prior     = dart_prior_file_set_temp[i]
    nc_fname_posterior = dart_posterior_file_set[i]

    # Remove the NetCDF if it already exists
    if os.path.isfile(nc_fname_prior):
        os.remove(nc_fname_prior)
    if os.path.isfile(nc_fname_posterior):
        os.remove(nc_fname_posterior)


    ###############################
    # Write the prior to NetCDF file
    ###############################
    with Dataset(nc_fname_prior, 'w') as root_nc_prior:
        # Add global attributes
        root_nc_prior.model_time          = model_time
        root_nc_prior.time_unit           = "day"
        root_nc_prior.assimilation_window = assim_window
        root_nc_prior.description         = "PFLOTRAN output/DART prior data"

        # Create the dimensions
        state_loc_d   = root_nc_prior.createDimension('state_location', nloc_state)
        state_time_d  = root_nc_prior.createDimension('state_time', None)
        para_loc_d   = root_nc_prior.createDimension('para_location', nloc_para)
        para_time_d   = root_nc_prior.createDimension('para_time', None)
        time_d   = root_nc_prior.createDimension('time', None)
        member_d = root_nc_prior.createDimension('member', 1)

        # Create the variables
        time_v             = root_nc_prior.createVariable('time', 'f8', ('time', ))
        assim_start_time_v = root_nc_prior.createVariable('assim_start_time', 'f8', ('time', ))
        assim_end_time_v   = root_nc_prior.createVariable('assim_end_time', 'f8', ('time', ))
        stime_v            = root_nc_prior.createVariable('state_time', 'f8', ('state_time', ))
        ptime_v            = root_nc_prior.createVariable('para_time', 'f8', ('para_time', ))
        member_v           = root_nc_prior.createVariable('member', 'i8', ('member', ))
        xloc_state_v       = root_nc_prior.createVariable('state_x_location', 'f8', ('state_location', ))
        yloc_state_v       = root_nc_prior.createVariable('state_y_location', 'f8', ('state_location', ))
        zloc_state_v       = root_nc_prior.createVariable('state_z_location', 'f8', ('state_location', ))
        xloc_para_v        = root_nc_prior.createVariable('para_x_location', 'f8', ('para_location', ))
        yloc_para_v        = root_nc_prior.createVariable('para_y_location', 'f8', ('para_location', ))
        zloc_para_v        = root_nc_prior.createVariable('para_z_location', 'f8', ('para_location', ))
        state_cell_ids_v   = root_nc_prior.createVariable('state_cell_ids', 'i8', ('state_location', ))
        para_cell_ids_v    = root_nc_prior.createVariable('para_cell_ids', 'i8', ('para_location', ))

        # Get the center of the assimilation window in unit day, as required by DART's read_model_time() subroutine time_v.units    = "day"
        time_v.units = "day"
        time_v[:]    = model_time
        time_v.calendar = 'none'
        assim_start_time_v.units = "day"
        assim_start_time_v[:]    = assim_start
        assim_start_time_v.calendar = 'none'
        assim_end_time_v.units = "day"
        assim_end_time_v[:]    = assim_end
        assim_end_time_v.calendar = 'none'
        member_v[:]     = ens

        # Get the observation time steps
        stime_v.units    = "day"
        # stime_v[:]       = time_vset_assim_day
        stime_v[:]       = time_state
        stime_v.calendar = 'none'

        # Get the parameter time steps
        ptime_v.units    = "day"
        ptime_v[:]       = time_para
        ptime_v.calendar = 'none'

        member_v.type, time_v.type = 'dimension_value', 'dimension_value'
        stime_v.type, ptime_v.type = 'dimension_value', 'dimension_value'
        xloc_state_v.units, xloc_state_v.type = 'm','dimension_value'
        yloc_state_v.units, yloc_state_v.type = 'm','dimension_value'
        zloc_state_v.units, zloc_state_v.type = 'm','dimension_value'

        # Write coordinates values
        # xloc_state_grid, yloc_state_grid, zloc_state_grid = np.meshgrid(x_loc_state, y_loc_state, z_loc_state)
        # xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = xloc_state_grid.flatten(), yloc_state_grid.flatten(), zloc_state_grid.flatten()
        # xloc_para_grid, yloc_para_grid, zloc_para_grid = np.meshgrid(x_loc_para, y_loc_para, z_loc_para)
        # xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = xloc_para_grid.flatten(), yloc_para_grid.flatten(), zloc_para_grid.flatten()
        xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = x_loc_state, y_loc_state, z_loc_state
        xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = x_loc_para, y_loc_para, z_loc_para

        # Write cell ids
        state_cell_ids_v[:] = cell_ids_state
        para_cell_ids_v[:] = cell_ids_para

        # Write the values to the variables
        # for varn in obs_set:
        for varn in obs_pflotran_set:
            # if dart_var_dict[varn] is None:
            if dart_state_dict[varn] is None:
                print( "%s is not available in PFLOTRAN HDF5 output" % varn)
                continue
            vargrp = root_nc_prior.createVariable(varn, 'f8', ('state_time', 'state_location'), fill_value=missing_value)
            vargrp.type = 'state_value'
            vargrp.unit = dart_state_dict[varn]["unit"]
            vargrp[:] = dart_state_dict[varn]["value"][i,:,:]

        # Write the values to the variables
        for varn in para_set:
            para_take_log = para_take_log_set[para_set.index(varn)]
            # if dart_var_dict[varn] is None:
            if dart_para_dict[varn] is None:
                print( "%s is not available in PFLOTRAN parameter HDF5 file" % varn)
                continue
            vargrp = root_nc_prior.createVariable(varn, 'f8', ('para_time', 'para_location'), fill_value=missing_value)
            vargrp.type = 'parameter_value'
            vargrp.unit = dart_para_dict[varn]["unit"]
            vargrp.history = "logged" if para_take_log else "not logged"
            vargrp[:] = dart_para_dict[varn]["value"][i,:,:]
            # print(i, np.sum(dart_para_dict[varn]["value"][i,:,:] != missing_value), varn)

    # root_nc_prior.close()

    ###############################
    # Construct NetCDF file for posterior
    ###############################
    # root_nc_posterior = Dataset(nc_fname_posterior, 'w')
    with Dataset(nc_fname_posterior, 'w') as root_nc_posterior:
        # Add global attributes
        root_nc_posterior.model_time          = model_time
        root_nc_posterior.time_unit           = "day"
        root_nc_posterior.assimilation_window = assim_window
        root_nc_posterior.description         = "PFLOTRAN output/DART prior data"

        # Create the dimensions
        state_loc_d   = root_nc_posterior.createDimension('state_location', nloc_state)
        state_time_d  = root_nc_posterior.createDimension('state_time', None)
        para_loc_d   = root_nc_posterior.createDimension('para_location', nloc_para)
        para_time_d   = root_nc_posterior.createDimension('para_time', None)
        time_d   = root_nc_posterior.createDimension('time', None)
        member_d = root_nc_posterior.createDimension('member', 1)

        # Create the variables
        time_v             = root_nc_posterior.createVariable('time', 'f8', ('time', ))
        assim_start_time_v = root_nc_posterior.createVariable('assim_start_time', 'f8', ('time', ))
        assim_end_time_v   = root_nc_posterior.createVariable('assim_end_time', 'f8', ('time', ))
        stime_v            = root_nc_posterior.createVariable('state_time', 'f8', ('state_time', ))
        ptime_v            = root_nc_posterior.createVariable('para_time', 'f8', ('para_time', ))
        member_v           = root_nc_posterior.createVariable('member', 'i8', ('member', ))
        xloc_state_v       = root_nc_posterior.createVariable('state_x_location', 'f8', ('state_location', ))
        yloc_state_v       = root_nc_posterior.createVariable('state_y_location', 'f8', ('state_location', ))
        zloc_state_v       = root_nc_posterior.createVariable('state_z_location', 'f8', ('state_location', ))
        xloc_para_v        = root_nc_posterior.createVariable('para_x_location', 'f8', ('para_location', ))
        yloc_para_v        = root_nc_posterior.createVariable('para_y_location', 'f8', ('para_location', ))
        zloc_para_v        = root_nc_posterior.createVariable('para_z_location', 'f8', ('para_location', ))
        state_cell_ids_v   = root_nc_prior.createVariable('state_cell_ids', 'i8', ('state_location', ))
        para_cell_ids_v    = root_nc_prior.createVariable('para_cell_ids', 'i8', ('para_location', ))

        # Get the center of the assimilation window in unit day, as required by DART's read_model_time() subroutine time_v.units    = "day"
        time_v.units = "day"
        time_v[:]       = model_time
        time_v.calendar = 'none'
        assim_start_time_v.units = "day"
        assim_start_time_v[:]    = assim_start
        assim_start_time_v.calendar = 'none'
        assim_end_time_v.units = "day"
        assim_end_time_v[:]    = assim_end
        assim_end_time_v.calendar = 'none'
        member_v[:]     = ens

        # Get the observation time steps
        stime_v.units    = "day"
        # stime_v[:]       = time_vset_assim_day
        stime_v[:]       = time_state
        stime_v.calendar = 'none'

        # Get the parameter time steps
        ptime_v.units    = "day"
        ptime_v[:]       = time_para
        ptime_v.calendar = 'none'

        member_v.type, time_v.type = 'dimension_value', 'dimension_value'
        stime_v.type, ptime_v.type = 'dimension_value', 'dimension_value'
        xloc_state_v.units, xloc_state_v.type = 'm','dimension_value'
        yloc_state_v.units, yloc_state_v.type = 'm','dimension_value'
        zloc_state_v.units, zloc_state_v.type = 'm','dimension_value'

        # Write coordinates values
        # xloc_state_grid, yloc_state_grid, zloc_state_grid = np.meshgrid(x_loc_state, y_loc_state, z_loc_state)
        # xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = xloc_state_grid.flatten(), yloc_state_grid.flatten(), zloc_state_grid.flatten()
        # xloc_para_grid, yloc_para_grid, zloc_para_grid = np.meshgrid(x_loc_para, y_loc_para, z_loc_para)
        # xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = xloc_para_grid.flatten(), yloc_para_grid.flatten(), zloc_para_grid.flatten()
        xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = x_loc_state, y_loc_state, z_loc_state
        xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = x_loc_para, y_loc_para, z_loc_para

        # Write cell ids
        state_cell_ids_v[:] = cell_ids_state
        para_cell_ids_v[:] = cell_ids_para

        # Write the values to the variables
        for varn in obs_pflotran_set:
            if dart_state_dict[varn] is None:
                print( "%s is not available in PFLOTRAN HDF5 output" % varn)
                continue
            vargrp = root_nc_posterior.createVariable(varn, 'f8', ('state_time', 'state_location'), fill_value=missing_value)
            vargrp.type = 'state_value'
            vargrp.unit = dart_state_dict[varn]["unit"]

        # Write the values to the variables
        for varn in para_set:
            para_take_log = para_take_log_set[para_set.index(varn)]
            if dart_para_dict[varn] is None:
                print( "%s is not available in PFLOTRAN parameter HDF5 file" % varn)
                continue
            vargrp = root_nc_posterior.createVariable(varn, 'f8', ('para_time', 'para_location'), fill_value=missing_value)
            vargrp.type = 'parameter_value'
            vargrp.history = "logged" if para_take_log else "not logged"
            vargrp.unit = dart_para_dict[varn]["unit"]

    # root_nc_posterior.close()

    ###############################
    # Copy the first ensemble to the prior template
    ###############################
    if i == 1:
        shutil.copyfile(nc_fname_prior, dart_prior_template)


###############################
# Get the prior when it is the first iteration at each time window
# or at each iteration if all iterations should be saved
###############################
# Get the file names of all ensembles for the prior file
dart_prior_file_set = [
    # re.sub(r"\[ENS\]", str(ens) + "_time" + str(model_time_ind), dart_prior_file)
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file)
    for ens in ens_set
]
dart_prior_file_set = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_prior_file_each)
    for dart_prior_file_each in dart_prior_file_set
]

if save_immediate_mda_result:
    dart_prior_file_set = [
        re.sub("[.]nc", "_iter"+str(iteration_step)+".nc", dart_prior_file_each)
        for dart_prior_file_each in dart_prior_file_set
    ]
    # Copy the temporary prior file to prior file
    for i in range(nens):
        subprocess.run("cd {0}; cp {1} {2}".format(dart_data_dir, dart_prior_file_set_temp[i], dart_prior_file_set[i]), shell=True, check=True)

else:
    if iteration_step == 1:
        # Copy the temporary prior file to prior file
        for i in range(nens):
            subprocess.run("cd {0}; cp {1} {2}".format(dart_data_dir, dart_prior_file_set_temp[i], dart_prior_file_set[i]), shell=True, check=True)
