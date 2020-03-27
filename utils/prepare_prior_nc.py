"""Generate the prior ensemble in NetCDF format by reading (1) the model output in HDF; (2) PFLOTRAN parameter.h5; (3) the list of variable names to be updated/analyzed/assimilated."""
# Note that this is for parsing structured grids.

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

dart_data_dir       = configs["other_dir_cfg"]["dart_data_dir"]
pflotran_out_file   = configs["file_cfg"]["pflotran_out_file"]
pflotran_para_file  = configs["file_cfg"]["pflotran_para_file"]
dart_prior_file     = configs["file_cfg"]["dart_prior_nc_file"]
dart_prior_file_temp = configs["file_cfg"]["dart_prior_nc_file_temp"]
dart_input_list     = configs["file_cfg"]["dart_input_list_file"]
dart_posterior_file = configs["file_cfg"]["dart_posterior_nc_file"]
dart_output_list    = configs["file_cfg"]["dart_output_list_file"]
dart_prior_template = configs["file_cfg"]["dart_prior_template_file"]
obs_set             = configs["obspara_set_cfg"]["obs_set"]
para_set            = configs["obspara_set_cfg"]["para_set"]
para_take_log       = configs["obspara_set_cfg"]["para_take_log"]
model_time          = float(configs["time_cfg"]["current_model_time"])  # days
model_time_list     = configs["time_cfg"]["model_time_list"]
ntimestep           = int(configs["da_cfg"]["ntimestep"])
nens                = configs["da_cfg"]["nens"]
spinup_time         = configs["time_cfg"]["spinup_length"]
para_homogeneous    = configs["obspara_set_cfg"]["para_homogeneous"]

# Convert the model_time_list to a list (model_time_list = 0 in the first model tim)
if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

# Get the start and end time of the current assimilation window
assim_window_fixed = configs["da_cfg"]["assim_window_fixed"]
if assim_window_fixed:
    assim_window = float(configs["da_cfg"]["assim_window_size"])  # days
else:
    assim_window_list = configs["da_cfg"]["assim_window_list"]  # days
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
if isinstance(obs_set, str):
    obs_set = [obs_set]
if isinstance(para_set, str):
    para_set = [para_set]
pflotran_var_set = obs_set + para_set

# Get some constants
one_sec        = 1. / 86400.  # one second in fractional days
ens_set        = np.arange(1, nens + 1)  # the set of ensembles
model_time_ind = len(model_time_list)  # the ith model step

# Parse the PFLOTRAN variables to be updated/analyzed/assimilated
p = re.compile('[A-Z_]+')
obs_set   = [p.search(v).group() for v in obs_set]
para_set  = [p.search(v).group() for v in para_set]
pflotran_var_set = [p.search(v).group() for v in pflotran_var_set]

# Get the file names of all ensembles for PFLOTRAN output
pflotran_out_file_set = [
    re.sub(r"\[ENS\]", str(ens), pflotran_out_file) for ens in ens_set
]

# Check the existences of these files
for f in pflotran_out_file_set:
    if not os.path.isfile(f):
        raise Exception("The PFLOTRAN output file %s does not exits!" % f)


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

###############################
# Save the list of prior.nc to dart_input_list
###############################
if os.path.isfile(dart_input_list):
    os.remove(dart_input_list)
with open(dart_input_list, "w") as f:
    for fname in dart_prior_file_set_temp:
        f.write(fname + "\n")

###############################
# Save the list of posterior.nc to dart_output_list
###############################
if os.path.isfile(dart_output_list):
    os.remove(dart_output_list)
with open(dart_output_list, "w") as f:
    for fname in dart_posterior_file_set:
        f.write(fname + "\n")

###############################
# Load the parameter.h5
###############################
f_para = h5py.File(pflotran_para_file, "r")
para_pflotran_set = f_para.keys()

###############################
# Initialize a dictionary for saving the state/parameter values
###############################
dart_var_dict = dict.fromkeys(pflotran_var_set)

###############################
# Convert the prior state/parameter information to NetCDF file and
# Create an empty NetCDF file for each posterior ensemble
###############################
for i in range(nens):

    ens = ens_set[i]

    # print("Converting state/parameter into NetCDF file for ensemble %d..." % ens)

    nc_fname_prior     = dart_prior_file_set_temp[i]
    nc_fname_posterior = dart_posterior_file_set[i]
    pf_fname           = pflotran_out_file_set[i]

    # Remove the NetCDF if it already exists
    if os.path.isfile(nc_fname_prior):
        os.remove(nc_fname_prior)
    root_nc_prior = Dataset(nc_fname_prior, 'w')
    if os.path.isfile(nc_fname_posterior):
        os.remove(nc_fname_posterior)
    root_nc_posterior = Dataset(nc_fname_posterior, 'w')

    ###############################
    # Read the model states from the model output in HDF
    # The following information is expected to be read from the model HDF output
    #   - the last time step/time unit
    #   - the grids information
    #   - the state parameter values required in pflotran_var_set
    ###############################
    # f_out = h5py.File(pf_fname, "r")
    with h5py.File(pf_fname, "r") as f_out:
        # Get the grids/coordinates of the domain
        coordinates = f_out["Coordinates"]
        # TODO: check whether PFLOTRAN generates outputs on the boundaries or the centers of cells.
        # # The following options are choosing the center of the grid
        # x_loc_state = coordinates['X [m]'][:]
        # y_loc_state = coordinates['Y [m]'][:]
        # z_loc_state = coordinates['Z [m]'][:]
        # x_loc_state = [(a + b) / 2 for a, b in zip(x_loc_state[:-1], x_loc_state[1:])]
        # y_loc_state = [(a + b) / 2 for a, b in zip(y_loc_state[:-1], y_loc_state[1:])]
        # z_loc_state = [(a + b) / 2 for a, b in zip(z_loc_state[:-1], z_loc_state[1:])]
        # The following choose the left boundaries
        x_loc_state = coordinates['X [m]'][:-1]
        y_loc_state = coordinates['Y [m]'][:-1]
        z_loc_state = coordinates['Z [m]'][:-1]
        nx_state, ny_state, nz_state = len(x_loc_state), len(y_loc_state), len(z_loc_state)

        # Get all the time steps
        time_set_o = np.array([t for t in list(f_out.keys()) if t.startswith("Time")])
        time_set   = np.array([t.split()[1:] for t in time_set_o])
        time_vset  = np.array([float(t[0]) for t in time_set])
        time_unit  = time_set[0][1]

        # Get the time steps within the assimilation window
        # print(time_vset, start_obs_sec, end_obs_sec)
        time_set_assim_ind = (time_vset > start_obs_sec) & (time_vset <= end_obs_sec)
        time_vset_assim    = time_vset[time_set_assim_ind]
        time_set_o_assim   = time_set_o[time_set_assim_ind]
        time_set_assim     = time_set[time_set_assim_ind]

        if time_unit in ["s", "sec", "second"]: # Convert from seconds to fractional days
            time_vset_assim_day  = time_vset_assim / 86400. 

        ntime_state = len(time_vset_assim)
        nloc_state  = nx_state*ny_state*nz_state

        # Initialize the dart_var_dict
        # for varn in pflotran_var_set:
        for varn in obs_set:
            dart_var_dict[varn] = {"value": np.zeros([ntime_state, nloc_state]), "unit": ""}

        # Get the state variable values required in obs_set
        for j in range(ntime_state):
            time_o           = time_set_o_assim[j]
            dataset          = f_out[time_o]
            pl_out_var_set_o = list(dataset.keys())
            # pl_out_var_dict  = dict()
            for v in pl_out_var_set_o:
                # Get the variable name and unit from the original variable name
                varinfo = v.split()
                if len(varinfo) == 1:
                    varn, varunit = varinfo[0].upper(), ''
                elif len(varinfo) == 2:
                    varn, varunit = varinfo[0].upper(), varinfo[1]
                else:
                    raise Exception('Invalid variable name %s!' % v)
                # Check if the variable is required by pflotran_var_set
                # if varn in pflotran_var_set:
                if varn in obs_set:
                    # dart_var_dict[varn]["value"].append(dataset[v][:]) 
                    # TODO: make sure flatten() follows the same order as np.meshgrid
                    dart_var_dict[varn]["value"][j,:] = dataset[v][:].flatten()
                    dart_var_dict[varn]["unit"] = varunit 

    ###############################
    # Read the parameter.h5
    # The following information is expected to be read from the model HDF output
    #   - the rest of the state parameter values required in pflotran_var_set
    # Note that parameter.h5 has a different structure compared with pflotran output
    # TODO a more strict structure is required. For now, we assume only one value is
    # needed or the whole domain is isotropic.
    ###############################
    if not all(elem in para_pflotran_set for elem in para_set):
        raise Exception("The parameter HDF file does not contain all the required parameters to be estimated!")

    if para_homogeneous:
        for varn in para_pflotran_set:
            if varn in pflotran_var_set:
                value = f_para[varn][:][i]
                # dart_var_dict[varn] = {"value": value * np.ones([ntime, nloc]), "unit": ""}
                if para_take_log:
                    dart_var_dict[varn] = {"value": np.log10(value), "unit": ""}
                else:
                    dart_var_dict[varn] = {"value": value, "unit": ""}
        x_loc_para, y_loc_para, z_loc_para = x_loc_state[int(nx_state/2)], y_loc_state[int(ny_state/2)], z_loc_state[int(nz_state/2)] 
        nx_para, ny_para, nz_para = 1, 1, 1
        time_para = [model_time]
    else:
        raise Exception("The heterongeous parameters setting has not been implemented yet!")
    
    nloc_para  = nx_para * ny_para * nz_para
    ntime_para = len(time_para)

    ###############################
    # Write the prior to NetCDF file
    ###############################
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
    time_v   = root_nc_prior.createVariable('time', 'f8', ('time', ))
    assim_start_time_v = root_nc_prior.createVariable('assim_start_time', 'f8', ('time', ))
    assim_end_time_v   = root_nc_prior.createVariable('assim_end_time', 'f8', ('time', ))
    stime_v  = root_nc_prior.createVariable('state_time', 'f8', ('state_time', ))
    ptime_v  = root_nc_prior.createVariable('para_time', 'f8', ('para_time', ))
    member_v = root_nc_prior.createVariable('member', 'i8', ('member', ))
    xloc_state_v = root_nc_prior.createVariable('state_x_location', 'f8', ('state_location', ))
    yloc_state_v = root_nc_prior.createVariable('state_y_location', 'f8', ('state_location', ))
    zloc_state_v = root_nc_prior.createVariable('state_z_location', 'f8', ('state_location', ))
    xloc_para_v = root_nc_prior.createVariable('para_x_location', 'f8', ('para_location', ))
    yloc_para_v = root_nc_prior.createVariable('para_y_location', 'f8', ('para_location', ))
    zloc_para_v = root_nc_prior.createVariable('para_z_location', 'f8', ('para_location', ))

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
    stime_v[:]       = time_vset_assim_day
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
    xloc_state_grid, yloc_state_grid, zloc_state_grid = np.meshgrid(x_loc_state, y_loc_state, z_loc_state)
    xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = xloc_state_grid.flatten(), yloc_state_grid.flatten(), zloc_state_grid.flatten() 
    xloc_para_grid, yloc_para_grid, zloc_para_grid = np.meshgrid(x_loc_para, y_loc_para, z_loc_para)
    xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = xloc_para_grid.flatten(), yloc_para_grid.flatten(), zloc_para_grid.flatten() 

    # Write the values to the variables
    for varn in obs_set:
        if dart_var_dict[varn] is None:
            print( "%s is not available in PFLOTRAN HDF5 output" % varn)
            continue
        vargrp = root_nc_prior.createVariable(varn, 'f8', ('state_time', 'state_location'))
        vargrp.type = 'state_value'
        vargrp.unit = dart_var_dict[varn]["unit"]
        vargrp[:] = dart_var_dict[varn]["value"]

    # Write the values to the variables
    for varn in para_set:
        if dart_var_dict[varn] is None:
            print( "%s is not available in PFLOTRAN parameter HDF5 file" % varn)
            continue
        vargrp = root_nc_prior.createVariable(varn, 'f8', ('para_time', 'para_location'))
        vargrp.type = 'parameter_value'
        vargrp.unit = dart_var_dict[varn]["unit"]
        vargrp.history = "logged" if para_take_log else "not logged"
        vargrp[:] = dart_var_dict[varn]["value"]

    root_nc_prior.close()

    ###############################
    # Construct NetCDF file for posterior
    ###############################
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
    time_v   = root_nc_posterior.createVariable('time', 'f8', ('time', ))
    assim_start_time_v = root_nc_posterior.createVariable('assim_start_time', 'f8', ('time', ))
    assim_end_time_v   = root_nc_posterior.createVariable('assim_end_time', 'f8', ('time', ))
    stime_v  = root_nc_posterior.createVariable('state_time', 'f8', ('state_time', ))
    ptime_v  = root_nc_posterior.createVariable('para_time', 'f8', ('para_time', ))
    member_v = root_nc_posterior.createVariable('member', 'i8', ('member', ))
    xloc_state_v = root_nc_posterior.createVariable('state_x_location', 'f8', ('state_location', ))
    yloc_state_v = root_nc_posterior.createVariable('state_y_location', 'f8', ('state_location', ))
    zloc_state_v = root_nc_posterior.createVariable('state_z_location', 'f8', ('state_location', ))
    xloc_para_v = root_nc_posterior.createVariable('para_x_location', 'f8', ('para_location', ))
    yloc_para_v = root_nc_posterior.createVariable('para_y_location', 'f8', ('para_location', ))
    zloc_para_v = root_nc_posterior.createVariable('para_z_location', 'f8', ('para_location', ))

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
    stime_v[:]       = time_vset_assim_day
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
    xloc_state_grid, yloc_state_grid, zloc_state_grid = np.meshgrid(x_loc_state, y_loc_state, z_loc_state)
    xloc_state_v[:], yloc_state_v[:], zloc_state_v[:] = xloc_state_grid.flatten(), yloc_state_grid.flatten(), zloc_state_grid.flatten() 
    xloc_para_grid, yloc_para_grid, zloc_para_grid = np.meshgrid(x_loc_para, y_loc_para, z_loc_para)
    xloc_para_v[:], yloc_para_v[:], zloc_para_v[:] = xloc_para_grid.flatten(), yloc_para_grid.flatten(), zloc_para_grid.flatten() 

    # Write the values to the variables
    for varn in obs_set:
        if dart_var_dict[varn] is None:
            print( "%s is not available in PFLOTRAN HDF5 output" % varn)
            continue
        vargrp = root_nc_posterior.createVariable(varn, 'f8', ('state_time', 'state_location'))
        vargrp.type = 'state_value'
        vargrp.unit = dart_var_dict[varn]["unit"]

    # Write the values to the variables
    for varn in para_set:
        if dart_var_dict[varn] is None:
            print( "%s is not available in PFLOTRAN parameter HDF5 file" % varn)
            continue
        vargrp = root_nc_posterior.createVariable(varn, 'f8', ('para_time', 'para_location'))
        vargrp.type = 'parameter_value'
        vargrp.history = "logged" if para_take_log else "not logged"
        vargrp.unit = dart_var_dict[varn]["unit"]

    root_nc_posterior.close()

    ###############################
    # Copy the first ensemble to the prior template
    ###############################
    if i == 1:
        shutil.copyfile(nc_fname_prior, dart_prior_template)

f_para.close()


###############################
# Get the prior when it is the first iteration at each window
###############################
enks_mda_iteration_step   = configs["da_cfg"]["enks_mda_iteration_step"]

if enks_mda_iteration_step == 1:
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

    # Copy the temporary prior file to prior file
    for i in range(nens):
        subprocess.run("cd {0}; cp {1} {2}".format(dart_data_dir, dart_prior_file_set_temp[i], dart_prior_file_set[i]), shell=True, check=True)