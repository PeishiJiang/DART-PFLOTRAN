"""Generate the prior ensemble in NetCDF format by reading (1) the model output in HDF; (2) PFLOTRAN parameter.h5; (3) the list of variable names to be updated/analyzed/assimilated."""

# Author: Peishi Jiang

import os
import re
import sys
import h5py
import shutil
import f90nml
import numpy as np
from math import ceil, log10
from netCDF4 import num2date, date2num, Dataset

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

pflotran_out_file   = configs["file_cfg"]["pflotran_out_file"]
pflotran_para_file  = configs["file_cfg"]["pflotran_para_file"]
dart_prior_file     = configs["file_cfg"]["dart_prior_nc_file"]
dart_input_list     = configs["file_cfg"]["dart_input_list_file"]
dart_posterior_file = configs["file_cfg"]["dart_posterior_nc_file"]
dart_output_list    = configs["file_cfg"]["dart_output_list_file"]
dart_prior_template = configs["file_cfg"]["dart_prior_template_file"]
obs_set             = configs["obspara_set_cfg"]["obs_set"]
para_set            = configs["obspara_set_cfg"]["para_set"]
model_time          = float(configs["time_cfg"]["current_model_time"])  # days
model_time_list     = configs["time_cfg"]["model_time_list"]
ntimestep           = int(configs["da_cfg"]["ntimestep"])
assim_window        = float(configs["da_cfg"]["assim_window_size"])  # days
nens                = configs["da_cfg"]["nens"]
spinup_time         = configs["time_cfg"]["spinup_length"]

# Get the start and end time of the current assimilation window
# start_obs    , end_obs     = model_time - assim_window / 2. + spinup_time, model_time + assim_window / 2. + spinup_time
start_obs    , end_obs     = model_time - assim_window / 2., model_time + assim_window / 2.
start_obs_sec, end_obs_sec = start_obs * 86400,              end_obs * 86400

# ndigit = np.ceil(np.log10(ntimestep), dtype=int)
ndigit_time = int(ceil(log10(ntimestep))) + 1
ndigit_ens = int(ceil(log10(nens))) + 1

# Get the list of all required PFLOTRAN variables
if isinstance(obs_set, str):
    obs_set = [obs_set]
if isinstance(para_set, str):
    para_set = [para_set]
pflotran_var_set = obs_set + para_set

# Convert the model_time_list to a list (model_time_list = 0 in the first model tim)
if not isinstance(model_time_list, list):
    model_time_list = [model_time_list]

# Get some constants
one_sec        = 1. / 86400.  # one second in fractional days
ens_set        = np.arange(1, nens + 1)  # the set of ensembles
model_time_ind = len(model_time_list)  # the ith model step

# Parse the PFLOTRAN variables to be updated/analyzed/assimilated
p = re.compile('[A-Z_]+')
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
# dart_prior_file_set, dart_posterior_file_set = [], []
dart_prior_file_set = [
    # re.sub(r"\[ENS\]", str(ens) + "_time" + str(model_time_ind), dart_prior_file)
    re.sub(r"\[ENS\]", str(ens).zfill(ndigit_ens), dart_prior_file)
    for ens in ens_set
]
dart_prior_file_set = [
    re.sub(r"\[TIME\]", str(model_time_ind).zfill(ndigit_time), dart_prior_file_each)
    for dart_prior_file_each in dart_prior_file_set
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
    for fname in dart_prior_file_set:
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

    nc_fname_prior     = dart_prior_file_set[i]
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
    # Read the model output in HDF
    # The following information is expected to be read from the model HDF output
    #   - the last time step/time unit
    #   - the grids information
    #   - the state parameter values required in pflotran_var_set
    ###############################
    f_out = h5py.File(pf_fname, "r")

    # Get the grids/coordinates of the domain
    coordinates = f_out["Coordinates"]
    # The following options are choosing the center of the grid
    # x_set = coordinates['X [m]'][:]
    # y_set = coordinates['Y [m]'][:]
    # z_set = coordinates['Z [m]'][:]
    # x_set = [(a + b) / 2 for a, b in zip(x_set[:-1], x_set[1:])]
    # y_set = [(a + b) / 2 for a, b in zip(y_set[:-1], y_set[1:])]
    # z_set = [(a + b) / 2 for a, b in zip(z_set[:-1], z_set[1:])]
    # Choose the left boundary of the grid
    x_set = coordinates['X [m]'][:-1]
    y_set = coordinates['Y [m]'][:-1]
    z_set = coordinates['Z [m]'][:-1]
    nx, ny, nz = len(x_set), len(y_set), len(z_set)

    ## Get the last time step
    # time_set_o = [t for t in list(f_out.keys()) if t.startswith("Time")]
    # time_set   = [t.split()[1:] for t in time_set_o]
    # time_vset  = [float(t[0]) for t in time_set]

    # last_time, last_time_ind = np.max(time_vset), np.argmax(time_vset)
    # time_unit, last_time_o   = time_set[last_time_ind][1], time_set_o[last_time_ind]

    # Get all the time steps
    time_set_o = np.array([t for t in list(f_out.keys()) if t.startswith("Time")])
    time_set   = np.array([t.split()[1:] for t in time_set_o])
    time_vset  = np.array([float(t[0]) for t in time_set])
    time_unit  = time_set[0][1]

    # Shift the time_vset by the model spinup time
    time_vset = time_vset - spinup_time * 86400

    # Get the time steps within the assimilation window
    # print(time_vset, start_obs_sec, end_obs_sec)
    time_set_assim_ind = (time_vset > start_obs_sec) & (time_vset <= end_obs_sec)
    time_vset_assim    = time_vset[time_set_assim_ind]
    time_set_o_assim   = time_set_o[time_set_assim_ind]
    time_set_assim     = time_set[time_set_assim_ind]

    if time_unit in ["s", "sec", "second"]: # Convert from seconds to fractional days
        time_vset_assim_day  = time_vset_assim / 86400. 

    ntime = len(time_vset_assim)

    # Initialize the dart_var_dict
    for varn in pflotran_var_set:
        dart_var_dict[varn] = {"value": np.zeros([1,nz,ntime,nx]), 
                               "unit": ""}

    # Get the state/parameter/variable values required in pflotran_var_set
    for j in range(ntime):
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
            # pl_out_var_dict[varn] = {"unit": varunit, "original_name": v}
            # Check if the variable is required by pflotran_var_set
            if varn in pflotran_var_set:
                # dart_var_dict[varn]["value"].append(dataset[v][:]) 
                # TODO: make sure the shapes between dataset and dart_var_dict[varn]["value"] are compatible.
                dart_var_dict[varn]["value"][0,:,j,0] = dataset[v][:] 
                dart_var_dict[varn]["unit"] = varunit 
                # if time_vset_assim[j] == 300 and ens == 1:
                #     print(time_vset_assim[j], varn, v)
                #     print(time_o)
                #     print(dataset[v][:])
    
    f_out.close()

    ###############################
    # Read the parameter.h5
    # The following information is expected to be read from the model HDF output
    #   - the rest of the state parameter values required in pflotran_var_set
    # Note that parameter.h5 has a different structure compared with pflotran output
    # TODO a more strict structure is required. For now, we assume only one value is
    # needed or the whole domain is isotropic.
    ###############################
    for varn in para_pflotran_set:
        if varn in pflotran_var_set:
            value = f_para[varn][:][i]
            dart_var_dict[varn] = {
                # "value": value * np.ones([nz, ny, nx]),
                "value": value * np.ones([1, nz, ntime, nx]),
                "unit": ""
            }

    ###############################
    # Write the prior to NetCDF file
    ###############################
    # Add global attributes
    root_nc_prior.model_time          = model_time
    root_nc_prior.time_unit           = "day"
    root_nc_prior.assimilation_window = assim_window
    root_nc_prior.description         = "PFLOTRAN output/DART prior data"

    # Create the dimensions
    xloc_d   = root_nc_prior.createDimension('x_location', nx)
    # TODO: This is a fake dimension
    # yloc_d   = root_nc_prior.createDimension('y_location', ny)
    yloc_d   = root_nc_prior.createDimension('y_location', ntime)
    # TODO: This is a fake dimension
    zloc_d   = root_nc_prior.createDimension('z_location', nz)
    time_d   = root_nc_prior.createDimension('time', None)
    member_d = root_nc_prior.createDimension('member', 1)

    # Create the variables
    time   = root_nc_prior.createVariable('time', 'f8', ('time', ))
    member = root_nc_prior.createVariable('member', 'i8', ('member', ))
    xloc   = root_nc_prior.createVariable('x_location', 'f8', ('x_location', ))
    yloc   = root_nc_prior.createVariable('y_location', 'f8', ('y_location', ))
    zloc   = root_nc_prior.createVariable('z_location', 'f8', ('z_location', ))

    # Get the center of the assimilation window in unit day, as required by DART's read_model_time() subroutine
    time.units    = "day"
    # time[:]       = model_time + (assim_window - one_sec) / 2.
    # time[:]       = model_time + (assim_window + one_sec) / 2.
    time[:]       = model_time
    time.calendar = 'none'
    member[:]     = ens

    member.type, time.type = 'dimension_value', 'dimension_value'
    yloc.units , yloc.type = 'm',               'dimension_value'
    zloc.units , zloc.type = 'm',               'dimension_value'
    xloc.units , xloc.type = 'm',               'dimension_value'

    # Write coordinates values
    # xloc[:], yloc[:], zloc[:] = x_set, y_set, z_set
    # TODO: This is a fake dimension
    xloc[:], zloc[:] = x_set, z_set
    # yloc[:] = time_vset_assim_day
    yloc[:] = time_vset_assim
    # TODO: This is a fake dimension

    # Write the values to the variables
    for varn in pflotran_var_set:
        if dart_var_dict[varn] is None:
            print(
                "%s is not available in both PFLOTRAN output and parameter.h5"
                % varn)
            continue
        vargrp = root_nc_prior.createVariable(
            varn, 'f8', ('time', 'z_location', 'y_location', 'x_location'))
        vargrp.type = 'observation_value'
        vargrp.unit = dart_var_dict[varn]["unit"]
        # print(dart_var_dict[varn]["value"].shape)
        # print(vargrp[:].shape)
        vargrp[:] = dart_var_dict[varn]["value"]

        # if ens == 1:
        #     ind300 = list(time_vset_assim).index(300)
        #     print("Let's check nc file...")
        #     print(vargrp[:].shape)
        #     print(dart_var_dict[varn]["value"].shape)
        #     print(dart_var_dict[varn]["value"][:,:,ind300,:])
        #     print(vargrp[:][:,:,ind300,:])

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
    xloc_d   = root_nc_posterior.createDimension('x_location', nx)
    # TODO: This is a fake dimension
    # yloc_d   = root_nc_posterior.createDimension('y_location', ny)
    yloc_d   = root_nc_posterior.createDimension('y_location', ntime)
    # TODO: This is a fake dimension
    zloc_d   = root_nc_posterior.createDimension('z_location', nz)
    time_d   = root_nc_posterior.createDimension('time', None)
    member_d = root_nc_posterior.createDimension('member', 1)

    # Create the variables
    time   = root_nc_posterior.createVariable('time', 'f8', ('time', ))
    member = root_nc_posterior.createVariable('member', 'i8', ('member', ))
    xloc   = root_nc_posterior.createVariable('x_location', 'f8', ('x_location', ))
    yloc   = root_nc_posterior.createVariable('y_location', 'f8', ('y_location', ))
    zloc   = root_nc_posterior.createVariable('z_location', 'f8', ('z_location', ))

    # Get the the center of the assimilation window in unit day, as required by DART's read_model_time() subroutine
    time.units = "day"
    # time[:]    = model_time + (assim_window - one_sec) / 2.
    time[:]       = model_time
    time.calendar = 'none'
    member[:]     = ens

    member.type, time.type = 'dimension_value', 'dimension_value'
    yloc.units , yloc.type = 'm',               'dimension_value'
    zloc.units , zloc.type = 'm',               'dimension_value'
    xloc.units , xloc.type = 'm',               'dimension_value'

    # Write coordinates values
    # xloc[:], yloc[:], zloc[:] = x_set, y_set, z_set
    # TODO: This is a fake dimension
    xloc[:], zloc[:] = x_set, z_set
    # yloc[:] = time_vset_assim_day
    yloc[:] = time_vset_assim
    # TODO: This is a fake dimension

    # Write the variables without values
    for varn in pflotran_var_set:
        if dart_var_dict[varn] is None:
            print(
                "%s is not available in both PFLOTRAN output and parameter.h5"
                % varn)
            continue
        vargrp = root_nc_posterior.createVariable(
            varn, 'f8', ('time', 'z_location', 'y_location', 'x_location'))
        vargrp.type = 'observation_value'
        vargrp.unit = dart_var_dict[varn]["unit"]

    root_nc_posterior.close()

    ###############################
    # Copy the first ensemble to the prior template
    ###############################
    if i == 1:
        shutil.copyfile(nc_fname_prior, dart_prior_template)

f_para.close()
