"""This file is used for parsing different PFLOTRAN parameter or output files."""
# Author: Peishi Jiang

import os
import re
import h5py
import f90nml
import glob
import json
import numpy as np
import pandas as pd
from scipy.stats import truncnorm


class pflotran_files:

    def __init__(self, config_nml):
        self.configs = f90nml.read(config_nml)
        # Parameters and observations configuration
        self.pflotran_out_file         = self.configs["file_cfg"]["pflotran_out_file"]
        self.pflotran_obs_file         = self.configs["file_cfg"]["pflotran_obs_file"]
        self.use_obs_tecfile_for_prior = self.configs["file_cfg"]["use_obs_tecfile_for_prior"]
        self.pflotran_para_file        = self.configs["file_cfg"]["pflotran_para_file"]
        self.para_reconditioned_cell_file = self.configs["file_cfg"]["pflotran_reconditioned_cell_file"]
        self.material_id_file          = self.configs["file_cfg"]["material_id_file"]
        self.model_time                = float(self.configs["time_cfg"]["current_model_time"])  # days
        self.model_time_list           = self.configs["time_cfg"]["model_time_list"]
        self.nens                      = self.configs["da_cfg"]["nens"]
        self.spinup_time               = self.configs["time_cfg"]["spinup_length"]

        self.obs_set                    = self.configs["obspara_set_cfg"]["obs_set"]
        self.obs_pflotran_set           = self.configs["obspara_set_cfg"]["obs_pflotran_set"]
        self.para_set                   = self.configs["obspara_set_cfg"]["para_set"]
        self.para_take_log_set          = self.configs["obspara_set_cfg"]["para_take_log_set"]
        self.para_homogeneous           = self.configs["obspara_set_cfg"]["para_homogeneous"]
        self.para_material_id_set       = self.configs["obspara_set_cfg"]["para_material_id_set"]
        self.para_hdf_dataset_name_set  = self.configs["obspara_set_cfg"]["para_hdf_dataset_name_set"]
        self.para_isotropic_set         = self.configs["obspara_set_cfg"]["para_isotropic_set"]
        self.para_anisotropic_ratio_set = self.configs["obspara_set_cfg"]["para_anisotropic_ratio_set"]
        self.para_min_set               = self.configs["obspara_set_cfg"]["para_min_set"]
        self.para_max_set               = self.configs["obspara_set_cfg"]["para_max_set"]
        self.para_mean_set              = self.configs["obspara_set_cfg"]["para_mean_set"]
        self.para_std_set               = self.configs["obspara_set_cfg"]["para_std_set"]
        # self.para_dist_set              = self.configs["obspara_set_cfg"]["para_dist_set"]
        # self.rescaled_set               = self.configs["obspara_set_cfg"]["para_prior_rescaled_set"]
        self.para_resampled_set         = self.configs["obspara_set_cfg"]["para_resampled_set"]
        self.para_sample_method_set     = self.configs["obspara_set_cfg"]["para_sample_method_set"]

        self.ens_set                   = np.arange(1, self.nens + 1)  # the set of ensembles
        if not isinstance(self.para_material_id_set, list):
            self.para_material_id_set = [self.para_material_id_set]
        if not isinstance(self.model_time_list, list):
            self.model_time_list = [self.model_time_list]


    ###############################
    # Read PFLOTRAN parameters
    ###############################
    def read_pflotran_parameter(self, missing_value):
        para_homogeneous = self.para_homogeneous
        if not para_homogeneous:
            return self.read_pflotran_heterogeneous_parameter_h5(missing_value=missing_value)
        else:
            return self.read_pflotran_homogeneous_parameter_h5(missing_value=missing_value)


    ###############################
    # Read PFLOTRAN heterogeneous parameters in HDF5
    ###############################
    def read_pflotran_heterogeneous_parameter_h5(self, missing_value):
        # Configurations
        nens                      = self.nens
        para_set                  = self.para_set
        para_material_id_set      = self.para_material_id_set
        para_hdf_dataset_name_set = self.para_hdf_dataset_name_set
        para_take_log_set         = self.para_take_log_set
        material_id_file          = self.material_id_file
        pflotran_para_file        = self.pflotran_para_file
        model_time                = self.model_time

        pflotran_para_dict = dict.fromkeys(para_set)

        # Obtain the material indices and the locations if heterogeneous parameters are required
        with h5py.File(material_id_file, "r") as h5file:
            # Get the locations
            # x_loc_all = h5file['Domain']['Vertices'][:][:,0]
            # y_loc_all = h5file['Domain']['Vertices'][:][:,1]
            # z_loc_all = h5file['Domain']['Vertices'][:][:,2]
            x_loc_all, y_loc_all, z_loc_all = \
                parse_unstructured_hdf5(h5file['Domain']['Cells'][:],h5file['Domain']['Vertices'][:])
            # Get the cell ids
            cell_ids = h5file['Materials']['Cell Ids'][:]
            # Get the material ids
            material_ids = h5file['Materials']['Material Ids'][:]
            # Get the locations of the parameters to be assimilated
            material_ids_para = material_ids[np.isin(material_ids, para_material_id_set)]
            cell_ids_para = cell_ids[np.isin(material_ids, para_material_id_set)]
            x_loc_para = x_loc_all[cell_ids_para-1]
            y_loc_para = y_loc_all[cell_ids_para-1]
            z_loc_para = z_loc_all[cell_ids_para-1]
            nloc_para = len(x_loc_para)

        # Load the prior ensemble of the parameter values
        with h5py.File(pflotran_para_file, "r") as h5file:
            # Get the cell ids
            # cell_ids_2 = list(h5file['Cell Ids'][:])
            cell_ids_2 = h5file['Cell Ids'][:]
            # Get the value for each parameter
            for i in range(len(para_set)):
                para_varn             = para_set[i]
                para_material_id      = para_material_id_set[i]
                para_hdf_dataset_name = para_hdf_dataset_name_set[i]
                para_take_log         = para_take_log_set[i]

                cell_ids_para_each = cell_ids_para[material_ids_para == para_material_id]
                cell_ids_para_each_index = cell_ids_2.searchsorted(cell_ids_para_each)
                # cell_ids_para_each_index = [cell_ids_2.index(cell) for cell in cell_ids_para_each]

                values = np.zeros([nens, 1, nloc_para]) + missing_value
                for j in range(nens):
                    if para_hdf_dataset_name+str(j+1) in h5file.keys():
                        para_values = h5file[para_hdf_dataset_name+str(j+1)][:]
                        para_values = para_values[cell_ids_para_each_index]
                    elif para_hdf_dataset_name+'X'+str(j+1) in h5file.keys():
                        para_values = h5file[para_hdf_dataset_name+'X'+str(j+1)][:]
                        para_values = para_values[cell_ids_para_each_index]
                    else:
                        raise Exception("Unknown parameter HDF dataset name: {}".format(para_hdf_dataset_name))
                    # Convert the logarithmic value if necessary
                    if para_take_log:
                        values[j, 0, material_ids_para==para_material_id] = np.log10(para_values)
                    else:
                        values[j, 0, material_ids_para==para_material_id] = para_values
                # Save it
                pflotran_para_dict[para_varn] = {"value": values, "unit":""}

        # Time for the parameter
        ntime_para = 1
        time_para  = [model_time]

        return pflotran_para_dict, ntime_para, time_para, nloc_para, x_loc_para, y_loc_para, z_loc_para, cell_ids_para


    ###############################
    # Read PFLOTRAN homogeneous parameters in HDF5
    ###############################
    def read_pflotran_homogeneous_parameter_h5(self, missing_value):
        # Configurations
        nens                      = self.nens
        para_set                  = self.para_set
        para_hdf_dataset_name_set = self.para_hdf_dataset_name_set
        para_take_log_set         = self.para_take_log_set
        pflotran_para_file        = self.pflotran_para_file
        model_time                = self.model_time

        pflotran_para_dict = dict.fromkeys(para_set)

        nloc_para = 1
        # Load the prior ensemble of the parameter values
        with h5py.File(pflotran_para_file, "r") as h5file:
            for i in range(len(para_set)):
                para_varn             = para_set[i]
                para_hdf_dataset_name = para_hdf_dataset_name_set[i]
                para_take_log         = para_take_log_set[i]

                values = np.zeros([nens, 0, nloc_para]) + missing_value
                if para_hdf_dataset_name in h5file.keys():
                    para_values = h5file[para_hdf_dataset_name][:]
                else:
                    raise Exception("Unknown parameter HDF dataset name: {}".format(para_hdf_dataset_name))
                # Conver the logarithmic value if necessary
                if para_take_log:
                    values[:, 0, i] = np.log10(para_values)
                else:
                    values[:, 0, i] = para_values
                # Save it
                pflotran_para_dict[para_varn] = {"value": values, "unit":""}

        # Time for the parameter
        ntime_para = 1
        time_para  = [model_time]

        # TODO: perhaps cell id is not 0
        cell_ids = 0

        return pflotran_para_dict, ntime_para, time_para, nloc_para, cell_ids


    ###############################
    # Read PFLOTRAN outputs
    ###############################
    def read_pflotran_output(self, start_time, end_time, missing_value):
        # Configuration
        use_obs_tecfile_for_prior = self.use_obs_tecfile_for_prior

        if use_obs_tecfile_for_prior:
            return self.read_pflotran_obs_output_tec(start_time, end_time, missing_value)
        else:
            return self.read_pflotran_output_h5(start_time, end_time, missing_value)

        # if pflotran_output_file.endswith('.hdf5'):
        #     return self.read_pflotran_output_h5(start_time, end_time, missing_value)
        # elif pflotran_output_file.endswith('.tec'):
        #     return self.read_pflotran_obs_output_tec(start_time, end_time, missing_value)
        # else:
        #     raise Exception("Unknown pflotran output file format: {}".format(pflotran_output_file))


    ###############################
    # Read PFLOTRAN outputs in HDF5
    ###############################
    def read_pflotran_output_h5(self, start_time, end_time, missing_value):
        # Configuration
        pflotran_out_file = self.pflotran_out_file
        material_id_file  = self.material_id_file
        obs_pflotran_set  = self.obs_pflotran_set
        ens_set           = self.ens_set
        nens              = self.nens

        # Get the file names of all ensembles for PFLOTRAN output
        pflotran_out_file_set = [
            re.sub(r"\[ENS\]", str(ens), pflotran_out_file) for ens in ens_set
        ]
        # Check the existences of these files
        for f in pflotran_out_file_set:
            if not os.path.isfile(f):
                raise Exception("The PFLOTRAN output file %s does not exits!" % f)

        pflotran_state_dict = dict.fromkeys(obs_pflotran_set)

        # Get the model states
        for i in range(nens):
            pf_fname = pflotran_out_file_set[i]

            # Parse the output state files
            with h5py.File(pf_fname, "r") as f_out:
                # Get the grids/coordinates of the domain
                if "Coordinates" in f_out.keys():
                    grid_type = "structured"
                    coordinates = f_out["Coordinates"]
                    # TODO: check whether PFLOTRAN generates outputs on the boundaries or the centers of cells.
                    # The following choose the left boundaries
                    x_loc_state = coordinates['X [m]'][:-1]
                    y_loc_state = coordinates['Y [m]'][:-1]
                    z_loc_state = coordinates['Z [m]'][:-1]
                    nx_state, ny_state, nz_state = len(x_loc_state), len(y_loc_state), len(z_loc_state)
                    nloc_state  = nx_state*ny_state*nz_state
                    x_loc_state, y_loc_state, z_loc_state = np.meshgrid(x_loc_state, y_loc_state, z_loc_state)
                    x_loc_state, y_loc_state, z_loc_state = x_loc_state.flatten(), y_loc_state.flatten(), z_loc_state.flatten()

                elif "Domain" in f_out.keys():
                    grid_type = "unstructured"
                    domain = f_out["Domain"]
                    x_loc_state = domain['XC']
                    y_loc_state = domain['YC']
                    z_loc_state = domain['ZC']
                    nloc_state  = len(x_loc_state)

                # Get all the time steps
                time_set_o = np.array([t for t in list(f_out.keys()) if t.startswith("Time")])
                time_set   = np.array([t.split()[1:] for t in time_set_o])
                time_vset  = np.array([float(t[0]) for t in time_set])
                time_unit  = time_set[0][1]

                # Get the time steps within the assimilation window
                # print(time_vset, start_obs_sec, end_obs_sec)
                time_set_assim_ind = (time_vset > start_time) & (time_vset <= end_time)
                time_vset_assim    = time_vset[time_set_assim_ind]
                time_set_o_assim   = time_set_o[time_set_assim_ind]
                # time_set_assim     = time_set[time_set_assim_ind]

                if time_unit in ["s", "sec", "second"]: # Convert from seconds to fractional days
                    time_vset_assim_day  = time_vset_assim / 86400.

                ntime_state = len(time_vset_assim)

                # Initialize the dart_var_dict
                # for varn in pflotran_var_set:
                if i == 0:
                    for varn in obs_pflotran_set:
                        pflotran_state_dict[varn] = {"value": np.zeros([nens, ntime_state, nloc_state]),
                                                     "unit": ""}

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
                        if varn in obs_pflotran_set:
                            # dart_var_dict[varn]["value"].append(dataset[v][:])
                            # TODO: make sure flatten() follows the same order as np.meshgrid
                            if grid_type == 'structured':
                                pflotran_state_dict[varn]["value"][i,j,:] = dataset[v][:].flatten()
                            elif grid_type == 'unstructured':
                                pflotran_state_dict[varn]["value"][i,j,:] = dataset[v][:]
                            pflotran_state_dict[varn]["unit"] = varunit

                    # Check whether all the states are obtained (do it once)
                    if i == 0:
                        for varn in obs_pflotran_set:
                            if pflotran_state_dict[varn] is None:
                                raise Exception("The following variable is not found in PFLOTRAN output: {}".format(varn))

        # Get cell ids
        if os.path.isfile(material_id_file):
            with h5py.File(material_id_file, "r") as h5file:
                cell_ids = list(h5file['Materials']['cell_ids'][:])
        else:
            cell_ids = np.arange(1, nloc_state+1)

        # Times steps for the model states
        time_state  = time_vset_assim_day
        ntime_state = len(time_state)

        return pflotran_state_dict, ntime_state, time_state, nloc_state, x_loc_state, y_loc_state, z_loc_state, cell_ids


    ###############################
    # Read PFLOTRAN outputs at observation points in HDF5
    ###############################
    def read_pflotran_obs_output_tec(self, start_time, end_time, missing_value):
        # Configuration
        pflotran_obs_file = self.pflotran_obs_file
        material_id_file  = self.material_id_file
        obs_pflotran_set  = self.obs_pflotran_set
        ens_set           = self.ens_set
        nens              = self.nens

        # Obtain the spatial locations if heterogeneous parameters are required
        with h5py.File(material_id_file, "r") as h5file:
            # cell_ids = list(h5file['Materials']['Cell Ids'][:])
            cell_ids = h5file['Materials']['Cell Ids'][:]
            x_loc_all, y_loc_all, z_loc_all = \
                parse_unstructured_hdf5(h5file['Domain']['Cells'][:],h5file['Domain']['Vertices'][:])

        # Get the file name format of all ensembles for PFLOTRAN output
        pflotran_obs_file_set = [
            re.sub(r"\[ENS\]", str(ens), pflotran_obs_file) for ens in ens_set
        ]
        pflotran_obs_file_set = [
            re.sub(r"\[ANY\]", "*", each_file) for each_file in pflotran_obs_file_set
        ]

        pflotran_state_dict = dict.fromkeys(obs_pflotran_set)

        # Read the data from each ensemble member
        for i in range(nens):
            pf_fname_fmt = pflotran_obs_file_set[i]

            # Get all the files with the same file name format
            pf_fname_list = glob.glob(pf_fname_fmt)

            # Check the existences of these files
            if len(pf_fname_list) == 0:
                raise Exception("Unknown PFLOTRAN output file name format: {}".format(pf_fname_fmt))

            # Read the tec files
            for j in range(len(pf_fname_list)):
                pf_fname = pf_fname_list[j]
                if i == 0:
                    tec_data_old, state_names, state_unit, time_state, time_unit = read_tec_file(pf_fname)
                else:
                    tec_data_old, state_names, state_unit, _, _ = read_tec_file(pf_fname)
                if j == 0:
                    tec_data = tec_data_old
                else:
                    tec_data = tec_data.append(tec_data_old)

            # Delete duplicates (it is possible that some wells are interpolated to the same cell)
            tec_data = tec_data.drop_duplicates()

            # Get the time and space
            if i == 0:
                # ------------ time -------------
                # Convert the time to the units of day
                if time_unit.lower() in ["h", "hr"]:
                    time_state = time_state / 24.
                    time_conversion = 24.
                elif time_unit.lower() in ["s", "sec"]:
                    time_state = time_state / 86400.
                    time_conversion = 86400.
                elif time_unit.lower() not in ["d", "day"]:
                    raise Exception("Unknown time unit: {}".format(time_unit))
                time_state = time_state[(time_state>=start_time) & (time_state<=end_time)]
                # Get the number of time steps
                ntime_state = len(time_state)
                # ------------ space -------------
                cell_ids_tec = tec_data["cell_id"].unique()
                # cell_ids_indices = [cell_ids.index(cell) for cell in cell_ids_tec]
                cell_ids_indices = cell_ids.searchsorted(cell_ids_tec)
                x_loc_state  = x_loc_all[cell_ids_indices]
                y_loc_state  = y_loc_all[cell_ids_indices]
                z_loc_state  = z_loc_all[cell_ids_indices]
                nloc_state   = len(cell_ids_tec)
            # print(time_state, start_time, end_time, i)

            # Get the data within the start and end time
            tec_data = tec_data[(tec_data['Time'] >= start_time*time_conversion) & (tec_data['Time'] <= end_time*time_conversion)]

            # Initialize pflotran_state_dict
            if i == 0:
                for varn in obs_pflotran_set:
                    # pflotran_state_dict[varn] = {"value": np.zeros([nens, nloc_state, ntime_state]) + missing_value, "unit":""}
                    pflotran_state_dict[varn] = {"value": np.zeros([nens, ntime_state, nloc_state]) + missing_value, "unit":""}

            # Check the shape
            if tec_data.shape[0] != ntime_state*nloc_state:
                print(os.getcwd())
                # tec_data.to_csv("tec_data_wrong.csv")
                raise Exception("The number of value {} does not match with the shape ({}, {})".format(tec_data.shape[0],nloc_state,ntime_state))

            # Get the required model states
            for varn in obs_pflotran_set:
                if varn in state_names:   # variable name exactly matches the required observaiton variable
                    tec_data_varn = tec_data[['Time', varn, 'cell_id']]
                    # TODO: reorganize the data
                    tec_data_varn_values = tec_data_varn[varn].values.reshape(nloc_state,ntime_state).T
                    # save it
                    pflotran_state_dict[varn]["value"][i,:,:] = tec_data_varn_values
                    pflotran_state_dict[varn]["unit"] = state_unit[state_names.index(varn)]
                elif varn == "WATER_LEVEL":
                    tec_data_varn = tec_data[['Time', 'LIQUID_PRESSURE', 'cell_id']]
                    # cell_id_varn_indices = [cell_ids.index(cell) for cell in tec_data_varn['cell_id']]
                    cell_id_varn_indices = cell_ids.searchsorted(tec_data_varn['cell_id'])
                    z_loc_varn = z_loc_all[cell_id_varn_indices]
                    # Define some constants
                    rho, g, patm = 1000., 9.8, 101325.
                    # convert liquid pressure to water level
                    tec_data_varn["WATER_LEVEL"] =(tec_data_varn['LIQUID_PRESSURE'] - patm)/(rho*g) + z_loc_varn
                    # TODO: reorganize the data
                    tec_data_varn_values = tec_data_varn["WATER_LEVEL"].values.reshape(nloc_state,ntime_state).T
                    # save it
                    pflotran_state_dict[varn]["value"][i,:,:] = tec_data_varn_values
                    # pflotran_state_dict[varn]["unit"] = state_unit[state_names.index('LIQUID_PRESSURE')]
                    pflotran_state_dict[varn]["unit"] = "m"
                elif varn == "GROUNDWATER_TRACER":
                    tec_data_varn = tec_data[['Time', 'TOTAL_TRACER_RIVER', "QLX", "QLY", 'cell_id']]
                    # TODO: convert total river tracer to groundwater tracer
                    tec_data_varn["TOTAL_TRACER_RIVER"][tec_data_varn["TOTAL_TRACER_RIVER"] > 1] = 1
                    tec_data_varn["GROUNDWATER_TRACER"] = 1 - tec_data_varn["TOTAL_TRACER_RIVER"]
                    # TODO: reorganize the data
                    tec_data_varn_values = tec_data_varn["GROUNDWATER_TRACER"].values.reshape(nloc_state,ntime_state).T
                    # save it
                    pflotran_state_dict[varn]["value"][i,:,:] = tec_data_varn_values
                    pflotran_state_dict[varn]["unit"] = state_unit[state_names.index('TOTAL_TRACER_RIVER')]
                else:
                    raise Exception("The following variable is not found in PFLOTRAN state file: {}".format(varn))

        return pflotran_state_dict, ntime_state, time_state, nloc_state, x_loc_state, y_loc_state, z_loc_state, cell_ids_tec


    ###############################
    # Update PFLOTRAN parameters
    ###############################
    def update_parameters_in_pflotran_parameter_file(self, para_prior_set, original_para_posterior_set):
        # Configurations
        nens                       = self.nens
        para_set                   = self.para_set
        para_material_id_set       = self.para_material_id_set
        para_hdf_dataset_name_set  = self.para_hdf_dataset_name_set
        para_take_log_set          = self.para_take_log_set
        material_id_file           = self.material_id_file
        pflotran_para_file         = self.pflotran_para_file
        model_time                 = self.model_time
        para_min_set               = self.para_min_set
        para_max_set               = self.para_max_set
        para_mean_set              = self.para_mean_set
        para_std_set               = self.para_std_set
        # para_dist_set              = self.para_dist_set
        para_sample_method_set     = self.para_sample_method_set
        # para_rescaled_set          = self.rescaled_set
        para_resampled_set         = self.para_resampled_set
        para_homogeneous           = self.para_homogeneous
        para_isotropic_set         = self.para_isotropic_set
        para_anisotropic_ratio_set = self.para_anisotropic_ratio_set
        model_time                 = self.model_time
        reconditioned_cellid_file  = self.para_reconditioned_cell_file
        # reconditioned_cellid_file_format  = self.para_reconditioned_cell_file

        # Get the MDA step
        try:
            update_obs_ens_posterior_now = self.configs["da_cfg"]["update_obs_ens_posterior_now"]
        except:
            update_obs_ens_posterior_now = False
        enks_mda_iteration_step        = self.configs["da_cfg"]["enks_mda_iteration_step"]

        # Update the parameters in PFLOTRAN parameter file
        para_posterior = dict.fromkeys(para_set)
        for i in range(len(para_set)):
            varn                         = para_set[i]
            para_prior_varn              = para_prior_set[varn]["value"]
            original_para_posterior_varn = original_para_posterior_set[varn]["value"]
            para_cell_ids_varn           = para_prior_set[varn]["cell_ids"]
            para_loc_set_varn            = para_prior_set[varn]["loc_set"]

            para_hdf_dataset_name        = para_hdf_dataset_name_set[i]
            para_material_id             = para_material_id_set[i]
            para_isotropic               = para_isotropic_set[i]
            para_anisotropic_ratio       = para_anisotropic_ratio_set[i]
            para_min                     = para_min_set[i]
            para_max                     = para_max_set[i]
            para_mean                    = para_mean_set[i]
            para_std                     = para_std_set[i]
            # para_dist                    = para_dist_set[i]
            para_sample_method           = para_sample_method_set[i]
            # para_rescaled                = para_rescaled_set[i]
            para_take_log                = para_take_log_set[i]
            para_resampled               = para_resampled_set[i]

            print("Updating {} in PFLOTRAN parameter file now ...".format(varn))

            # Correct the original posterior first
            if para_sample_method.lower() != "rescale":
                if para_min != -99999:
                    original_para_posterior_varn[original_para_posterior_varn < para_min] = para_min
                if para_max != 99999:
                    original_para_posterior_varn[original_para_posterior_varn > para_max] = para_max

            # When the parameters are homogeneous
            if para_homogeneous:
                # Get the posterior to be updated from the original posterior and prior
                # resample the prior if
                # (1) it is required and
                # (2) it is not the time for updating observation ensemble posterior
                # (3) it is not during ES-MDA iteration
                if para_resampled and not update_obs_ens_posterior_now and enks_mda_iteration_step == 1:
                    posterior_varn = update_para_prior_from_original_posterior_homogeneous(
                        para_prior_varn, original_para_posterior_varn,
                        para_sample_method, para_min, para_max, para_mean, para_std, para_take_log)
                else:
                    posterior_varn = original_para_posterior_varn

                # Convert from the logarithmic form
                if para_take_log:
                    posterior_varn = np.power(10, posterior_varn, dtype=float)

                # Assign the posterior to pflotran parameter file
                with h5py.File(pflotran_para_file, "r+") as h5file:
                    if para_hdf_dataset_name in h5file.keys():
                        h5file[para_hdf_dataset_name][:] = posterior_varn
                    else:
                        raise Exception("Unknown parameter HDF dataset name: {}".format(para_hdf_dataset_name))

            else:
                # # Obtain the material indices and the locations if heterogeneous parameters are required
                # with h5py.File(material_id_file, "r") as h5file:
                #     # Get the locations
                #     x_loc_all, y_loc_all, z_loc_all = \
                #         parse_unstructured_hdf5(h5file['Domain']['Cells'][:],h5file['Domain']['Vertices'][:])
                #     # Get the cell ids
                #     cell_ids = h5file['Materials']['Cell Ids'][:]
                #     # Get the material ids
                #     material_ids = h5file['Materials']['Material Ids'][:]
                #     # Get the cell ids of the parameter
                #     para_cell_ids_varn = cell_ids[material_ids == para_material_id]

                # Get the posterior to be updated from the original posterior and prior
                # resample the prior if
                # (1) it is required and
                # (2) it is not the time for updating observation ensemble posterior
                # (3) it is not during ES-MDA iteration
                # print(enks_mda_iteration_step)
                # raise Exception("Stop")
                # if para_resampled and not update_obs_ens_posterior_now:
                if para_resampled and not update_obs_ens_posterior_now and enks_mda_iteration_step == 1:
                    # reconditioned_cellid_file = re.sub(r"\[PARA\]", varn, reconditioned_cellid_file_format)
                    posterior_varn = update_para_prior_from_original_posterior_heterogeneous(
                        para_prior_varn, original_para_posterior_varn, varn,
                        para_loc_set_varn, para_cell_ids_varn, reconditioned_cellid_file, model_time,
                        para_sample_method, para_min, para_max, para_mean, para_std, para_take_log)
                else:
                    print("No sampling is required for updating parameters at this stage ...")
                    posterior_varn = original_para_posterior_varn

                # Quality control
                print("How many invalid values do we got here ...")
                print(para_min, para_max, posterior_varn.max(), posterior_varn.min())
                print(np.sum((posterior_varn<=para_min)|(posterior_varn>=para_max)))
                if np.sum((posterior_varn<=para_min)|(posterior_varn>=para_max)) != 0:
                    print("Correct those invalid values")
                    # posterior_varn[posterior_varn<=para_min] = para_min
                    # posterior_varn[posterior_varn>=para_max] = para_max
                    posterior_varn[posterior_varn<=para_min] = original_para_posterior_varn.mean()
                    posterior_varn[posterior_varn>=para_max] = original_para_posterior_varn.mean()
                # if np.sum(posterior_varn == np.inf) != 0:

                # Convert from the logarithmic form
                posterior_varn_o = np.copy(posterior_varn)
                if para_take_log:
                    posterior_varn = np.power(10, posterior_varn, dtype=float)
                
                # print("How many improper values do we got here ...")
                # print(np.sum((posterior_varn<=1e-30)|(posterior_varn>1)))
                # if np.sum((posterior_varn<=1e-30)|(posterior_varn>1)) != 0:
                # if np.sum(posterior_varn == np.inf) != 0:
                    # print(posterior_varn_o[(posterior_varn<=1e-30)|(posterior_varn>1)])
                    # print(original_para_posterior_varn[(posterior_varn<=1e-30)|(posterior_varn>1)])
                    # print(original_para_posterior_varn.max(), original_para_posterior_varn.min())
                    # raise Exception('Stop for check!')

                # Assign the posterior to pflotran parameter file
                with h5py.File(pflotran_para_file, "r+") as h5file:
                    # Get the cell ids
                    # cell_ids_2 = list(h5file['Cell Ids'][:])
                    cell_ids_2 = h5file['Cell Ids'][:]
                    # cell_ids_para_index = [cell_ids_2.index(cell) for cell in cell_ids_para]
                    cell_ids_para_index = cell_ids_2.searchsorted(para_cell_ids_varn)
                    cell_ids_para_index = list(cell_ids_para_index)
                    # Update
                    for j in range(nens):
                        if para_hdf_dataset_name+str(j+1) in h5file.keys():
                            # h5file[para_hdf_dataset_name+str(j+1)][cell_ids_para_index] = posterior_varn[j,:]
                            para_from_h5file = h5file[para_hdf_dataset_name+str(j+1)][:]
                            para_from_h5file[cell_ids_para_index] = posterior_varn[j,:]
                            h5file[para_hdf_dataset_name+str(j+1)][:] = para_from_h5file
                        elif para_hdf_dataset_name+'X'+str(j+1) in h5file.keys():
                            # h5file[para_hdf_dataset_name+'X'+str(j+1)][cell_ids_para_index] = posterior_varn[j,:]
                            # h5file[para_hdf_dataset_name+'Y'+str(j+1)][cell_ids_para_index] = posterior_varn[j,:]
                            # if para_isotropic:
                            #     h5file[para_hdf_dataset_name+'Z'+str(j+1)][cell_ids_para_index] = posterior_varn[j,:]
                            # else:
                            #     h5file[para_hdf_dataset_name+'Z'+str(j+1)][cell_ids_para_index] = posterior_varn[j,:] * para_anisotropic_ratio
                            # X
                            para_from_h5file = h5file[para_hdf_dataset_name+'X'+str(j+1)][:]
                            para_from_h5file[cell_ids_para_index] = posterior_varn[j,:]
                            h5file[para_hdf_dataset_name+'X'+str(j+1)][:] = para_from_h5file
                            # Y
                            para_from_h5file = h5file[para_hdf_dataset_name+'Y'+str(j+1)][:]
                            para_from_h5file[cell_ids_para_index] = posterior_varn[j,:]
                            h5file[para_hdf_dataset_name+'Y'+str(j+1)][:] = para_from_h5file
                            # Z
                            para_from_h5file = h5file[para_hdf_dataset_name+'Z'+str(j+1)][:]
                            if para_isotropic:
                                para_from_h5file[cell_ids_para_index] = posterior_varn[j,:]
                            else:
                                para_from_h5file[cell_ids_para_index] = posterior_varn[j,:] * para_anisotropic_ratio
                            h5file[para_hdf_dataset_name+'Z'+str(j+1)][:] = para_from_h5file
                        else:
                            raise Exception("Unknown parameter HDF dataset name: {}".format(para_hdf_dataset_name))
                
        # raise Exception('Stop for check!')



###############################
# Function for parsing unstructured domain in HDF5
###############################
def parse_unstructured_hdf5(cells, vertices):
    ncell, _ = cells.shape
    x_loc_all = [vertices[cells[i,1:]-1,0].mean() for i in range(ncell)]
    y_loc_all = [vertices[cells[i,1:]-1,1].mean() for i in range(ncell)]
    z_loc_all = [vertices[cells[i,1:]-1,2].mean() for i in range(ncell)]

    x_loc_all = np.array(x_loc_all)
    y_loc_all = np.array(y_loc_all)
    z_loc_all = np.array(z_loc_all)

    return x_loc_all, y_loc_all, z_loc_all


###############################
# Function for read a tec file
###############################
def read_tec_file(f_tec_file):
    # Get the data and head
    data_all = np.loadtxt(f_tec_file, skiprows=1)
    data, time = data_all[:,1:], data_all[:,0]
    time = np.expand_dims(time, axis=1)
    columns = pd.read_csv(f_tec_file, nrows=0).keys()
    var_columns = columns[1:]

    # Get the number of datapoints
    N, _ = data.shape

    # Split the columns
    var_keys = [c.split("Well")[0] for c in var_columns]
    well_levels = ["Well"+c.split("Well")[1] for c in var_columns]
    # well_levels = [w.strip() for w in well_levels]
    well_levels_unique = set(well_levels)

    # Get the variable names and units
    var_keys_unique = list(np.unique(var_keys))
    var_names = [key.split("[")[0].split() for key in var_keys_unique]
    var_names = ["_".join(key).upper() for key in var_names]
    # var_names = [key.split("[")[1].split() for key in var_keys_unique]
    var_units = [re.match(r"^.*\[(.*)\].*$", key).group(1) if "[" in key else "" for key in var_keys_unique]

    time_unit = re.match(r"^.*\[(.*)\].*$", columns[0]).group(1)

    # Create the keys
    keys = ["Time"] + var_names + \
           ["northing", "easting", "elevation", "cell_id"]
           # ["well_name", "northing", "easting", "elevation", "cell_id"]

    # Create an empty dataframe
    tec_data = pd.DataFrame(columns=keys)

    # For each level
    for level in well_levels_unique:

        # Get the corresponding column indices in data
        indices = [i for i, x in enumerate(well_levels) if x == level]

        # Get the data of the current level
        data_level = data[:,indices]

        # Combine data_level with time
        # data_level = np.hstack((time, data_level, np.empty([N,5])))
        data_level = np.hstack((time, data_level, np.empty([N,4])))

        # Put data_level to a dataframe
        pd_level = pd.DataFrame(data_level, columns=keys)

        # Get the northing, easting, elevation,
        # cell id, well name, and well level
        well_name, cell_id, east, north, ele = level.split(" ")
        # well_name = level_info.split("_")[1]
        # well_level = int(level_info.split("_")[-1])
        cell_id = float(cell_id[1:-1])
        east, north, ele = float(east[1:]), float(north), float(ele[:-1])

        # Creat new columns recording northing, easting, elevation,
        # cell id, well name, and well level
        # pd_level["well_name"] = well_name
        # pd_level["well_level"] = well_level
        pd_level["cell_id"] = cell_id
        pd_level["easting"] = east
        pd_level["northing"] = north
        pd_level["elevation"] = ele

#         print(pd_level.head())
        # Append the data of the current level to the previous level
        tec_data = tec_data.append(pd_level)

    return tec_data, var_names, var_units, time, time_unit


##############################################################
# Function for getting the posterior to be updated from the original posterior and prior
# (homogeneous parameter)
##############################################################
def update_para_prior_from_original_posterior_homogeneous(
    para_prior, original_para_posterior,
    para_sample_method, para_min, para_max, para_mean, para_std, para_take_log):

    # Get the number of ensemble member
    nens = len(original_para_posterior)

    # Get the mean and std from the original posterior
    para_mean_post = np.mean(original_para_posterior)
    para_std_post  = np.std(original_para_posterior)

    print("Sampling method -- {} -- is used ...".format(para_sample_method))
    # If from rescaling the original posterior parameter
    if para_sample_method.lower() == "rescale":
        posterior = (original_para_posterior - para_mean_post) / para_std_post * para_std + para_mean_post

    # Sampled from normal distribution
    elif para_sample_method.lower() == "normal":
        posterior = np.random.normal(para_mean_post, para_std, nens)

    # Sampled from lognormal distribution
    elif para_sample_method.lower() == "lognormal":
        if para_take_log:
            posterior  = np.random.normal(para_mean_post, para_std, nens)
        else:
            logmean = np.exp(para_mean_post + para_std**2 / 2.)
            logstd  = np.exp(2 * para_mean_post + para_std**2) * (np.exp(para_std**2) - 1)
            posterior = np.random.lognormal(logmean, logstd)

    # Sampled from truncated normal distribution
    elif para_sample_method.lower() == "truncated_normal":
        posterior = truncnorm.rvs(para_min, para_max, loc=para_mean_post, scale=para_std, size=nens)

    # Sampled from uniformation distribution
    elif para_sample_method.lower() == "uniform":
        posterior = np.random.uniform(para_min, para_max, nens)

    else:
        raise Exception("Unknown resampling method: {}".format(para_sample_method))

    # Discard any value exceeds the range of [para_min, para_max] except when the para_sample_method is rescale
    if para_sample_method.lower() != "rescale":
        if para_min != -99999:
            posterior[posterior < para_min] = para_min
        if para_max != 99999:
            posterior[posterior > para_max] = para_max

    return posterior


##############################################################
# Function for getting the posterior to be updated from the original posterior and prior
# (heterogeneous parameter)
##############################################################
def update_para_prior_from_original_posterior_heterogeneous(
    para_prior, original_para_posterior, para_name,
    para_loc_set, para_cell_ids, reconditioned_cellid_file, model_time,
    para_sample_method, para_min, para_max, para_mean, para_std, para_take_log):

    # Get the number of ensemble member
    nens, nloc = original_para_posterior.shape

    posterior = np.zeros([nens, nloc])

    # If it is sampled from a given distribution, then this is done by each individual location
    if para_sample_method.lower() in ["rescale","normal","lognormal","truncated_normal","uniform"]:
        print("Sampling method -- {} -- is used ...".format(para_sample_method))
        for i in range(nloc):
            # Get the mean and std from the original posterior
            para_mean_post = np.mean(original_para_posterior[:,i])
            para_std_post  = np.std(original_para_posterior[:,i])

            # If from rescaling the original posterior parameter
            if para_sample_method.lower() == "rescale":
                posterior[:,i] = (original_para_posterior[:,i] - para_mean_post) / para_std_post * para_std + para_mean_post

            # Sampled from normal distribution
            elif para_sample_method.lower() == "normal":
                posterior[:,i] = np.random.normal(para_mean_post, para_std, nens)

            # Sampled from lognormal distribution
            elif para_sample_method.lower() == "lognormal":
                if para_take_log:
                    posterior[:,i]  = np.random.normal(para_mean_post, para_std, nens)
                else:
                    logmean = np.exp(para_mean_post + para_std**2 / 2.)
                    logstd  = np.exp(2 * para_mean_post + para_std**2) * (np.exp(para_std**2) - 1)
                    posterior[:,i] = np.random.lognormal(logmean, logstd)

            # Sampled from truncated normal distribution
            elif para_sample_method.lower() == "truncated_normal":
                posterior[:,i] = truncnorm.rvs(para_min, para_max, loc=para_mean_post, scale=para_std, size=nens)

            # Sampled from uniformation distribution
            elif para_sample_method.lower() == "uniform":
                posterior[:,i] = np.random.uniform(para_min, para_max, nens)

    # If it is sampled from sequential gaussian simulation or kriging
    elif para_sample_method.lower() in ["sgs","kriging"]:
        print("Sampling method -- {} -- is used ...".format(para_sample_method))
        # Look for the new reconditioned points
        new_reconditioned_cell_ids = get_reconditioned_points(para_prior, original_para_posterior, para_cell_ids)
        new_reconditioned_cell_ids = new_reconditioned_cell_ids.tolist() # Convert it to python list
        # Add the new reconditioned points to the existing pool
        if os.path.exists(reconditioned_cellid_file):
            # reconditioned_cell_ids = np.loadtxt(reconditioned_cellid_file, dtype=int)
            # reconditioned_cell_ids = np.concatenate([reconditioned_cell_ids, new_reconditioned_cell_ids])
            with open(reconditioned_cellid_file, "r+") as f:
                all_reconditioned_cellids = json.load(f)
                if para_name not in all_reconditioned_cellids.keys():
                    reconditioned_cell_ids = new_reconditioned_cell_ids
                    # all_reconditioned_cellids[para_name] = {}
                    all_reconditioned_cellids[para_name] = {"initial":[], 
                                                         str(model_time):reconditioned_cell_ids,
                                                         "total":list(reconditioned_cell_ids)}
                else:
                    # Get the total ids
                    previous_total_reconditioned_cell_ids = all_reconditioned_cellids[para_name]["total"]
                    reconditioned_cell_ids = list(set(previous_total_reconditioned_cell_ids+new_reconditioned_cell_ids))
                    # Update the all_reconditioned_cellids
                    all_reconditioned_cellids[para_name][str(model_time)] = new_reconditioned_cell_ids
                    all_reconditioned_cellids[para_name]["total"] = reconditioned_cell_ids
                    # all_reconditioned_cellids[para_name] = {str(model_time):new_reconditioned_cell_ids,
                    #                                     "total":updated_total_reconditioned_cell_ids}
                f.seek(0)
                json.dump(all_reconditioned_cellids, f)
                f.truncate()
        else:
            reconditioned_cell_ids = new_reconditioned_cell_ids
            with open(reconditioned_cellid_file, "w") as f:
                all_reconditioned_cellids = {para_name: {"initial":[], "total":list(reconditioned_cell_ids)}}
                json.dump(all_reconditioned_cellids, f)
        # reconditioned_cell_ids = np.unique(reconditioned_cell_ids)
        reconditioned_cell_ids = np.array(reconditioned_cell_ids)  # Convert it back to numpy array
        # Check whether the reconditioned_cell_ids are within para_cell_ids
        if np.sum(np.isin(reconditioned_cell_ids, para_cell_ids)) != reconditioned_cell_ids.size:
            raise Exception("Only {} of {} reconditioned cells are within parameter cells".format(np.sum(np.isin(reconditioned_cell_ids, para_cell_ids)), reconditioned_cell_ids.size))
        print("The total number of conditioned point: {}".format(len(reconditioned_cell_ids)))
        print("Saving the conditioned points to file {} ...".format(reconditioned_cellid_file))
        # np.savetxt(reconditioned_cellid_file, reconditioned_cell_ids, fmt="%d")
        # Generate the ensemble from the conditioning simulation
        posterior = conditional_simulation(original_para_posterior, para_loc_set, para_cell_ids, reconditioned_cell_ids, para_min, para_max)

    else:
        raise Exception("Unknown resampling method: {}".format(para_sample_method))

    # Discard any value exceeds the range of [para_min, para_max] except when the para_sample_method is rescale
    if para_sample_method.lower() != "rescale":
        if para_min != -99999 and para_min != None:
            posterior[posterior < para_min] = para_min
        if para_max != 99999 and para_max != None:
            posterior[posterior > para_max] = para_max
    print(posterior.max(axis=1), posterior.min(axis=1))

    return posterior


##############################################################
# Function for choosing the reconditioning points
##############################################################
def get_reconditioned_points(prior, posterior, cell_ids):
    print("Selecting the reconditioned points ...")
    from scipy.stats import ks_2samp
    nens, nloc = prior.shape

    reconditioned_pool_first = []
    reconditioned_pool_second = []

    # The following should be configurable in the future
    ks_p_thres = 0.05
    ks_v_thres = 0.05

    # Compute KS statistics and p-value
    ks_p_set, ks_v_set = np.zeros(nloc), np.zeros(nloc)
    for i in range(nloc):
        prior_loc, posterior_loc = prior[:,i], posterior[:,i]
        ks_result = ks_2samp(prior_loc, posterior_loc)
        ks_v_set[i], ks_p_set[i] = ks_result.statistic, ks_result.pvalue

    # First round, pickup whose KS p-value is below the threshold
    reconditioned_pool_first = np.where(ks_p_set < ks_p_thres)[0]
    ks_v_set_filtered_first  = ks_v_set[ks_p_set < ks_p_thres]
    print("The number of conditioned points at the first round: {}".format(len(reconditioned_pool_first)))

    # Second round, pickup the points whose KS v-values are in the top ks_v_thres in the filtered result from the first round
    n_pool_first  = len(ks_v_set_filtered_first)
    n_pool_second = int(np.floor(n_pool_first * ks_v_thres))
    if n_pool_second == 0:
        print("The number of new points added to the reconditioning pool is zero!")
        reconditioned_pool_second = []
    else:
        reconditioned_pool_second = reconditioned_pool_first[np.argsort(ks_v_set_filtered_first)[-n_pool_second:]]
        print("The number of conditioned points at the second round: {}".format(len(reconditioned_pool_second)))

    # cell_ids = cell_ids.astype(int)
    reconditioned_cell_ids = cell_ids[reconditioned_pool_second]

    return reconditioned_cell_ids


#TODO
##############################################################
# Function for conducting Kriging for each realization
##############################################################
def krige_per_realization(i, conditioned, estimated, original_posterior):
    from rpy2 import robjects
    from rpy2.robjects import r, pandas2ri
    from rpy2.robjects import Formula, Environment
    from rpy2.robjects.vectors import IntVector, FloatVector, StrVector
    from rpy2.robjects.lib import grid, grdevices
    from rpy2.robjects.packages import importr, data
    import rpy2.robjects.packages as rpackages
    # from rpy2.rinterface import RRuntimeError

    # Activate and load the required R packages
    pandas2ri.activate()
    packageNames = ('ggplot2', 'gstat', 'sp','lazyeval', 'automap', 'stats')
    utils = importr('utils')
    utils.chooseCRANmirror(ind=1)
    packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
    # print(packnames_to_install)
    if len(packnames_to_install) > 0:
        utils.install_packages(StrVector(packnames_to_install))
    r_base  = importr('base')
    gstat   = importr('gstat')
    sp      = importr('sp')
    datasets = importr('datasets')
    stats = importr('stats')
    automap = importr('automap')
    print('a')

    conditioned["para"] = original_posterior
    # Convert to R objects
    conditioned_r = pandas2ri.py2rpy(conditioned)
    estimated_r   = pandas2ri.py2rpy(estimated)

    est_grids = 'est.grids'+str(i)
    obs_values = 'obs.values'+str(i)
    varig_m   = 'm'+str(i)
    est_values = 'sim_df'+str(i)

    robjects.globalenv[est_grids] = estimated_r
    robjects.globalenv[obs_values] = conditioned_r

    # Conditional simulation
    robjects.r("coordinates("+obs_values+") = ~x+y+z")
    robjects.r("coordinates("+est_grids+") = ~x+y+z")
    robjects.r("{} <- autofitVariogram(para~1,{},model=c('Sph','Exp','Gau','Ste'))".format(varig_m,obs_values))
    robjects.r("{} <- krige(para~1, {}, {}, model={}[['var_model']], nmax=20)".format(
        est_values, obs_values, est_grids, varig_m
    ))
    robjects.r("{} <- as.data.frame({})".format(est_values, est_values))
    sim_df = robjects.globalenv[est_values]

    # Remove variables
    robjects.r("rm("+est_grids+")")
    robjects.r("rm("+obs_values+")")
    # robjects.r("rm("+est_values+")")
    del conditioned, estimated
    print('d')

    # print(sim_df.shape)
    # print(sim_df.head())
    return i, sim_df["var1.pred"]


##############################################################
# Function for conducting conditional simulation
##############################################################
kriging_results = []
def call_back_from_poolapply(results):
    kriging_results.append(results)

def conditional_simulation(original_posterior, para_loc_set, para_cell_ids, reconditioned_cell_ids, para_min, para_max):

    nens, nloc = original_posterior.shape
    posterior  = np.zeros([nens, nloc])

    if len(reconditioned_cell_ids) == 0:
        print("There is no reconditioned points for conducting conditional simulation, so")
        print("use the original posterior generated from DART directly.")
        return original_posterior
        # raise Exception("There is no reconditioned points for conducting conditional simulation")
    else:
        from multiprocessing import Pool


    # Get the conditioned value and locations
    conditioned_loc_ind = np.searchsorted(para_cell_ids, reconditioned_cell_ids)
    conditioned_loc_set = para_loc_set[:, conditioned_loc_ind]
    conditioned         = pd.DataFrame(conditioned_loc_set.T, columns=["x","y","z"])

    # Get the locations to be estimated
    # estimated_loc_ind   = np.array(set(range(1,nloc+1))-set(conditioned_loc_ind))
    estimated_loc_ind   = np.setxor1d(np.arange(0,nloc), conditioned_loc_ind)
    estimated_loc_set   = para_loc_set[:, estimated_loc_ind]
    estimated           = pd.DataFrame(estimated_loc_set.T, columns=["x","y","z"])
    
    # conditioned["para"] = original_posterior[i, conditioned_loc_ind]
    # Conduct the kriging parallelly using multiprocessing Pool
    # with Pool() as p:
        # kriging_results = p.map(krig_per_realization, range(nens))
    global kriging_results
    kriging_results = []
    p = Pool()
    back = []
    for i in range(nens):
        b = p.apply_async(krige_per_realization, 
                          args=(i,conditioned.copy(),estimated.copy(),original_posterior[i, conditioned_loc_ind]), 
                          callback=call_back_from_poolapply)
        back.append(b)
        # p.apply(krige_per_realization, args=(i,conditioned.copy(),estimated.copy(),original_posterior[i]), callback=call_back_from_poolapply)
    p.close()
    # p.terminate()
    p.join()

    # Check whether each subprocess is successful
    for b in back:
        if not b.successful():
            print(b.get())
    
    # Re-sort kriging_results
    # if len(kriging_results) != nens:
    #     print("Something Wrong !!!")
    #     raise Exception("The number of kriging results {} does not equal to the realization number {}".format(len(kriging_results), nens))
    kriging_results_sorted = sorted(kriging_results, key=lambda tup: tup[0])
    kriging_results_sorted = [r[1] for r in kriging_results_sorted]
    print(os.cpu_count())
    # print(kriging_results)

    for i in range(nens):
        posterior[i,conditioned_loc_ind] = original_posterior[i,conditioned_loc_ind]
        posterior[i,estimated_loc_ind]   = kriging_results_sorted[i]
        # Conduct quality control
        # This is because kriging might result in some unreasonable values.
        # To avoid that, we convert those outliers to the mean of the conditional points
        improper_ind =(posterior[i,:]<=para_min)|(posterior[i,:]>=para_max)|(np.isnan(posterior[i,:]))
        posterior[i, improper_ind] = original_posterior[i,conditioned_loc_ind].mean()
        print(original_posterior[i,conditioned_loc_ind].mean(), posterior[i,:].mean())

    # print(posterior[:5,:10])
    # print(np.sum(posterior==0))

    # # Conduct the simulation for each realization using Kriging or one-time SGS
    # # TODO: to enable parallel computing to speed it up
    # for i in range(nens):
    #     conditioned["para"] = original_posterior[i, conditioned_loc_ind]
    #     # print(conditioned)
    #     # print(estimated)
    #     # Convert to R objects
    #     conditioned_r = pandas2ri.py2rpy(conditioned)
    #     estimated_r   = pandas2ri.py2rpy(estimated)

    #     robjects.globalenv['est.grids'] = estimated_r
    #     robjects.globalenv['obs.values'] = conditioned_r

    #     # Conditional simulation
    #     robjects.r("coordinates(obs.values) = ~x+y+z")
    #     robjects.r("coordinates(est.grids) = ~x+y+z")
    #     robjects.r("m <- autofitVariogram(para~1,obs.values,model=c('Sph','Exp','Gau','Ste'))")
    #     robjects.r("test.sim <- krige(para~1, obs.values, est.grids, model=m[['var_model']], nmax=100 )")
    #     robjects.r("test.sim_df <- as.data.frame(test.sim)")
    #     sim_df = robjects.globalenv['test.sim_df']
    #     # print(sim_df.shape)
    #     # print(sim_df.head())

    #     posterior[i,conditioned_loc_ind] = original_posterior[i,conditioned_loc_ind]
    #     posterior[i,estimated_loc_ind] = sim_df["var1.pred"]
        # raise Exception("stop")

    return posterior


    
###############################
# Conversion from PFLOTRAN output to the required state variables
# LIQUID_PRESSURE --> WATER_LEVEL
# Total Tracer --> TRACER_GROUNDWATER
###############################
