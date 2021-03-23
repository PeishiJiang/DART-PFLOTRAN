"""Plot the data assimilation results."""

# Author: Peishi Jiang

import os
import re
import sys
import f90nml
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import spatial
from datetime import datetime
from math import ceil, log10
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
from matplotlib import gridspec

# TODO: A more sophisticated OOP architecture will be required in the future to avoid the redundant coding.

calendar     = 'gregorian'
tunits_begin = 'days since '

def plot_compare_multiple_daresults(var_name,         # the variable name to be compared
                                    settings_list,    # the list of settings/scenario names
                                    config_file_list, # the list of config.nml files to be compared
                                    ax,               # the axis for plotting
                                    true_used):       # the true values
                                    # true_file_name):  # the file containing the true values
    """This function is used for compare multiple data assimilation results. Therefore, it has to have the same variables and the same
       assimilation time length & time window."""
    n_settings = len(config_file_list)

    # Get the true values
    # true_set        = pd.read_csv(true_file_name)
    # true_set_raw_time  = true_set.iloc[:, 0].values
    # # true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
    # true               = true_set.iloc[:, 1].values
    # true               = true[2:, :]
    # true_set_used_ind  = (true_set_time >= model_start_time) & (true_set_time <= model_end_time)
    # true_set_used      = true[true_set_used_ind]

    # Initialize the dictionary for ensemble errors
    error_pd = pd.DataFrame(columns=['type', 'error', 'setting'])
    # errors_prior     = dict.fromkeys(settings_list)
    # errors_posterior = dict.fromkeys(settings_list)

    # Get the simulated values
    for i in range(n_settings):
        # Read the namelist file
        configs = f90nml.read(config_file_list[i])

        # Get the file names for prior and posterior
        dart_prior_all_ens_file     = configs["file_cfg"]["dart_prior_nc_all_ens_file"]
        dart_posterior_all_ens_file = configs["file_cfg"]["dart_posterior_nc_all_ens_file"]

        # Read the prior and posterior data
        root_prior     = Dataset(dart_prior_all_ens_file, 'r')
        root_posterior = Dataset(dart_posterior_all_ens_file, 'r')

        # Get the spatial averaged values
        prior          , posterior           = root_prior.variables[var_name][:], root_posterior.variables[var_name][:]
        prior_space_ave, posterior_space_ave = np.mean(prior, axis=(2, 3, 4)),    np.mean(posterior, axis=(2, 3, 4))

        nens, ntime = prior_space_ave.shape
        # true_used = true[:nens, :]

        errors_prior     = prior_space_ave - np.tile(true_used, [nens, 1])
        errors_posterior = posterior_space_ave - np.tile(true_used, [nens, 1])
        # errors_prior, errors_posterior = errors_prior.flatten(), errors_posterior.flatten()

        errors_prior_pd = pd.DataFrame({'type': ["prior"] * nens * ntime,
                                        'error': errors_prior.flatten(),
                                        'setting': [settings_list[i]] * nens * ntime })
        errors_posterior_pd = pd.DataFrame({'type': ["posterior"] * nens * ntime,
                                            'error': errors_posterior.flatten(),
                                            'setting': [settings_list[i]] * nens * ntime })

        error_pd = error_pd.append(errors_prior_pd)
        error_pd = error_pd.append(errors_posterior_pd)

        # errors_prior[settings_list[i]]     = prior_space_ave - np.tile(true, [nens, 1])
        # errors_posterior[settings_list[i]] = posterior_space_ave - np.tile(true, [nens, 1])

    # Plot
    sns.violinplot(x='setting', y='error', hue='type', data=error_pd, split=True, ax=ax)
    ax.axhline(y=0, color='k', linestyle='--')


class DaResults(object):
    def __init__(self, config_nml):
        """Initialization: Read in the required configurations from config_nml."""
        # Parse the configuration in Fortran namelist
        if not os.path.isfile(config_nml):
            raise Exception(
                "The configuration file of this application does not exist under the path: %s!"
                % config_nml)
        self.config_nml = config_nml
        self.configs    = f90nml.read(config_nml)

        # Variables
        self.obs_var_set  = self.configs["obspara_set_cfg"]["obs_set"]
        self.para_var_set = self.configs["obspara_set_cfg"]["para_set"]
        if isinstance(self.obs_var_set, str):
            self.obs_var_set = [self.obs_var_set]
        if isinstance(self.para_var_set, str):
            self.para_var_set = [self.para_var_set]
        self.pflotran_var_set = self.obs_var_set + self.para_var_set
        self.nvar_state = len(self.obs_var_set)
        self.nvar_para  = len(self.para_var_set)
        self.nvar = len(self.pflotran_var_set)

        # Model time, assimilation start time, number of ensemble
        self.assim_start_time = self.configs["time_cfg"]["assim_start"]
        self.model_start_time = self.configs["time_cfg"]["model_start"]
        self.model_time       = float(self.configs["time_cfg"]["current_model_time"])  # days
        self.model_time_list  = self.configs["time_cfg"]["model_time_list"]
        self.model_time_list  = [self.model_time_list] if not isinstance(
            self.model_time_list, list) else self.model_time_list
        # self.model_time_list  = [t+model_time_offset for t in self.model_time_list]
        self.exceeds_obs_time = self.configs["time_cfg"]["exceeds_obs_time"]
        # self.assim_window     = float(self.configs["da_cfg"]["assim_window_size"])
        self.nens             = self.configs["da_cfg"]["nens"]
        self.model_ntime      = len(self.model_time_list)
        self.model_ntimestep  = int(self.configs["da_cfg"]["ntimestep"])

        # Convert the model time to time units
        self.tunits                = tunits_begin + self.model_start_time
        self.model_time_dates_list = num2date(self.model_time_list, self.tunits, calendar)

        # Observation file in NetCDF
        self.obs_nc = self.configs["file_cfg"]["obs_nc_file"]
        
        # Whether immediate mda results are saved
        self.has_immediate_mda_results = self.configs["file_cfg"]["save_immediate_mda_result"]
        self.total_mda_iterations = self.configs["da_cfg"]["enks_mda_total_iterations"]

        if not self.exceeds_obs_time:
            warnings.warn("The data assimilation is not completed!")

        # # Get the model spatial domain
        # self.dart_prior_template = self.configs["file_cfg"]["dart_prior_template_file"]
        # root_template = Dataset(self.dart_prior_template, 'r')
        # x_loc = root_template.variables['x_location'][:]
        # y_loc = root_template.variables['y_location'][:]
        # z_loc = root_template.variables['z_location'][:]
        # self.model_loc_set = np.array([x_loc, y_loc, z_loc]).T
        # print(self.model_loc_set.shape)
        # self.model_nloc = self.model_loc_set.shape[0]
        # root_template.close()


    def setup(self):
        """Read in the observation, prior, and posterior data"""
        # Get the parameters
        obs_nc           = self.obs_nc
        obs_var_set      = self.obs_var_set
        para_var_set     = self.para_var_set
        pflotran_var_set = self.pflotran_var_set

        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        nvar, nens, niter     = self.nvar, self.nens, self.total_mda_iterations
        has_immediate_mda_results = self.has_immediate_mda_results

        # DART prior and posterior files
        self.dart_posterior_all_ens_file = self.configs["file_cfg"]["dart_posterior_nc_all_ens_file"]
        self.dart_prior_all_ens_file     = self.configs["file_cfg"]["dart_prior_nc_all_ens_file"]
        dart_prior_all_ens_file, dart_posterior_all_ens_file = self.dart_prior_all_ens_file, self.dart_posterior_all_ens_file

        prior     = dict.fromkeys(["para","state"])
        posterior = dict.fromkeys(["para","state"])

        # Get the spatial-temporal domains for both model parameters and states
        with Dataset(dart_prior_all_ens_file, 'r') as root_prior:
            # Get the list of assimilation start and end time
            if has_immediate_mda_results:
                assim_start_set = root_prior.variables['assim_start_time'][:][0,0,:]
                assim_end_set   = root_prior.variables['assim_end_time'][:][0,0,:]
            else:
                assim_start_set = root_prior.variables['assim_start_time'][:][0,:]
                assim_end_set   = root_prior.variables['assim_end_time'][:][0,:]

            # Get the numbers of time steps and model grids
            ntime_state = root_prior.dimensions['state_time'].size
            ntime_para  = root_prior.dimensions['para_time'].size
            nloc_state = root_prior.dimensions['state_location'].size
            nloc_para  = root_prior.dimensions['para_location'].size

            # The temporal domain
            state_time_set = root_prior.variables['state_time'][:]
            para_time_set = root_prior.variables['para_time'][:]
            # ntime_state, ntime_para = len(state_time_set), len(para_time_set)

            # The spatial domain
            if has_immediate_mda_results:
                xloc_set_state = root_prior.variables['state_x_location'][:][0,0,:]
                yloc_set_state = root_prior.variables['state_y_location'][:][0,0,:]
                zloc_set_state = root_prior.variables['state_z_location'][:][0,0,:]
                xloc_set_para = root_prior.variables['para_x_location'][:][0,0,:]
                yloc_set_para = root_prior.variables['para_y_location'][:][0,0,:]
                zloc_set_para = root_prior.variables['para_z_location'][:][0,0,:]
            else:
                xloc_set_state = root_prior.variables['state_x_location'][:][0,:]
                yloc_set_state = root_prior.variables['state_y_location'][:][0,:]
                zloc_set_state = root_prior.variables['state_z_location'][:][0,:]
                xloc_set_para = root_prior.variables['para_x_location'][:][0,:]
                yloc_set_para = root_prior.variables['para_y_location'][:][0,:]
                zloc_set_para = root_prior.variables['para_z_location'][:][0,:]
            nxloc_state, nyloc_state, nzloc_state = len(xloc_set_state), len(yloc_set_state), len(zloc_set_state)
            nxloc_para, nyloc_para, nzloc_para = len(xloc_set_para), len(yloc_set_para), len(zloc_set_para)
            # nloc_state = nxloc_state * nyloc_state * nzloc_state
            # nloc_para = nxloc_para * nyloc_para * nzloc_para

        # Sort the time
        state_time_arg_sort = state_time_set.argsort()
        para_time_arg_sort  = para_time_set.argsort()
        state_time_set.sort(); para_time_set.sort()

        # Read in the prior data
        if has_immediate_mda_results:
            prior['state'] = np.zeros([nvar_state, niter, nens, ntime_state, nloc_state])
            prior['para'] = np.zeros([nvar_para, niter, nens, ntime_para, nloc_para])
        else:
            prior['state'] = np.zeros([nvar_state, nens, ntime_state, nloc_state])
            prior['para'] = np.zeros([nvar_para, nens, ntime_para, nloc_para])
        
        with Dataset(dart_prior_all_ens_file, 'r') as root_prior:
            # State variables
            for i in range(nvar_state):
                varn      = obs_var_set[i]
                prior_var = root_prior.variables[varn][:]
                if has_immediate_mda_results:
                    prior["state"][i, :, :, :, :] = prior_var[:, :,state_time_arg_sort,:]
                else:
                    prior["state"][i, :, :, :] = prior_var[:,state_time_arg_sort,:]
            # Parameter variables
            for i in range(nvar_para):
                varn      = para_var_set[i]
                prior_var = root_prior.variables[varn][:]
                if has_immediate_mda_results:
                    prior["para"][i, :, :, :, :] = prior_var[:, :,para_time_arg_sort,:]
                else:
                    prior["para"][i, :, :, :] = prior_var[:,para_time_arg_sort,:]

        # Read in the posterior data
        if has_immediate_mda_results:
            posterior['state'] = np.zeros([nvar_state, niter, nens, ntime_state, nloc_state])
            posterior['para'] = np.zeros([nvar_para, niter, nens, ntime_para, nloc_para])
        else:
            posterior['state'] = np.zeros([nvar_state, nens, ntime_state, nloc_state])
            posterior['para'] = np.zeros([nvar_para, nens, ntime_para, nloc_para])
        with Dataset(dart_posterior_all_ens_file, 'r') as root_posterior:
            # State variables
            for i in range(nvar_state):
                varn          = obs_var_set[i]
                posterior_var = root_posterior.variables[varn][:]
                if has_immediate_mda_results:
                    posterior["state"][i, :, :, :, :] = posterior_var[:, :,state_time_arg_sort,:]
                else:
                    posterior["state"][i, :, :, :] = posterior_var[:,state_time_arg_sort,:]
            # Parameter variables
            for i in range(nvar_para):
                varn          = para_var_set[i]
                posterior_var = root_posterior.variables[varn][:]
                if has_immediate_mda_results:
                    posterior["para"][i, :, :, :, :] = posterior_var[:, :,para_time_arg_sort,:]
                else:
                    posterior["para"][i, :, :, :] = posterior_var[:,para_time_arg_sort,:]

        # Read in the observation data
        obs_value_set      = dict.fromkeys(obs_var_set)
        obs_value_set_used = dict.fromkeys(obs_var_set)
        with Dataset(obs_nc, 'r') as root_obs:
            # Time in observation
            obs_time_set = root_obs.variables['time'][:]
            obs_used_ind = (obs_time_set >= state_time_set.min()) & (obs_time_set <= state_time_set.max())
            obs_time_set_used = obs_time_set[obs_used_ind]
            # Locations in observation
            xloc = root_obs.variables['x_location']
            yloc = root_obs.variables['y_location']
            zloc = root_obs.variables['z_location']
            obs_loc = np.array([xloc, yloc, zloc]).T
            # Values in observation
            for i in range(len(obs_var_set)):
                obs_var = obs_var_set[i]
                obs_value_set[obs_var] = root_obs.variables[obs_var][:]
                obs_value_set_used[obs_var] = obs_value_set[obs_var][:, obs_used_ind]
                
        self.prior, self.posterior = prior, posterior
        self.assim_start_set, self.assim_end_set = assim_start_set, assim_end_set
        self.state_time_set, self.para_time_set = state_time_set, para_time_set
        self.ntime_state, self.ntime_para = ntime_state, ntime_para
        self.xloc_set_state, self.yloc_set_state, self.zloc_set_state = xloc_set_state, yloc_set_state, zloc_set_state
        self.nxloc_state, self.nyloc_state, self.nzloc_state = nxloc_state, nyloc_state, nzloc_state
        self.xloc_set_para, self.yloc_set_para, self.zloc_set_para = xloc_set_para, yloc_set_para, zloc_set_para
        self.nxloc_para, self.nyloc_para, self.nzloc_para = nxloc_para, nyloc_para, nzloc_para
        self.nloc_state, self.nloc_para = nloc_state, nloc_para
        self.obs_value_set_used = obs_value_set_used
        self.obs_time_set_used = obs_time_set_used
        self.obs_loc_set = obs_loc
        


    def plot_spatial_average(self, axes, ylim=None, plot_averaged_obs=False):
        """Plot the spatial averaged DA results"""
        nvar, nens = self.nvar, self.nens
        prior, posterior      = self.prior, self.posterior

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        has_immediate_mda_results = self.has_immediate_mda_results

        # The shape of axes has to be (nvar,2)
        if axes.shape != (nvar,2):
            raise Exception("The shape {} of the axes is incorrect".format(axes.shape))

        # Plot states
        for i in range(nvar_state):
            varn = obs_var_set[i]
            if has_immediate_mda_results:
                prior_var, posterior_var = prior["state"][i][0], posterior["state"][i][-1]
            else:
                prior_var, posterior_var = prior["state"][i], posterior["state"][i]

            # Plot the prior
            ax1 = axes[i, 0]
            for j in range(nens):
                prior_var_ens_mean = np.mean(prior_var[j, :, :], axis=(1))
                line1, = ax1.plot(state_time_set, prior_var_ens_mean, color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            prior_var_mean = np.mean(prior_var[:, :, :], axis=(0, 2))
            line2, = ax1.plot(state_time_set, prior_var_mean, color='red', linewidth=1, label='mean')
            if varn in obs_var_set and plot_averaged_obs:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                line3, = ax1.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=1, label='obs')

            # Plot the posterior
            ax2 = axes[i, 1]
            for j in range(nens):
                posterior_var_ens_mean = np.mean(posterior_var[j, :, :], axis=(1))
                ax2.plot(state_time_set, posterior_var_ens_mean, color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            posterior_var_mean = np.mean(posterior_var[:, :, :], axis=(0, 2))
            ax2.plot(state_time_set, posterior_var_mean, color='red', linewidth=1, label='mean')
            if varn in obs_var_set and plot_averaged_obs:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                ax2.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=1, label='obs')
            ax2.set_ylim(ylim)

            # Plot the variable name
            ax1.set_ylabel(varn)

        # Plot parameters
        for i in range(nvar_para):
            varn = para_var_set[i]
            if has_immediate_mda_results:
                prior_var, posterior_var = prior["para"][i][0], posterior["para"][i][-1]
            else:
                prior_var, posterior_var = prior["para"][i], posterior["para"][i]

            # Plot the prior
            ax1 = axes[nvar_state+i, 0]
            for j in range(nens):
                prior_var_ens_mean = np.mean(prior_var[j, :, :], axis=(1))
                line1, = ax1.plot(para_time_set, prior_var_ens_mean, color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            prior_var_mean = np.mean(prior_var[:, :, :], axis=(0, 2))
            line2, = ax1.plot(para_time_set, prior_var_mean, color='red', linewidth=1, label='mean')

            # Plot the posterior
            ax2 = axes[nvar_state+i, 1]
            for j in range(nens):
                posterior_var_ens_mean = np.mean(posterior_var[j, :, :], axis=(1))
                ax2.plot(para_time_set, posterior_var_ens_mean, color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            posterior_var_mean = np.mean(posterior_var[:, :, :], axis=(0, 2))
            ax2.plot(para_time_set, posterior_var_mean, color='red', linewidth=1, label='mean')
            ax2.set_ylim(ylim)

            # Plot the variable name
            ax1.set_ylabel(varn)

        for i in range(2):
            ax = axes[-1, i]
            ax.set_xlabel('Time (day)')
        axes[0, 0].set_title("Prior")
        axes[0, 1].set_title("Posterior")

        # Plot the legends
        if plot_averaged_obs:
            plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
                        frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))
        else:
            plt.legend((line1, line2), ('ensemble', 'mean'),
                        frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))


    def plot_obs_at_point(self, obs_name, axes, obs_loc_ind=0, plot_time_offset=0,
                        #   figsize=None, constrained_layout=True,
                          vmin=None, vmax=None, ylim=None, option='both'):
        """Plot the temporal evolution of DA results for the observation data along one dimension"""
        nvar, nens = self.nvar, self.nens
        prior_state, posterior_state = self.prior["state"], self.posterior["state"]

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        has_immediate_mda_results = self.has_immediate_mda_results

        xloc_set_state, yloc_set_state, zloc_set_state = self.xloc_set_state, self.yloc_set_state, self.zloc_set_state

        model_loc_set = np.array([xloc_set_state, yloc_set_state, zloc_set_state]).T

        if obs_name not in obs_var_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the observation point and values
        obs_value = obs_value_set_used[obs_name][obs_loc_ind, :]
        obs_loc   = obs_loc_set[obs_loc_ind]

        # Find the closest location in the model grids for this observation
        dist, model_loc_ind = spatial.KDTree(model_loc_set).query(obs_loc)
        model_loc = model_loc_set[model_loc_ind]
        print("The observation location of interest is {}".format(obs_loc))
        print("The corresponding model grid is {}".format(model_loc))

        # Get the index of observation
        obs_var_ind = obs_var_set.index(obs_name)

        ##############################
        # Plot the temporal evolution of the ensemble, the mean, and the observation
        # at the given observed location
        # Get the locations of observation and the corresponding index of DA results
        if option == 'both':
            # Plot the prior
            ax1 = axes[0]
            for j in range(nens):
                if has_immediate_mda_results:
                    prior_var_ens = prior_state[obs_var_ind, 0, j, :, model_loc_ind]
                else:
                    prior_var_ens = prior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax1.plot(state_time_set[plot_time_offset:], prior_var_ens[plot_time_offset:], 
                                  color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                prior_var_mean = np.mean(prior_state[obs_var_ind, 0, :, :, model_loc_ind], axis=(0))
            else:
                prior_var_mean = np.mean(prior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax1.plot(state_time_set[plot_time_offset:], prior_var_mean[plot_time_offset:], 
                              color='red', linewidth=1, label='mean')
            line3, = ax1.plot(obs_time_set_used[plot_time_offset:], obs_value[plot_time_offset:], 
                              color='black', linewidth=1, label='obs')
            ax1.set_ylim(ylim)

            # Plot the posterior
            ax2 = axes[1]
            for j in range(nens):
                if has_immediate_mda_results:
                    posterior_var_ens = posterior_state[obs_var_ind, -1, j, :, model_loc_ind]
                else:
                    posterior_var_ens = posterior_state[obs_var_ind, j, :, model_loc_ind]
                ax2.plot(state_time_set[plot_time_offset:], posterior_var_ens[plot_time_offset:], 
                         color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, -1, :, :, model_loc_ind], axis=(0))
            else:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            ax2.plot(state_time_set[plot_time_offset:], posterior_var_mean[plot_time_offset:], 
                     color='red', linewidth=1, label='mean')
            ax2.plot(obs_time_set_used[plot_time_offset:], obs_value[plot_time_offset:], 
                     color='black', linewidth=1, label='obs')
            ax2.set_ylim(ylim)
        
        elif option == 'prior':
            ax1 = axes[0]
            for j in range(nens):
                if has_immediate_mda_results:
                    prior_var_ens = prior_state[obs_var_ind, 0, j, :, model_loc_ind]
                else:
                    prior_var_ens = prior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax1.plot(state_time_set[plot_time_offset:], prior_var_ens[plot_time_offset:], 
                                  color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                prior_var_mean = np.mean(prior_state[obs_var_ind, 0, :, :, model_loc_ind], axis=(0))
            else:
                prior_var_mean = np.mean(prior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax1.plot(state_time_set[plot_time_offset:], prior_var_mean[plot_time_offset:], 
                              color='red', linewidth=1, label='mean')
            line3, = ax1.plot(obs_time_set_used[plot_time_offset:], obs_value[plot_time_offset:], 
                              color='black', linewidth=1, label='obs')
            ax1.set_ylim(ylim)
            
        elif option == 'posterior':
            ax2 = axes[0]
            for j in range(nens):
                if has_immediate_mda_results:
                    posterior_var_ens = posterior_state[obs_var_ind, -1, j, :, model_loc_ind]
                else:
                    posterior_var_ens = posterior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax2.plot(state_time_set[plot_time_offset:], posterior_var_ens[plot_time_offset:], 
                         color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, -1, :, :, model_loc_ind], axis=(0))
            else:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax2.plot(state_time_set[plot_time_offset:], posterior_var_mean[plot_time_offset:], 
                     color='red', linewidth=1, label='mean')
            line3, = ax2.plot(obs_time_set_used[plot_time_offset:], obs_value[plot_time_offset:], 
                     color='black', linewidth=1, label='obs')
            ax2.set_ylim(ylim)

        return line1, line2, line3

    
    def plot_1to1_obs_at_point(self, obs_name, axes, obs_loc_ind=0, plot_time_offset=0,
                        #   figsize=None, constrained_layout=True,
                          vmin=None, vmax=None, ylim=None, option='prior'):
        """Plot the temporal evolution of DA results for the observation data along one dimension"""
        nvar, nens = self.nvar, self.nens
        prior_state, posterior_state = self.prior["state"], self.posterior["state"]

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        has_immediate_mda_results = self.has_immediate_mda_results

        xloc_set_state, yloc_set_state, zloc_set_state = self.xloc_set_state, self.yloc_set_state, self.zloc_set_state

        model_loc_set = np.array([xloc_set_state, yloc_set_state, zloc_set_state]).T

        if obs_name not in obs_var_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the observation point and values
        obs_value = obs_value_set_used[obs_name][obs_loc_ind, :]
        obs_loc   = obs_loc_set[obs_loc_ind]

        # Find the closest location in the model grids for this observation
        dist, model_loc_ind = spatial.KDTree(model_loc_set).query(obs_loc)
        model_loc = model_loc_set[model_loc_ind]
        print("The observation location of interest is {}".format(obs_loc))
        print("The corresponding model grid is {}".format(model_loc))

        # Get the index of observation
        obs_var_ind = obs_var_set.index(obs_name)

        ##############################
        # Plot the temporal evolution of the ensemble, the mean, and the observation
        # at the given observed location
        # Get the locations of observation and the corresponding index of DA results
        if option == 'both':
            # Plot the prior
            ax1 = axes[0]
            ax1.plot(ylim, ylim, 'k--')
            for j in range(nens):
                if has_immediate_mda_results:
                    prior_var_ens = prior_state[obs_var_ind, 0, j, :, model_loc_ind]
                else:
                    prior_var_ens = prior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax1.plot(obs_value[plot_time_offset:], prior_var_ens[plot_time_offset:], '.',
                                  color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                prior_var_mean = np.mean(prior_state[obs_var_ind, 0, :, :, model_loc_ind], axis=(0))
            else:
                prior_var_mean = np.mean(prior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax1.plot(obs_value[plot_time_offset:], prior_var_mean[plot_time_offset:], '.',
                              color='red', linewidth=1, label='mean')
            ax1.set_xlim(ylim); ax1.set_ylim(ylim)

            # Plot the posterior
            ax2 = axes[1]
            ax2.plot(ylim, ylim, 'k--')
            for j in range(nens):
                if has_immediate_mda_results:
                    posterior_var_ens = posterior_state[obs_var_ind, -1, j, :, model_loc_ind]
                else:
                    posterior_var_ens = posterior_state[obs_var_ind, j, :, model_loc_ind]
#                 ax2.plot(state_time_set[plot_time_offset:], posterior_var_ens[plot_time_offset:], 
                ax2.plot(obs_value[plot_time_offset:], posterior_var_ens[plot_time_offset:], '.',
                         color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, -1, :, :, model_loc_ind], axis=(0))
            else:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
#             ax2.plot(state_time_set[plot_time_offset:], posterior_var_mean[plot_time_offset:],
            ax2.plot(obs_value[plot_time_offset:], posterior_var_mean[plot_time_offset:], '.',
                     color='red', linewidth=1, label='mean')
            ax2.set_xlim(ylim); ax2.set_ylim(ylim)
        
        elif option == 'prior':
            ax1 = axes[0]
            ax1.plot(ylim, ylim, 'k--')
            for j in range(nens):
                if has_immediate_mda_results:
                    prior_var_ens = prior_state[obs_var_ind, 0, j, :, model_loc_ind]
                else:
                    prior_var_ens = prior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax1.plot(obs_value[plot_time_offset:], prior_var_ens[plot_time_offset:], '.',
                                  color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                prior_var_mean = np.mean(prior_state[obs_var_ind, 0, :, :, model_loc_ind], axis=(0))
            else:
                prior_var_mean = np.mean(prior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax1.plot(obs_value[plot_time_offset:], prior_var_mean[plot_time_offset:], '.',
                              color='red', linewidth=1, label='mean')
            ax1.set_xlim(ylim); ax1.set_ylim(ylim)
            
        elif option == 'posterior':
            ax2 = axes[0]
            ax2.plot(ylim, ylim, 'k--')
            for j in range(nens):
                if has_immediate_mda_results:
                    posterior_var_ens = posterior_state[obs_var_ind, -1, j, :, model_loc_ind]
                else:
                    posterior_var_ens = posterior_state[obs_var_ind, j, :, model_loc_ind]
                line1, = ax2.plot(obs_value[plot_time_offset:], posterior_var_ens[plot_time_offset:], '.',
                         color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            if has_immediate_mda_results:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, -1, :, :, model_loc_ind], axis=(0))
            else:
                posterior_var_mean = np.mean(posterior_state[obs_var_ind, :, :, model_loc_ind], axis=(0))
            line2, = ax2.plot(obs_value[plot_time_offset:], posterior_var_mean[plot_time_offset:], '.',
                     color='red', linewidth=1, label='mean')
            ax2.set_xlim(ylim); ax2.set_ylim(ylim)

        return line1, line2

    
    def compute_obs_at_point_diff(self, obs_name, obs_loc_ind=0, method='mae'):
        """Compute the bias of the updated observation variable against the true values"""
        nvar, nens = self.nvar, self.nens
        prior_state, posterior_state = self.prior["state"], self.posterior["state"]

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        has_immediate_mda_results = self.has_immediate_mda_results

        xloc_set_state, yloc_set_state, zloc_set_state = self.xloc_set_state, self.yloc_set_state, self.zloc_set_state

        model_loc_set = np.array([xloc_set_state, yloc_set_state, zloc_set_state]).T

        if obs_name not in obs_var_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the observation point and values
        obs_value = obs_value_set_used[obs_name][obs_loc_ind, :]
        obs_loc   = obs_loc_set[obs_loc_ind]

        # Find the closest location in the model grids for this observation
        dist, model_loc_ind = spatial.KDTree(model_loc_set).query(obs_loc)
        model_loc = model_loc_set[model_loc_ind]
        print("The observation location of interest is {}".format(obs_loc))
        print("The corresponding model grid is {}".format(model_loc))

        # Get the index of observation
        obs_var_ind = obs_var_set.index(obs_name)

        # Get the difference between the estimated and the true for each realization
        if has_immediate_mda_results:
            prior_ens     = prior_state[obs_var_ind, 0, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, -1, :, :, model_loc_ind]
        else:
            prior_ens     = prior_state[obs_var_ind, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, :, :, model_loc_ind]
        diff_prior = prior_ens[:, :] - obs_value[:]
        diff_post = posterior_ens[:, :] - obs_value[:]
        if method.lower() == 'mae':
            diff_prior = np.mean(np.abs(diff_prior), axis=0)
            diff_post = np.mean(np.abs(diff_post), axis=0)
        elif method.lower() == 'bias':
            diff_prior = np.mean(diff_prior, axis=0)
            diff_post = np.mean(diff_post, axis=0)
        
        return diff_prior, diff_post


    def compute_obs_at_point_sd(self, obs_name, obs_loc_ind=0):
        """Compute the std of the updated observation variable against the true values"""
        nvar, nens = self.nvar, self.nens
        prior_state, posterior_state = self.prior["state"], self.posterior["state"]

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        has_immediate_mda_results = self.has_immediate_mda_results

        xloc_set_state, yloc_set_state, zloc_set_state = self.xloc_set_state, self.yloc_set_state, self.zloc_set_state

        model_loc_set = np.array([xloc_set_state, yloc_set_state, zloc_set_state]).T

        if obs_name not in obs_var_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the observation point and values
        obs_value = obs_value_set_used[obs_name][obs_loc_ind, :]
        obs_loc   = obs_loc_set[obs_loc_ind]

        # Find the closest location in the model grids for this observation
        dist, model_loc_ind = spatial.KDTree(model_loc_set).query(obs_loc)
        model_loc = model_loc_set[model_loc_ind]
        print("The observation location of interest is {}".format(obs_loc))
        print("The corresponding model grid is {}".format(model_loc))

        # Get the index of observation
        obs_var_ind = obs_var_set.index(obs_name)

        # Get the difference between the estimated and the true for each realization
        if has_immediate_mda_results:
            prior_ens     = prior_state[obs_var_ind, 0, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, -1, :, :, model_loc_ind]
        else:
            prior_ens     = prior_state[obs_var_ind, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, :, :, model_loc_ind]
        std_prior = np.std(prior_ens, axis=0)
        std_post  = np.std(posterior_ens, axis=0)
        
        return std_prior, std_post
        
    
    def compare_obs_at_point_bias(self, obs_name, axes, obs_loc_ind=0, constrained_layout=True, 
                                  plot_time_offset=0, ylim=None, xlim=None):
        """Plot the temporal evolution of the bias and MAE of the updated observation variables against the true value."""
        nvar, nens = self.nvar, self.nens
        prior_state, posterior_state = self.prior["state"], self.posterior["state"]

        ntime_state, ntime_para = self.ntime_state, self.ntime_para
        nvar_state, nvar_para = self.nvar_state, self.nvar_para
        obs_var_set           = self.obs_var_set
        para_var_set          = self.para_var_set
        state_time_set        = self.state_time_set
        para_time_set         = self.para_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        has_immediate_mda_results = self.has_immediate_mda_results

        xloc_set_state, yloc_set_state, zloc_set_state = self.xloc_set_state, self.yloc_set_state, self.zloc_set_state

        model_loc_set = np.array([xloc_set_state, yloc_set_state, zloc_set_state]).T

        if obs_name not in obs_var_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the observation point and values
        obs_value = obs_value_set_used[obs_name][obs_loc_ind, :]
        obs_loc   = obs_loc_set[obs_loc_ind]

        # Find the closest location in the model grids for this observation
        dist, model_loc_ind = spatial.KDTree(model_loc_set).query(obs_loc)
        model_loc = model_loc_set[model_loc_ind]
        print("The observation location of interest is {}".format(obs_loc))
        print("The corresponding model grid is {}".format(model_loc))

        # Get the index of observation
        obs_var_ind = obs_var_set.index(obs_name)

        # Get the difference between the estimated and the true for each realization
        if has_immediate_mda_results:
            prior_ens     = prior_state[obs_var_ind, 0, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, -1, :, :, model_loc_ind]
        else:
            prior_ens     = prior_state[obs_var_ind, :, :, model_loc_ind]
            posterior_ens = posterior_state[obs_var_ind, :, :, model_loc_ind]
        diff_prior = prior_ens[:, :] - obs_value[:]
        bias_prior = np.mean(diff_prior, axis=0)
        mae_prior  = np.mean(np.abs(diff_prior), axis=0)
        diff_post = posterior_ens[:, :] - obs_value[:]
        bias_post = np.mean(diff_post, axis=0)
        mae_post  = np.mean(np.abs(diff_post), axis=0)

        # Plot the bias and MAE
        # Prior
        ax = axes[0]
        line1, = ax.plot(state_time_set[plot_time_offset:], bias_prior[plot_time_offset:], color='blue', linewidth=1, label='BIAS')
        line2, = ax.plot(state_time_set[plot_time_offset:], mae_prior[plot_time_offset:], color='red', linestyle='--', linewidth=1, label='MAE')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{|BIAS|}: %.3f;$ $\mu_{MAE}: %.3f$" % (np.mean(np.abs(bias_prior[plot_time_offset:])), np.mean(mae_prior[plot_time_offset:])))
        # Plot the labels and titles
        # ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # Posterior
        ax = axes[1]
        line1, = ax.plot(state_time_set[plot_time_offset:], bias_post[plot_time_offset:], color='blue', linewidth=1, label='BIAS')
        line2, = ax.plot(state_time_set[plot_time_offset:], mae_post[plot_time_offset:], color='red', linestyle='--', linewidth=1, label='MAE')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{|BIAS|}: %.3f;$ $\mu_{MAE}: %.3f$" % (np.mean(np.abs(bias_post[plot_time_offset:])), np.mean(mae_post[plot_time_offset:]))) 
        # Plot the labels and titles
        # ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        return line1, line2 


    def compare_univar_spatial_average(self, var_name, true_file_name, axes, constrained_layout=True, 
                                       model_time_offset=0, plot_time_offset=0, ylim=None, xlim=None, option='both'):
        """Plot the temporal evolution of a spatial averaged analyzed variable against the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values
        true_set_used_ind  = (true_set_time >= para_time_set[0]) & (true_set_time <= para_time_set[-1])
        true_set_time_used = true_set_time[true_set_used_ind]
        true_set_used      = true[true_set_used_ind]
        # # TODO: fix the lag in plotting the true values
        # print(true_set_used_ind)
        # print(np.where(true_set_used_ind)[0]-12)
        # true_set_used      = true[np.where(true_set_used_ind)[0]-24]
        # TODO: fix the true flux
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (para_time_set[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, 0, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # If plotting both prior and posterior
        if option == 'both':
            # Plot the prior
            ax1 = axes[0]
            for j in range(nens):
                prior_ens = analyzed_prior_ens[j, :]
                line1, = ax1.plot(para_time_set[plot_time_offset:], prior_ens[plot_time_offset:], color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            prior_mean = np.mean(analyzed_prior_ens, axis=(0))
            line2, = ax1.plot(para_time_set[plot_time_offset:], prior_mean[plot_time_offset:], color='red', linewidth=1, label='mean')
            # line3, = ax1.plot(true_set_time_used[plot_time_offset:], true_set_used[plot_time_offset:], 
            line3, = ax1.plot(para_time_set[plot_time_offset:], true_set_used_ave[plot_time_offset:], color='black', linewidth=1, label='obs')

            # Plot the posterior
            ax2 = axes[1]
            for j in range(nens):
                posterior_ens = analyzed_posterior_ens[j, :]
                line1, = ax2.plot(para_time_set[plot_time_offset:], posterior_ens[plot_time_offset:], color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            posterior_mean = np.mean(analyzed_posterior_ens, axis=(0))
            line2, = ax2.plot(para_time_set[plot_time_offset:], posterior_mean[plot_time_offset:], color='red', linewidth=1, label='mean')
            line3, = ax2.plot(para_time_set[plot_time_offset:], true_set_used_ave[plot_time_offset:], color='black', linewidth=1, label='obs')

            # # Plot the legends
            # plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
            #            frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))

            # Plot the labels and titles
            # ax1.set_title("Prior ({})".format(var_name))
            # ax2.set_title("Posterior ({})".format(var_name))
            ax1.set_xlabel("Time (day)")
            ax2.set_xlabel("Time (day)")
            ax1.set_ylim(ylim)
            ax2.set_ylim(ylim)
            ax1.set_xlim(xlim)
            ax2.set_xlim(xlim)
        
        # If plotting only prior
        elif option == 'prior':
            ax1 = axes[0]
            for j in range(nens):
                prior_ens = analyzed_prior_ens[j, :]
                line1, = ax1.plot(para_time_set[plot_time_offset:], prior_ens[plot_time_offset:], color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            prior_mean = np.mean(analyzed_prior_ens, axis=(0))
            line2, = ax1.plot(para_time_set[plot_time_offset:], prior_mean[plot_time_offset:], color='red', linewidth=1, label='mean')
            # line3, = ax1.plot(true_set_time_used[plot_time_offset:], true_set_used[plot_time_offset:], 
            line3, = ax1.plot(para_time_set[plot_time_offset:], true_set_used_ave[plot_time_offset:], color='black', linewidth=1, label='obs')

            # Plot the labels and titles
            ax1.set_xlabel("Time (day)")
            ax1.set_ylim(ylim)
            ax1.set_xlim(xlim)
        
        # If plotting only posterior
        elif option == 'posterior':
            # Plot the posterior
            ax2 = axes[0]
            for j in range(nens):
                posterior_ens = analyzed_posterior_ens[j, :]
                line1, = ax2.plot(para_time_set[plot_time_offset:], posterior_ens[plot_time_offset:], color='grey', linewidth=0.5, linestyle=':', label='ensemble')
            posterior_mean = np.mean(analyzed_posterior_ens, axis=(0))
            line2, = ax2.plot(para_time_set[plot_time_offset:], posterior_mean[plot_time_offset:], color='red', linewidth=1, label='mean')
            line3, = ax2.plot(para_time_set[plot_time_offset:], true_set_used_ave[plot_time_offset:], color='black', linewidth=1, label='obs')

            # Plot the labels and titles
            ax2.set_xlabel("Time (day)")
            ax2.set_ylim(ylim)
            ax2.set_xlim(xlim)

        return line1, line2, line3

    
    def compare_univar_1to1_spatial_average(self, var_name, true_file_name, axes, constrained_layout=True, 
                                           model_time_offset=0, plot_time_offset=0, ylim=None, xlim=None, option='both'):
        """Plot the temporal evolution of a spatial averaged analyzed variable against the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values
        true_set_used_ind  = (true_set_time >= para_time_set[0]) & (true_set_time <= para_time_set[-1])
        true_set_time_used = true_set_time[true_set_used_ind]
        true_set_used      = true[true_set_used_ind]
        # # TODO: fix the lag in plotting the true values
        # print(true_set_used_ind)
        # print(np.where(true_set_used_ind)[0]-12)
        # true_set_used      = true[np.where(true_set_used_ind)[0]-24]
        # TODO: fix the true flux
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (para_time_set[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, 0, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # If plotting both prior and posterior
        if option == 'both':
            # Plot the prior
            ax1 = axes[0]
            ax1.plot(xlim, ylim, 'k--')
            for j in range(nens):
                prior_ens = analyzed_prior_ens[j, :]
                line1, = ax1.plot(true_set_used_ave[plot_time_offset:], prior_ens[plot_time_offset:], '.', color='grey', label='ensemble')
            prior_mean = np.mean(analyzed_prior_ens, axis=(0))
            line2, = ax1.plot(true_set_used_ave[plot_time_offset:], prior_mean[plot_time_offset:], '.', color='red', label='mean')
            ax1.set_xlim(xlim); ax1.set_ylim(ylim)

            # Plot the posterior
            ax2 = axes[1]
            ax2.plot(xlim, ylim, 'k--')
            for j in range(nens):
                posterior_ens = analyzed_posterior_ens[j, :]
                line1, = ax2.plot(true_set_used_ave[plot_time_offset:], posterior_ens[plot_time_offset:], '.', color='grey', label='ensemble')
            posterior_mean = np.mean(analyzed_posterior_ens, axis=(0))
            line2, = ax2.plot(true_set_used_ave[plot_time_offset:], posterior_mean[plot_time_offset:], '.', color='red', label='mean')
            ax2.set_xlim(xlim); ax2.set_ylim(ylim)
        
        # If plotting only prior
        elif option == 'prior':
            ax1 = axes[0]
            ax1.plot(xlim, ylim, 'k--')
            for j in range(nens):
                prior_ens = analyzed_prior_ens[j, :]
                line1, = ax1.plot(true_set_used_ave[plot_time_offset:], prior_ens[plot_time_offset:], '.', color='grey', label='ensemble')
            prior_mean = np.mean(analyzed_prior_ens, axis=(0))
            line2, = ax1.plot(true_set_used_ave[plot_time_offset:], prior_mean[plot_time_offset:], '.', color='red', label='mean')
            ax1.set_xlim(xlim); ax1.set_ylim(ylim)
        
        # If plotting only posterior
        elif option == 'posterior':
            ax2 = axes[0]
            ax2.plot(xlim, ylim, 'k--')
            for j in range(nens):
                posterior_ens = analyzed_posterior_ens[j, :]
                line1, = ax2.plot(true_set_used_ave[plot_time_offset:], posterior_ens[plot_time_offset:], '.', color='grey', label='ensemble')
            posterior_mean = np.mean(analyzed_posterior_ens, axis=(0))
            line2, = ax2.plot(true_set_used_ave[plot_time_offset:], posterior_mean[plot_time_offset:], '.', color='red', label='mean')
            ax2.set_xlim(xlim); ax2.set_ylim(ylim)

        return line1, line2
    


    def compare_univar_spatial_average_diff(self, var_name, true_file_name, ax, plot_time_offset=0,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the mean and std of the difference between a spatial averaged analyzed variable
           and the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, 0, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # Get the difference between the estimated and the true for each realization
        diff      = analyzed_posterior_ens[:, :] - true_set_used_ave[:]
        diff_mean = np.mean(diff, axis=0)
        diff_std  = np.std(diff, axis=0)
        upper, bottom = diff_mean+diff_std, diff_mean-diff_std

        # Plot the difference
        ax.plot(para_time_set[plot_time_offset:], diff_mean[plot_time_offset:], color='blue', linewidth=0.5, label='mean of the difference')
        ax.fill_between(para_time_set[plot_time_offset:], bottom[plot_time_offset:], upper[plot_time_offset:], color='blue', alpha=0.1)
        ax.axhline(y=0, color='black', linestyle='--')

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        # ax.set_ylabel("ensemble mean - the true")
        ax.set_ylabel(unit)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return diff_mean, analyzed_posterior_ens.mean(axis=0), true_set_used_ave, true, true_set_time, true_set_dates

    
    def compute_univar_spatial_average_diff(self, var_name, true_file_name, model_time_offset=0., method='mae'):
        """Compute the bias between a spatial averaged analyzed variable and the true values from other source
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, 0, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))


        # Get the difference between the estimated and the true for each realization
        diff_prior = analyzed_prior_ens[:, :] - true_set_used_ave[:]
        diff_post = analyzed_posterior_ens[:, :] - true_set_used_ave[:]
        if method.lower() == 'bias':
            diff_prior = np.mean(diff_prior, axis=0)
            diff_post = np.mean(diff_post, axis=0)
        elif method.lower() == 'mae':
            diff_prior = np.mean(np.abs(diff_prior), axis=0)
            diff_post = np.mean(np.abs(diff_post), axis=0)
        
        return diff_prior, diff_post

        
    def compare_univar_spatial_average_bias(self, var_name, true_file_name, axes, plot_time_offset=0,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the bias between a spatial averaged analyzed variable
           and the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, 0, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens     = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))


        # Get the difference between the estimated and the true for each realization
        diff_prior = analyzed_prior_ens[:, :] - true_set_used_ave[:]
        bias_prior = np.mean(diff_prior, axis=0)
        mae_prior  = np.mean(np.abs(diff_prior), axis=0)
        diff_post = analyzed_posterior_ens[:, :] - true_set_used_ave[:]
        bias_post = np.mean(diff_post, axis=0)
        mae_post  = np.mean(np.abs(diff_post), axis=0)

        # Plot the bias and MAE
        # Prior
        ax = axes[0]
        line1, = ax.plot(para_time_set[plot_time_offset:], bias_prior[plot_time_offset:], color='blue', linewidth=1, label='BIAS')
        line2, = ax.plot(para_time_set[plot_time_offset:], mae_prior[plot_time_offset:], color='red', linestyle='--', linewidth=1, label='MAE')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{|BIAS|}: %.3f;$ $\mu_{MAE}: %.3f$" % (np.mean(np.abs(bias_prior[plot_time_offset:])), np.mean(mae_prior[plot_time_offset:])))
        # Plot the labels and titles
        # ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # Posterior
        ax = axes[1]
        line1, = ax.plot(para_time_set[plot_time_offset:], bias_post[plot_time_offset:], color='blue', linewidth=1, label='BIAS')
        line2, = ax.plot(para_time_set[plot_time_offset:], mae_post[plot_time_offset:], color='red', linestyle='--', linewidth=1, label='MAE')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{|BIAS|}: %.3f;$ $\mu_{MAE}: %.3f$" % (np.mean(np.abs(bias_post[plot_time_offset:])), np.mean(mae_post[plot_time_offset:]))) 
        # Plot the labels and titles
        # ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return line1, line2


    def compare_univar_spatial_average_variance(self, var_name, true_file_name, ax, plot_time_offset=None,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the variance of a spatial averaged analyzed variable.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # Plot the variance
        variance = np.var(analyzed_posterior_ens, axis=0)
        line1, = ax.plot(para_time_set[plot_time_offset:], variance[plot_time_offset:], color='black', linewidth=1, label='VAR')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{VAR}: %.3f;$" % (np.mean(variance[plot_time_offset:])))

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return line1


    def compute_univar_spatial_average_sd(self, var_name, model_time_offset=0.):
        """Compute the std of a spatial averaged analyzed variable.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_prior_ens = np.mean(prior_para[var_ind, -1, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_prior_ens = np.mean(prior_para[var_ind, :, :, :], axis=(2))
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # Plot the std
        std_prior = np.std(analyzed_prior_ens, axis=0)
        std_post = np.std(analyzed_posterior_ens, axis=0)
        return std_prior, std_post
    

    def compare_univar_spatial_average_sd(self, var_name, true_file_name, ax, plot_time_offset=None,
                                           model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the std of a spatial averaged analyzed variable.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        nvar, nens = self.nvar, self.nens
        prior_para, posterior_para = self.prior["para"], self.posterior["para"]

        ntime_para    = self.ntime_para
        nvar_para     = self.nvar_para
        para_var_set  = self.para_var_set
        para_time_set = self.para_time_set
        assim_start_set = self.assim_start_set
        assim_end_set = self.assim_end_set
        tunits        = self.tunits
        model_start_str = self.model_start_time
        has_immediate_mda_results = self.has_immediate_mda_results

        if var_name not in para_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(model_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_start_set[i]+model_time_offset)) &
                                          (true_set_time <= (assim_end_set[i]+model_time_offset))])
                            #  if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(para_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = para_var_set.index(var_name)
        if has_immediate_mda_results:
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, -1, :, :, :], axis=(2))
        else:
            analyzed_posterior_ens = np.mean(posterior_para[var_ind, :, :, :], axis=(2))

        # Plot the variance
        variance = np.std(analyzed_posterior_ens, axis=0)
        line1, = ax.plot(para_time_set[plot_time_offset:], variance[plot_time_offset:], color='black', linewidth=1, label='SD')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{SD}: %.3f;$" % (np.mean(variance[plot_time_offset:])))

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return line1