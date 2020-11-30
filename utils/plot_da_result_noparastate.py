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

# TODO: A more sophisticated OOP architecture will be required in the future to avoid the many-lines coding.

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
        self.assim_window     = float(self.configs["da_cfg"]["assim_window_size"])
        self.nens             = self.configs["da_cfg"]["nens"]
        self.model_ntime      = len(self.model_time_list)
        self.model_ntimestep  = int(self.configs["da_cfg"]["ntimestep"])
        # self.model_start_time = self.model_time_list[0]
        # self.model_end_time   = self.model_time_list[-1]
        # self.model_end_time   = self.model_time_list[-1] + self.assim_window
        # self.ndigit_time = int(ceil(log10(self.ntimestep))) + 1
        # self.ndigit_ens = int(ceil(log10(self.nens))) + 1

        # Convert the model time to time units
        self.tunits                = tunits_begin + self.assim_start_time
        self.model_time_dates_list = num2date(self.model_time_list, self.tunits, calendar)

        # Observation file in NetCDF
        self.obs_nc = self.configs["file_cfg"]["obs_nc_file"]

        if not self.exceeds_obs_time:
            warnings.warn("The data assimilation is not completed!")

        # Get the model spatial domain
        self.dart_prior_template = self.configs["file_cfg"]["dart_prior_template_file"]
        root_template = Dataset(self.dart_prior_template, 'r')
        x_loc    = root_template.variables['x_location'][:]
        y_loc    = root_template.variables['y_location'][:]
        z_loc    = root_template.variables['z_location'][:]
        self.model_loc_set = np.array([x_loc, y_loc, z_loc]).T
        print(self.model_loc_set.shape)
        self.model_nloc = self.model_loc_set.shape[0]
        root_template.close()


    def setup(self):
        """Read in the observation, prior, and posterior data"""
        # Get the parameters
        obs_nc           = self.obs_nc
        obs_var_set      = self.obs_var_set
        para_var_set     = self.para_var_set

        nvar_obs, nvar_para = len(obs_var_set), len(para_var_set)
        nvar, nens          = self.nvar, self.nens
        nloc                = self.model_nloc

        pflotran_var_set = self.pflotran_var_set
        # model_start_time = self.model_start_time
        # model_end_time   = self.model_end_time

        # DART prior and posterior files
        self.dart_posterior_all_ens_file = self.configs["file_cfg"]["dart_posterior_nc_all_ens_file"]
        self.dart_prior_all_ens_file     = self.configs["file_cfg"]["dart_prior_nc_all_ens_file"]
        dart_prior_all_ens_file, dart_posterior_all_ens_file = self.dart_prior_all_ens_file, self.dart_posterior_all_ens_file

        # prior     = dict.fromkeys(["para","state"])
        # posterior = dict.fromkeys(["para","state"])

        # Get the assimilated time
        with Dataset(dart_prior_all_ens_file, 'r') as root_prior:
            assim_time = root_prior.variables['state_time'][:]
            assim_ntime = len(assim_time)

        # Read in the prior data
        prior     = np.zeros([nvar, nens, assim_ntime, nloc])
        with Dataset(dart_prior_all_ens_file, 'r') as root_prior:
            for k in range(nvar):
                varn              = pflotran_var_set[k]
                prior[k, :, :, :] = root_prior.variables[varn][:]

        # Read in the posterior data
        posterior = np.zeros([nvar, nens, assim_ntime, nloc])
        with Dataset(dart_posterior_all_ens_file, 'r') as root_posterior:
            for k in range(nvar):
                varn                  = pflotran_var_set[k]
                posterior[k, :, :, :] = root_posterior.variables[varn][:]

        # TODO: Get the true time data
        # Read in the observation data
        obs_value_set      = dict.fromkeys(obs_var_set)
        obs_value_set_used = dict.fromkeys(obs_var_set)
        with Dataset(obs_nc, 'r') as root_obs:
            # Time in observation
            obs_time_set = root_obs.variables['time'][:]
            obs_used_ind = (obs_time_set >= assim_time[0]) & (obs_time_set <= assim_time[-1])
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
        self.assim_time_set, self.assim_ntime = assim_time, assim_ntime
        self.obs_value_set_used = obs_value_set_used
        self.obs_time_set_used = obs_time_set_used
        self.obs_loc_set = obs_loc


    def plot_spatial_average(self, axes, ylim=None, plot_averaged_obs=False):
        """Plot the spatial averaged DA results"""
        # Get the parameter
        nvar, nens, ntime = self.nvar, self.nens, self.assim_ntime
        prior, posterior = self.prior, self.posterior

        obs_var_set           = self.obs_var_set
        pflotran_var_set      = self.pflotran_var_set
        assim_time_set        = self.assim_time_set
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used

        # Plot
        for i in range(nvar):
            varn = pflotran_var_set[i]
            # Plot the prior
            ax1 = axes[i, 0]
            for j in range(nens):
                prior_ens_mean = np.mean(prior[i, j, :, :], axis=(1))
                line1, = ax1.plot(assim_time_set, prior_ens_mean, color='grey',
                                  linewidth=0.5, linestyle=':', label='ensemble')
            prior_mean = np.mean(prior[i, :, :, :], axis=(0, 2))
            line2, = ax1.plot(assim_time_set, prior_mean, color='red',
                              linewidth=1, label='mean')
            if varn in obs_var_set and plot_averaged_obs:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                line3, = ax1.plot(obs_time_set_used, obs_used_mean,
                                  color='black', linewidth=1, label='obs')

            # Plot the posterior
            ax2 = axes[i, 1]
            for j in range(nens):
                posterior_ens_mean = np.mean(posterior[i, j, :, :], axis=(1))
                ax2.plot(assim_time_set, posterior_ens_mean, color='grey',
                         linewidth=0.5, linestyle=':', label='ensemble')
            posterior_mean = np.mean(posterior[i, :, :, :], axis=(0, 2))
            ax2.plot(assim_time_set, posterior_mean, color='red',
                     linewidth=1, label='mean')
            if varn in obs_var_set and plot_averaged_obs:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                ax2.plot(obs_time_set_used, obs_used_mean, color='black',
                         linewidth=1, label='obs')
            ax2.set_ylim(ylim)

        # Plot the labels and titles
        for i in range(nvar):
            ax = axes[i, 0]
            ax.set_ylabel(pflotran_var_set[i])
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


    def plot_obs_at_point(self, axes, obs_name, obs_loc_ind=0,
                        #   figsize=None, constrained_layout=True,
                          vmin=None, vmax=None, ylim=None):
        """Plot the temporal evolution of DA results for the observation data along one dimension"""
        # Get the parameter
        nens , ntime     = self.nens,  self.assim_ntime
        prior, posterior = self.prior, self.posterior

        obs_var_set           = self.obs_var_set
        pflotran_var_set      = self.pflotran_var_set
        assim_time_set        = self.assim_time_set
        # model_time_list       = self.model_time_list
        # model_time_dates_list = self.model_time_dates_list
        obs_time_set_used     = self.obs_time_set_used
        obs_value_set_used    = self.obs_value_set_used
        obs_loc_set           = self.obs_loc_set
        model_loc_set         = self.model_loc_set

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
        obs_var_ind = pflotran_var_set.index(obs_name)

        ##############################
        # Plot the temporal evolution of the ensemble, the mean, and the observation
        # at the given observed location
        # Get the locations of observation and the corresponding index of DA results
        # Plot the prior
        ax1 = axes[0]
        for j in range(nens):
            prior_ens = prior[obs_var_ind, j, :, model_loc_ind]
            line1, = ax1.plot(assim_time_set, prior_ens, color='grey', linewidth=0.5,
                              linestyle=':', label='ensemble')
        prior_mean = np.mean(prior[obs_var_ind, :, :, model_loc_ind], axis=(0))
        line2, = ax1.plot(assim_time_set, prior_mean, color='red',
                          linewidth=1, label='mean')
        line3, = ax1.plot(obs_time_set_used, obs_value, color='black', linewidth=1, label='obs')
        ax1.set_ylim(ylim)

        # Plot the posterior
        ax2 = axes[1]
        for j in range(nens):
            posterior_ens = posterior[obs_var_ind, j, :, model_loc_ind]
            ax2.plot(assim_time_set, posterior_ens, color='grey', linewidth=0.5,
                     linestyle=':', label='ensemble')
        posterior_mean = np.mean(posterior[obs_var_ind, :, :, model_loc_ind], axis=(0))
        ax2.plot(assim_time_set, posterior_mean, color='red',
                 linewidth=1, label='mean')
        ax2.plot(obs_time_set_used, obs_value, color='black', linewidth=1, label='obs')
        ax2.set_ylim(ylim)

        return line1, line2, line3

        # Plot the legends
        # plt.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.9, -0.3))
        # plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
        #            frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))

        # # Plot the labels and titles
        # axes[0, 0].set_title("Prior")
        # axes[0, 1].set_title("Posterior")
        # axes[-1, 0].set_xlabel("Time (day)")
        # axes[-1, 1].set_xlabel("Time (day)")
        # axes[0, 0].set_ylabel("Dimension %s (m)" % (dim_str, ))
        # for i in range(nloc):
        #     axes[i + 1, 0].set_ylabel("Dimension %s: \n %.2f (m)" %
        #                               (dim_str, obs_loc_set[dim, i]))

    # TODO
    def compare_univar_spatial_average(self, var_name, true_file_name, axes, constrained_layout=True, 
                                       model_time_offset=0, plot_time_offset=0, ylim=None, xlim=None):
        """Plot the temporal evolution of a spatial averaged analyzed variable against the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        # Get the parameter
        nens , ntime     = self.nens,  self.assim_ntime
        prior, posterior = self.prior, self.posterior

        model_start_str       = self.model_start_time
        pflotran_var_set      = self.pflotran_var_set
        assim_time_set        = self.assim_time_set
        # model_time_list       = self.model_time_list
        # model_time_dates_list = self.model_time_dates_list
        # model_start_time      = self.model_start_time
        # model_end_time        = self.model_end_time
        tunits                = self.tunits

        if var_name not in pflotran_var_set:
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
        true_set_used_ind  = (true_set_time >= assim_time_set[0]) & (true_set_time <= assim_time_set[-1])
        true_set_time_used = true_set_time[true_set_used_ind]
        true_set_used      = true[true_set_used_ind]
        # # TODO: fix the lag in plotting the true values
        # print(true_set_used_ind)
        # print(np.where(true_set_used_ind)[0]-12)
        # true_set_used      = true[np.where(true_set_used_ind)[0]-24]
        true_set_used_ave = [np.mean(true[(true_set_time >  (assim_time_set[i-1]+model_time_offset)) &
                                          (true_set_time <= (assim_time_set[i]+model_time_offset))])
                             if i != 0 else np.mean(true[(true_set_time <= (assim_time_set[i]+model_time_offset))])
                             for i in range(len(assim_time_set))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = pflotran_var_set.index(var_name)
        analyzed_prior_ens     = np.mean(prior[var_ind, :, :, :], axis=(2))
        analyzed_posterior_ens = np.mean(posterior[var_ind, :, :, :], axis=(2))

        # Plot the prior
        ax1 = axes[0]
        for j in range(nens):
            prior_ens = analyzed_prior_ens[j, :]
            line1, = ax1.plot(assim_time_set[plot_time_offset:], prior_ens[plot_time_offset:],
                              color='grey', linewidth=0.5, linestyle=':', label='ensemble')
        prior_mean = np.mean(analyzed_prior_ens, axis=(0))
        line2, = ax1.plot(assim_time_set[plot_time_offset:], prior_mean[plot_time_offset:], 
                          color='red', linewidth=1, label='mean')
        # line3, = ax1.plot(true_set_time_used[plot_time_offset:], true_set_used[plot_time_offset:], 
        line3, = ax1.plot(assim_time_set[plot_time_offset:], true_set_used_ave[plot_time_offset:], 
                          color='black', linewidth=1, label='obs')

        # Plot the posterior
        ax2 = axes[1]
        for j in range(nens):
            posterior_ens = analyzed_posterior_ens[j, :]
            line1, = ax2.plot(assim_time_set[plot_time_offset:], posterior_ens[plot_time_offset:],
                              color='grey', linewidth=0.5, linestyle=':', label='ensemble')
        posterior_mean = np.mean(analyzed_posterior_ens, axis=(0))
        line2, = ax2.plot(assim_time_set[plot_time_offset:], posterior_mean[plot_time_offset:], 
                          color='red', linewidth=1, label='mean')
        line3, = ax2.plot(assim_time_set[plot_time_offset:], 
                          true_set_used_ave[plot_time_offset:], color='black', linewidth=1, label='obs')

        # Plot the legends
        plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
                   frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))

        # Plot the labels and titles
        # ax1.set_title("Prior ({})".format(var_name))
        # ax2.set_title("Posterior ({})".format(var_name))
        ax1.set_xlabel("Time (day)")
        ax2.set_xlabel("Time (day)")
        ax1.set_ylim(ylim)
        ax2.set_ylim(ylim)
        ax1.set_xlim(xlim)
        ax2.set_xlim(xlim)


    # TODO
    def compare_univar_spatial_average_diff(self, var_name, true_file_name, ax, plot_time_offset=0,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the mean and std of the difference between a spatial averaged analyzed variable
           and the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        # Get the parameter
        nens , ntime     = self.nens,  self.ntime
        prior, posterior = self.prior, self.posterior

        assim_start_str       = self.assim_start_time
        pflotran_var_set      = self.pflotran_var_set
        model_time_list       = self.model_time_list
        model_time_dates_list = self.model_time_dates_list
        model_start_time      = self.model_start_time
        model_end_time        = self.model_end_time
        tunits                = self.tunits

        if var_name not in pflotran_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values
        # true_set_used_ind  = (true_set_time >= model_start_time) & (true_set_time <= model_end_time)
        # true_set_time_used = true_set_time[true_set_used_ind]
        # true_set_used      = true[true_set_used_ind]
        # # TODO: fix the lag in plotting the true values
        # print(true_set_used_ind)
        # print(np.where(true_set_used_ind)[0]-12)
        # true_set_used      = true[np.where(true_set_used_ind)[0]-24]

        # Compute the temporal averaged true values
        # true_set_used_ave = [true[true_set_time<=model_time_list[i]][-1]
        #                      for i in range(len(model_time_list))]
        # true_set_used_ave = [np.mean(true[(true_set_time>model_time_list[i-1]-model_time_offset) & (true_set_time<=model_time_list[i]-model_time_offset)])
        #                      if i != 0 else np.mean(true[(true_set_time<=model_time_list[i])])
        #                      for i in range(len(model_time_list))]
        true_set_used_ave = [np.mean(true[(true_set_time >  (model_time_list[i-1]+model_time_offset)) &
                                          (true_set_time <= (model_time_list[i]+model_time_offset))])
                             if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(model_time_list))]
        # true_set_used_ave = [true[(true_set_time >  (model_time_list[i-1]+model_time_offset)) &
        #                           (true_set_time <= (model_time_list[i]+model_time_offset))][-1]
        #                      if i != 0 else true[(true_set_time <= (model_time_list[i]+model_time_offset))][-1]
        #                      for i in range(len(model_time_list))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = pflotran_var_set.index(var_name)
        analyzed_prior_ens     = np.mean(prior[var_ind, :, :, :, :, :], axis=(2, 3, 4))
        analyzed_posterior_ens = np.mean(posterior[var_ind, :, :, :, :, :], axis=(2, 3, 4))

        # Get the difference between the estimated and the true for each realization
        diff      = analyzed_posterior_ens[:, :] - true_set_used_ave[:]
        diff_mean = np.mean(diff, axis=0)
        # diff_mean = np.mean(analyzed_posterior_ens, axis=0) - true_set_used_ave
        diff_std  = np.std(diff, axis=0)
        upper, bottom = diff_mean+diff_std, diff_mean-diff_std

        # Plot the difference
        ax.plot(model_time_list[plot_time_offset:], diff_mean[plot_time_offset:], color='blue', linewidth=0.5, label='mean of the difference')
        ax.fill_between(model_time_list[plot_time_offset:], bottom[plot_time_offset:], upper[plot_time_offset:], color='blue', alpha=0.1)
        ax.axhline(y=0, color='black', linestyle='--')

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        # ax.set_ylabel("ensemble mean - the true")
        ax.set_ylabel(unit)
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return diff_mean, analyzed_posterior_ens.mean(axis=0), true_set_used_ave, true, true_set_time, true_set_dates


    # TODO
    def compare_univar_spatial_average_bias(self, var_name, true_file_name, ax, plot_time_offset=0,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the bias between a spatial averaged analyzed variable
           and the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        # Get the parameter
        nens , ntime     = self.nens,  self.ntime
        prior, posterior = self.prior, self.posterior

        assim_start_str       = self.assim_start_time
        pflotran_var_set      = self.pflotran_var_set
        model_time_list       = self.model_time_list
        model_time_dates_list = self.model_time_dates_list
        model_start_time      = self.model_start_time
        model_end_time        = self.model_end_time
        tunits                = self.tunits

        if var_name not in pflotran_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (model_time_list[i-1]+model_time_offset)) &
                                          (true_set_time <= (model_time_list[i]+model_time_offset))])
                             if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(model_time_list))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = pflotran_var_set.index(var_name)
        analyzed_prior_ens     = np.mean(prior[var_ind, :, :, :, :, :], axis=(2, 3, 4))
        analyzed_posterior_ens = np.mean(posterior[var_ind, :, :, :, :, :], axis=(2, 3, 4))

        # Get the difference between the estimated and the true for each realization
        diff = analyzed_posterior_ens[:, :] - true_set_used_ave[:]
        bias = np.mean(diff, axis=0)
        mae  = np.mean(np.abs(diff), axis=0)
        # mae  = np.abs(bias)

        # Plot the bias
        line1, = ax.plot(model_time_list[plot_time_offset:], bias[plot_time_offset:], color='blue', linewidth=1, label='BIAS')
        line2, = ax.plot(model_time_list[plot_time_offset:], mae[plot_time_offset:], color='red', linestyle='--', linewidth=1, label='MAE')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{BIAS}: %.3f;$ $\mu_{MAE}: %.3f$" % (np.mean(bias[plot_time_offset:]), np.mean(mae[plot_time_offset:])))

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return line1, line2


    # TODO
    def compare_univar_spatial_average_variance(self, var_name, true_file_name, ax, plot_time_offset=None,
                                            model_time_offset=0., constrained_layout=True, ylim=None, xlim=None, unit=''):
        """Plot the temporal evolution of the bias between a spatial averaged analyzed variable
           and the true values from other source.
           Note that the true_file_name has to be a csv file with two columns (the first for time and the second for the values)
        """
        # Get the parameter
        nens , ntime     = self.nens,  self.ntime
        prior, posterior = self.prior, self.posterior

        assim_start_str       = self.assim_start_time
        pflotran_var_set      = self.pflotran_var_set
        model_time_list       = self.model_time_list
        model_time_dates_list = self.model_time_dates_list
        model_start_time      = self.model_start_time
        model_end_time        = self.model_end_time
        tunits                = self.tunits

        if var_name not in pflotran_var_set:
            raise Exception('Unknown analyzed variable name %s' % var_name)

        # Get the reference, start and end dates
        ref_time = datetime.strptime(assim_start_str, "%Y-%m-%d %H:%M:%S")
        # model_start_date, model_end_date = model_time_dates_list[0], model_time_dates_list[-1]

        # Read the true value from file_name
        true_set           = pd.read_csv(true_file_name)
        true_set_raw_time  = true_set.iloc[:, 0].values
        true_set_dates     = [datetime.strptime(t, '%m/%d/%Y %H:%M') for t in true_set_raw_time]
        dates_ref          = [t-ref_time for t in true_set_dates]
        true_set_time      = np.array([t.days+float(t.seconds)/86400. for t in dates_ref])
        true               = true_set.iloc[:, 1].values

        # Compute the temporal averaged true values
        true_set_used_ave = [np.mean(true[(true_set_time >  (model_time_list[i-1]+model_time_offset)) &
                                          (true_set_time <= (model_time_list[i]+model_time_offset))])
                             if i != 0 else np.mean(true[(true_set_time <= (model_time_list[i]+model_time_offset))])
                             for i in range(len(model_time_list))]
        true_set_used_ave = np.array(true_set_used_ave)

        # Get the spatially averaged analyzed variable (prior and posterior)
        var_ind                = pflotran_var_set.index(var_name)
        analyzed_posterior_ens = np.mean(posterior[var_ind, :, :, :, :, :], axis=(2, 3, 4))

        # Plot the variance
        variance = np.var(analyzed_posterior_ens, axis=0)
        line1, = ax.plot(model_time_list[plot_time_offset:], variance[plot_time_offset:], color='blue', linewidth=1, label='VAR')
        ax.axhline(y=0, color='black', linestyle='--')
        ax.set_title("$\mu_{VAR}: %.3f;$" % (np.mean(variance[plot_time_offset:])))

        # Plot the labels and titles
        ax.set_xlabel("Time (day)")
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # return
        return line1