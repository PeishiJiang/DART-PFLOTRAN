"""Plot the data assimilation results."""

# Author: Peishi Jiang

import os
import re
import sys
import f90nml
import warnings
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import gridspec

# TODO: A more sophisticated OOP architecture will be required in the future to avoid the many-lines coding.


class DaResults(object):

    def __init__(self, config_nml):
        """Initialization: Read in the required configurations from config_nml."""
        # Parse the configuration in Fortran namelist
        if not os.path.isfile(config_nml):
            raise Exception(
                "The configuration file of this application does not exist under the path: %s!" % config_nml)
        self.config_nml = config_nml
        self.configs = f90nml.read(config_nml)

        # DART prior and posterior files
        self.dart_posterior_file = self.configs["file_cfg"]["dart_posterior_nc_file"]
        self.dart_prior_file     = self.configs["file_cfg"]["dart_prior_nc_file"]
        self.dart_prior_template = self.configs["file_cfg"]["dart_prior_template_file"]

        # Variables
        self.obs_set             = self.configs["obspara_set_cfg"]["obs_set"]
        self.para_set            = self.configs["obspara_set_cfg"]["para_set"]
        if isinstance(self.obs_set, str):
            self.obs_set = [self.obs_set]
        if isinstance(self.para_set, str):
            self.para_set = [self.para_set]
        self.pflotran_var_set  = self.obs_set + self.para_set
        self.nvar              = len(self.pflotran_var_set)

        # Model time, assimilation start time, number of ensemble
        self.assim_start_time    = self.configs["time_cfg"]["assim_start"]
        self.model_time          = float(self.configs["time_cfg"]["current_model_time"])   # days
        self.model_time_list     = self.configs["time_cfg"]["model_time_list"]
        self.model_time_list     = [self.model_time_list] if not isinstance(self.model_time_list, list) else self.model_time_list
        self.exceeds_obs_time    = self.configs["time_cfg"]["exceeds_obs_time"]
        self.assim_window        = float(self.configs["da_cfg"]["assim_window_size"])
        self.nens                = self.configs["da_cfg"]["nens"]
        self.ntime               = len(self.model_time_list)
        self.model_start_time    = self.model_time_list[0]
        self.model_end_time      = self.model_time_list[-1]
        # self.model_end_time      = self.model_time_list[-1] + self.assim_window

        # Observation file in NetCDF
        self.obs_nc = self.configs["file_cfg"]["obs_nc_file"]

        if not self.exceeds_obs_time:
            warnings.warn("The data assimilation is not completed!")

        # Get the model spatial domain
        root_template = Dataset(self.dart_prior_template, 'r')
        self.x_loc = root_template.variables['x_location'][:]
        self.y_loc = root_template.variables['y_location'][:]
        self.z_loc = root_template.variables['z_location'][:]
        self.nx, self.ny, self.nz = len(self.x_loc), len(self.y_loc), len(self.z_loc)
        root_template.close()

    def setup(self):
        """Read in the observation, prior, and posterior data"""
        # Get the parameters
        obs_nc            = self.obs_nc
        nvar, nens, ntime = self.nvar, self.nens, self.ntime
        nx, ny, nz        = self.nx, self.ny, self.nz
        obs_set           = self.obs_set
        pflotran_var_set  = self.pflotran_var_set
        model_start_time  = self.model_start_time
        model_end_time    = self.model_end_time
        dart_prior_file, dart_posterior_file = self.dart_prior_file, self.dart_posterior_file

        # Read in the prior and posterior data
        prior     = np.zeros([nvar, nens, ntime, nz, ny, nx])
        posterior = np.zeros([nvar, nens, ntime, nz, ny, nx])
        for i in range(nens):
            ens = i + 1
            for j in range(ntime):
                model_time_ind = j + 1
                # Prior data
                prior_nc_file = re.sub(r"\[ENS\]", str(ens)+"_time"+str(model_time_ind),
                                       dart_prior_file)
                root_prior = Dataset(prior_nc_file, 'r')
                for k in range(nvar):
                    varn = pflotran_var_set[k]
                    prior_var = root_prior.variables[varn][:]
                    # print(prior_var.shape, nz, ny, nx)
                    prior[k, i, j, :, :, :] = prior_var
                root_prior.close()
                # Posterior data
                posterior_nc_file = re.sub(r"\[ENS\]", str(ens)+"_time"+str(model_time_ind),
                                           dart_posterior_file)
                root_posterior = Dataset(posterior_nc_file, 'r')
                for k in range(nvar):
                    varn = pflotran_var_set[k]
                    posterior_var = root_posterior.variables[varn][:]
                    posterior[k, i, j, :, :, :] = posterior_var
                root_posterior.close()

        # TODO: Get the true time data
        # Read in the observation data
        obs_value_set      = dict.fromkeys(obs_set)
        obs_value_set_used = dict.fromkeys(obs_set)
        root_obs = Dataset(obs_nc, 'r')
        # Time in observation
        obs_time_set = root_obs.variables['time'][:]
        obs_used_ind = (obs_time_set >= model_start_time) & (obs_time_set <= model_end_time)
        obs_time_set_used = obs_time_set[obs_used_ind]
        # Locations in observation
        xloc = root_obs.variables['x_location']
        yloc = root_obs.variables['y_location']
        zloc = root_obs.variables['z_location']
        obs_loc = np.array([xloc, yloc, zloc])
        # Values in observation
        for i in range(len(obs_set)):
            obs_var = obs_set[i]
            obs_value_set[obs_var]      = root_obs.variables[obs_var][:]
            obs_value_set_used[obs_var] = obs_value_set[obs_var][:, obs_used_ind]
        root_obs.close()

        self.prior, self.posterior = prior, posterior
        self.obs_value_set_used    = obs_value_set_used
        self.obs_time_set_used     = obs_time_set_used
        self.obs_loc_set           = obs_loc

    def plot_spatial_average(self, figsize=None, constrained_layout=True, ylim=None):
        """Plot the spatial averaged DA results"""
        # Get the parameter
        nvar, nens, ntime  = self.nvar, self.nens, self.ntime
        obs_set            = self.obs_set
        pflotran_var_set   = self.pflotran_var_set
        model_time_list    = self.model_time_list
        prior, posterior   = self.prior, self.posterior
        obs_time_set_used  = self.obs_time_set_used
        obs_value_set_used = self.obs_value_set_used

        # Plot
        if figsize is not None:
            fig = plt.figure(num=1, dpi=150, figsize=figsize, constrained_layout=constrained_layout)
        else:
            fig = plt.figure(num=1, dpi=150, figsize=(12, 6 * nvar), constrained_layout=constrained_layout)
        gs  = gridspec.GridSpec(nvar, 2, width_ratios=[1,1])

        # Define axes array
        axes = np.empty([nvar, 2], dtype=object)
        for i in range(nvar):
            varn = pflotran_var_set[i]
            # Plot the prior
            ax1 = plt.subplot(gs[i, 0])
            for j in range(nens):
                prior_ens_mean = np.mean(prior[i, j, :, :, :, :], axis=(1, 2, 3))
                line1, = ax1.plot(model_time_list, prior_ens_mean, color='grey',
                        linewidth=0.5, linestyle=':', label='ensemble')
            prior_mean = np.mean(prior[i, :, :, :, :, :], axis=(0, 2, 3, 4))
            line2, = ax1.plot(model_time_list, prior_mean, color='red', linewidth=2, label='mean')
            if varn in obs_set:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                line3, = ax1.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=2, label='obs')

            # Plot the posterior
            ax2 = plt.subplot(gs[i, 1], sharey=ax1, sharex=ax1)
            # ax2 = plt.subplot(gs[i, 1])
            for j in range(nens):
                posterior_ens_mean = np.mean(posterior[i, j, :, :, :, :], axis=(1, 2, 3))
                ax2.plot(model_time_list, posterior_ens_mean, color='grey',
                        linewidth=0.5, linestyle=':', label='ensemble')
            posterior_mean = np.mean(posterior[i, :, :, :, :, :], axis=(0, 2, 3, 4))
            ax2.plot(model_time_list, posterior_mean, color='red', linewidth=2, label='mean')
            if varn in obs_set:
                obs_used_mean = np.mean(obs_value_set_used[varn], axis=0)
                ax2.plot(obs_time_set_used, obs_used_mean, color='black', linewidth=2, label='obs')
            ax2.set_ylim(ylim)

            axes[i, :] = [ax1, ax2]

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
        plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
                   frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.2))

    def plot_oned_obs(self, obs_name, dim_str='z', figsize=None, constrained_layout=True,
                      vmin=None, vmax=None, ylim=None):
        """Plot the temporal evolution of DA results for the observation data along one dimension"""
        # Get the parameter
        nens, ntime        = self.nens, self.ntime
        obs_set            = self.obs_set
        pflotran_var_set   = self.pflotran_var_set
        model_time_list    = self.model_time_list
        prior, posterior   = self.prior, self.posterior
        obs_time_set_used  = self.obs_time_set_used
        obs_value_set_used = self.obs_value_set_used
        obs_loc_set        = self.obs_loc_set
        xloc_set, yloc_set, zloc_set = self.x_loc, self.y_loc, self.z_loc

        if obs_name not in obs_set:
            raise Exception('Unknown observation variable name %s' % obs_name)

        # Get the locations associated with this observation variable
        _, nloc = obs_loc_set.shape
        # nvar_obs = len(obs_set)

        # Plot
        if figsize is not None:
            fig = plt.figure(num=1, dpi=150, figsize=figsize, constrained_layout=constrained_layout)
        else:
            fig = plt.figure(num=1, dpi=150, figsize=(12, 6 * (nloc+1)), constrained_layout=constrained_layout)
        gs  = gridspec.GridSpec(nloc+2, 2, width_ratios=[1, 1])

        # Define axes array
        axes = np.empty([nloc+1, 2], dtype=object)

        # Get the index of observation
        obs_ind = pflotran_var_set.index(obs_name)

        # Get the averaged dimensions
        dim        = 0 if dim_str == 'x' else 1 if dim_str == 'y' else 2
        ave_dims   = [2, 1] if dim_str == 'z' else [2, 0] if dim_str == 'y' else [1, 0]
        ave_dims_d = [e+2 for e in ave_dims]

        ##############################
        # Plot the temporal evolution of the mean along the given dimension
        dim_set = xloc_set if dim_str == 'x' else yloc_set if dim_str == 'y' else zloc_set
        timev, locv = np.meshgrid(model_time_list, dim_set)
        # max_m, min_m = np.max([prior.max(), posterior.max()]), np.min([prior.min(), posterior.min()])
        # Plot the prior
        ax1 = plt.subplot(gs[:2, 0])
        prior_mean = np.mean(prior[obs_ind, :, :, :, :, :], axis=(0,)+tuple(ave_dims_d))
        cs = ax1.contourf(timev, locv, prior_mean.T, vmin=vmin, vmax=vmax)
        for i in range(nloc):
            ax1.axhline(y=obs_loc_set[dim, i], color='k', ls='--')
        ax1.set_yticks(obs_loc_set[dim, :])
        plt.colorbar(cs, ax=ax1, orientation='horizontal', fraction=0.1)

        # Plot the posterior
        ax2 = plt.subplot(gs[:2, 1], sharey=ax1, sharex=ax1)
        posterior_mean = np.mean(posterior[obs_ind, :, :, :, :, :], axis=(0,)+tuple(ave_dims_d))
        cs = ax2.contourf(timev, locv, posterior_mean.T, vmin=vmin, vmax=vmax)
        for i in range(nloc):
            ax2.axhline(y=obs_loc_set[dim, i], color='k', ls='--')
        ax2.set_yticks(obs_loc_set[dim, :])
        plt.colorbar(cs, ax=ax2, orientation='horizontal', fraction=0.1)

        axes[0, :] = [ax1, ax2]

        ##############################
        # Plot the temporal evolution of the ensemble, the mean, and the observation
        # at each observed location
        for i in range(nloc):
            # Get the locations of observation and the corresponding index of DA results
            obs_xloc, obs_yloc, obs_zloc = obs_loc_set[:, i]
            # print(np.where(xloc_set == obs_xloc))
            # print(zloc_set, obs_zloc)
            # xind = np.where(xloc_set == obs_xloc)[0]
            # yind = np.where(yloc_set == obs_yloc)[0]
            # zind = np.where(zloc_set == obs_zloc)[0]
            xind = np.argmin(np.abs(xloc_set - obs_xloc))
            yind = np.argmin(np.abs(yloc_set - obs_yloc))
            zind = np.argmin(np.abs(zloc_set - obs_zloc))
            # print(obs_xloc, obs_yloc, obs_zloc, xind, yind, zind)

            # Plot the prior
            ax1 = plt.subplot(gs[i+2, 0], sharex=ax2)
            for j in range(nens):
                # prior_ens_mean = np.mean(prior[obs_ind, j, :, xind, yind, zind])
                prior_ens = prior[obs_ind, j, :, zind, yind, xind]
                line1, = ax1.plot(model_time_list, prior_ens, color='grey',
                         linewidth=0.5, linestyle=':', label='ensemble')
            prior_mean = np.mean(prior[obs_ind, :, :, zind, yind, xind], axis=(0))
            line2, = ax1.plot(model_time_list, prior_mean, color='red', linewidth=2, label='mean')
            obs_used = obs_value_set_used[obs_name][i, :]
            line3, = ax1.plot(obs_time_set_used, obs_used, color='black', linewidth=2, label='obs')
            # for k in range(ntime):
            #     ax1.axvline(x=model_time_list[k], color='k', ls='--')

            # Plot the posterior
            ax2 = plt.subplot(gs[i+2, 1], sharey=ax1, sharex=ax1)
            for j in range(nens):
                posterior_ens = posterior[obs_ind, j, :, zind, yind, xind]
                ax2.plot(model_time_list, posterior_ens, color='grey',
                         linewidth=0.5, linestyle=':', label='ensemble')
            posterior_mean = np.mean(posterior[obs_ind, :, :, zind, yind, xind], axis=(0))
            ax2.plot(model_time_list, posterior_mean, color='red', linewidth=2, label='mean')
            obs_used = obs_value_set_used[obs_name][i, :]
            ax2.plot(obs_time_set_used, obs_used, color='black', linewidth=2, label='obs')
            # for k in range(ntime):
            #     ax2.axvline(x=model_time_list[k], color='k', ls='--')
            ax2.set_ylim(ylim)

            axes[i+1, :] = [ax1, ax2]

        # Plot the legends
        # plt.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.9, -0.3))
        plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
                   frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5))

        # Plot the labels and titles
        axes[0, 0].set_title("Prior")
        axes[0, 1].set_title("Posterior")
        axes[-1, 0].set_xlabel("Time (day)")
        axes[-1, 1].set_xlabel("Time (day)")
        axes[0, 0].set_ylabel("Dimension %s (m)" % (dim_str,))
        for i in range(nloc):
            axes[i+1, 0].set_ylabel("Dimension %s: \n %.2f (m)" % (dim_str, obs_loc_set[dim, i]))

