# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'manhattan/models/pflotran/applications/workflow_ipynb'))
	print(os.getcwd())
except:
	pass

# %% [markdown]
# # Overview
# The **objective** of this notebook is to visualize the data assimilation results generated from [the DART_PFLOTRAN_Integrate notebook](./DART_PFLOTRAN_Integrate.ipynb):
# 1. [x] [Configuration](#parameter): define directories, file locations, and other parameters
# - [x] [Visualize spatial average](#plot_spatial_average): plot the time evolution of spatial-averaged prior and posterior results
# - [x] [Visualize time evolution along one axis](#plot_along_zaxis):
# %% [markdown]
# <a id='parameter'></a>
# # Configuration

# %%
import os
import sys
import pickle
import f90nml
import subprocess
import numpy as np
from math import floor
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib import gridspec

# %% [markdown]
# ****************
# **Define the locations of application folder, DART-PFLOTRAN interface folder, and the configuation file of interest**

# %%
# Main directory names
# app_dir_name = "1dthermal_test_1month_4mda_v2"
# app_par_dir  = "/mnt/4tbb/peishi/1dthermal_simulation"
# app_par_dir  = "/home/jian449/DART/manhattan/models/pflotran/applications"
# dart_dir     = "/home/jian449/DART/manhattan"
# dart_pf_dir  = "/home/jian449/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
# fig_dir      = "/home/jian449/figures"
app_par_dir  = "/Users/jian449/Codes/DART/manhattan/models/pflotran/applications"
dart_dir     = "/Users/jian449/Codes/DART/manhattan"
dart_pf_dir  = "/Users/jian449/Codes/DART/manhattan/models/pflotran"     # The dart pflotran utitlity folder name
fig_dir      = "/Users/jian449/OneDrive - PNNL/Documents/Publications/2020DART_PFLOTRAN/simulation_plots"

# The list of application folders
# app_dir_name_list = ['1dthermal_test_1month_7200s_1mda', '1dthermal_test_1month_7200s_2mda',
#                      '1dthermal_test_1month_7200s_3mda', '1dthermal_test_1month_7200s_4mda',
#                      '1dthermal_test_1month_14400s_1mda', '1dthermal_test_1month_14400s_2mda',
#                      '1dthermal_test_1month_14400s_3mda', '1dthermal_test_1month_14400s_4mda']
# app_dir_name_list = ['1dthermal_test_1month_1800s_1mda', '1dthermal_test_1month_1800s_2mda',
app_dir_name_list = ['1dthermal_parastate_1month_1mda']
# app_dir_name_list = ['1dthermal_parastate_1month_3600s_2mda', '1dthermal_parastate_1month_3600s_3mda',
#                       '1dthermal_parastate_1month_3600s_4mda']

# Plotting time offset
plot_time_offset = 0 # for two days of 5min resolution data

for app_dir_name in app_dir_name_list:
    app_dir     = os.path.join(app_par_dir, app_dir_name)          # The application folder name
    # Get the configuration file
    config_file = os.path.join(app_dir, "work/config.nml")

    print("Generate plots from application folder: {}".format(app_dir))

    # %% [markdown]
    # ****************
    # **Change the locations of the all the files saved in the original config_file according to the new defined locations if necessary (this is used for the case that the application folder is moved/copied from another location)**

    # %%
    change_path_name_files = os.path.join(dart_pf_dir, "utils/change_file_paths_in_confignml.py")
    subprocess.run("python {} {} {} {}".format(change_path_name_files, config_file, dart_pf_dir, app_dir),
                shell=True, check=True)

    # %% [markdown]
    # ****************
    # **Import the visualization functions**

    # %%
    sys.path.append(dart_pf_dir)
    from utils.plot_da_result import DaResults
    dar = DaResults(config_file)
    # dar.obs_nc = os.path.join(app_dir, "pflotran_input/obs_true_pflotran_clipped_tafterspinup.nc") # Put the true data instead of observation data
    dar.setup()
    nvar = dar.nvar

    # %%
    # fig = plt.figure(num=1,dpi=150, figsize=(10,5), constrained_layout=True)
    # gs = gridspec.GridSpec(nvar, 2, width_ratios=[1, 1])
    # axes = np.empty([nvar, 2], dtype=object)
    # for i in range(nvar):
    #     axes[i, 0] = plt.subplot(gs[i, 0])
    #     axes[i, 1] = plt.subplot(gs[i, 1], sharey=axes[i, 0], sharex=axes[i, 0])
    # dar.plot_spatial_average(axes)
    # plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux.png"), dpi=150)
    # plt.close()
    # %%
    # Plot assimilated temperature profile
    nloc_obs = len(dar.obs_loc_set)
    fig = plt.figure(num=1,dpi=150,figsize=(12,8),constrained_layout=True)
    gs = gridspec.GridSpec(nloc_obs, 2, width_ratios=[1, 1])
    for i in range(nloc_obs):
        axes = []
        axes.append(plt.subplot(gs[i,0]));axes.append(plt.subplot(gs[i,1]))
        line1,line2,line3 = dar.plot_obs_at_point('TEMPERATURE', axes, obs_loc_ind=0, vmin=4, vmax=10, ylim=[4, 10])
        if i == 0: axes[0].set_title('Prior'); axes[1].set_title('Posterior')
    plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'),
           frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.0, -0.5));
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_temperature.png"), dpi=150)
    plt.close()

    # %%
    flux_file = os.path.join(app_dir, "pflotran_input/flux_1d_new.csv")

    # %%
    # Flow flux
    fig = plt.figure(num=1,dpi=150, figsize=(5,4), constrained_layout=False)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], figure=fig)
    axes = np.empty(2, dtype=object)
    axes[0] = plt.subplot(gs[0, 0])
    axes[1] = plt.subplot(gs[1, 0], sharey=axes[0], sharex=axes[0])
    line1,line2,line3 = dar.compare_univar_spatial_average(var_name='FLOW_FLUX', true_file_name=flux_file, axes=axes, ylim=[-4, 4])
    plt.subplots_adjust(wspace=0.2)
    plt.legend((line1, line2, line3), ('ensemble', 'mean', 'obs'), frameon=False, ncol=3, loc="center", bbox_to_anchor=(0.5, -0.5))
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux.png"), dpi=150)
    plt.close()

    # %%
    # Plot the difference between true and updated (flow flux)
    _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    dar.compare_univar_spatial_average_diff(var_name='FLOW_FLUX', true_file_name=flux_file, ax=ax,
                                            model_time_offset=0., ylim=[-1.8, 1.8], xlim=None, plot_time_offset=plot_time_offset)
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_diff.png"), dpi=150, bbox_inches='tight')
    plt.close()

    # %%
    # Plot the bias and mae between true and updated (flow flux)
    # _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    fig = plt.figure(num=1,dpi=150, figsize=(5,4), constrained_layout=True)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], figure=fig)
    axes = np.empty(2, dtype=object)
    axes[0] = plt.subplot(gs[0, 0])
    axes[1] = plt.subplot(gs[1, 0], sharey=axes[0], sharex=axes[0])
    line1, line2 = dar.compare_univar_spatial_average_bias(var_name='FLOW_FLUX', true_file_name=flux_file, axes=axes, plot_time_offset=plot_time_offset,
                                            model_time_offset=0., ylim=[-1.8, 1.8])
    axes[1].set_xlabel('Time (day)')
    plt.legend((line1, line2), ('BIAS', 'MAE'), frameon=False, ncol=2, loc="center", bbox_to_anchor=(0.5, -0.5))
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_posterior_bias.png"), dpi=150, bbox_inches='tight')
    plt.close()

    # %%
    # Plot the variance of the updated (flow flux)
    _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    dar.compare_univar_spatial_average_variance(var_name='FLOW_FLUX', true_file_name=flux_file, ax=ax, plot_time_offset=plot_time_offset,
                                            model_time_offset=0., ylim=[-0.01, 0.3])
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_posterior_variance.png"), dpi=150, bbox_inches='tight')
    plt.close()
