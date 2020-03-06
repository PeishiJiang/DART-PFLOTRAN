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
app_dir_name_list = ['1dthermal_test_1month_obserror005_abs_3600s_1mda']
# app_dir_name_list = ['1dthermal_test_1month_obserror005_abs_3600s_1mda', '1dthermal_test_1month_obserror005_abs_3600s_2mda',
                     # '1dthermal_test_1month_obserror005_abs_3600s_3mda', '1dthermal_test_1month_obserror005_abs_3600s_4mda']

# Plotting time offset
plot_time_offset = 12 # for two days of 5min resolution data

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
    from utils.plot_da_result import DaResults, plot_compare_multiple_daresults
    dar = DaResults(config_file)
    dar.setup(from_concatenated=0)
    nvar = dar.nvar

    # %% [markdown]
    # <a id='plot_spatial_average'></a>
    # # Visualize spatial average results

    # %%
    # fig = plt.figure(num=1,dpi=150, figsize=(10,5), constrained_layout=True)
    # gs = gridspec.GridSpec(nvar, 2, width_ratios=[1, 1])
    # axes = np.empty([nvar, 2], dtype=object)
    # for i in range(nvar):
    #     axes[i, 0] = plt.subplot(gs[i, 0])
    #     axes[i, 1] = plt.subplot(gs[i, 1], sharey=axes[i, 0], sharex=axes[i, 0])
    # dar.plot_spatial_average(axes)

    # %% [markdown]
    # <a id='plot_along_zaxis'></a>
    # # Visualize temporal evolution along one dimension

    # %%
    dar.plot_oned_obs('TEMPERATURE', figsize=(10,10), dim_str='z', vmin=4, vmax=10, ylim=[4, 10])
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_temperature.png"), dpi=150)

    # %% [markdown]
    # # Compare the true with analyzed results (one variable, spatially averaged)

    # %%
    flux_file = os.path.join(app_dir, "pflotran_input/flux_1d_new.csv")
    # fig = plt.figure(num=1,dpi=150, figsize=(10,2))
    # gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], figure=fig)
    # axes = np.empty(2, dtype=object)
    # axes[0] = plt.subplot(gs[0, 0])
    # axes[1] = plt.subplot(gs[0, 1], sharey=axes[0], sharex=axes[0])
    fig, axes = plt.subplots(1, 2, figsize=(10,2))
    dar.compare_univar_spatial_average(var_name='FLOW_FLUX', true_file_name=flux_file, axes=axes, 
                                       model_time_offset=1800./86400., plot_time_offset=plot_time_offset, ylim=[-4, 4], xlim=[-2,32])
    # dar.compare_univar_spatial_average(var_name='FLOW_FLUX', true_file_name=flux_file, axes=axes, ylim=[-10, 10])
    # axes[1].set_xlim([8,9.5])
    plt.subplots_adjust(wspace=0.2)
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux.png"), dpi=150, bbox_inches='tight')

    # %% [markdown]
    # # Compare the true with analyzed results (one variable, spatially averaged),  Cont'd

    # %%
    _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    diff_mean, posterior_mean, true_mean, true, true_time_set, true_set_dates = dar.compare_univar_spatial_average_diff(var_name='FLOW_FLUX', true_file_name=flux_file, ax=ax,
                                            # model_time_offset=1800./86400., ylimm=[-4, 4], ylimstd=[0.0, 1], unit='m/s')
                                            model_time_offset=1800./86400., ylim=[-1.8, 1.8], xlim=[-2,32], unit='m/s', plot_time_offset=plot_time_offset)
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_diff.png"), dpi=150, bbox_inches='tight')

    _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    dar.compare_univar_spatial_average_bias(var_name='FLOW_FLUX', true_file_name=flux_file, ax=ax, plot_time_offset=plot_time_offset,
                                            model_time_offset=1800./86400., ylim=[-1.8, 1.8], xlim=[-2,32], unit='m/s')
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_posterior_bias.png"), dpi=150, bbox_inches='tight')

    _, ax = plt.subplots(1, 1, figsize=(4,2), dpi=150)
    dar.compare_univar_spatial_average_variance(var_name='FLOW_FLUX', true_file_name=flux_file, ax=ax, plot_time_offset=plot_time_offset,
                                            model_time_offset=1800./86400., ylim=[-0.01, 0.3], xlim=[-2,32], unit='m/s')
    plt.savefig(os.path.join(fig_dir, app_dir_name+"_flowflux_posterior_variance.png"), dpi=150, bbox_inches='tight')
