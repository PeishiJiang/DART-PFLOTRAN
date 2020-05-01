"""Prepare PFLOTRAN's ensemble parameters based on configuration file. This is done before the model spinup."""

# Author: Peishi Jiang

import os
import sys
import h5py
import f90nml
import subprocess
from scipy.stats import truncnorm
import numpy as np

from parse_pflotran_files_utils import pflotran_files

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_para_file = configs["file_cfg"]["pflotran_para_file"]
pflotran_para_backup_file = configs["file_cfg"]["pflotran_para_backup_file"]
use_para_backup    = False if configs["da_cfg"]["use_para_initial_at_nth_window"] < 0 else True
use_default_para_initial = configs["obspara_set_cfg"]["use_default_para_initial"]

para_set           = configs["obspara_set_cfg"]["para_set"]
para_min_set       = configs["obspara_set_cfg"]["para_min_set"]
para_max_set       = configs["obspara_set_cfg"]["para_max_set"]
para_mean_set      = configs["obspara_set_cfg"]["para_mean_set"]
para_std_set       = configs["obspara_set_cfg"]["para_std_set"]
# para_dist_set      = configs["obspara_set_cfg"]["para_dist_set"]
para_sample_method_set = configs["obspara_set_cfg"]["para_sample_method_set"]
para_resampled_set = configs["obspara_set_cfg"]["para_resampled_set"]
nens               = configs["da_cfg"]["nens"]

if not isinstance(para_set, list):
    para_set      = [para_set]
    para_min_set  = [para_min_set]
    para_max_set  = [para_max_set]
    para_mean_set = [para_mean_set]
    para_std_set  = [para_std_set]
    para_dist_set = [para_dist_set]

if not isinstance(para_resampled_set, list):
    para_resampled_set = [para_resampled_set]


###############################
# Generate ensemble parameters
###############################
if use_default_para_initial:
    print("Use the defaul parameter file as the initial")
    if not os.path.isfile(pflotran_para_file):
        raise Exception("No detault parameter initial file exists!")

else:
# TODO: For now, a single value for each variable is assumed, it should be more generic in the future by considering the 3D case.
    raise Exception("Generating the initial parameters needs to be implemented!")
    # if os.path.isfile(pflotran_para_file):
    #     os.remove(pflotran_para_file)

    # h5file = h5py.File(pflotran_para_file, 'w')

    # for i in range(len(para_set)):
    #     varn = para_set[i]
    #     dist = para_dist_set[i]

    #     mean, std  = para_mean_set[i], para_std_set[i]
    #     maxv, minv = para_max_set[i],  para_min_set[i]

    #     # Generate the ensemble
    #     if dist.lower() == 'normal':
    #         values = np.random.normal(mean, std, nens)
    #         # Exclude those values outside of [minv, maxv]
    #         if minv != -99999:
    #             values[values < minv] = minv
    #         if maxv != 99999:
    #             values[values > maxv] = maxv

    #     elif dist.lower() == 'lognormal':
    #         # logmean = np.exp(mean + std**2 / 2.)
    #         # logstd  = np.exp(2 * mean + std**2) * (np.exp(std**2) - 1)
    #         # values  = np.random.lognormal(logmean, logstd)
    #         logvalues  = np.random.normal(mean, std, nens)
    #         # Exclude those values outside of [minv, maxv]
    #         if minv != -99999:
    #             logvalues[logvalues < minv] = minv
    #         if maxv != 99999:
    #             logvalues[logvalues > maxv] = maxv
    #         values = np.power(10, logvalues)

    #     elif dist.lower() == 'truncated_normal':
    #         values = truncnorm.rvs(minv, maxv, loc=mean, scale=std, size=nens)
    #         # Exclude those values outside of [minv, maxv]
    #         if minv != -99999:
    #             values[values < minv] = minv
    #         if maxv != 99999:
    #             values[values > maxv] = maxv

    #     elif dist.lower() == 'uniform':
    #         values = np.random.uniform(minv, maxv, nens)
    #         # Exclude those values outside of [minv, maxv]
    #         if minv != -99999:
    #             values[values < minv] = minv
    #         if maxv != 99999:
    #             values[values > maxv] = maxv

    #     elif dist.lower() == 'test':
    #         values = np.linspace(minv, maxv, nens)
    #         # Exclude those values outside of [minv, maxv]
    #         if minv != -99999:
    #             values[values < minv] = minv
    #         if maxv != 99999:
    #             values[values > maxv] = maxv

    #     else:
    #         raise Exception("unknown distribution %s" % dist)


    #     h5dset = h5file.create_dataset(varn, data=values)

    # h5file.close()


###############################
# Need to backup the intial parameters?
###############################
if use_para_backup:
   subprocess.run("cp {} {}".format(pflotran_para_file, pflotran_para_backup_file), shell=True, check=True) 