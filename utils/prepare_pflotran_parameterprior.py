"""Prepare PFLOTRAN's ensemble parameters based on configuration file. This is done before the model spinup."""

# Author: Peishi Jiang

import os
import sys
import h5py
import f90nml
from scipy.stats import truncnorm
import numpy as np

# TODO: For now, a single value for each variable is assumed, it should be more generic in the future by considering the 3D case.

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml_file = sys.argv[1]
configs         = f90nml.read(config_nml_file)

pflotran_para_file  = configs["file_cfg"]["pflotran_para_file"]
para_set            = configs["obspara_set_cfg"]["para_set"]
para_min_set        = configs["obspara_set_cfg"]["para_min_set"]
para_max_set        = configs["obspara_set_cfg"]["para_max_set"]
para_mean_set       = configs["obspara_set_cfg"]["para_mean_set"]
para_std_set        = configs["obspara_set_cfg"]["para_std_set"]
para_dist_set       = configs["obspara_set_cfg"]["para_dist_set"]
nens                = configs["da_cfg"]["nens"]

if not isinstance(para_set, list):
    para_set      = [para_set]
    para_min_set  = [para_min_set]
    para_max_set  = [para_max_set]
    para_mean_set = [para_mean_set]
    para_std_set  = [para_std_set]
    para_dist_set = [para_dist_set]


###############################
# Generate ensemble parameters
###############################
if os.path.isfile(pflotran_para_file):
    os.remove(pflotran_para_file)

h5file = h5py.File(pflotran_para_file, 'w')

for i in range(len(para_set)):
    varn = para_set[i]
    dist = para_dist_set[i]
    mean, std  = para_mean_set[i], para_std_set[i]
    maxv, minv = para_max_set[i], para_min_set[i]

    # Generate the ensemble
    if dist.lower() == 'normal':
        values = np.random.normal(mean, std, nens)
    elif dist.lower() == 'lognormal':
        logmean = np.exp(mean + std**2 / 2.)
        logstd  = np.exp(2*mean + std**2) * (np.exp(std**2) - 1)
        values = np.random.lognormal(logmean, logstd)
    elif dist.lower() == 'truncated_normal':
        values = truncnorm.rv(minv, maxv, loc=mean, scale=std, size=nens)
    elif dist.lower() == 'uniform':
        values = np.random.uniform(minv, maxv, nens)
    else:
        raise Exception("unknown distribution %s" % dist)

    h5dset = h5file.create_dataset(varn, data=values)

h5file.close()

