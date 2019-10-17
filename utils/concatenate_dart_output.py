"""This python file is used for concatenate dart output"""

import os
import re
import sys
import f90nml
import subprocess
from math import ceil, log10
import numpy as np

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs = f90nml.read(config_nml)

dart_prior_file         = configs["file_cfg"]["dart_prior_nc_file"]
dart_posterior_file     = configs["file_cfg"]["dart_posterior_nc_file"]
dart_prior_all_file     = configs["file_cfg"]["dart_prior_nc_all_file"]
dart_posterior_all_file = configs["file_cfg"]["dart_posterior_nc_all_file"]
ntimestep               = int(configs["da_cfg"]["ntimestep"])
model_time_list         = configs["time_cfg"]["model_time_list"]
nens                    = configs["da_cfg"]["nens"]

# ndigit = np.ceil(np.log10(ntimestep), dtype=int)
ndigit = int(ceil(log10(ntimestep)))
# Use the real number of time steps here
ntimestep = len(model_time_list)


###############################
# Concatenate all the time steps per ensemble
###############################
for i in range(nens):
    ens = i + 1

    # Get the first prior and posterior files
    dart_prior_first     = re.sub(r"\[ENS\]", str(ens), dart_prior_file)
    dart_prior_first     = re.sub(r"\[TIME\]", str(1).zfill(ndigit), dart_prior_first)
    dart_posterior_first = re.sub(r"\[ENS\]", str(ens), dart_posterior_file)
    dart_posterior_first = re.sub(r"\[TIME\]", str(1).zfill(ndigit), dart_posterior_first)

    dart_prior_all     = re.sub(r"\[ENS\]", str(ens), dart_prior_all_file)
    dart_posterior_all = re.sub(r"\[ENS\]", str(ens), dart_posterior_all_file)

    # Concatanate the prior data
    subprocess.run("ncrcat -n {},{},1 {} {}".format(ntimestep, ndigit,
                                                    dart_prior_first,
                                                    dart_prior_all),
                   shell=True,
                   check=True)

    # Concatanate the posterior data
    subprocess.run("ncrcat -n {},{},1 {} {}".format(ntimestep, ndigit,
                                                    dart_posterior_first,
                                                    dart_posterior_all),
                   shell=True,
                   check=True)
