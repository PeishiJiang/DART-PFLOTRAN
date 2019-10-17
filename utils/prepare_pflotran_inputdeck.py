"""Generate the PFLOTRAN input file (i.e., PFLOTRAN.in) before and during the assimilation."""

# Author: Peishi Jiang

import os
from os.path import dirname
import sys
import copy
import h5py
import f90nml
import numpy as np

# Let's define the observation and state vector model first
# Codes are modified from Kewei's 1D thermal model
# Credit: https://github.com/pnnl-sbrsfa/DA-HEF-paper/blob/master/src/util.py
#

########################################################
# Parameter settings
########################################################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

# Read them from user input
pflotran_in_file = configs["file_cfg"]["pflotran_in_file"]
spinup_length    = float(configs["time_cfg"]["spinup_length"])  #unit: day,   spinup time
assim_window     = float(configs["da_cfg"]["assim_window_size"])  #unit: day, assimilation window
spinup_done      = bool(configs["time_cfg"]["is_spinup_done"])

pflotran_checkpoint_file = configs["file_cfg"]["pflotran_checkpoint_file"]
# pflotran_checkpoint_file = os.path.basename(pflotran_checkpoint_file)

para_set = configs["obspara_set_cfg"]["para_set"]
if not isinstance(para_set, list):
    para_set = [para_set]

#TODO should it be read from the observation file?
# from temperature.csv
# therm_loc = [-0.01, -0.05, -0.65]  # unit:m, location of thermistor, negative means below the riverbed
# Configure model domain and PFLOTRAN running environment
hz = 0.64  # unit: m, height of the 1-D column

# Get the application folder
pflotran_in_path = dirname(os.path.abspath(pflotran_in_file))
# print(os.path.abspath(os.pardir))
# print(os.getcwd())
# print(dirname(os.path.abspath(__file__)))

model_run_time = (spinup_length + assim_window) * 86400.

########################################################
# Obtain the time information
########################################################
obs_data = np.loadtxt(os.path.join(pflotran_in_path, 'obs_data.dat'), skiprows=1)

########################################################
# Write the PFLOTRAN input card file
########################################################
# Read the input card template
with open(os.path.join(pflotran_in_path, '1dthermal_template.in'), 'r') as f:
    pflotranin = f.readlines()

# Write the new input card
with open(pflotran_in_file, 'w') as f:
    for i, s in enumerate(pflotranin):
        if 'NXYZ' in s:
            pflotranin[i]     = "  NXYZ 1 1 {}".format(int(hz * 100)) + "\n"
            pflotranin[i + 2] = "    0.d0 0.d0 {}".format(-hz) + "d0" + "\n"
        if 'REGION all' in s and 'COORDINATES' in pflotranin[i + 1]:
            pflotranin[i + 2] = "    0.d0 0.d0 {}".format(-hz) + "d0" + "\n"
        if 'REGION bottom' in s and "FACE" in pflotranin[i + 1]:
            pflotranin[i + 3] = "    0.d0 0.d0 {}".format(-hz) + "d0" + "\n"
            pflotranin[i + 4] = "    1.d0 1.d0 {}".format(-hz) + "d0" + "\n"
        if "FLOW_CONDITION flow_top" in s and "TYPE" in pflotranin[i + 1]:
            pflotranin[i + 2] = "    FLUX NEUMANN" + "\n"
            pflotranin[i + 5] = "\n"
            pflotranin[i + 6] = "  FLUX DBASE_VALUE FLOW_FLUX m/day" + "\n"
        if 'FLOW_CONDITION initial' in s and "TYPE" in pflotranin[i + 1]:
            # pflotranin[i+6] = "  TEMPERATURE " + str(np.mean(obs.value[0, :])) + "d0" + "\n"
            pflotranin[i + 6] = "  TEMPERATURE " + str(np.mean(
                obs_data[0, :])) + "d0" + "\n"
        if 'THERMAL_CONDUCTIVITY_WET' in s:
            if 'THERMAL_CONDUCTIVITY' in para_set:
                pflotranin[i] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE THERMAL_CONDUCTIVITY" + "\n"
        if 'POROSITY' in s:
            if 'POROSITY' in para_set:
                pflotranin[i] = "  POROSITY DBASE_VALUE POROSITY" + "\n"
        # if "FILENAME 1dthermal" in s:
        #     pflotranin[i-1] = "  RESTART"+"\n"
        #     pflotranin[i] = "    FILENAME " + pflotran_checkpoint_file + " \n"
        #     pflotranin[i+1] = "    REALIZATION_DEPENDENT"+"\n"
        #     pflotranin[i+2] = "    RESET_TO_TIME_ZERO /" + "\n"
        if "FINAL_TIME" in s:
            pflotranin[i] = "  FINAL_TIME {} sec".format(model_run_time) + "\n"

    f.writelines(pflotranin)

print("Finished generating the input card for PFLOTRAN...")

# ########################################################
# # Modify and save the configuration namelist
# ########################################################
# print("Finished generating the DBASE for PFLOTRAN...")
# configs["time_cfg"]["is_spinup_done"] = not spinup_done
# configs.write(config_nml, force=True)
#
