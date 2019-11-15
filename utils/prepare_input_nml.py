"""Generate the input.nml file by reading (1) the model information (variables, assimilation time step); (2) the grid information; (3) filter configuration."""

# Author: Peishi Jiang

import os
import sys
import pickle
import f90nml

###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs    = f90nml.read(config_nml)

input_nml_file          = configs["file_cfg"]["input_nml_file"]
input_nml_template_file = configs["file_cfg"]["input_nml_template_file"]
input_nml_dict          = configs["inputnml_cfg"]
# configure_dict_file     = sys.argv[3]

# Get the location of DART and application paths
directories  = configs["main_dir_cfg"]
APP_DIR      = directories["app_dir"]
DART_DIR     = directories["dart_dir"]
APP_WORK_DIR = configs["other_dir_cfg"]["app_work_dir"]

if not os.path.isfile(input_nml_template_file):
    raise Exception("The file %s does not exists!" % input_nml_template_file)

###############################
# Read the namelist information
###############################
# Get the input namelists
added_namelist = input_nml_dict.keys()

###############################
# Read the namelist information from
# input_nml_dict
###############################
# Read in the template file
nml = f90nml.read(input_nml_template_file)

# Remove the current one
if os.path.isfile(input_nml_file):
    os.remove(input_nml_file)

# Modify the template file
for name in added_namelist:
    added_nml_ele = input_nml_dict[name]
    if name in nml.keys():
        for key in added_nml_ele:
            nml[name][key] = added_nml_ele[key]
    else:
        nml[name] = added_nml_ele

# Replace the paths of some files
for name in nml.keys():
    nml_ele = nml[name]
    for item, value in nml_ele.items():
        if not isinstance(value, str):
            continue
        if "[app_dir]" in value:
            nml_ele[item] = value.replace("[app_dir]", APP_DIR)
        elif "[dart_dir]" in value:
            nml_ele[item] = value.replace("[dart_dir]", DART_DIR)
        # Get the relative directory w.r.t. the application work directory
        if DART_DIR in nml_ele[item]:
            nml_ele[item] = os.path.relpath(nml_ele[item], APP_WORK_DIR)
    nml[name] = nml_ele

# Now, let's write
nml.write(input_nml_file)

print("Finished generating the input namelist file...")
