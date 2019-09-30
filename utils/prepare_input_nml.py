"""Generate the input.nml file by reading (1) the model information (variables, assimilation time step); (2) the grid information; (3) filter configuration."""

# Author: Peishi Jiang

import os
import sys
import pickle
import f90nml

###############################
# Parameters
###############################
input_nml_file          = sys.argv[1]
input_nml_template_file = sys.argv[2]
configure_dict_file     = sys.argv[3]

if not os.path.isfile(configure_dict_file):
    raise Exception("The file %s does not exists!" % configure_dict_file)

# dir_path = os.path.dirname(os.path.realpath(input_nml_dict_file))
# input_nml_template_file = os.path.join(dir_path, "input.nml.template")
if not os.path.isfile(input_nml_template_file):
    raise Exception("The file %s does not exists!" % input_nml_template_file)


###############################
# Read the namelist information from
# configure_pickle
###############################
with open(configure_dict_file, "rb") as f:
    configure = pickle.load(f)

# Get the input namelists
input_nml_dict = configure["inputnml"]
added_namelist = input_nml_dict.keys()

# Get the location of DART and application paths
directories = configure["directories"]
APP_DIR     = directories["APP_DIR"]
DART_DIR    = directories["DART_DIR"]


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
        if not isinstance(value,str):
            continue
        if "[app_dir]" in value:
            nml_ele[item] = value.replace("[app_dir]", APP_DIR)
        elif "[dart_dir]" in value:
            nml_ele[item] = value.replace("[dart_dir]", DART_DIR)
    nml[name] = nml_ele

# Now, let's write
nml.write(input_nml_file)

print("Finished generating the input namelist file...")
