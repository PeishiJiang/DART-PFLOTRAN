"""Generate the input.nml file by reading (1) the model information (variables, assimilation time step); (2) the grid information; (3) filter configuration."""

# Author: Peishi Jiang

import os
import sys
import pickle
import f90nml

###############################
# Parameters
###############################
input_nml_file = sys.argv[1]
input_nml_dict_file = sys.argv[2]

if not os.path.isfile(input_nml_dict_file):
    raise Exception("The file %s does not exists!" % input_nml_dict_file)

dir_path = os.path.dirname(os.path.realpath(input_nml_dict_file))
input_nml_template_file = os.path.join(dir_path, "input.nml.template")
if not os.path.isfile(input_nml_template_file):
    raise Exception("The file %s does not exists!" % input_nml_template_file)


###############################
# Read the namelist information from
# input_nml_dict
###############################
with open(input_nml_dict_file, "rb") as f:
    input_nml_dict = pickle.load(f)
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

# Now, let's write
nml.write(input_nml_file)

os.remove(input_nml_dict_file)

print("Finished generating the input namelist file...")
