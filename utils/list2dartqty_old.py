"""Add new PFLOTRAN variable(s) to DART and create a corresponding obs_def_mod.f90 file."""

# Author: Peishi Jiang

import re
import sys
import shutil
import numpy as np

########################
# Parameters
########################
obs_type_file = sys.argv[1]
obs_kind_file = sys.argv[2]
obs_qty_type_ind_file = sys.argv[3]
pflotran_set  = sys.argv[4:]

obs_kind_template = '../obs_kind/DEFAULT_obs_kind_mod_template.F90'

# Parse the PFLOTRAN variables and get DART variable quantities
p = re.compile('[A-Z_]+')
pflotran_set = [p.search(v).group() for v in pflotran_set]
dart_set     = ['QTY_PFLOTRAN_'+v for v in pflotran_set]
dart_ind_set = []

########################
# Add new PFLOTRAN variable(s) to DART obs_kind_file
########################
with open(obs_kind_template, 'r') as f:
    data = f.readlines()

n_line = len(data)

p2 = re.compile('[1-9]+')
with open(obs_kind_file, 'w') as f:
    for i in range(n_line):
        line = data[i]

        # Add the variable quantity index number
        if line.startswith('integer, parameter :: max_defined_quantities'):
            # Get the current number of quantities
            num_qty = int(p2.search(line).group())
            # Add PFLOTRAN variable quantities from here
            f.write("! PFLOTRAN -- Added by Peishi\n")
            f.write("integer, parameter, public :: &\n")
            for j in range(len(dart_set)):
                num_qty = num_qty + 1
                if j == len(dart_set)-1:
                    f.write("  "+dart_set[j]+"  = "+str(num_qty)+"\n")
                else:
                    f.write("  "+dart_set[j]+"  = "+str(num_qty)+", &\n")
                dart_ind_set.append(num_qty)
            # Add the line for max_defined_quantities
            f.write("\n")
            f.write("\n")
            f.write("integer, parameter :: max_defined_quantities = "+str(num_qty)+"\n")
        # Add the variable quantity definition
        elif line.startswith("! count here, then output below"):
            f.write("! PFLOTRAN -- Added by Peishi\n")
            for j in range(len(dart_set)):
                dt_ind = dart_ind_set[j]
                dt_var = dart_set[j]
                f.write("obs_kind_names("+str(dt_ind)+")"+" = obs_kind_type("+dt_var+", '"+dt_var+"')\n")
            f.write("\n")
            f.write(line)
        else:
            f.write(line)

print("The added DART variable quantities names and indices are ...")
print(pflotran_set)
print(dart_set)
print(dart_ind_set)


########################
# Create the obs_qty_type_ind_file
# This file is only used by the users for reference
########################
with open(obs_qty_type_ind_file, 'w') as f:
    f.write(','.join(pflotran_set)+'\n')
    f.write(','.join(dart_set)+'\n')
    f.write(','.join([str(ind) for ind in dart_ind_set])+"\n")


########################
# Create the obs_type_file
########################
with open(obs_type_file, 'w') as f:
    f.write('! BEGIN DART PREPROCESS KIND LIST\n')
    for i in range(len(dart_set)):
        f.write("!"+pflotran_set[i]+",  "+dart_set[i]+", COMMON_CODE\n")
    f.write('! END DART PREPROCESS KIND LIST\n')
