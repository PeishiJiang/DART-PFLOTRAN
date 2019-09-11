"""Add new PFLOTRAN variable(s) to DART and create a corresponding obs_def_mod.f90 file."""

# Author: Peishi Jiang

import os
import re
import sys
import shutil
import numpy as np

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

########################
# Parameters
########################
obs_type_file = sys.argv[1]
obs_kind_file = sys.argv[2]
pflotran_set  = sys.argv[3:]

# obs_kind_template = '../obs_kind/DEFAULT_obs_kind_mod_template.F90'

# Parse the PFLOTRAN variables and get DART variable quantities
p = re.compile('[A-Z_]+')
pflotran_set = [p.search(v).group() for v in pflotran_set]
dart_set     = ['QTY_PFLOTRAN_'+v for v in pflotran_set]
dart_ind_set = []

########################
# Add new PFLOTRAN variable(s) to DART obs_kind_file
########################
with open(obs_kind_file, 'r') as f:
    data = f.readlines()

# Remove the current one
os.remove(obs_kind_file)

n_line = len(data)
existing_dart_set = []
dart_set_unique   = np.copy(dart_set)

count = False
finished_add = False

p2 = re.compile('[1-9]+')
line = ''
with open(obs_kind_file, 'w') as f:
    for i in range(n_line):
        line = data[i]

        # Start to count the avaiable PFLOTRAN variables
        if not finished_add and line.startswith("! PFLOTRAN -- Added by Peishi"):
            count = True
            f.write(line)
            continue

        # Skip the first line after "! PFLOTRAN -- Added by Peishi"
        if count and line.startswith("integer, parameter, public ::"):
            f.write(line)
            continue
        # Get the current available PFLOTRAN variables
        elif count and not line.startswith("integer, parameter, public ::") and not line == "\n":
            if ", &\n" not in line:
        # TODO change the line of te last current PFLOTRAN variable, i.e., adding ', &\n'
                line = line[:-1] + ", &\n"
            dart_var = line.split()[0]
            existing_dart_set.append(dart_var)
            last_dartvar = dart_var
            num_qty = int(p2.search(line).group())
            last_num_qty = num_qty
            f.write(line)
            continue
        # End of getting the current PFLOTRAN variables in DART and do the following
        elif count and line == "\n":
        # Get the unique variables in our list
            dart_set_inters = intersection(dart_set, existing_dart_set)
            dart_set_unique = list(set(dart_set)-set(dart_set_inters))
        # Add these unique variables to our DEFAULT_obs_kind_mod.f90 file
            for j in range(len(dart_set_unique)):
                num_qty = num_qty + 1
                if j == len(dart_set_unique)-1:
                    f.write("  "+dart_set_unique[j]+"  = "+str(num_qty)+"\n")
                else:
                    f.write("  "+dart_set_unique[j]+"  = "+str(num_qty)+", &\n")
                dart_ind_set.append(num_qty)
            count = False
            f.write("\n")
            f.write("\n")
            f.write("integer, parameter :: max_defined_quantities = "+str(num_qty)+"\n")
            finished_add = True
            count = False
            # print(last_dartvar)
            continue

        # Skip this line since it has already been written
        if finished_add and line.startswith("integer, parameter :: max_defined_quantities"):
            continue

        # Add the variable quantity definition
        if finished_add and line.startswith("obs_kind_names("+str(last_num_qty)+")"):
            # print(line)
            f.write(line)
            for j in range(len(dart_set_unique)):
                dt_ind = dart_ind_set[j]
                dt_var = dart_set_unique[j]
                f.write("obs_kind_names("+str(dt_ind)+")"+" = obs_kind_type("+dt_var+", '"+dt_var+"')\n")
            f.write("\n")
            continue

        # For other cases, just write the line
        f.write(line)


print("The added DART variable quantities names and indices are ...")
print(dart_set_unique)
print(dart_ind_set)


########################
# Create the obs_type_file
########################
with open(obs_type_file, 'w') as f:
    f.write('! BEGIN DART PREPROCESS KIND LIST\n')
    for i in range(len(dart_set)):
        f.write("!"+pflotran_set[i]+",  "+dart_set[i]+", COMMON_CODE\n")
    f.write('! END DART PREPROCESS KIND LIST\n')
