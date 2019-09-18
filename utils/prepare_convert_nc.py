"""Prepare the convert_nc.f90 file by reading the variable names in the observation .nc file."""

# Author: Peishi Jiang

import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import num2date, date2num, Dataset
from datetime import datetime, timedelta

###############################
# Parameters
###############################
obs_nc       = sys.argv[1]
convert_file = '../utils/convert_nc.f90'
convert_template_file = '../utils/convert_nc_template.f90'

###############################
# Read the obs_nc file
# and get the list of observation variables
###############################
root_nc = Dataset(obs_nc, 'r')
var_set = []
for v in root_nc.variables.keys():
    vtype = root_nc.variables[v].type
    if vtype == 'observation_value':
        var_set.append(v)

# print(var_set)

###############################
# Modify the convert_nc_template.f90
# by including these observation
# variables
###############################
# Read in the template file
with open(convert_template_file, 'r') as f:
    data = f.readlines()

# Remove the current one
os.remove(convert_file)

n_line = len(data)

added = 0
var_abbr_set     = [v.lower()+"_val" for v in var_set]
var_abbr_def_set = [v.lower()+"_val(:,:)" for v in var_set]
var_miss_set     = [v.lower()+"_miss" for v in var_set]
var_qc_set       = ["qc_"+v.lower() for v in var_set]
var_qc_def_set   = ["qc_"+v.lower()+"(:,:)" for v in var_set]
nvar             = len(var_set)
# Now, let's write
with open(convert_file, 'w') as f:
    for i in range(n_line):
        line = data[i]
        # Add 1 -- definition for observation variables
        if line.startswith("!$TODO Add here 1"):
            f.write("use obs_kind_mod, only: "+",".join(var_set)+"\n")
            added = 1
            continue
        elif added == 1 and not line.startswith("!$END"):
            continue
        elif added == 1 and line.startswith("!$END"):
            added = 0
            continue

        # Add 2 -- definition for observation variables values
        if line.startswith("!$TODO Add here 2"):
            f.write("real(r8), allocatable :: "+",".join(var_abbr_def_set)+"\n")
            added = 2
            continue
        elif added == 2 and not line.startswith("!$END"):
            continue
        elif added == 2 and line.startswith("!$END"):
            added = 0
            continue

        # Add 3 -- definition for observation variables missing values
        if line.startswith("!$TODO Add here 3"):
            f.write("real(r8) :: "+",".join(var_miss_set)+"\n")
            added = 3
            continue
        elif added == 3 and not line.startswith("!$END"):
            continue
        elif added == 3 and line.startswith("!$END"):
            added = 0
            continue

        # Add 4 -- definition for quality control value
        if line.startswith("!$TODO Add here 4"):
            f.write("integer, allocatable :: "+",".join(var_qc_def_set)+"\n")
            added = 4
            continue
        elif added == 4 and not line.startswith("!$END"):
            continue
        elif added == 4 and line.startswith("!$END"):
            added = 0
            continue

        # Add 5 -- definition for number of variable
        if line.startswith("!$TODO Add here 5"):
            f.write("integer, parameter :: nvar="+str(nvar)+"\n")
            added = 5
            continue
        elif added == 5 and not line.startswith("!$END"):
            continue
        elif added == 5 and line.startswith("!$END"):
            added = 0
            continue

        # Add 6 -- allocate the memory for observational variables
        if line.startswith("!$TODO Add here 6"):
            for j in range(nvar):
                f.write("allocate("+var_abbr_set[j]+"(ntime,nloc))\n")
                f.write("allocate("+var_qc_set[j]+"(ntime,nloc))\n")
            added = 6
            continue
        elif added == 6 and not line.startswith("!$END"):
            continue
        elif added == 6 and line.startswith("!$END"):
            added = 0
            continue

        # Add 7 -- definition for number of variable
        if line.startswith("!$TODO Add here 7"):
            for j in range(nvar):
                f.write("call getvar_real_2d(ncid, '"+var_set[j]+"',"+var_abbr_set[j]+","+var_miss_set[j]+")\n")
            added = 7
            continue
        elif added == 7 and not line.startswith("!$END"):
            continue
        elif added == 7 and line.startswith("!$END"):
            added = 0
            continue

        # Add 8 -- getting quality control value for each variable
        if line.startswith("!$TODO Add here 8"):
            f.write("! Define or get the quality control value for each observation variable\n")
            f.write("if (use_input_qc) then\n")
            for j in range(nvar):
                f.write("call getvar_int_2d(ncid, '"+var_set[j]+"QCR"+"', "+var_qc_set[j]+")\n")
            f.write("else\n")
            for j in range(nvar):
                f.write(var_qc_set[j]+" = 0\n")
            f.write("endif\n")
            added = 8
            continue
        elif added == 8 and not line.startswith("!$END"):
            continue
        elif added == 8 and line.startswith("!$END"):
            added = 0
            continue

        # Add 9 -- add observation value
        if line.startswith("!$TODO Add here 9"):
            f.write("! Add each observation value here\n")
            f.write("if ( &\n")
            for j in range(nvar-1):
                f.write("  "+var_abbr_set[j]+"(n,k) /= "+var_miss_set[j]+" .and. "+var_qc_set[j]+"(n,k) == 0 .and. &\n")
            f.write("  "+var_abbr_set[-1]+"(n,k) /= "+var_miss_set[-1]+" .and. "+var_qc_set[-1]+"(n,k) == 0) then\n")
            for j in range(nvar):
                f.write("   call create_3d_obs(xloc(k), yloc(k), zloc(k), 0, "+var_abbr_set[j]+"(n,k), "+var_set[j]+", oerr, oday, osec, qc, obs)\n")
                f.write("   call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)\n")
            f.write("endif\n")
            added = 9
            continue
        elif added == 9 and not line.startswith("!$END"):
            continue
        elif added == 9 and line.startswith("!$END"):
            added = 0
            continue

        # Else, just write the line
        f.write(line)
