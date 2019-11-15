"""Change the file paths in config namelist files given new locations of DART-PFLOTRAN and APPLICATION folders.
   This is very useful when we analyze or visualize the results from other machines with different locations."""

# Author: Peishi Jiang

import os
import sys
import f90nml
from os.path import dirname as up

###############################
# Get the new locations of DART-PFLOTRAN and applications folder
###############################
dart_pflotran_dir = ''
application_dir   = ''

if dart_pflotran_dir == '' and application_dir == '':
    dart_pflotran_dir = sys.argv[2]
    application_dir   = sys.argv[3]

dart_dir = up(up(up(dart_pflotran_dir)))


###############################
# Parameters
###############################
# Parse the configuration in Fortran namelist
config_nml = sys.argv[1]
configs = f90nml.read(config_nml)


###############################
# Get the old locations of DART-PFLOTRAN and applications folder
###############################
main_dir_cfg          = configs["main_dir_cfg"]
dart_pflotran_dir_old = main_dir_cfg["dart_pf_dir"]
dart_dir_old          = main_dir_cfg["dart_dir"]
application_dir_old   = main_dir_cfg["app_dir"]


###############################
# Change the locations in the other_dir_cfg and file_cfg
###############################
other_dir_cfg = configs["other_dir_cfg"]
file_cfg      = configs["file_cfg"]

for index in other_dir_cfg:
    loc = other_dir_cfg[index]
    if not isinstance(loc, str):
        continue
    if application_dir_old in loc:
        configs["other_dir_cfg"][index] = loc.replace(application_dir_old, application_dir)
    elif dart_pflotran_dir_old in loc:
        configs["other_dir_cfg"][index] = loc.replace(dart_pflotran_dir_old, dart_pflotran_dir)

for index in file_cfg:
    loc = file_cfg[index]
    if not isinstance(loc, str):
        continue
    if application_dir_old in loc:
        configs["file_cfg"][index] = loc.replace(application_dir_old, application_dir)
    elif dart_pflotran_dir_old in loc:
        configs["file_cfg"][index] = loc.replace(dart_pflotran_dir_old, dart_pflotran_dir)


###############################
# Save the new config namelist file
###############################
configs["main_dir_cfg"]["dart_pf_dir"] = dart_pflotran_dir
configs["main_dir_cfg"]["dart_dir"]    = dart_dir
configs["main_dir_cfg"]["app_dir"]     = application_dir

configs.write(config_nml, force=True)