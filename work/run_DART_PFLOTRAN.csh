#!/bin/tcsh
# This is a C-shell script for conducting data assimilation on PFLOTRAN by using DART framework
#
# Author: Peishi Jiang

# TODO Enabled the EnKS-MDA capability

set INPUT_NML  = $1  # the input namelist file required by DART filter
set CONFIG_NML = $2  # the configuration files

##########################################
# Check whether the input namelist exists
##########################################
if ( ! -e $INPUT_NML ) then
   echo "ERROR - input.nml does not exist in local directory."
   echo "ERROR - input.nml needed to determine several settings for this script."
   exit 1
endif

if ( ! -e $CONFIG_NML ) then
   echo "ERROR - config.nml does not exist in local directory."
   exit 1
endif


##########################################
# Define several utility functions
##########################################
# Set MATH alias - takes an arithmetic assignment statemen  as argument, e.g., newvar = var1 + var2
# Separate all items and operators in the expression with blanks
alias MATH 'echo "\!:1-$" | bc -l'
# Define a function for getting item value in either input.nml or config.nml
alias GET_ITEM 'grep \!:1 \!:2 | sed -e "s#.*=\(\)#\1#" -e "s/['\'']//g" -e "s/[[:space:]]*//g"'
#alias GET_ITEM 'set \!:1 = `grep \!:2 \!:3 | sed -e "s#.*=\(\)#\1#" -e "s/['\'']//g" -e "s/[[:space:]]*//g"`'
#alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'


##########################################
# Define/Get variables
##########################################
# Get the application work directory
set APP_WORK_DIR = `GET_ITEM app_work_dir $CONFIG_NML`   # the location of the executable MPI

# Get the MPI information
@ NCORE_DA = `GET_ITEM ncore_da $CONFIG_NML`   # the number of cores used in DA
@ NCORE_PF = `GET_ITEM ncore_pf $CONFIG_NML`   # the number of cores used in PFLOTRAN simulation
set MPI_RUN = `GET_ITEM mpi_exe $CONFIG_NML`   # the location of the executable MPI

# Get the time information
set MODEL_TIME    = `GET_ITEM current_model_time $CONFIG_NML`  # the current model time (day)
set LAST_OBS_TIME = `GET_ITEM last_obs_time_size $CONFIG_NML`  # the last observation time (day)
set ASSIM_WINDOW  = `GET_ITEM assim_window_size $CONFIG_NML`   # the assimilation window (day)

# Get the locations of a bunch of files
set PFLOTRAN_SH           = `GET_ITEM pflotran_sh_file $CONFIG_NML`           # shell script for running PFLOTRAN
set FILTER_EXE            = `GET_ITEM filter_exe $CONFIG_NML`                 # the executable filter file
set CONVERT_NC_EXE        = `GET_ITEM convert_nc_exe $CONFIG_NML`             # the executable for generating DART obs
set PREP_PFLOTRAN_INPUT   = `GET_ITEM prep_pflotran_input_file $CONFIG_NML`   # python script for preparing PFLOTRAN input files
set PREP_PRIOR_NC         = `GET_ITEM prep_prior_nc $CONFIG_NML`  # python script for converting PFLOTRAN HDF 5 output to DART NetCDF prior data
set UPDATE_CONFIGNML_TIME = `GET_ITEM update_confignml_time_file $CONFIG_NML` # python script for updating the time data in config.nml
# TODO
set UPDATE_PFLOTRAN_INPUT = `GET_ITEM update_pflotran_input_file $CONFIG_NML` # python script for updating PFLOTRAN input files


##########################################
# Go to the application working directory
##########################################
cd $APP_WORK_DIR


##########################################
# Data assimilation workflow starts here!
##########################################
# Continue the assimilation if MODEL_TIME is smaller than LAST_OBS_TIME
@ WITHIN_OBS_TIME = `MATH $MODEL_TIME < $LAST_OBS_TIME`
while ($WITHIN_OBS_TIME)
  ##########################################
  # Step 1 -- Data Assimilation
  ##########################################
  if ($NCORE_DA == 1) then
    $FILTER_EXE
  # TODO parallel usage of DART filter
  else
    $MPI_RUN -N $NCORE_DA $FILTER_EXE
  endif

  ##########################################
  # TODO
  # Step 2 -- Update PFLOTRAN input files
  # based on DART posterior output
  # It includes PFLOTRAN.in and parameter_prior.h5 files
  ##########################################
  python $UPDATE_PFLOTRAN_INPUT

  ##########################################
  # Step 3 -- Conduct PFLOTRAN forward simulation
  # by using PFLOTRAN.sh file
  ##########################################
  $PFLOTRAN_SH $CONFIG_NML

  ##########################################
  # Step 4 -- Convert the PFLOTRAN output to NetCDF format
  ##########################################
  python $PREP_PRIOR_NC $CONFIG_NML

  ##########################################
  # Step 5 -- Update the model time and observation start/end time
  # for the next assimilation window in the input namelist file
  ##########################################
  python $UPDATE_CONFIGNML_TIME $CONFIG_NML $INPUT_NML
  set MODEL_TIME = `GET_ITEM current_model_time $CONFIG_NML`

  ##########################################
  # Step 6 -- Generate the data observation
  # in the current assimilation window
  ##########################################
  $CONVERT_NC_EXE

end
