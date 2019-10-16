#!/bin/tcsh
# This is a C-shell script for conducting data assimilation on PFLOTRAN by using DART framework
#
# Author: Peishi Jiang

# TODO Enabled the EnKS-MDA capability
# TODO If there is no observation in one time window, continue running the model and move to the next time step without the data assimilation

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
alias PRINT 'echo "\!:1-$"'
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
set MPI_RUN = `GET_ITEM mpi_exe_da $CONFIG_NML`   # the location of the executable MPI
#set MPI_RUN = mpirun   # the location of the executable MPI

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
set UPDATE_PFLOTRAN_INPUT = `GET_ITEM update_pflotran_input_file $CONFIG_NML` # python script for updating PFLOTRAN input files


##########################################
# Go to the application working directory
##########################################
cd $APP_WORK_DIR


##########################################
# Conduct the data assimilation at time t0
##########################################
echo ""
echo ""
echo "------------------------------------------------------------"
echo "Start the first assimilation at the current model time $MODEL_TIME [day] ..."
echo ""
if ($NCORE_DA == 1) then
  $FILTER_EXE  || exit 1
else
#  $MPI_RUN -n $NCORE_DA -nolocal $FILTER_EXE || exit 1
  $MPI_RUN -n $NCORE_DA $FILTER_EXE || exit 1
  wait
endif

##########################################
# Update the model time and observation start/end time
##########################################
echo ""
echo ""
echo "------------------------------------------------------------"
echo "Move the time forward ..."
echo ""
python $UPDATE_CONFIGNML_TIME $CONFIG_NML $INPUT_NML  || exit 5
set MODEL_TIME = `GET_ITEM current_model_time $CONFIG_NML`  || exit 6
set EXCEEDS_OBS_TIME = `GET_ITEM exceeds_obs_time $CONFIG_NML`  || exit 6


##########################################
# Data assimilation workflow starts here!
##########################################
# Continue the assimilation if MODEL_TIME is smaller than LAST_OBS_TIME
#@ EXCEEDS_OBS_TIME = `MATH $MODEL_TIME >= $LAST_OBS_TIME`
#@ EXCEEDS_OBS_TIME = `echo "$MODEL_TIME >= $LAST_OBS_TIME" | bc -l`
while ($EXCEEDS_OBS_TIME == ".false.")

#  echo "---"
#  echo $MODEL_TIME
#  echo $LAST_OBS_TIME
#  echo $EXCEEDS_OBS_TIME
#  echo "---"

  ##########################################
  # Step 1 -- Update PFLOTRAN input files
  # based on DART posterior output
  # It includes PFLOTRAN.in and parameter_prior.h5 files
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Update PFLOTRAN input files ..."
  echo ""
  python $UPDATE_PFLOTRAN_INPUT  $CONFIG_NML || exit 2

  ##########################################
  # Step 2 -- Conduct PFLOTRAN forward simulation
  # by using PFLOTRAN.sh file
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Conduct the ensemble forward simulation for PFLOTRAN ..."
  echo ""
  $PFLOTRAN_SH $CONFIG_NML  || exit 3

  ##########################################
  # Step 3 -- Convert the PFLOTRAN output to NetCDF format
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Convert the PFLOTRAN HDF output to DART prior NetCDF data ..."
  echo ""
  python $PREP_PRIOR_NC $CONFIG_NML  || exit 4

  ##########################################
  # Step 4 -- Generate the DART data observation
  # in the current assimilation window
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Generate the DART data observations in the current assimilation window ..."
  echo ""
  $CONVERT_NC_EXE  || exit 7

  ##########################################
  # Step 5 -- Data Assimilation
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Start the assimilation at the current model time $MODEL_TIME [day] ..."
  echo ""
  if ($NCORE_DA == 1) then
    $FILTER_EXE  || exit 1
  else
    $MPI_RUN -n $NCORE_DA $FILTER_EXE || exit 1
    wait
  endif

  ##########################################
  # Step 6 -- Update the model time and observation start/end time
  # for the next assimilation window in the input namelist file
  ##########################################
  echo ""
  echo ""
  echo "------------------------------------------------------------"
  echo "Move the time forward ..."
  echo ""
  python $UPDATE_CONFIGNML_TIME $CONFIG_NML $INPUT_NML  || exit 5
  set MODEL_TIME = `GET_ITEM current_model_time $CONFIG_NML`  || exit 6
  set EXCEEDS_OBS_TIME = `GET_ITEM exceeds_obs_time $CONFIG_NML`  || exit 6

end

echo ""
echo ""
echo "------------------------------------------------------------"
echo "The entire data assimilation completes!"
echo ""
