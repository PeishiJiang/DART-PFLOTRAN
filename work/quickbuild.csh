#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: quickbuild.csh 13137 2019-04-30 15:48:00Z nancy@ucar.edu $

#----------------------------------------------------------------------
# compile all programs in the current directory with a mkmf_xxx file.
#
# usage: [ -mpi | -nompi ]
#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so this MUST be run first.
#----------------------------------------------------------------------

set app_work_dir  = $1   # The application work folder
#set mpisetting    = $2   # With or without mpi (i.e., -mpi or -nompi)
set dart_work_dir = $PWD # The DART-PFLOTRAN work folder

echo "---------------------------------------------------------------"
echo "Removing *.o *.mod files"
\rm -f preprocess *.o *.mod Makefile .cppdefs
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "PFLOTRAN"

@ n = 1

echo
echo
echo "---------------------------------------------------------------"
echo "${MODEL} build number ${n} is preprocess"

csh  mkmf_preprocess
make || exit $n

mv preprocess ${app_work_dir} || exit $n
cd ${app_work_dir}
./preprocess || exit $n

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------

cd ${dart_work_dir}
foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   switch ( $TARGET )
   case mkmf_preprocess:
      breaksw
   default:
      @ n = $n + 1
      echo
      echo "---------------------------------------------------"
      echo "${MODEL} build number ${n} is ${PROG}"
      \rm -f ${PROG}
      #csh $TARGET || exit $n
      csh $TARGET
      make        || exit $n
      mv ${PROG} ${app_work_dir} || exit 6
      breaksw
   endsw
end

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

if ( $#argv == 2 && "$2" == "-mpi" ) then
  echo "Success: All single task DART programs compiled."
  echo "Script now compiling MPI parallel versions of the DART programs."
else if ( $#argv == 2 && "$2" == "-nompi" ) then
  echo "Success: All single task DART programs compiled."
  echo "Script is exiting without building the MPI version of the DART programs."
  exit 0
else
  echo ""
  echo "Success: All DART programs compiled."
  echo "Script is exiting before building the MPI version of the DART programs."
  echo "Run the quickbuild.csh script with a -mpi argument or"
  echo "edit the quickbuild.csh script and remove the exit line"
  echo "to compile with MPI to run in parallel on multiple cpus."
  echo ""
  exit 0
endif

#----------------------------------------------------------------------
# to enable an MPI parallel version of filter for this model,
# call this script with the -mpi argument, or if you are going to build
# with MPI all the time, remove or comment out the entire section above.
#----------------------------------------------------------------------

cd ${app_work_dir}
\rm -f filter
cd ${dart_work_dir}

@ n = $n + 1
echo
echo "---------------------------------------------------"
echo "build number $n is mkmf_filter"
csh   mkmf_filter -mpi
make

mv filter ${app_work_dir} || exit 6

if ($status != 0) then
   echo
   echo "If this died in mpi_utilities_mod, see code comment"
   echo "in mpi_utilities_mod.f90 starting with 'BUILD TIP' "
   echo
   exit $n
endif

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs

#echo
#echo 'time to run filter here:'
#echo ' for lsf run "bsub < runme_filter"'
#echo ' for pbs run "qsub runme_filter"'
#echo ' for lam-mpi run "lamboot" once, then "runme_filter"'
#echo ' for mpich run "mpd" once, then "runme_filter"'

exit 0

# <next few lines under version control, do not edit>
# $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/models/template/work/quickbuild.csh $
# $Revision: 13137 $
# $Date: 2019-04-30 08:48:00 -0700 (Tue, 30 Apr 2019) $

