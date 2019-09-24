#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: quickbuild.csh 2019-09-10 peishi.jiang@pnnl.gov $
#
# compile all converter programs

#----------------------------------------------------------------------
# 'preprocess' is a program that culls the appropriate sections of the
# observation module for the observations types in 'input.nml'; the
# resulting source file is used by all the remaining programs,
# so this MUST be run first.
#----------------------------------------------------------------------

set app_work_dir = $1    # The application work folder
set dart_work_dir = $PWD # The DART-PFLOTRAN work folder

set nonomatch
echo "---------------------------------------------------------------"
echo "Removing *.o *.mod files"
\rm -f convert_nc preprocess *.o *.mod Makefile
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "PFLOTRAN"

echo
echo
echo "---------------------------------------------------------------"
echo "${MODEL}: preprocessing DART-PFLOTRAN generic variable quantities"

csh  mkmf_preprocess
make || exit 2

mv preprocess ${app_work_dir} exit || 6
cd ${app_work_dir}
./preprocess || exit 99

echo "---------------------------------------------------------------"
echo "${MODEL}: generating the DART observation file converter"
cd ${dart_work_dir}
csh  mkmf_convert_nc
make || exit 4

mv convert_nc ${app_work_dir} || exit 5

##----------------------------------------------------------------------
## Build all the single-threaded targets
##----------------------------------------------------------------------


#foreach TARGET ( mkmf_* )

   #set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   #switch ( $TARGET )
   #case mkmf_preprocess:
      #breaksw
   #default:
      #@ n = $n + 1
      #echo
      #echo "---------------------------------------------------"
      #echo "${MODEL} build number ${n} is ${PROG}"
      #\rm -f ${PROG}
      #csh $TARGET || exit $n
      #make        || exit $n
      #breaksw
   #endsw
#end

\rm -f *.o *.mod input.nml*_default Makefile .cppdefs
#\rm -f ../utils/convert_nc.f90

echo "Success: All ${MODEL} programs compiled."

exit 0
