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

set nonomatch
echo "---------------------------------------------------------------"
echo "Removing *.o *.mod files"
\rm -f preprocess *.o *.mod Makefile
\rm -f ../../../obs_def/obs_def_mod.f90
\rm -f ../../../obs_kind/obs_kind_mod.f90

set MODEL = "NetCDF converters"

@ n = 1

echo
echo
echo "---------------------------------------------------------------"
echo "${MODEL} build number ${n} is preprocess"

csh  mkmf_preprocess || exit 1
make || exit 2

./preprocess || exit 99

echo "---------------------------------------------------------------"
echo "Generating the DART observation file converter"
csh  mkmf_convert_nc  || exit 3
make || exit 4

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

echo "Success: All ${MODEL} programs compiled."

exit 0
