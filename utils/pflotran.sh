#!/bin/bash -l

config_file=$1

exeprg_str=`grep -A 50 exe_cfg $config_file | grep pflotran_exe`
mpirun_str=`grep -A 50 exe_cfg $config_file | grep mpi_exe`
ncore_str=`grep -A 50 exe_cfg $config_file | grep ncore_pf`
pflotranin_str=`grep -A 50 file_cfg $config_file | grep pflotran_in_file`
indir_str=`grep -A 50 other_dir_cfg $config_file | grep pflotran_in_dir`
outdir_str=`grep -A 50 other_dir_cfg $config_file | grep pflotran_out_dir`
nreaz_str=`grep -A 50 da_cfg $config_file | grep nens`

exeprg=`echo $exeprg_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
pflotranin=`echo $pflotranin_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
indir=`echo $indir_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
outdir=`echo $outdir_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
nreaz=`echo $nreaz_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
mpirun=`echo $mpirun_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`
ncore=`echo $ncore_str | sed -e 's#.*=\(\)#\1#' -e "s/[']//g" -e 's/[[:space:]]*//g'`

#exeprg=$1
#pflotranin=$2
#indir=$3
#outdir=$4
#nreaz=$5
#mpirun=$6
#ncore=$7

#echo $pflotranin
#echo $nreaz

if [ $ncore -eq 1 ]
then
    $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1 -screen_output off
else
    $mpirun -n $ncore $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1 -screen_output off
fi

# Move the PFLOTRAN outputs files to output folder
cd $indir
mv pflotran*.h5 $outdir
mv pflotran*.chk $outdir
mv pflotran*.out $outdir

#$mpirun -n $ncore $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1 -screen_output off
#$exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups $ncore -screen_output off ;
wait
echo "Finished running PFLOTRAN..."
