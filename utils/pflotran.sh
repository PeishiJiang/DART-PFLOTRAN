#!/bin/bash -l
exeprg=$1
pflotranin=$2
dir=$3
nreaz=$4
mpirun=$5
ncore=$6

#name="1dthermal"
#dir="pflotran_results/"

#echo $exeprg
echo $pflotranin
#echo $dir
#echo $nreaz
#echo $mpirun
#echo $ncore
cd $dir
#$mpirun -n $ncore $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1
$mpirun -n $ncore $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1 -screen_output off
#$exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups $ncore -screen_output off ;
wait
echo "Finished running PFLOTRAN..."
