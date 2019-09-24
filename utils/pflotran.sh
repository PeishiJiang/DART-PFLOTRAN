#!/bin/bash -l
exeprg=$1
pflotranin=$2
indir=$3
outdir=$4
nreaz=$5
mpirun=$6
ncore=$7

#name="1dthermal"
#dir="pflotran_results/"

if [ $ncore -eq 1 ]
then
    $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1
else
    $mpirun -n $ncore $exeprg -pflotranin $pflotranin -stochastic -num_realizations $nreaz -num_groups 1
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
