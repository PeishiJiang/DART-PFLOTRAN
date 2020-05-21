#!/bin/bash
#SBATCH -A m1800
#SBATCH -q debug
#SBATCH -N 32
#SBATCH -t 00:30:00
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH -J DTPF_300A
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=peishi.jiang@pnnl.gov

python workflow_300A_richards_2months.py 
wait