#!/bin/bash 
#SBATCH -A m1800 
#SBATCH -q regular 
#SBATCH -N 96 
#SBATCH -t 48:00:00 
#SBATCH -L SCRATCH 
#SBATCH -C haswell 
#SBATCH -J DTPF_300A 
#SBATCH --mail-type=begin,end,fail 
#SBATCH --mail-user=peishi.jiang@pnnl.gov 

python workflow.py 
wait