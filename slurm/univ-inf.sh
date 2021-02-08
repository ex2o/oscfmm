#!/bin/bash
#SBATCH --job-name=univinf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G    # memory (MB)
#SBATCH --time=0-00:05      # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out  # STDOUT
#SBATCH -e slurm.%N.%j.err  # STDERR

module load R

Rscript /home/uqdfrye1/univ-inf/Analysis.R