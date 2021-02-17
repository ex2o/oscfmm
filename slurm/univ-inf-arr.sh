#!/bin/bash
#SBATCH --array=1-2
#SBATCH --job-name=univinf_arr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G       # memory (MB)
#SBATCH --time=0-00:10          # time (D-HH:MM)
#SBATCH -o slurm_arr_%A_%a.out  # STDOUT
#SBATCH -e slurm_arr_%A_%a.err  # STDERR

module load R

Rscript /home/uqdfrye1/univ-inf/Simulation.R