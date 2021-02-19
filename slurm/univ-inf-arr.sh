#!/bin/bash
#SBATCH --array=1-40
#SBATCH --job-name=univinf_arr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1G        # memory (MB)
#SBATCH --time=0-15:30          # time (D-HH:MM)
#SBATCH -o slurm_arr_%A_%a.out  # STDOUT
#SBATCH -e slurm_arr_%A_%a.err  # STDERR

module load R

Rscript /home/uqdfrye1/univ-inf/Simulation.R