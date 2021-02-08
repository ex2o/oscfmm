This folder contains shell scripts for High Performance Computer jobs using SLURM.

This is probably not useful to you unless you are running the simulation on a HPC cluster.

Jobs are submitted with

```
sbatch my-job.sh
```

The array job [univ-inf-arr.sh](univ-inf-arr.sh) will have `#SBATCH --array=1-n` where n  must be equal to `length(NN)` from whichever config list is in use.

