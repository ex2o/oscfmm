This folder contains shell scripts for High Performance Computer jobs using SLURM.

This is probably not useful to you unless you are running the simulation on a HPC cluster.

Jobs are submitted with

```
sbatch my-job.sh
```

The array job [univ-inf-arr.sh](univ-inf-arr.sh) will have `#SBATCH --array=1-n` where n  must be equal to `length(x)` where `x` is the variable whose name is specified in `config$slurm_array_col`. 

## Checks before submitting array SLURM job

1. Check that `--array=1-n` has correct `n` in [univ-inf-arr.sh](univ-inf-arr.sh) (and check the corresponding variable referenced by `config$slurm_array_col`)
2. Set `config$slurm_array=T` 
3. Check the `max_elapsed` in config suits the job time in [univ-inf-arr.sh](univ-inf-arr.sh)
4. Check that `config$parallel=T`
5. Check that `use_rcpp` is set to the desired value

