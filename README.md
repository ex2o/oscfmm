# oscfmm
Order Selection with Confidence for Finite Mixture Models

### Details

The script [Analysis.R](scripts/Analysis.R) will call the script [Helpers.R](scripts/Helpers.R) from the working directory.

Parallelisation can be turned off by setting `parallel = F` in [Analysis.R](scripts/Analysis.R). If left on, the number of cores will be detected automatically, leaving one core for the OS.

Parallelisation is handled in the function `repeat_sim`. It works with any simulation function, producing `reps` independent repetitions of the simulation function `sim` with parameters given by the `params` list. For example, choosing `sim <- sim_NN` and setting `NN_range = c(100,200,300)` in the `params` list will simulate the order selection procedure for sample sizes 100, 200 and 300 for `reps` number of times each.

Currently, a new `sim` function needs to be written if you want to vary any parameter other than `NN`. This will likely change in an update.