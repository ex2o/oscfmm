# oscfmm

Scripts and simulation results and for the paper [Nguyen, H., Fryer, D., & Mclachlan, G. (2021). Order selection with confidence for finite mixture models](https://arxiv.org/abs/2103.10640).

### Details

The script [Simulation.R](scripts/Simulation.R) implements the simulations. The script [Descriptives.R](scripts/Descriptives.R) loads and combines the simulation results, transforms them from a nested list to a simpler (`data.frame`) format, saves this as [formatted_results.Rds](results/formatted_results.Rds), then produces some basic descriptives and tables of interest. 
