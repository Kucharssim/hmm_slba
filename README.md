# Hidden Markov model of evidence accumulation

This repository is associated with the article

> Kucharský Š., Tran, N.-H., Veldkamp, K., Raijmakers, M.E.J., & Visser, I. (submitted) Hidden Markov Models of Ecidence Accumulation in Speeded Decision Tasks. Preprint at *PsyArXiv*.

### Dependencies (setting up correct version of Stan)

The models were written in Stan and require `CmdStan` version above `2.24.0`. The code was run using `cmdstanr` package in R. To set up the computational environment correctly in the analysis scripts, the correct version of `CmdStan` is loaded using the following command:

```
set_cmdstan_path(readRDS("path_to_cmdstan.Rds"))
```

The object `path_to_cmdstan.Rds` actually contains just a string that specifies the folder where the appropriate version of CmdStan is installed:

 ```
 > readRDS("path_to_cmdstan.Rds")[1] "~/.cmdstan/cmdstan-2.24.0-rc1/"
 ```
 
 For it to work on a computer with different location to `CmdStan`, run the following command in R:
 
 ```
 path_to_cmdstan <- "my/path/to/cmdstan/installation"
 saveRDS(path_to_cmdstan, "path_to_cmdstan.Rds")
 ```

### Structure of this repository

1. [data/](data/) folder contains cleaned data from Dutilh, et al. (2010) that are reanalyzed in the article.
2. [stan/](stan/) folder contains stan models and scripts:
	* [helpers/](stan/helpers/) contain some functions that are used in the stan models
	* [hmm/](stan/hmm/) contain manually written backward, forward, and forward-backward algorithm that was used to check the use of the new `hmm_marginal_lpdf()` function implemented in CmdStan version 2.24.0.
	* [later/](stan/later/) contains functions that implement the simplified LBA model in Stan
	* Additionally, the folder contains models like `hmm_normal.stan` that were used to check the correctness of use of the new HMM functions in Stan language (see also script [hmm_normal.R](scripts/hmm_normal.R)).
3. [scripts/](scripts) folder contains the R scripts used to produce the output presented in the article.
4. [saves/](saves/) folder contains some R objects saved during the project for convenience.
5. [R/](R/) folder contains some R functions for convenience.
6. [figures/](figures) folder contains all figures generated in this project.


