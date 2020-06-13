# adr: Advantage Doubly Robust (ADR) Estimator for Learning When-to-Treat Policies

This repository implements ADR for learning when-to-treat policies, as proposed by Nie, Brunskill and Wager (2019). 

### Authors
This package is written and maintained by Xinkun Nie (xinkun@stanford.edu).

### Reproducibility:
The code has the option to run on a multi-core cluster or on an individual machine. You can set the flag `cluster=T` if the former in `main.R`. 
There are also a few additional flags in `main.R` you can set for running with different options, with comments/documentations in the script.

Here we include instructions on how to reproduce results in the paper. 

First, make new directories for experiments.
```bash
mkdir results
mkdir plots
mkdir logs
```
The following two commands start running two bashscripts that write to disk Monte-Carlo evaluated oracle values for all policies specified in a policy class, for the two simulation setups respectively. 
`start_oracle_1.sh` is a bash script runs the oracle simulations for the multiple-treatment setup specified in the paper, and `start_oracle_2.sh` is a bash script runs the oracle simulations for the binary-treatment setup specified in the paper  

```bash
./start_oracle_1.sh
./start_oracle_2.sh
```
Running the above two commands will write to the directory `results_oracle`. Note that the results of the runs are included in the folder `results_oracle` so you can skip these two commands if you choose to.

Next run the following two commands to actually start the experiments, which read oracle values from the `results_oracle` directory.
`results_oracle`
```bash
./start_1.sh
./start_2.sh
```

Finally, to produce tables and make plots, run the following files in the following order `combine_results.R`, `parse_results.R`, `make_plot.R`. 
You might need to make small modifications to adjust for which setup's results is being compiled and plotted, etc.

To produce plots that compare Q-Opt decision boundaries with ADR, turn `log_trajs=T` in `main.R`, run a few repeats, and then run `plot_policy.R`.

### References
Xinkun Nie, Emma Brunskill and Stefan Wager.
<b>Learning When-to-Treat Policies.</b>
2019.
[<a href="https://arxiv.org/abs/1905.09751">arxiv</a>]
