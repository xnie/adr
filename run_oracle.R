# pre-computes oracle values (via monte carlo rollouts) for each policy in the policy class
rm(list = ls())
options(error = recover)
start_time = Sys.time()
library(data.table);
library(lubridate);
library(dplyr);
library(splines)
source("utils.R")
cluster = F
results_dir = "results_oracle"
if (cluster) {
  args=(commandArgs(T))
  setup= as.character(args[1])
} else {
  setup ="simu2"
}

if (setup == "simu1") {

  source("simu1.R")
  outcome_interest = "lifetime"
  treat_once = T
  n_oracle_eval = 20000
  H=10

  setup_params = list()
  setup_params$z_prob = c(0.3, 0.3, 0.4)

  if (cluster) {
  setup_params$obs_noise = as.numeric(args[2])
  } else {
    setup_params$obs_noise = 0
    REPIDX = 1
  }

  params = get_policy_params()

  fnm = paste0(paste(paste0(results_dir, "/", setup), setup_params$obs_noise, sep="-"), ".csv")
  print(fnm)

} else if (setup == "simu2")  {

  source("simu2.R")
  treat_once = F
  outcome_interest = "other"
  n_oracle_eval = 20000
  H = 10

  setup_params = list()
  if (cluster) {
    setup_params$obs_noise = as.numeric(args[2])
    setup_params$beta = as.numeric(args[3])
    setup_params$sigma = as.numeric(args[4])
  } else {
    setup_params$obs_noise = 0.5
    setup_params$beta = 5
    setup_params$sigma = 1
    REPIDX = 1
  }

  params = get_policy_params()

  fnm = paste0(paste(paste0(results_dir, "/", setup), setup_params$obs_noise, setup_params$beta, setup_params$sigma, sep="-"), ".csv")
  print(fnm)

}
adv = apply(params, 1, function(pp) {
    print(pp)
    oracle = get_oracle_eval(n_oracle_eval, H, pp, setup_params, outcome_interest = outcome_interest, treat_once = treat_once)
    adv_oracle = oracle$V_pi - oracle$V_pi0
    c(oracle$V_pi,
      oracle$V_pi0,
      adv_oracle)
})

adv = matrix(unlist(adv), 3, dim(params)[1])
oracle_pp_idx = which(adv[3,] == max(adv[3,]))
if (length(oracle_pp_idx) > 1) {
  oracle_pp_idx = oracle_pp_idx[1]
}

oracle_pp = params[oracle_pp_idx,]
if (!is.null(dim(oracle_pp))) {
  if (dim(oracle_pp)[1] > 1){
    oracle_pp = pp[1,]
  }
}

oracle_eval = get_oracle_eval(n_oracle_eval, H, oracle_pp, setup_params, outcome_interest = outcome_interest, treat_once = treat_once)
oracle_best_value = oracle_eval$V_pi - oracle_eval$V_pi0
adv[1,oracle_pp_idx] = oracle_eval$V_pi
adv[2,oracle_pp_idx] = oracle_eval$V_pi0
adv[3,oracle_pp_idx] = oracle_best_value
adv = cbind(adv,c(oracle_pp_idx, oracle_eval$V_pi, oracle_best_value))
write.table(adv,file=fnm, row.names = FALSE)
