rm(list = ls())
options(error = recover)
start_time = Sys.time()
library(data.table);
library(lubridate);
library(dplyr);
library(splines)
source("utils.R")
source("adr.R")
source("ipw.R")
source("qlearning.R")
source("aipw.R")

cluster = F
oracle_from_file = T # using pre-computed oracle/monte carlo rollout values for each policy in the policy class. Requires running run_oracle.R first
results_dir = "results"
plots_dir = "plots"
results_oracle_dir = "results_oracle"
log_trajs = F # log trajectories into RData files to produce state evolution plots later that compare Fitted-Q decision boundaries vs ADR
if (cluster) { # if running on a cluster
  args=(commandArgs(T))
  setup= as.character(args[1])
} else {
  setup ="simu2" # simu1 corresponds to multiple-treatment setup (the second setup in the paper), and simu2 corresponds to the binary treatment setup (the first setup in the paper)
}

mu_learner="grf"
e_learner="grf"
d_learner="grf"
c_learner = "grf"
q_learner = "grf"
run_adr = T
run_ipw = T
run_aipw = F
run_q_learning = F # run Q-eval for each policy in the policy class
run_q_opt = T # run Q-Opt to find the best policy (fitted-Q iteration)
run_q_aipw_oracle = T # run AIPW and Q-eval for evaluating the oracle-found best policy

if (setup == "simu1") {
  source("simu1.R")
  e_rct=T
  outcome_interest = "lifetime"
  Q_max = T # Q_max is whether we take the max in fitted-Q (vs. take the min). e.g. are we maximizing or minimizing the outcome
  treat_once = T
  if (cluster) {
    n = as.numeric(args[2])
  } else {
    n = 10000
  }
  n_oracle_eval = 20000
  H=10

  setup_params = list()
  setup_params$z_prob = c(0.3, 0.3, 0.4)

  if (cluster) {
  setup_params$obs_noise = as.numeric(args[3])
  REPIDX = as.numeric(args[4])
  } else {
    setup_params$obs_noise = 0
   REPIDX = 1
  }

  params = get_policy_params()

  nfold=5
  foldid = get_foldid(n, nfold)

  fnm = paste0(paste(paste0(results_dir, "/", setup), n, setup_params$obs_noise, REPIDX, sep="-"), ".csv")
  print(fnm)
  if (oracle_from_file) {
    oracle_fnm = paste0(paste(paste0(results_oracle_dir, "/", setup), setup_params$obs_noise, sep="-"), ".csv")
  }
  data_res = get_data(n=n, H=H, setup_params, outcome_interest = outcome_interest, e_rct=e_rct, treat_once=treat_once)

} else if (setup == "simu2")  {

  source("simu2.R")
  e_rct=F
  treat_once = F
  outcome_interest = "other"
  Q_max = T
  if (cluster) {
    n = as.numeric(args[2])
  } else {
    n = 1000
  }
  n_oracle_eval = 20000
  H = 10
  setup_params = list()
  if (cluster) {
    setup_params$obs_noise = as.numeric(args[3])
    setup_params$beta = as.numeric(args[4])
    setup_params$sigma = as.numeric(args[5])
    REPIDX = as.numeric(args[6])
  } else {
    setup_params$obs_noise = 0
    setup_params$beta = 0.5
    setup_params$sigma = 1
    REPIDX = 1
  }

  if (log_trajs) {
    cmp_fnm = paste(paste0(plots_dir, "/plot"), n, setup_params$obs_noise, setup_params$beta, setup_params$sigma, "cmp", REPIDX, sep="-")
  }

  params = get_policy_params()[1:3] # TODO change

  nfold=5
  foldid = get_foldid(n, nfold)
  fnm = paste0(paste(paste0(results_dir, "/", setup), n, setup_params$obs_noise, setup_params$beta, setup_params$sigma, REPIDX, sep="-"), ".csv")
  print(fnm)
  oracle_fnm = paste0(paste(paste0(results_oracle_dir, "/", setup), setup_params$obs_noise, setup_params$beta, setup_params$sigma, sep="-"), ".csv")
  data_res = get_data(n=n, H=H, setup_params, outcome_interest = outcome_interest, e_rct=e_rct, treat_once=treat_once)

}

data_static = data_res$data_static
data_dynamic = data_res$data_dynamic
censor_prob = data_res$censor_prob
n = dim(data_static)[1]


if (run_q_opt) {
 if (setup == "simu1") {
     normal_idx= which(data_static$patient_type == 2)
     data_static_normal = data_static[normal_idx, -c("patient_type"), drop=F]

     q_opt_res = qlearning(data_static_normal,
                           data_dynamic[data_dynamic$id %in% normal_idx,,drop=F],
                           H,
                           foldid[normal_idx],
                           treat_policy=NULL,
                           learner=q_learner,
                           q_opt = T,
                           Q_max = Q_max,
                           outcome_interest = outcome_interest)
   } else {
      q_opt_res = qlearning(data_static,
                            data_dynamic,
                            H,
                            foldid,
                            treat_policy=NULL,
                            learner=q_learner,
                            q_opt = T,
                            Q_max = Q_max,
                            outcome_interest = outcome_interest)
  }
  Q_hat_opt = q_opt_res$Q_hat
  Q_regs_opt = q_opt_res$Q_reg_list

  q_opt_eval_res = q_opt_eval(Q_regs_opt, n_oracle_eval, H, max(data_static$treatment_option), setup_params, outcome_interest, setup, treat_once, q_learner, Q_max)
  q_opt_adv = q_opt_eval_res$V_pi - q_opt_eval_res$V_pi0
  if (setup == "simu2" & log_trajs) {
    q_Z = q_opt_eval_res$Z
    q_A_time = q_opt_eval_res$A_time
    eval_Z0 = q_opt_eval_res$Z0
    eval_X0 = q_opt_eval_res$X0
  }
} else {
  q_opt_adv = 0
}
if (e_rct) {
  e_hat = get_e_hat_rct(n, H, data_static)
} else {
  e_hat = get_est(data_static, data_dynamic, foldid, H, regress_type="e", learner=e_learner)
}


traj_prob = get_traj_prob(e_hat, data_static$treat_start_time, data_static$treatment_option, H) # n by H+1
traj_prob_matrix = get_traj_prob_matrix(e_hat, data_static$treat_start_time, data_static$treatment_option, H) # n by H+1

if (run_adr) {
  traj_after_prob_matrix = get_traj_after_prob_matrix(e_hat, data_static$treat_start_time,H) # n by H+1
  if (setup == "simu1") {
     die_idx = which(data_static$patient_type == 0)
     live_idx = which(data_static$patient_type == 1)
     normal_idx= which(data_static$patient_type == 2)
     data_static_normal = data_static[normal_idx, -c("patient_type"), drop=F]
     e_hat_normal = e_hat[normal_idx,,,drop=F]

     mu_now = get_est(data_static_normal, data_dynamic[data_dynamic$id %in% normal_idx,,drop=F], foldid[normal_idx], H, regress_type="mu_now", learner=mu_learner, outcome_interest=outcome_interest)
     mu_next_internal = get_est(data_static_normal, data_dynamic[data_dynamic$id %in% normal_idx,,drop=F], foldid[normal_idx], H, regress_type="mu_next_internal", learner=mu_learner, outcome_interest=outcome_interest, e_hat=e_hat_normal)
     d_hat = get_est(data_static_normal, data_dynamic[data_dynamic$id %in% normal_idx,,drop=F], foldid[normal_idx], H, regress_type="d", learner=d_learner, outcome_interest = outcome_interest)
     mu_next = get_mu_next(data_static_normal, mu_next_internal, d_hat, H, outcome_interest=outcome_interest, mu_combined=F)
     mu_final = get_mu_final(die_idx, live_idx, normal_idx, mu_now, mu_next, mu_next_internal, H)
     mu_now = mu_final$mu_now
     mu_next = mu_final$mu_next
     mu_next_internal = mu_final$mu_next_internal

  } else {
    mu_now = get_est(data_static, data_dynamic, foldid, H, regress_type="mu_now", learner=mu_learner, outcome_interest=outcome_interest)
    mu_next_internal = get_est(data_static, data_dynamic, foldid, H, regress_type="mu_next_internal", learner=mu_learner, outcome_interest=outcome_interest, e_hat=e_hat)
    d_hat = get_est(data_static, data_dynamic, foldid, H, regress_type="d", learner=d_learner, outcome_interest = outcome_interest)
    mu_next = get_mu_next(data_static, mu_next_internal, d_hat, H, outcome_interest=outcome_interest, mu_combined=F)
  }

  adr_ret= adr(A_time = data_static$treat_start_time,
                   A_choice = data_static$treatment_option,
                   lifetime = data_static$total_months,
                   R = data_static$outcome,
                   mu_now = mu_now,
                   mu_next = mu_next,
                   mu_next_internal = mu_next_internal,
                   censor_prob = censor_prob,
                   traj_prob = traj_prob,
                   traj_after_prob_matrix = traj_after_prob_matrix)
}

num_choices = max(data_static$treatment_option)
if (run_aipw | run_q_learning | run_q_aipw_oracle) {
      pi0_policy = list(month = rep(H+1, n),
                        choice = rep(1,n))

     if (setup == "simu1") {
       normal_idx= which(data_static$patient_type == 2)
       die_idx = which(data_static$patient_type == 0)
       live_idx = which(data_static$patient_type == 1)
       data_static_normal = data_static[normal_idx, -c("patient_type"), drop=F]
       pi0_policy_normal = list(month = pi0_policy$month[normal_idx],
                                choice = pi0_policy$choice[normal_idx])

       q_opt_res = qlearning(data_static_normal,
                             data_dynamic[data_dynamic$id %in% normal_idx,,drop=F],
                             H,
                             foldid[normal_idx],
                             treat_policy=pi0_policy_normal,
                             learner=q_learner,
                             q_opt = F,
                             Q_max = Q_max,
                             outcome_interest = outcome_interest)
      Q_hat_pi0 = array(NA, dim=c(n, H, num_choices+1))
      Q_hat_pi0[die_idx,,] = 1
      Q_hat_pi0[live_idx,,] = H+1
      Q_hat_pi0[normal_idx,,] = q_opt_res$Q_hat

      V_hat_pi0 = matrix(0, n, H)
      V_hat_pi0[normal_idx,] = get_V_hat(Q_hat_pi0[normal_idx,,,drop=F], data_static$treat_start_time[normal_idx], data_static$treatment_option[normal_idx], pi0_policy$month[normal_idx], pi0_policy$choice[normal_idx], data_static$total_months[normal_idx], outcome_interest)
      V_hat_pi0[die_idx,] = 1
      V_hat_pi0[live_idx,] = H+1

     } else {
      q_opt_res_pi0 = qlearning(data_static,
                            data_dynamic,
                            H,
                            foldid,
                            treat_policy=pi0_policy,
                            learner=q_learner,
                            q_opt = F,
                            Q_max = Q_max,
                            outcome_interest = outcome_interest)
      Q_hat_pi0 = q_opt_res_pi0$Q_hat
      V_hat_pi0 = get_V_hat(Q_hat_pi0, data_static$treat_start_time, data_static$treatment_option, pi0_policy$month, pi0_policy$choice, data_static$total_months, outcome_interest)
     }
      if (run_aipw | run_q_aipw_oracle) {
        aipw_pi0 = aipw(A_time= data_static$treat_start_time,
                        A_choice= data_static$treatment_option,
                        policy_time = rep(H+1, n),
                        policy_choice = rep(1, n),
                        R = data_static$outcome,
                        censor_prob = censor_prob,
                        traj_prob = traj_prob,
                        traj_prob_matrix = traj_prob_matrix,
                        Q_hat = Q_hat_pi0,
                        V_hat = V_hat_pi0,
                        weighted=F)

        waipw_pi0 = aipw(A_time= data_static$treat_start_time,
                       A_choice= data_static$treatment_option,
                       policy_time = rep(H+1, n),
                       policy_choice = rep(1, n),
                       R = data_static$outcome,
                       censor_prob = censor_prob,
                       traj_prob = traj_prob,
                       traj_prob_matrix = traj_prob_matrix,
                       Q_hat = Q_hat_pi0,
                       V_hat = V_hat_pi0,
                       weighted=T)
      }
      if (run_q_learning | run_q_aipw_oracle) {
        q_value_pi0 = mean(V_hat_pi0[,1])
      }
}
if (oracle_from_file) {
  oracle_matrix = matrix(as.vector(unlist(read.table(oracle_fnm, header=T))), 3,dim(params)[1]+1)
}

adv = sapply(1:dim(params)[1], function(pp_count) {
    pp = params[pp_count,]

    if (setup == "simu2" | setup == "simu1") {
    if (oracle_from_file) {
      adv_oracle = oracle_matrix[3,pp_count]
    } else {
      oracle = get_oracle_eval(n_oracle_eval, H, pp, setup_params, outcome_interest = outcome_interest, treat_once = treat_once)
      adv_oracle = oracle$V_pi - oracle$V_pi0
    }
   }

   treat_policy = get_treat_policy(data_static, data_dynamic, pp, H)


    if (run_adr) {
      adr_res = get_adr_value(adr_ret, treat_policy, data_static$total_months)

      adv_adr = adr_res$adr
      adv_adr_weighted = adr_res$adr_weighted

      adv_adr_reg = adr_res$adr_reg
      adv_adr_reg_weighted = adr_res$adr_reg_weighted

    } else{
      adv_adr = 0
      adv_adr_weighted = 0

      adv_adr_reg = 0
      adv_adr_reg_weighted = 0
    }

    if (run_ipw) {
      ipw_pi = ipw(A_time = data_static$treat_start_time,
                        A_choice = data_static$treatment_option,
                        policy_time = treat_policy$month,
                        policy_choice = treat_policy$choice,
                        R = data_static$outcome,
                        censor_prob = censor_prob,
                        traj_prob = traj_prob,
                        weighted=F)

      ipw_pi0 = ipw(A_time = data_static$treat_start_time,
                        A_choice = data_static$treatment_option,
                        policy_time = rep(H+1, n),
                        policy_choice = rep(1, n),
                        R = data_static$outcome,
                        censor_prob = censor_prob,
                        traj_prob = traj_prob,
                        weighted=F)

      adv_ipw = ipw_pi - ipw_pi0

      wipw_pi = ipw(A_time = data_static$treat_start_time,
                        A_choice = data_static$treatment_option,
                        policy_time = treat_policy$month,
                        policy_choice = treat_policy$choice,
                        R = data_static$outcome,
                        censor_prob = censor_prob,
                        traj_prob = traj_prob,
                        weighted=T)

      wipw_pi0 = ipw(A_time = data_static$treat_start_time,
                        A_choice = data_static$treatment_option,
                        policy_time = rep(H+1, n),
                        policy_choice = rep(1, n),
                        R = data_static$outcome,
                        censor_prob = censor_prob,
                        traj_prob = traj_prob,
                        weighted=T)

      adv_ipw_weighted = wipw_pi - wipw_pi0

    } else{
      adv_ipw = 0
      adv_ipw_weighted = 0
    }
   if (run_aipw | run_q_learning) {
      if (setup == "simu1") {
       treat_policy_normal = list(month = treat_policy$month[normal_idx],
                                  choice = treat_policy$month[normal_idx])
       q_opt_res = qlearning(data_static_normal,
                             data_dynamic[data_dynamic$id %in% normal_idx,,drop=F],
                             H,
                             foldid[normal_idx],
                             treat_policy=treat_policy_normal,
                             learner=q_learner,
                             q_opt = F,
                             Q_max = Q_max,
                             outcome_interest = outcome_interest)
      Q_hat = array(NA, dim=c(n, H, num_choices+1))
      Q_hat[die_idx,,] = 1
      Q_hat[live_idx,,] = H+1
      Q_hat[normal_idx,,] = q_opt_res$Q_hat

      V_hat = matrix(0, n, H)
      V_hat[normal_idx,] = get_V_hat(Q_hat[normal_idx,,,drop=F], data_static$treat_start_time[normal_idx], data_static$treatment_option[normal_idx], treat_policy$month[normal_idx], treat_policy$choice[normal_idx], data_static$total_months[normal_idx], outcome_interest)
      V_hat[die_idx,] = 1
      V_hat[live_idx,] = H+1

     } else {
      q_opt_res = qlearning(data_static,
                            data_dynamic,
                            H,
                            foldid,
                            treat_policy=treat_policy,
                            learner=q_learner,
                            q_opt = F,
                            Q_max = Q_max,
                            outcome_interest = outcome_interest)
      Q_hat = q_opt_res$Q_hat
      V_hat = get_V_hat(Q_hat, data_static$treat_start_time, data_static$treatment_option, treat_policy$month, treat_policy$choice, data_static$total_months, outcome_interest)
     }
   }

      if (run_aipw) {
        aipw_pi = aipw(A_time= data_static$treat_start_time,
                       A_choice= data_static$treatment_option,
                       policy_time = treat_policy$month,
                       policy_choice = treat_policy$choice,
                       R = data_static$outcome,
                       censor_prob = censor_prob,
                       traj_prob = traj_prob,
                       traj_prob_matrix = traj_prob_matrix,
                       Q_hat = Q_hat,
                       V_hat = V_hat,
                       weighted=F)

        waipw_pi = aipw(A_time= data_static$treat_start_time,
                       A_choice= data_static$treatment_option,
                       policy_time = treat_policy$month,
                       policy_choice = treat_policy$choice,
                       R = data_static$outcome,
                       censor_prob = censor_prob,
                       traj_prob = traj_prob,
                       traj_prob_matrix = traj_prob_matrix,
                       Q_hat = Q_hat,
                       V_hat = V_hat,
                       weighted=T)


        adv_aipw = aipw_pi - aipw_pi0
        adv_aipw_weighted = waipw_pi - waipw_pi0
     } else {
        adv_aipw = 0
        adv_aipw_weighted = 0
     }
   if (run_q_learning) {
     q_value_pi = mean(V_hat[,1])
     adv_q = q_value_pi - q_value_pi0
   } else {
     adv_q = 0
   }

    if (setup == "simu2" | setup == "simu1") {
      c(adv_adr, adv_adr_reg, adv_adr_weighted, adv_adr_reg_weighted,
        adv_ipw, adv_ipw_weighted,
        adv_aipw, adv_aipw_weighted,
        adv_q,
        adv_oracle)
    } else {
      c(adv_adr, adv_adr_reg, adv_adr_weighted, adv_adr_reg_weighted,
        adv_ipw, adv_ipw_weighted,
        adv_aipw, adv_aipw_weighted,
        adv_q)
    }
})

if (setup == "simu2" | setup=="simu1") {
  adv = matrix(unlist(adv), 10, dim(params)[1])
  if (oracle_from_file) {
    oracle_pp_idx = oracle_matrix[1, dim(oracle_matrix)[2]]
    oracle_pp = params[oracle_pp_idx,]
    if (!is.null(dim(oracle_pp))) {
      if (dim(oracle_pp)[1] > 1){
        oracle_pp = pp[1,]
      }
    }
    oracle_best_value = oracle_matrix[3, dim(oracle_matrix)[2]]
  } else {
    oracle_pp_idx = which(adv[10,] == max(adv[10,]))
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
     adv[10,which.max(adv[10,])] = oracle_best_value # an independent run for oracle to avoid winner's curse
  }
   oracle_treat_policy = get_treat_policy(data_static, data_dynamic, oracle_pp, H)
   if (log_trajs & setup == "simu2") {
       wadr_pp = params[which(adv[3,] == max(adv[3,])),]
       if (!is.null(dim(wadr_pp))) {
         if (dim(wadr_pp)[1] > 1){
             wadr_pp = wadr_pp[1,]
         }
       }
       wipw_pp = params[which(adv[6,] == max(adv[6,])),]
       if (!is.null(dim(wipw_pp))) {
         if (dim(wipw_pp)[1] > 1){
             wipw_pp = wipw_pp[1,]
         }
       }
       adr_eval = get_oracle_eval(n = n_oracle_eval, H = H, param = wadr_pp, setup_params = setup_params, outcome_interest=outcome_interest, treat_once=treat_once, X0=eval_X0, Z0=eval_Z0)
       ipw_eval = get_oracle_eval(n = n_oracle_eval, H = H, param = wipw_pp, setup_params = setup_params, outcome_interest=outcome_interest, treat_once=treat_once, X0=eval_X0, Z0=eval_Z0)
       t_limit = max(max(q_A_time), max(adr_eval$policy$month), max(ipw_eval$policy$month))
       save(list = ls(all.names = TRUE), file = paste0(cmp_fnm,".RData"), envir = environment())
    }
   if (run_q_aipw_oracle) {
    if (setup == "simu1") {
       oracle_treat_policy_normal = list(month=oracle_treat_policy$month[normal_idx],
                                         choice=oracle_treat_policy$choice[normal_idx])
       q_opt_res_oracle = qlearning(data_static_normal,
                             data_dynamic[data_dynamic$id %in% normal_idx,,drop=F],
                             H,
                             foldid[normal_idx],
                             treat_policy=oracle_treat_policy_normal,
                             learner=q_learner,
                             q_opt = F,
                             Q_max = Q_max,
                             outcome_interest = outcome_interest)
      Q_hat_oracle = array(NA, dim=c(n, H, num_choices+1))
      Q_hat_oracle[die_idx,,] = 1
      Q_hat_oracle[live_idx,,] = H+1
      Q_hat_oracle[normal_idx,,] = q_opt_res_oracle$Q_hat

      V_hat_oracle = matrix(0, n, H)
      V_hat_oracle[normal_idx,] = get_V_hat(Q_hat_oracle[normal_idx,,,drop=F], data_static$treat_start_time[normal_idx], data_static$treatment_option[normal_idx], oracle_treat_policy$month[normal_idx], oracle_treat_policy$choice[normal_idx], data_static$total_months[normal_idx], outcome_interest)
      V_hat_oracle[die_idx,] = 1
      V_hat_oracle[live_idx,] = H+1

     } else {
      q_opt_res_oracle = qlearning(data_static,
                            data_dynamic,
                            H,
                            foldid,
                            treat_policy=oracle_treat_policy,
                            learner=q_learner,
                            q_opt = F,
                            Q_max = Q_max,
                            outcome_interest = outcome_interest)
      Q_hat_oracle = q_opt_res_oracle$Q_hat
      V_hat_oracle = get_V_hat(Q_hat_oracle, data_static$treat_start_time, data_static$treatment_option, oracle_treat_policy$month, oracle_treat_policy$choice, data_static$total_months, outcome_interest)
     }

        aipw_pi_oracle = aipw(A_time= data_static$treat_start_time,
                              A_choice= data_static$treatment_option,
                              policy_time = oracle_treat_policy$month,
                              policy_choice = oracle_treat_policy$choice,
                              R = data_static$outcome,
                              censor_prob = censor_prob,
                              traj_prob = traj_prob,
                              traj_prob_matrix = traj_prob_matrix,
                              Q_hat = Q_hat_oracle,
                              V_hat = V_hat_oracle,
                              weighted=F)

        waipw_pi_oracle = aipw(A_time= data_static$treat_start_time,
                               A_choice= data_static$treatment_option,
                               policy_time = oracle_treat_policy$month,
                               policy_choice = oracle_treat_policy$choice,
                               R = data_static$outcome,
                               censor_prob = censor_prob,
                               traj_prob = traj_prob,
                               traj_prob_matrix = traj_prob_matrix,
                               Q_hat = Q_hat_oracle,
                               V_hat = V_hat_oracle,
                               weighted=T)


       adv_aipw_oracle = aipw_pi_oracle - aipw_pi0
       adv_waipw_oracle = waipw_pi_oracle - waipw_pi0
       q_value_pi_oracle = mean(V_hat_oracle[,1])
       adv_q_oracle = q_value_pi_oracle - q_value_pi0


} else {
  adv_aipw_oracle = 0
  adv_waipw_oracle = 0
  adv_q_oracle = 0
}


   result = c(adr_best_value = adv[10,which.max(adv[1,])],
     adr_reg_best_value = adv[10,which.max(adv[2,])],

     adr_weighted_best_value = adv[10,which.max(adv[3,])],
     adr_reg_weighted_best_value = adv[10,which.max(adv[4,])],

     ipw_best_value = adv[10,which.max(adv[5,])],
     ipw_weighted_best_value = adv[10,which.max(adv[6,])],

     aipw_best_value = adv[10,which.max(adv[7,])],
     aipw_weighted_best_value = adv[10,which.max(adv[8,])],

     q_best_value = adv[10,which.max(adv[9,])],

     oracle=oracle_best_value,
     fitted_q = q_opt_adv,

     adr_mse = mean((adv[1,] - adv[10,])^2),
     adr_reg_mse = mean((adv[2,] - adv[10,])^2),

     adr_weighted_mse = mean((adv[3,] - adv[10,])^2),
     adr_reg_weighted_mse = mean((adv[4,] - adv[10,])^2),

     ipw_mse = mean((adv[5,] - adv[10,])^2),
     ipw_weighted_mse = mean((adv[6,] - adv[10,])^2),

     aipw_mse = mean((adv[7,] - adv[10,])^2),
     aipw_weighted_mse = mean((adv[8,] - adv[10,])^2),

     q_mse = mean((adv[9,] - adv[10,])^2),

     adv_adr_oracle = adv[1, oracle_pp_idx],
     adv_adr_reg_oracle = adv[2, oracle_pp_idx],
     adv_wadr_oracle = adv[3, oracle_pp_idx],
     adv_wadr_reg_oracle = adv[4, oracle_pp_idx],
     adv_ipw_oracle = adv[5, oracle_pp_idx],
     adv_wipw_oracle = adv[6, oracle_pp_idx],
     adv_aipw_oracle = adv_aipw_oracle,
     adv_waipw_oracle= adv_waipw_oracle,
     adv_q_oracle= adv_q_oracle)

   write.table(result, fnm, row.names=FALSE)

} else {
  stop("setup needs to be simu1 or simu2")
}

end_time = Sys.time()
print(end_time - start_time)
