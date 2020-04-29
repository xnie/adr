qlearning = function(data_static,
                     data_dynamic,
                     H,
                     foldid,
                     treat_policy=NULL,
                     learner="grf",
                     q_opt = F,
                     Q_max=F,
                     outcome_interest) {

  num_choices = max(data_static$treatment_option)
  n = dim(data_static)[1]
  if (q_opt) {
    foldid = rep(1,n)
  } else {
    if (is.null(treat_policy)) {
      stop("need to supply treat policy if not learning optimal q function")
    }
  }

  Q_hat = array(NA, dim=c(n, H, num_choices+1))
  Q_reg_list = vector("list", H)

  if (num_choices > 1) {
    treat_option = model.matrix(~ factor(c(1:num_choices)) + 0)
    treat_option_0 = rep(0,dim(treat_option)[2])
    treat_option = rbind(treat_option_0, treat_option)
  } else {
    treat_option = matrix(c(0,1),2,1)
  }

  data_treat_option = matrix(0, n, dim(treat_option)[2])
  data_treat_option[which(data_static$treat_start_time <= H),] = treat_option[c(1+data_static$treatment_option[which(data_static$treat_start_time <= H)]),]

  for (t in H:1) {
    if (t==H) {
      Q_next = data_static$outcome
    } else {
      Q_next = rep(0, n)
      Q_next[which(data_static$treat_start_time <= t)] = Q_hat[cbind(which(data_static$treat_start_time <= t), t+1, 1 + data_static$treatment_option[which(data_static$treat_start_time<=t)])]
      if (q_opt) {
        # fitted Q to learn best policy
        if (Q_max) {
          Q_next[which(data_static$treat_start_time > t)] = apply(Q_hat[,t+1,], 1, function(x) max(x))[which(data_static$treat_start_time > t)]
        } else {
          Q_next[which(data_static$treat_start_time > t)] = apply(Q_hat[,t+1,], 1, function(x) min(x))[which(data_static$treat_start_time > t)]
        }
      } else {
        # Q-learning for single policy evaluation
        Q_next[which(data_static$treat_start_time > t & treat_policy$month > t+1)] = Q_hat[cbind(which(data_static$treat_start_time > t & treat_policy$month> t+1),t+1,1)]
        Q_next[which(data_static$treat_start_time > t & treat_policy$month <= t+1)] = Q_hat[cbind(which(data_static$treat_start_time > t & treat_policy$month <= t+1), t+1, 1 + treat_policy$choice[which(data_static$treat_start_time > t & treat_policy$month <= t+1)])]
      }
    }
   reg= get_regmat_Q(data_static, data_dynamic, data_treat_option, Q_next, foldid, H, t, learner, q_opt)

   Q_reg_list[[t]] = reg$reg_folds
   df = reg$df
   for (ff in unique(foldid)){
     cols <- grep("z.*", names(df), value=T)

     # treated up to t-1
     index_treated = which(data_static$treat_start_time < t & data_static$id %in% df$id[which(df$f==ff)])
     if (length(index_treated) > 0) {
       newx_treated = data.frame(df[which(df$id %in% data_static$id[index_treated]),cols]) %>% make_matrix
       if (dim(newx_treated)[1] > 0) {
         Q_hat[cbind(index_treated,t,data_static$treatment_option[index_treated]+1)] = pred_func(reg$reg_folds[[ff]], newx_treated, regress_type = "q", learner=learner, num_choices = num_choices)$pred
       }
     }

     # untreated up to t-1
     index_untreated = which(data_static$treat_start_time >= t & data_static$id %in% df$id[which(df$f==ff)])
     if (length(index_untreated) > 0) {
       for (kk in c(0:num_choices)) {
         newx_untreated = data.frame(df[which(df$id %in% data_static$id[index_untreated]),cols]) %>% make_matrix
         index_treatment = dim(newx_untreated)[2] - dim(treat_option)[2]
         newx_untreated[, index_treatment:(dim(newx_untreated)[2]-1)] = matrix(treat_option[kk+1,,drop=F],nrow=length(index_untreated),ncol=num_choices,byrow=TRUE)
         if (kk == 0) {
           newx_untreated[, dim(newx_untreated)[2]] = 0
         } else {
           newx_untreated[, dim(newx_untreated)[2]] = 1
         }
         Q_hat[cbind(index_untreated,t,kk+1)] = pred_func(reg$reg_folds[[ff]], newx_untreated, regress_type = "q", learner=learner, num_choices = num_choices)$pred
       }
     }
   }

    if (outcome_interest == "survival") {
      Q_hat[which(data_static$total_months < t), t, ] = 1
    } else if (outcome_interest == "lifetime") {
      lifetime = data_static$total_months
      Q_hat[which(data_static$total_months < t), t, ] = lifetime[which(data_static$total_months < t)]
    }
  }
  return(list(Q_hat = Q_hat,
              Q_reg_list = Q_reg_list))
}

get_regmat_Q = function(data_static, data_dynamic, data_treat_option, Q_next, foldid, H, t, learner, q_opt=F){
  num_choices = max(data_static$treatment_option)
  data_dynamic_t = data_dynamic[data_dynamic$month == t,,drop=F]
  data_t = merge(data_static[data_static$id %in% data_dynamic_t$id,-c("treat_start_time", "total_months", "outcome", "treatment_option", "censor"),drop=F], data_dynamic_t[,-c("evertreat","month")], by="id")
  Q_next_t = Q_next[data_static$id %in% data_dynamic_t$id]
  foldid_t = foldid[data_static$id %in% data_dynamic_t$id]
  data_static_t = data_static[data_static$id %in% data_dynamic_t$id,,drop=F]
  data_treat_option_t = data_treat_option[data_static$id %in% data_dynamic_t$id,,drop=F]

  index_not_treated = which(data_static_t$treat_start_time > t)
  index_treated = which(data_static_t$treat_start_time <= t)
  treat_duration = as.vector(t+1 - data_static_t$treat_start_time[index_treated])

  treat_history_t = matrix(NA, dim(data_static_t)[1], num_choices+1)
  treat_history_t[index_not_treated,] = 0
  treat_history_t[index_treated,] = cbind(data_treat_option_t[index_treated,,drop=F], treat_duration)
  z = cbind(as.matrix(data_t[,-c("id")]), treat_history_t)
  df = data.frame(z=z,
                  y=Q_next_t,
                  f=foldid_t,
                  id = data_t$id)

  num_z_covariates = dim(z)[2]
  colnames(df)[1:num_z_covariates] <- sapply(1:num_z_covariates, function(i) paste0("z.",i))
  reg_folds = vector("list", max(foldid))
  if (q_opt) {
    reg_folds[[1]] = reg_func(df, 2, regress_type = "q", learner=learner, num_choices = num_choices) # 2 is a dummy variable for the fold
  } else {
    for (f in unique(foldid)){
      reg_folds[[f]] = reg_func(df, f, regress_type = "q", learner=learner, num_choices = num_choices)
    }
  }
  return(list(reg_folds = reg_folds,
             df = df))
}

q_opt_eval = function(Q_regs, n, H, num_choices, setup_params, outcome_interest, setup, treat_once=F, learner="grf", Q_max=F) {
  init = get_init(n, setup_params)

  A_time = rep(H+1, n)
  A_choice = rep(1, n)
  X = get_X_init(init, A_time, A_choice, n, H, setup_params)
  Z = get_Z(X, setup_params)
  X0 = X
  Z0 = Z

  if (length(dim(Z)) == 2) {
    dim(Z) = c(n, H+1, 1) # make Z of dimension n by H+1 by # state features (in this case, 1)
  }


  Q_hat = array(NA, dim=c(n, H, num_choices + 1))

  for (tt in 1:H) {
    lifetime = get_lifetime(Z = Z, A_time = A_time, A_choice = A_choice, setup_params = setup_params)
    R = get_reward(lifetime, Z, A_choice, A_time, H, setup_params, outcome_interest = outcome_interest)
    if (tt == 1) {
      V_pi0 = mean(R)
    }
    data_formatted = get_data_formatted(Z, A_time, A_choice, lifetime, R)
    data_static = data_formatted$data_static
    data_dynamic = data_formatted$data_dynamic
    data_dynamic_t = data_dynamic[data_dynamic$month == tt,,drop=F]
    data_t = merge(data_static[data_static$id %in% data_dynamic_t$id,-c("treat_start_time", "total_months", "outcome", "treatment_option", "censor"),drop=F], data_dynamic_t[,-c("evertreat","month")], by="id")
    if (setup == "simu1") {
      data_t = data_t[,-c("patient_type")]
    }
    for (kk in 0:num_choices) {
      reg = Q_regs[[tt]][[1]]
      treat_history_t = matrix(NA, n, num_choices+1)
      if (num_choices > 1) {
        treat_option = model.matrix(~ factor(c(1:num_choices)) + 0)
        treat_option_0 = rep(0,dim(treat_option)[2])
        treat_option = rbind(treat_option_0, treat_option)
      } else {
        treat_option = matrix(c(0,1),2,1)
      }

      if (tt==1) {
        if (kk == 0) {
          treat_duration = 0
        } else {
          treat_duration = 1
        }
        treat_history_t = matrix(cbind(treat_option[kk+1,,drop=F],treat_duration), nrow=n,ncol=num_choices + 1,byrow=TRUE)
        Q_hat[, tt, kk+1] = pred_func(reg, cbind(as.matrix(data_t[,-c("id")]), treat_history_t), regress_type="q", learner=learner, num_choices = num_choices)$pred

      } else {
        if (kk == 0) {
          index_treated = which(A_time < H+1 & data_static$id %in% data_dynamic_t$id)
          index_untreated = which(A_time == H+1 & data_static$id %in% data_dynamic_t$id)
          if (length(index_untreated) > 0) {
            id_untreated = data_static$id[index_untreated]
            z = data_t[which(data_t$id %in% id_untreated), -c("id"), drop=F]
            treat_history_t = matrix(cbind(treat_option[kk+1,,drop=F],0), nrow=n, ncol= num_choices+1, byrow=TRUE)
            Q_hat[index_untreated, tt, kk+1] = pred_func(reg, cbind(z, treat_history_t[index_untreated,,drop=F]), regress_type="q", learner=learner, num_choices = num_choices)$pred
          }
          Q_hat[index_treated, tt, kk+1] = NA
        } else {
          index_untreated = which(A_time == H+1 & data_static$id %in% data_dynamic_t$id)
          if (length(index_untreated) > 0) {
            id_untreated = data_static$id[index_untreated]
            z = data_t[which(data_t$id %in% id_untreated), -c("id")]
            treat_history_t = matrix(cbind(treat_option[kk+1,,drop=F],1), nrow=n,ncol= num_choices+1, byrow=TRUE)
            Q_hat[index_untreated, tt, kk+1] = pred_func(reg, cbind(z, treat_history_t[index_untreated,,drop=F]), regress_type="q", learner=learner, num_choices = num_choices)$pred
          }

          index_treated = which(A_time < H+1 & A_choice == kk & data_static$id %in% data_dynamic_t$id)
          if (length(index_treated) > 0) {
            id_treated = data_static$id[index_treated]
            z = data_t[which(data_t$id %in% id_treated), -c("id")]
            treat_history_t = cbind(matrix(rep(treat_option[kk+1,,drop=F],each=n,nrow=n), n, num_choices), tt-A_time+1)
            Q_hat[index_treated, tt, kk+1] = pred_func(reg, cbind(z, treat_history_t[index_treated,, drop=F]), regress_type="q", learner=learner, num_choices = num_choices)$pred
          }

          index_treated_other = which(A_time < H+1 & A_choice != kk & data_static$id %in% data_dynamic_t$id)
          Q_hat[index_treated_other, tt, kk+1] = NA
        }
      }
    }
    if (Q_max) {
      A_opt_t = argmax(Q_hat[,tt,], rows = TRUE) # need to -1
    } else {
      A_opt_t = argmin(Q_hat[,tt,], rows = TRUE)  # ignoring the NAs because we only care about changing actions for trajs that still haven't started treating
    }
    index_treatment = which(A_time == H+1 & A_opt_t>1)
    index_same = which(!(A_time == H+1 & A_opt_t>1))

    if (length(index_treatment) > 0) {
      A_time[index_treatment] = tt
      A_choice[index_treatment] = unlist(A_opt_t[index_treatment]) - 1

      X_new = get_X_traj(X, A_time, A_choice, H, setup_params, tt, treat_once=treat_once) # only up to the t+1 th column is valid

      old_Z = Z
      Z = get_Z(X_new, setup_params)
      Z_new = postprocess_Z(old_Z[index_treatment,,,drop=F], Z[index_treatment,,,drop=F], A_time[index_treatment])
      Z[index_treatment,,] = Z_new
      Z[index_same,,] = old_Z[index_same,,,drop=F]

      X = X_new
    }
  }

  lifetime = get_lifetime(Z = Z, A_time = A_time, A_choice = A_choice, setup_params = setup_params)
  R = get_reward(lifetime, Z, A_choice, A_time, H, setup_params, outcome_interest = outcome_interest)

  return(list(A_choice = A_choice,
              A_time = A_time,
              V_pi = mean(R),
              V_pi0 = V_pi0,
              Z = Z,
              X0 = X0,
              Z0 = Z0))
}
