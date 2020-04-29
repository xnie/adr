options(error = recover)
library(abind)
library(ramify)
library(magrittr)
library(scales)

get_data = function(n, H, setup_params,outcome_interest="survival", e_rct=F, treat_once=F){
  init = get_init(n, setup_params)

  A0_time = rep(H+1, n)
  A0_choice = rep(1, n)
  X0 = get_X_init(init, A0_time, A0_choice, n, H, setup_params)
  Z0 = get_Z(X0, setup_params)
  action = get_A_sim(Z0, e_rct)
  A_time = action$A_time
  A_choice = action$A_choice

  X= get_X_traj(X0, A_time, A_choice, H, setup_params, 1, treat_once=treat_once)

  Z_ = get_Z(X, setup_params)

  if (length(dim(Z_)) == 2) {
    dim(Z_) = c(n, H+1, 1) # make Z of dimension n by H+1 by # state features (in this case, 1)
  }

  lifetime = get_lifetime(Z = Z_, A_time = A_time, A_choice = A_choice, setup_params = setup_params)
  action_post = postprocess_A(A_time = A_time, A_choice = A_choice, lifetime=lifetime, H=H)
  A_time = action_post$A_time
  A_choice = action_post$A_choice
  R = get_reward(lifetime, Z_, A_choice, A_time, H, setup_params, outcome_interest = outcome_interest)


  data_formatted = get_data_formatted(Z_, A_time, A_choice, lifetime, R)
  censor_prob = get_censor_prob(data_formatted$data_static, data_formatted$data_dynamic)

  return(list(data_static = data.table(data_formatted$data_static),
              data_dynamic = data.table(data_formatted$data_dynamic),
              censor_prob = censor_prob))

}
get_X_traj = function(X, A_time, A_choice, H, setup_params, start, treat_once = F){
  # start is the time index for which X[1:start] would be fixed, and we only roll out starting from X[,start+1]
  A = 1 - t(sapply(A_time, function(tt) as.numeric(1:H < tt)))
  for (t in (start+1):(H+1)) {
    if (sum(as.numeric(A_time<t)) > 0) {
      if (treat_once) {
        X[which(A_time<t), t, ] = transition(adrop(X[which(A_time<t), t-1, , drop=FALSE], drop=2), A_time[which(A_time<t)] + 1 == t, A_choice[which(A_time<t)], t, H, setup_params)
      } else {
        X[which(A_time<t), t, ] = transition(adrop(X[which(A_time<t), t-1, , drop=FALSE], drop=2), A[which(A_time<t),t-1], A_choice[which(A_time<t)], t, H, setup_params)
      }
    }
  }
  return(X)
}

postprocess_Z = function(Z0, Z, policy_time){
  for (t in 1:(H+1)) {
    if (sum(as.numeric(policy_time>=t)) > 0) {
      idx = which(policy_time >= t)
      Z[idx,t,] = Z0[idx,t,,drop=F]
    }
  }
  return(Z)
}


get_oracle_eval = function(n, H, param, setup_params, outcome_interest="survival", treat_once=F, X0=NULL, Z0=NULL) {
  A0_time = rep(H+1, n)
  A0_choice = rep(1, n)
  if (is.null(X0) | is.null(Z0)) {
    init0 = get_init(n, setup_params)

    X0 = get_X_init(init0, A0_time, A0_choice, n, H, setup_params)
    Z0 = get_Z(X0, setup_params)

    if (length(dim(Z0)) == 2) {
      dim(Z0) = c(n, H+1, 1) # make Z of dimension n by H+1 by # state features (in this case, 1)
    }
  }

  lifetime0 = get_lifetime(Z = Z0, A_time = A0_time, A_choice = A0_choice, setup_params = setup_params)
  R0 = get_reward(lifetime0, Z0, A0_choice, A0_time, H, setup_params, outcome_interest = outcome_interest)

  data_formatted0 = get_data_formatted(Z0, A0_time, A0_choice, lifetime0, R0)
  policy = get_treat_policy(data_formatted0$data_static, data_formatted0$data_dynamic, param, H)

  X= get_X_traj(X0, policy$month, policy$choice, H, setup_params, 1, treat_once=treat_once)

  Z = get_Z(X, setup_params)
  Z = postprocess_Z(Z0, Z, policy$month)

  lifetime = get_lifetime(Z = Z, A_time =policy$month, A_choice = policy$choice, setup_params = setup_params)
  R = get_reward(lifetime, Z, policy$choice, policy$month, H, setup_params, outcome_interest = outcome_interest)

  return(list(V_pi = mean(R),
              V_pi0 = mean(R0),
              lifetime0 = lifetime0,
              X0 = X0,
              Z0 = Z0,
              X = X,
              Z = Z,
              policy=policy))
}
get_X_init = function(init, A_time, A_choice, n, H, setup_params) {
  traj <- vector("list", H)
  traj[[1]] = init
  curr = init
  n = dim(init)[1]
  A = 1 - t(sapply(A_time, function(tt) as.numeric(1:H < tt)))
  for (t in 2:(H+1)) {
    new = transition(curr, A[,t-1], A_choice, t, H, setup_params)
    traj[[t]] = new
    curr = new
  }
  traj.concat = abind(traj, along=3)
  traj.concat = aperm(traj.concat, c(1,3,2)) # dimension n by H+1 by # state features (in this case, 1)
  return(traj.concat)
}

get_e_hat_rct = function(n, H, data_static) {
  num_choices = max(data_static$treatment_option)
  e_hat = matrix(unlist(lapply(1:H, function(t){
    rct_prob = rep(1 / num_choices / (H - t + 2), n)
    rct_prob
  })), n, H)


  for (t in 1:H) {
    e_hat[data_static$total_months < t, t] = 0
    e_hat[data_static$treat_start_time < t, t] = NA
  }
  replicate(num_choices, e_hat)
}

get_e_p_hat_rct = function(n, H, data_static) {
  num_choices = max(data_static$treatment_option)
  e_p_hat = matrix(unlist(lapply(1:H, function(t){
    rct_prob = rep(1 / num_choices / (H - t + 1), n)
    rct_prob
  })), n, H)

  e_p_hat = replicate(num_choices, e_p_hat)
  for (c in 1:num_choices) {
    if (c==1) {
      e_p_hat[,H,c] = 1
    } else {
      e_p_hat[,H,c] = 0
    }
  }
  for (t in 1:H) {
    e_p_hat[data_static$total_months <= t, t, ] = NA
    e_p_hat[data_static$treat_start_time <= t, t, ] = NA
  }
  e_p_hat
}
postprocess_A = function(A_time, A_choice, lifetime, H) {
  #A_time[A_time > lifetime] = lifetime[A_time > lifetime] + 1
  A_time[A_time > lifetime] = H+1
  A_choice[A_time > lifetime] = 1
  return(list(A_time = A_time,
              A_choice = A_choice))
}

