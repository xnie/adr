aipw = function(A_time,
                A_choice,
                policy_time,
                policy_choice,
                R,
                censor_prob,
                traj_prob,
                traj_prob_matrix,
                Q_hat,
                V_hat,
                weighted=F) {
  # TODO: currently censoring is NOT implemented with AIPW

  ipw_term =  ipw(A_time,
                      A_choice,
                      policy_time,
                      policy_choice,
                      R,
                      censor_prob,
                      traj_prob,
                      weighted=weighted)
  n = dim(Q_hat)[1]
  H = dim(Q_hat)[2]

  V_weight_vector = rep(1,H)
  Q_weight_vector = rep(1,H)
  V_terms_vector = rep(1,H)
  Q_terms_vector = rep(1,H)

  V_weight_pi = rep(1,n)

  for (t in 1:H) {

    index_untreated =  which(t < A_time & t < policy_time)
    index_treated =  which(A_time == policy_time & t >= A_time & A_choice == policy_choice) # This assumes if never treated (A_time or policy_time = H+1), then the corresponding treatment choice is 1.

    Q_weight_pi_untreated = 1 / traj_prob_matrix[index_untreated, t]
    Q_weight_pi_treated = 1 / traj_prob_matrix[index_treated,t]

    Q_terms_vector[[t]] = sum(Q_hat[index_untreated, t, 1] * Q_weight_pi_untreated) + sum(Q_hat[cbind(index_treated, t, A_choice[index_treated]+1)] * Q_weight_pi_treated)
    V_terms_vector[[t]] = sum(V_hat[, t] * V_weight_pi)

    V_weight_vector[[t]] = sum(V_weight_pi)
    Q_weight_vector[[t]] = sum(Q_weight_pi_untreated) + sum(Q_weight_pi_treated)

    V_weight_pi = rep(0, n)
    V_weight_pi[index_untreated] = Q_weight_pi_untreated
    V_weight_pi[index_treated] = Q_weight_pi_treated
  }

  if (weighted) {
    res = ipw_term - sum(Q_terms_vector / Q_weight_vector) + sum(V_terms_vector / V_weight_vector)
  } else {
    res = ipw_term - sum(Q_terms_vector) / n + sum(V_terms_vector) / n
  }
  return(res)
}

get_V_hat = function(Q_hat, A_time, A_choice, policy_time, policy_choice, lifetime, outcome_interest){
  n = dim(Q_hat)[1]
  H = dim(Q_hat)[2]
  num_choices = dim(Q_hat)[3]-1
  V_hat = matrix(0, n, H)
  for (tt in 1:H) {
    V_hat[which(A_time < tt),tt] = Q_hat[cbind(which(A_time<tt), tt, A_choice[which(A_time<tt)]+1)]
    V_hat[which(A_time >= tt & policy_time <= tt),tt] = Q_hat[cbind(which(A_time >= tt & policy_time <= tt), tt, policy_choice[which(A_time >= tt & policy_time <= tt)] + 1)]
    V_hat[which(A_time >= tt & policy_time > tt),tt] = Q_hat[which(A_time >= tt & policy_time > tt), tt, 1]
    if (outcome_interest == "survival") {
      V_hat[which(lifetime < tt), tt ] = 1
    } else if (outcome_interest == "lifetime") {
      V_hat[which(lifetime < tt), tt ] = lifetime[which(lifetime < tt)]
    }
  }
  V_hat
}
