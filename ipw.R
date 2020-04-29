ipw = function(A_time,
                   A_choice,
                   policy_time,
                   policy_choice,
                   R,
                   censor_prob,
                   traj_prob,
                   weighted=T) {
  if (is.null(traj_prob)) {
    stop("must supply traj_prob to ipw")
  }
  n = length(R)

  weight_pi = as.numeric(A_time == policy_time & A_choice == policy_choice) / traj_prob / censor_prob

  ipw_terms_pi = R * weight_pi

  if (weighted) {
    if (sum(weight_pi) == 0) {
      result = 0
    } else {
      result = sum(ipw_terms_pi) / (sum(weight_pi))
    }
    return(result)
  }
  else {
    return(mean(ipw_terms_pi))
  }
}
