adr_hiv = function(A_time,
                   A_choice,
                   lifetime,
                   R,
                   mu_now,
                   mu_next,
                   mu_next_internal,
                   censor_prob,
                   traj_prob,
                   traj_after_prob_matrix){
  n = dim(mu_now)[1]
  H = dim(mu_now)[2]
  num_choices = max(A_choice)

  list_args <- Vectorize( function(a,b) c( as.list(a), as.list(b) ),
                         SIMPLIFY = FALSE)

  make_args_mtx <- function( alist ) {
    Reduce(function(x, y) outer(x, y, list_args), alist)
  }

  multi.outer <- function(f, ... ) {
    args <- make_args_mtx(list(...))
    apply(args, 1:length(dim(args)), function(a) do.call(f, a[[1]] ) )
  }

  Gamma_reg_func <- function(ii, tt, kk) {
    indicator = as.numeric((A_time[ii] >= tt) & (lifetime[ii] >= tt))
    if (indicator > 0) {
      if (tt == H) {
        indicator / traj_after_prob_matrix[ii, tt]  * (mu_now[ii, tt, kk] - mu_next[ii, tt, 1])
      } else{
        indicator / traj_after_prob_matrix[ii, tt]  * (mu_now[ii, tt, kk] - mu_next[ii, tt, kk])
      }
    } else {
      0
    }
  }

  Gamma_ipw_pi_func <- function(ii, tt, kk) {
    indicator = as.numeric((A_time[ii] == tt) & (A_choice[ii] == kk) & (lifetime[ii] >= tt))
    if (indicator > 0) {
      indicator / traj_prob[ii]  * (R[ii] - mu_now[ii, tt, kk])
    }
    else {
      0
    }
  }

  Gamma_ipw_pi0_func <- function(ii, tt, kk) {
    if (tt == H) {
      kk = 1
    }
    indicator = as.numeric((A_time[ii] == tt+1) & (A_choice[ii] == kk) & (lifetime[ii] >= tt+1))
    if (indicator > 0) {
      indicator / traj_prob[ii] * (R[ii] - mu_next_internal[ii, tt, kk])
    } else {
      0
    }
  }

  Gamma_ipw_pi_R_only_func <- function(ii, tt,kk) {
    indicator = as.numeric((A_time[ii] == tt )& (A_choice[ii]==kk) & (lifetime[ii] >= tt))
    if (indicator > 0) {
      indicator / traj_prob[ii] * R[ii]
    }
    else {
      0
    }
  }

  Gamma_ipw_pi0_R_only_func <- function(ii, tt, kk) {
    if (tt == H) {
      kk = 1
    }
    indicator = as.numeric(A_time[ii] == tt+1 & A_choice[ii]==kk  & lifetime[ii] >= tt+1)
    if (indicator > 0) {
      indicator / traj_prob[ii] * R[ii]
    } else {
      0
    }
  }

  weight_reg_func <- function(ii, tt, kk) {
    indicator = as.numeric(A_time[ii] >= tt)
    if (indicator  > 0) {
      indicator / traj_after_prob_matrix[ii, tt]
    }
    else {
      0
    }
  }

  weight_ipw_pi_func <- function(ii, tt, kk) {
    indicator = as.numeric((A_time[ii] == tt )& (A_choice[ii] == kk) & (lifetime[ii] >= tt))
    if (indicator > 0) {
      indicator / traj_prob[ii]
    }
    else {
      0
    }
  }

  weight_ipw_pi0_survived_func <- function(ii, tt, kk) {
    if (tt== H) {
      kk = 1
    }
    indicator = as.numeric((A_time[ii] == tt+1) & (A_choice[ii] == kk) & (lifetime[ii] >= tt+1))
    if (indicator > 0) {
      indicator / traj_prob[ii]
    } else {
      0
    }
  }

  weight_ipw_pi0_death_func <- function(ii, tt, kk) {
    indicator = as.numeric((A_time[ii] > tt) & (lifetime[ii] == tt))
    if (indicator > 0) {
      indicator / (traj_after_prob_matrix[ii, tt+1])
    } else {
      0
    }
  }

  Gamma_reg <- multi.outer(Gamma_reg_func, 1:n, 1:H, 1:num_choices)
  Gamma_ipw_pi <- multi.outer(Gamma_ipw_pi_func, 1:n, 1:H, 1:num_choices)
  Gamma_ipw_pi0 <- multi.outer(Gamma_ipw_pi0_func, 1:n, 1:H, 1:num_choices)
  Gamma_ipw_pi_R_only <- multi.outer(Gamma_ipw_pi_R_only_func, 1:n, 1:H, 1:num_choices)
  Gamma_ipw_pi0_R_only <- multi.outer(Gamma_ipw_pi0_R_only_func, 1:n, 1:H, 1:num_choices)

  weight_reg <- multi.outer(weight_reg_func, 1:n, 1:H, 1:num_choices)
  weight_ipw_pi <- multi.outer(weight_ipw_pi_func, 1:n, 1:H, 1:num_choices)
  weight_ipw_pi0_survived <- multi.outer(weight_ipw_pi0_survived_func, 1:n, 1:H, 1:num_choices)
  weight_ipw_pi0_death <- multi.outer(weight_ipw_pi0_death_func, 1:n, 1:H, 1:num_choices)

  Gamma_reg <- sweep(Gamma_reg, MARGIN=1, 1 / censor_prob, `*`)
  Gamma_ipw_pi <- sweep(Gamma_ipw_pi, MARGIN=1, 1 / censor_prob, `*`)
  Gamma_ipw_pi0 <- sweep(Gamma_ipw_pi0, MARGIN=1, 1 / censor_prob, `*`)
  Gamma_ipw_pi_R_only <- sweep(Gamma_ipw_pi_R_only, MARGIN=1, 1 / censor_prob, `*`)
  Gamma_ipw_pi0_R_only <- sweep(Gamma_ipw_pi0_R_only, MARGIN=1, 1 / censor_prob, `*`)

  weight_reg <- sweep(weight_reg, MARGIN=1, 1 / censor_prob, `*`)
  weight_ipw_pi <- sweep(weight_ipw_pi, MARGIN=1, 1 / censor_prob, `*`)
  weight_ipw_pi0_survived <- sweep(weight_ipw_pi0_survived, MARGIN=1, 1 / censor_prob, `*`)
  weight_ipw_pi0_death <- sweep(weight_ipw_pi0_death, MARGIN=1, 1 / censor_prob, `*`)

  return(list(Gamma_reg = Gamma_reg,
              Gamma_ipw_pi = Gamma_ipw_pi,
              Gamma_ipw_pi0 = Gamma_ipw_pi0,
              Gamma_ipw_pi_R_only = Gamma_ipw_pi_R_only,
              Gamma_ipw_pi0_R_only = Gamma_ipw_pi0_R_only,
              weight_reg = weight_reg,
              weight_ipw_pi = weight_ipw_pi,
              weight_ipw_pi0_survived = weight_ipw_pi0_survived,
              weight_ipw_pi0_death = weight_ipw_pi0_death))
}

get_adr_value = function(adr_ret, treat_policy, lifetime) {
  Gamma_reg = adr_ret$Gamma_reg
  Gamma_ipw_pi = adr_ret$Gamma_ipw_pi
  Gamma_ipw_pi0 = adr_ret$Gamma_ipw_pi0

  Gamma_ipw_pi_R_only = adr_ret$Gamma_ipw_pi_R_only
  Gamma_ipw_pi0_R_only = adr_ret$Gamma_ipw_pi0_R_only

  weight_reg = adr_ret$weight_reg
  weight_ipw_pi = adr_ret$weight_ipw_pi
  weight_ipw_pi0_survived = adr_ret$weight_ipw_pi0_survived
  weight_ipw_pi0_death = adr_ret$weight_ipw_pi0_death

  n = dim(Gamma_reg)[1]

  terms_reg = sapply(1:H, function(t) {
      sum((t >= treat_policy$month) * Gamma_reg[cbind(1:n, t, treat_policy$choice)])
  })

  terms_ipw_pi = sapply(1:H, function(t) {
      sum((t >= treat_policy$month) * Gamma_ipw_pi[cbind(1:n, t, treat_policy$choice)])
  })

  terms_ipw_pi0 = sapply(1:H, function(t) {
      sum((t >= treat_policy$month) * Gamma_ipw_pi0[cbind(1:n, t, treat_policy$choice)])
  })

  terms_ipw_pi_R_only = sapply(1:H, function(t) {
      sum((t >= treat_policy$month) * Gamma_ipw_pi_R_only[cbind(1:n, t, treat_policy$choice)])
  })

  terms_ipw_pi0_R_only = sapply(1:H, function(t) {
      sum((t >= treat_policy$month) * Gamma_ipw_pi0_R_only[cbind(1:n, t, treat_policy$choice)])
  })

  adr_reg = sum(terms_reg) / n
  adr_ipw_pi = sum(terms_ipw_pi) / n
  adr_ipw_pi0 = sum(terms_ipw_pi0) / n
  adv_adr = adr_reg + adr_ipw_pi - adr_ipw_pi0

  adr_ipw_pi_R_only = sum(terms_ipw_pi_R_only) / n
  adr_ipw_pi0_R_only = sum(terms_ipw_pi0_R_only) / n
  adv_ipw_R_only = adr_ipw_pi_R_only - adr_ipw_pi0_R_only

  weight_reg_t = sapply(1:H, function(t) {
      s = sum(weight_reg[cbind(1:n, t, treat_policy$choice)])
      if (s == 0) {
        s = 1
      }
      s
  })


  weight_ipw_pi_t= sapply(1:H, function(t) {
      s = sum((t >= treat_policy$month & lifetime >= t) * weight_ipw_pi[cbind(1:n, t, treat_policy$choice)]) +
            sum((!(t >= treat_policy$month & lifetime >= t)) * weight_reg[cbind(1:n, t, treat_policy$choice)])
      if (s==0) {
        s = 1
      }
      s
  })

  weight_ipw_pi0_t= sapply(1:H, function(t) {
      s = sum((t >= treat_policy$month & lifetime >= t) * weight_ipw_pi0_survived[cbind(1:n, t, treat_policy$choice)]) +
            sum((t >= treat_policy$month & lifetime == t) * weight_ipw_pi0_death[cbind(1:n, t, treat_policy$choice)]) +
            sum((!(t >= treat_policy$month & lifetime >= t)) * weight_reg[cbind(1:n, t, treat_policy$choice)])
      if (s==0) {
        s = 1
      }
      s
  })

  adr_reg_weighted = sum(terms_reg / weight_reg_t)
  adr_ipw_pi_weighted = sum(terms_ipw_pi / weight_ipw_pi_t)
  adr_ipw_pi0_weighted = sum(terms_ipw_pi0 / weight_ipw_pi0_t)

  adv_adr_weighted = adr_reg_weighted + adr_ipw_pi_weighted - adr_ipw_pi0_weighted

  adr_ipw_pi_R_only_weighted = sum(terms_ipw_pi_R_only / weight_ipw_pi_t)
  adr_ipw_pi0_R_only_weighted = sum(terms_ipw_pi0_R_only / weight_ipw_pi0_t)

  adv_ipw_R_only_weighted = adr_ipw_pi_R_only_weighted - adr_ipw_pi0_R_only_weighted

  return(list(adr = adv_adr,
              adr_weighted = adv_adr_weighted,

              adr_reg = adr_reg,
              adr_reg_weighted = adr_reg_weighted,

              adr_ipw = adv_ipw_R_only,
              adr_ipw_weighted = adv_ipw_R_only_weighted,

              mean_weight_reg = mean(weight_reg), # for debugging purpose
              mean_weight_pi = mean(weight_ipw_pi_t/n),
              mean_weight_pi0 = mean(weight_ipw_pi0_t/n)))
}
