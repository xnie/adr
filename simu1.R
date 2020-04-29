library(abind)
library(ramify)
library(magrittr)
source("simu_common.R")

# params
NOACTION=0
NONINV=1
INV=2
x_max = 10
y_max = 16
z_max = 2

get_policy_params = function(){
  params_treat_t = cbind(0,0,1, c(1:11))
  params_noninv = cbind(params_treat_t, matrix(rep(c(0,0,0,1),each=11),nrow=11))
  params_inv = cbind(params_treat_t, matrix(rep(c(0,0,1,0),each=11),nrow=11))

  params_treat = expand.grid(c(seq(0.2, 1, 0.5), seq(1,5,2)), 1, c(0,1), seq(1,10,2))
  colnames(params_treat) = names(params_treat_t)
  params_treat = as.matrix(params_treat)
  params_triangle_inv = params_treat
  params_triangle_inv[,1] = params_treat[,1]/2
  params_triangle = cbind(params_treat, params_triangle_inv)
  colnames(params_triangle) = names(params_inv)
  params_triangle = as.matrix(params_triangle)

  params_cross = expand.grid(c(seq(0.2, 1, 0.5), seq(1,5,2)), 1, c(0,1), seq(1,10,2), -0.5, 1, c(0,1), seq(-5, 5, 3))
  colnames(params_cross) = names(params_triangle)
  params_cross = as.matrix(params_cross)

  params = rbind(params_noninv, params_inv, params_triangle, params_cross)
  params
}

get_A_sim = function(Z, e_rct=F){
  n = dim(Z)[1]
  H = dim(Z)[2] - 1

  A_choice = sample.int(2, n, replace=TRUE)
  A_time = sample.int(H+1, n, replace=TRUE)
  A_choice[A_time==H+1] = 1

  return(list(A_time = A_time,
              A_choice = A_choice))
}

transition = function(prev, A, A_choice, t, H, setup_params){
  n = dim(prev)[1]

  new_x = rep(0,n)
  new_y = rep(0,n)
  new_z = prev[,4]

  obs_idx = which(A == 0)
  noninv_idx = which((A_choice == NONINV) & (A == 1))
  inv_idx = which((A_choice == INV) & (A == 1))

  random_step = matrix(rnorm(n * 2, sd=0.5), n, 2)

  new_x[obs_idx] = abs(prev[obs_idx,1] + random_step[obs_idx, 1])
  new_y[obs_idx] = abs(prev[obs_idx,2] + (0.5 * prev[obs_idx, 1] + random_step[obs_idx, 2]))

  new_x[noninv_idx] = abs(prev[noninv_idx,1] + random_step[noninv_idx, 1])
  new_y[noninv_idx] = abs(prev[noninv_idx,2] / 2 +  random_step[noninv_idx, 2])

  new_x[inv_idx] = (abs(pmax(prev[inv_idx,1]^2, 1.5 * prev[inv_idx,1]) + random_step[inv_idx, 1] - prev[inv_idx,1])) + prev[inv_idx, 1]
  new_y[inv_idx] = rep(0, length(inv_idx))

  z_0_idx = which(prev[,4] == 0) # always die
  z_1_idx = which(prev[,4] == 1) # always live
  z_2_idx = which(prev[,4] == 2)

  new_alive_low = as.numeric(rbinom(n, 1, exp(-0.02 * prev[,2])))
  new_alive_high = as.numeric(rbinom(n, 1, exp(-0.06 * prev[,2])))

  new_alive = new_alive_high
  new_alive[prev[,2] >= 14] = 0
  new_alive[prev[,2] <= 5] = new_alive_low[prev[,2] <= 5 ]
  new_alive[prev[,3] == 0] = 0

  new_alive[z_0_idx] = 0
  new_alive[z_1_idx] = 1

  new_x[new_alive == 0] = 0
  new_y[new_alive == 0] = 0

  new_x = pmin(new_x, x_max)
  new_y = pmin(new_y, y_max)

  new = matrix(c(new_x, new_y, new_alive, new_z), nrow = n)
  return(new)
}


get_reward = function(lifetime, Z, A_choice, A_time, H, setup_params, outcome_interest = "survival"){
  return(lifetime)
}

get_Z = function(X, setup_params) {

  obs_noise = setup_params$obs_noise

  n = dim(X)[1]
  H = dim(X)[2] - 1
  Z = X
  Z[,,1:2] = Z[,,1:2] + obs_noise * array(rnorm(n * (H+1) * 2),c(n, H+1, 2))
  Z[,,1] = pmin(x_max, pmax(0, Z[,,1]))
  Z[,,2] = pmin(y_max, pmax(0, Z[,,2]))
  return(Z)

}

get_init = function(n, setup_params) {
  z_prob = setup_params$z_prob

  init_x =  rexp(n, rate = 1)
  init_y = 0.5* rexp(n, rate = 3)
  init_alive = rep(1, n)
  init_z = rep(0,n)
  z = rmultinom(n, size = 1, prob = z_prob)
  init_z[which(z[1,] == 1)] = 0
  init_z[which(z[2,] == 1)] = 1
  init_z[which(z[3,] == 1)] = 2
  init = matrix(c(init_x, init_y, init_alive, init_z), nrow = n)
  return(init)

}

get_lifetime = function(Z, A_time, A_choice, setup_params) {
  rowSums(Z[,,3])
}

get_data_formatted = function(Z, A_time, A_choice, lifetime, R) {
  n = dim(Z)[1]
  H = dim(Z)[2]-1
  dynamic_max_t = pmin(lifetime, H)
  data_static = data.frame(id = c(1:n), treatment_option = A_choice, total_months = lifetime, patient_type = Z[,1,4], treat_start_time = A_time, outcome = R, censor=0)
  data_dynamic = Reduce(rbind, lapply(1:n, function(id){
    data.frame(id = id,
               month = c(1:dynamic_max_t[[id]]),
               evertreat= as.numeric(A_time[[id]] <= c(1:dynamic_max_t[[id]])),
               covar1 = Z[id, 1:dynamic_max_t[[id]], 1],
               covar2 = Z[id, 1:dynamic_max_t[[id]], 2])
}))
  return(list(data_static = data.table(data_static),
              data_dynamic = data.table(data_dynamic)))

}

get_censor_prob = function(data_static, data_dynamic) {
  n = dim(data_static)[1]
  return(rep(1,n))
}

get_treat_policy = function(data_static, data_dynamic, param, H) {

  n = dim(data_static)[1]

  policy_choice = rep(NOACTION, n)
  policy_time = rep(H+1, n)

  noninv_time = data_dynamic%>%
                      group_by(id) %>%
                      filter(month == month[min(which((param[[1]] * covar1 + param[[2]] * covar2 + param[[3]] * month >= param[[4]]) &
                                                      (param[[5]] * covar1 + param[[6]] * covar2 + param[[7]] * month < param[[8]])))]) %>%
                      select(id, month)%>%
                      ungroup()

  inv_time = data_dynamic%>%
                      group_by(id) %>%
                      filter(month == month[min(which((param[[1]] * covar1 + param[[2]] * covar2 + param[[3]] * month >= param[[4]]) &
                                                      (param[[5]] * covar1 + param[[6]] * covar2 + param[[7]] * month >= param[[8]])))]) %>%
                      select(id, month)%>%
                      ungroup()

  noninv_time = process_treat_time(data_static, noninv_time)
  inv_time = process_treat_time(data_static, inv_time)

  choice_noninv_index = noninv_time$month <= inv_time$month
  treat_time = inv_time
  treat_time$month[choice_noninv_index] = noninv_time$month[choice_noninv_index]
  treat_time$choice = rep(INV, n)
  treat_time$choice[choice_noninv_index] = NONINV

  return(treat_time)

}

process_treat_time = function(data_static, treat_time){
    never_treat_id = setdiff(data_static$id, treat_time$id)
    if (length(never_treat_id) > 0) {
      never_treat = cbind(never_treat_id, H+1)
      colnames(never_treat) = c("id", "month")
      treat_time_all = rbind(as.data.frame(treat_time),  data.frame(never_treat))
    } else {
      treat_time_all = as.data.frame(treat_time)
    }
    treat_time_all = treat_time_all[match(data_static$id, treat_time_all$id),]
    treat_time_all
}

get_mu_final = function(die_idx, live_idx, normal_idx, mu_now, mu_next, mu_next_internal, H) {
  n = length(die_idx) + length(live_idx) + length(normal_idx)
  num_choices = dim(mu_now)[3]
  mu_now_final = array(0, dim=c(n, H, num_choices))
  mu_next_internal_final = array(0, dim=c(n, H, num_choices))
  mu_next_final = array(0, dim=c(n, H, num_choices))

  mu_now_final[normal_idx,,] = mu_now
  mu_next_final[normal_idx,,] = mu_next
  mu_next_internal_final[normal_idx,,] = mu_next_internal

  mu_now_final[die_idx,,] = array(1, dim=c(length(die_idx), H, num_choices))
  mu_now_final[live_idx,,] = array(H+1, dim=c(length(live_idx), H, num_choices))

  mu_next_final[die_idx,,] = array(1, dim=c(length(die_idx), H, num_choices))
  mu_next_final[live_idx,,] = array(H+1, dim=c(length(live_idx), H, num_choices))

  mu_next_internal_final[die_idx,,] = array(1, dim=c(length(die_idx), H, num_choices))
  mu_next_internal_final[live_idx,,] = array(H+1, dim=c(length(live_idx), H, num_choices))

  return(list(mu_now = mu_now_final,
              mu_next = mu_next_final,
              mu_next_internal = mu_next_internal_final))
}

