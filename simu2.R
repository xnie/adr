library(abind)
library(ramify)
library(magrittr)
library(latex2exp)
library(scales)
library(RColorBrewer)
source("simu_common.R")

pal = brewer.pal(n = 8, name = "Dark2")
get_policy_params = function() {
  params = cbind(0, -1, c(1:11)) # always treat, and never treat corresponding to the 1st and 11th policy here
  params = rbind(params, cbind(1,0, seq(-.5, 5, 0.5)))
  params_0 = cbind(1, expand.grid(c(-1/4, -1/3, -1/2, -1, -2, -3, -4), seq(1, 15, 1)))
  colnames(params) = names(params_0)
  params = rbind(params, params_0)
  params
}

get_A_sim = function(Z, e_rct=F){
  n = dim(Z)[1]
  H = dim(Z)[2] - 1

  #A_choice = sample.int(2, n, replace=TRUE)
  #A_time = sample.int(H+1, n, replace=TRUE)
  #A_choice[A_time==H+1] = 1
#
  #return(list(A_time = A_time,
  #            A_choice = A_choice))
   n = dim(Z)[1]
   H = dim(Z)[2] - 1
   stop_matrix = matrix(unlist(lapply(1:H, function(t) {
     stop_prob = 1/(1+exp(-(Z[,t,1]-1.5)) + exp(-(t-3)))
     stop = rbinom(n, 1, stop_prob)
     stop
   }
   )), n, H)
   A_time = unlist(apply(stop_matrix, 1, function(row) {
     nonzero = which(row!=0)
     if (length(nonzero) == 0) {
       H+1
     } else{
       head(which(row!=0),1)
     }
   }))
   A_choice = rep(1,n)
   return(list(A_time = A_time,
               A_choice = A_choice))
}

transition = function(prev, A, A_choice, t, H, setup_params){

  sigma = setup_params$sigma
  delta = sigma / sqrt(2 * H) * rnorm(dim(prev)[1])

  new = prev + delta + 1 / (1+exp(0.3 * prev[,1])) * (1-A)
  new[prev[,1] < -0.5,] = prev[prev < -0.5,]

  new
}


get_reward = function(lifetime, Z, A_choice, A_time, H, setup_params, outcome_interest = "survival"){
  beta = setup_params$beta
  success = as.numeric(Z[,H+1,1] > 0)
  cost = (A_time - 1) / H
  return(beta * success - cost)
}


get_Z = function(X, setup_params) {

  obs_noise = setup_params$obs_noise
  n = dim(X)[1]
  H = dim(X)[2] - 1
  Z = X + obs_noise * array(rnorm(n*(H+1)), dim = c(n, H+1, 1))

  return(Z)

}

get_init = function(n, setup_params) {

  sigma = setup_params$sigma
  init = sigma * rnorm(n)
  return(matrix(init, n, 1))

}


get_lifetime = function(Z, A_time, A_choice, setup_params) {
  n = dim(Z)[1]
  H = dim(Z)[2] - 1
  rep(H+1, n)
}


get_data_formatted = function(Z, A_time, A_choice, lifetime, R) {
  n = dim(Z)[1]
  H = dim(Z)[2]-1
  dynamic_max_t = pmin(lifetime, H)
  data_static = data.frame(id = c(1:n), treatment_option = A_choice, total_months = lifetime, treat_start_time = A_time, outcome = R, censor=0)
  data_dynamic = Reduce(rbind, lapply(1:n, function(id){
    data.frame(id = id,
               month = c(1:dynamic_max_t[[id]]),
               evertreat= as.numeric(A_time[[id]] <= c(1:dynamic_max_t[[id]])),
               indicator = Z[id, 1:dynamic_max_t[[id]], 1])
}))
  return(list(data_static = data.table(data_static),
              data_dynamic = data.table(data_dynamic)))

}

get_censor_prob = function(data_static, data_dynamic) {
  n = dim(data_static)[1]
  return(rep(1,n))
}

get_treat_policy = function(data_static, data_dynamic, param, H) {
  treat_time = data_dynamic%>%
                      group_by(id) %>%
                      filter(month == month[min(which(param[[1]] * indicator >= param[[3]] + param[[2]] * month))]) %>%
                      select(id, month)%>%
                      ungroup()
  never_treat_id = setdiff(data_static$id, treat_time$id)
  if (length(never_treat_id) > 0) {
    never_treat = cbind(never_treat_id, H+1)
    colnames(never_treat) = c("id", "month")
    treat_time_all = rbind(as.data.frame(treat_time),  data.frame(never_treat))
  } else {
    treat_time_all = as.data.frame(treat_time)
  }
  treat_time_all = treat_time_all[match(data_static$id, treat_time_all$id),]
  treat_time_all$choice = 1
  treat_time_all$choice[treat_time_all$month == H+1] = 1
  return(treat_time_all)
}
