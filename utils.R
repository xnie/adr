library(mgcv)

get_mu_next = function(data_static, mu_next_internal, d_hat, H, outcome_interest="survival"){
  n = dim(data_static)[1]
  lifetime = data_static$total_months
  num_choices = max(data_static$treatment_option)
  est = array(NA, dim=c(n, H, num_choices))
  if (outcome_interest == "survival") {
    mu_next_raw_dead = matrix(1,n,H)
  } else if (outcome_interest == "lifetime") {
    mu_next_raw_dead = matrix(rep(1:H,each=n),nrow=n)
  } else if (outcome_interest == "other"){
    mu_next_raw_dead = matrix(0, n, H)
  } else {
    mu_next_raw_dead = NA # TODO: implement this for other outcomes
  }
  for (j in 1:H) {
    for (c in 1:num_choices) {
      est[,j,c] = d_hat[,j] * mu_next_raw_dead[,j] + (1-d_hat[,j]) * mu_next_internal[,j,c]
    }
    if (outcome_interest == "survival") {
      est[which(lifetime < j), j, ] = 1
    } else if (outcome_interest == "lifetime") {
      est[which(lifetime < j), j, ] = lifetime[which(lifetime<j)]
    } else {
      est[which(lifetime < j), j, ] = NA # TODO: implement this for other outcomes
    }

    est[which(data_static$treat_start_time < j), j,] = NA
  }
  return(est)
}

get_est= function(data_static, data_dynamic, foldid, H, regress_type="e", learner="grf", outcome_interest="survival", e_hat=NULL){
  if (regress_type == "mu_next_internal") {
    if (is.null(e_hat)) {
      stop("needs to supply e_hat for learning mu_next_internal")
    }
  }
  n = dim(data_static)[1]
  nfold = length(unique(foldid))
  num_choices = max(data_static$treatment_option)
  lifetime = data_static$total_months

  if (regress_type == "e" | regress_type == "mu_now" | regress_type == "mu_next_internal"){
    est = array(NA, dim=c(n, H, num_choices))
  }
  if (regress_type == "c" | regress_type == "d"){
    if (outcome_interest == "other" & regress_type == "d") {
      return(matrix(0, n, H))
    }
    est = matrix(NA, n, H)
  }
  regmat_res = regmat_func(data_static, data_dynamic, foldid, H, regress_type=regress_type, e_hat=e_hat)
  if ((regress_type == "mu_next_internal" & (outcome_interest == "survival" | outcome_interest == "lifetime"))){
    regmat_t_all = lapply(1:(H-1), function(t){
      regmat_res[[t]]$df
    })
  } else {
    regmat_t_all = lapply(1:H, function(t){
      regmat_res[[t]]$df
    })
  }

  regmat = Reduce(rbind,
                  regmat_t_all[!vapply(regmat_t_all, is.null, logical(1))])


  for (ff in 1:nfold) {
    print(ff)
    reg_joint_time = reg_func(regmat, ff, regress_type=regress_type, learner=learner, num_choices=num_choices)
    for (tt in 1:H) {
      if (tt == H) {
        if (regress_type == "mu_next_internal" & (outcome_interest == "survival" | outcome_interest == "lifetime")) {
          next
        }
      } else {
        reg = reg_joint_time
      }
      idx = (foldid == ff) & (data_static$id %in% regmat_res[[tt]]$data_t$id)
      if (sum(idx) == 0) {
        next
      }
      newx = data.frame(regmat_res[[tt]]$data_t[regmat_res[[tt]]$foldid_t == ff, -c("id")]) %>% make_matrix
      colnames(newx)[1:dim(newx)[2]] <- sapply(1:dim(newx)[2], function(i) paste0("z.",i))
      pred_res = pred_func(reg, newx, regress_type=regress_type, learner = learner, num_choices=num_choices)
      if (regress_type == "e" | regress_type == "mu_now" | regress_type == "mu_next_internal") {
        est[idx, tt, ] = pred_res$pred
      } else if (regress_type == "c" | regress_type == "d") {
        est[idx, tt] = pred_res$pred
      }
    }
  }

if (regress_type == "e") {

        for (j in 1:H) {
          est[which(data_static$total_months < j), j, ] = 0
          est[which(data_static$treat_start_time < j), j, ] = NA
      }
  } else if (regress_type == "d") {
        for (j in 1:H) {
          est[which(data_static$total_months < j), j] = 1
          est[which(data_static$treat_start_time < j), j] = NA
        }
  } else if (regress_type == "mu_now") {
        for (j in 1:H) {
          if (outcome_interest == "survival") {
            est[which(data_static$total_months < j), j, ] = 1
          } else if (outcome_interest == "lifetime") {
            est[which(data_static$total_months < j), j, ] = lifetime[which(data_static$total_months < j)]
          }
          est[which(data_static$treat_start_time < j), j, ] = NA
        }
  } else if (regress_type == "mu_next_internal") {
        for (j in 1:H) {
          if (outcome_interest == "survival") {
            if (j == H) {
              for (c in 1:num_choices) {
                if (c!=1) {
                  est[, j, c] = NA # TODO: change for other outcomes
                } else {
                  est[, j, c] = 0 # TODO: change for other outcomes
                }
              }
            }
            est[which(data_static$total_months < j), j, ] = 1
          } else if (outcome_interest == "lifetime") {
            if (j == H) {
              for (c in 1:num_choices) {
                if (c!=1) {
                  est[, j, c] = NA # TODO: change for other outcomes
                } else {
                  est[, j, c] = H+1 # TODO: change for other outcomes
                }
              }
            }
            est[which(data_static$total_months < j), j, ] = lifetime[which(data_static$total_months < j)]
          }
          est[which(data_static$treat_start_time < j), j,  ] = NA
        }
  }
  return(est)
}

pred_func = function(reg, data, regress_type="e", learner = "grf", num_choices=1) {
  if (is.null(reg)) {
    pred = rep(NA, dim(data)[1])
    return(list(pred=pred))
  } else {
    if (regress_type == "e") {
      preds = sapply(0:num_choices, function(c){
          if (learner == "grf") {
            pred = predict(reg[[c+1]], data)$predictions
          } else {
            stop("only learner grf is supported currently")
          }
        })
      preds = sweep(preds, 1, rowSums(preds), `/`)
      preds = preds[,-1] # drop the predicting no action column because of normalization
      return(list(pred=preds))
    }
    if (regress_type == "mu_now" | regress_type == "mu_next_internal"){
       preds = sapply(1:num_choices, function(c){
          if (learner == "grf") {
            pred = predict(reg[[c]], data)$predictions
          } else {
            stop("only learner grf is supported currently")
          }
        })
      return(list(pred=preds))
    }

    if (regress_type == "c" | regress_type == "d" | regress_type == "q") {
      if (learner == "grf") {
            pred = predict(reg, data)$predictions
      } else {
        stop("only learner grf is supported currently")
      }
      return(list(pred=as.vector(pred)))
    }
  }
}

make_matrix = function(x) stats::model.matrix(~.-1, x)

train_func = function(x,y,f,w, regress_type="e", learner="grf", weight=NULL){
  newfolds <- unique(f)
  custom_folds <- vector("list", length(newfolds))
  i = 1
  for( id in newfolds){
    custom_folds[[i]] <- which(f == id)
    i = i+1
  }
  if (regress_type == "e" | regress_type == "c" | regress_type == "d" | regress_type == "mu_now" | regress_type == "mu_next_internal" | regress_type == "q") {
    if (learner == "grf") {
      if (regress_type == "mu_next_internal") {
        if (is.null(weight)) {
          stop("regress type mu_next_interal needs to provide regression weight")
        }
        fit = grf::regression_forest(x, y, num.threads = 1, num.trees = 500, sample.weights = weight)
      } else {
        fit = grf::regression_forest(x, y, num.threads = 1, num.trees = 500)
      }
    } else {
      stop("only learner grf is supported currently")
    }
  }
  return(fit)
}

reg_func = function(regmat, fold, regress_type = "e", learner="grf", num_choices=1) {
  data = regmat[regmat$f != fold,]

  cols <- grep("z.*", names(data), value=T)
  x = data.frame(data[, cols]) %>% make_matrix


  f = data[,"f"]
  if (max(f) > fold) {
    f[f==max(f)] = fold # reorder foldid so it starts with 1, 2, 3,...
  }
  if (regress_type=="e"){
    model = lapply(0:num_choices, function(choice) {
      w = NULL
      y = data[,"w"] == choice
      model_temp = train_func(x,y,f,w, regress_type=regress_type, learner=learner)
    })
  }
  if (regress_type == "mu_now" | regress_type == "mu_next_internal"){
    model = lapply(1:num_choices, function(choice) {
        relevant = data[,"c"] == choice
        w = NULL
        y = data[,'y']
        if (regress_type == "mu_now") {
          weight = NULL
        } else {
          weight = 1/data[,"e"]
        }
        model_temp = train_func(x[relevant,,drop=F],y[relevant],f[relevant],w, regress_type=regress_type, learner=learner, weight=weight[relevant])
    })
  } else if (regress_type == "c" | regress_type == "d") {
    w = NULL
    y = data[,"w"]

    model = train_func(x,y,f,w, regress_type=regress_type, learner=learner)
  } else if (regress_type == "q") {
      w = NULL
      y = data[,'y']
      model = train_func(x,y,f,w, regress_type=regress_type, learner=learner)
  }
  return(model)
}

regmat_func = function(data_static, data_dynamic, foldid, H, regress_type="e", e_hat=NULL) {
  A_time = data_static$treat_start_time
  A_choice = data_static$treatment_option
  lifetime = data_static$total_months
  outcome = data_static$outcome
  if (regress_type == "c") {
    censor = data_static$censor
  }
  regmat_t_all = lapply(1:H, function(t) {
    data_dynamic_t = data_dynamic[data_dynamic$month == t,,drop=F]
    if (regress_type == "c") {
      data_t = merge(data_static[data_static$id %in% data_dynamic_t$id,-c("treat_start_time", "total_months", "outcome", "treatment_option", "censor"),drop=F], data_dynamic_t[,-c("evertreat")], by="id")
      censor_t = censor[data_static$id %in% data_dynamic_t$id]
    } else {
      data_t = merge(data_static[data_static$id %in% data_dynamic_t$id,-c("treat_start_time", "total_months", "outcome", "treatment_option", "censor"),drop=F], data_dynamic_t[,-c("evertreat")], by="id")
    }

    A_time_t = A_time[data_static$id %in% data_dynamic_t$id]
    A_choice_t = A_choice[data_static$id %in% data_dynamic_t$id]
    lifetime_t = lifetime[data_static$id %in% data_dynamic_t$id]
    foldid_t = foldid[data_static$id %in% data_dynamic_t$id]
    outcome_t = outcome[data_static$id %in% data_dynamic_t$id]
    e_hat_t = e_hat[data_static$id %in% data_dynamic_t$id,,,drop=F]
    if (regress_type == "e") {
        relevant = (t <= A_time_t) & (t <= lifetime_t)
        A_time_yet = A_time_t[relevant] == t
    } else if (regress_type == "mu_now") {
          relevant = (A_time_t == t) & (lifetime_t >= t)
          outcome_now = outcome_t[relevant]
    } else if (regress_type == "mu_next_internal") {
          relevant = (A_time_t == t+1) & (lifetime_t > t)
          outcome_next = outcome_t[relevant]
    } else if (regress_type == "d") {
         relevant = (t < A_time_t) & (t <= lifetime_t)
         D_time_yet = (lifetime_t[relevant] == t)
    } else if (regress_type == "c") {
      relevant =  t <= lifetime_t
      censor_yet = (lifetime_t[relevant] == t) & (censor_t[relevant] == 1)
    }
    if (all(relevant==FALSE)) {
      df = NULL
    } else{ # TODO: add in history
      if (regress_type == "e") {
        df = data.frame(
          z=data_t[relevant, -c("id"), with=F],
          w=as.numeric(A_time_yet) * A_choice_t[relevant],
          f=foldid_t[relevant])
      } else if (regress_type == "mu_now") {
        df = data.frame(
          z=data_t[relevant, -c("id"), with=F],
          c=A_choice_t[relevant],
          y=as.numeric(outcome_now),
          f=foldid_t[relevant])
      } else if (regress_type == "mu_next_internal") {
        if (t==H) {
          e = rep(1, sum(relevant))
          e[which(A_choice_t[relevant]!=1)] = 0
        } else {
          e=e_hat_t[cbind(which(relevant), t+1, A_choice_t[relevant])]
        }
        if (any(is.na(e))) {
          browser()
        }
        df = data.frame(
          z=data_t[relevant, -c("id"), with=F],
          c=A_choice_t[relevant],
          y=as.numeric(outcome_next),
          e=e,
          f=foldid_t[relevant])
      } else if (regress_type == "d") {
        df = data.frame(
          z=data_t[relevant, -c("id"), with=F],
          w=as.numeric(D_time_yet),
          f=foldid_t[relevant])
      } else if (regress_type == "c") {
        df = data.frame(
          z=data_t[relevant, -c("id"), with=F],
          w=as.numeric(censor_yet),
          f=foldid_t[relevant])
      }
      num_z_covariates = dim(data_t)[2] - 1
      colnames(df)[1:num_z_covariates] <- sapply(1:num_z_covariates, function(i) paste0("z.",i))
    }
    list(df = df,
         foldid_t = foldid_t,
         data_t = data_t)
  })
  return(regmat_t_all)
}

get_traj_prob = function(e_hat, A_time, A_choice, H) {
  # e_hat is P(W=1 | history so far)
  # e_hat is n by H
  # return probability of traj at time t if we were to start treating at time t.
  # except if lifetime < t, then return never treating probability.
  # when this function is called, we always look at the traj probability of min(start_treat_time, lifetime)

  n = dim(e_hat)[1]
  traj_prob = rep(0,n)
  for (t in 1:(H+1)) {
    index = which(A_time==t)
    if (t == 1) {
      traj_prob[index] = e_hat[cbind(index, t, A_choice[index])]
    } else if (t == H+1){
      traj_prob[index] = apply((1-rowSums(e_hat[index,1:(t-1), , drop=F], dims=2)),1,prod)
    } else {
      traj_prob[index] = apply((1-rowSums(e_hat[index,1:(t-1), , drop=F], dims=2)),1,prod) * e_hat[cbind(index, t, A_choice[index])]
    }
  }
  traj_prob
}

get_traj_prob_matrix = function(e_hat, A_time, A_choice, H){
  n = dim(e_hat)[1]
  traj_prob = matrix(NA, n, H)

  for (t in 1:H) {
    untreat_index = which(A_time>t)
    traj_prob[untreat_index, t] = apply((1-rowSums(e_hat[untreat_index, 1:t, , drop=F], dims=2)),1,prod)

    treat_index = which(A_time==t)
    traj_prob[treat_index, t] = apply((1-rowSums(e_hat[treat_index,1:(t-1), , drop=F], dims=2)),1,prod) * e_hat[cbind(treat_index, t, A_choice[treat_index])]

    treated_index = which(A_time < t)
    traj_prob[treated_index, t] = traj_prob[treated_index, t-1]
  }
  traj_prob
}

get_traj_after_prob_matrix = function(e_hat, A_time,  H){
  n = dim(e_hat)[1]
  traj_prob = matrix(NA, n, H+1)
  for (t in 1:(H+1)) {
    if (t==1) {
      traj_prob[,t] = rep(1,n)
    } else{
      index = which(A_time>=t)
      traj_prob[index,t] = apply((1-rowSums(e_hat[index,1:(t-1), , drop=F], dims=2)),1,prod)
    }
  }
  traj_prob
}


get_foldid = function(n, nfold) {
  if (n - nfold * floor(n/nfold) > 0) {
    foldid = sample(c(rep(1:nfold, n/nfold), 1:(n - nfold * floor(n/nfold))))
  } else {
    foldid = sample((rep(1:nfold, n/nfold)))
  }
  return(foldid)
}
