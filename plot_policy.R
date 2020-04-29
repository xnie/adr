rm(list=ls())
library(data.table)
library(scales)
library(ggplot2)
library(ggpubr)
library(latex2exp)
get_scatterpoints = function(Z, A_time) {
   regmat_t_all = lapply(1:H, function(t) {
   num_covariates = dim(Z)[3]
     relevant = A_time >= t
     if (all(relevant==FALSE)) {
       NULL
     } else {
       f = data.frame(
         z=Z[relevant, t, 1:num_covariates, drop=F],
         t=t,
         action=(A_time==t)[relevant])
       colnames(f)[1:num_covariates] = sapply(1:num_covariates, function(i) paste0("z",i))
       f
     }
   })
  Reduce(rbind,
         regmat_t_all[!vapply(regmat_t_all, is.null, logical(1))])

}

plot_policy_all = function(adr_eval_Z,
                           adr_eval_A_time,
                           fq_Z,
                           fq_A_time,
                           H,
                           t_limit,
                           pp,
                           wadr_pp,
                           wipw_pp,
                           fnm) {
  Z = adr_eval_Z
  A_time = adr_eval_A_time
  n = dim(Z)[1]
  scatterpoints = get_scatterpoints(Z, A_time)
  if (length(unique(scatterpoints$action)) > 1) {
    scatterpoints$treatment = factor(scatterpoints$action, labels=c("on", "off"))
  } else {
    scatterpoints$treatment = factor(scatterpoints$action, labels=c("on"))
  }


   d_lines <- data.frame(int = c(wadr_pp[[3]]/wadr_pp[[1]], pp[[3]]/pp[[1]]),
                        sl = c(wadr_pp[[2]]/wadr_pp[[1]], pp[[2]]/pp[[1]]),
                        lty = c("ADR", "oracle"))

  # extract the linetypes from the data.frame (with names to map the values correclty)
  ltys <- as.character(d_lines$lty)
  names(ltys) <- as.character(d_lines$lty)

  off_color = hue_pal()(4)[1]
  on_color = hue_pal()(4)[3]
  plot_base = ggplot(scatterpoints, aes(x=t, y=z1, shape=treatment, color=treatment, size=treatment)) +
        geom_point() +
        xlim(1, t_limit) +
        ylab(TeX('S_t')) +
        theme_grey(base_size = 22) +
        scale_shape_manual(values=c(17,4)) +
        scale_size_manual(values=c(3,4))+
        scale_color_manual(values=c(on_color, off_color))
 plot_base = plot_base + geom_abline(data = d_lines,
              aes(intercept = int, slope = sl, lty = lty ), size=1.5) +
 scale_linetype_manual(name = "policy", values = c("ADR"="solid","oracle"="dotdash","IPW"="dotted")) +
 ggtitle("ADR") +
 theme(plot.title = element_text(hjust = 0.5), plot.margin=unit(c(1,2.8,1,1),"cm"))

  fq_scatterpoints = get_scatterpoints(fq_Z, fq_A_time)

  if (length(unique(fq_scatterpoints$action)) > 1) {
    fq_scatterpoints$treatment = factor(fq_scatterpoints$action, labels=c("on", "off"))
  } else {
    fq_scatterpoints$treatment = factor(fq_scatterpoints$action, labels=c("on"))
  }

  fq_plot= ggplot(fq_scatterpoints, aes(x=t, y=z1, shape=treatment, color=treatment, size=treatment)) +
        geom_point() +
        xlim(1, t_limit) +
        ylab(TeX('S_t')) +
        theme_grey(base_size = 22) +
        scale_shape_manual(values=c(17,4)) +
        scale_size_manual(values=c(3,4))+
        scale_color_manual(values=c(on_color, off_color)) +
        theme(legend.position="none") +
        ggtitle("Q-Opt") +
        theme(plot.title = element_text(hjust = 0.5), plot.margin=unit(c(1,1,1,2.8),"cm"))

  ggarrange(plot_base, fq_plot, ncol=2, nrow=1, common.legend = TRUE, legend="right")
  ggsave(paste0(fnm, ".pdf"), width = 16, height = 9, dpi = 320)
}

H=10
t_limit=5
#load("plots/plot-1000-0.5-5-1-cmp-10.RData", env=env1)
for (n in c(5000)) {
  for (obsnoise in c(0.5)) {
  	for (beta in c(5)){
  	    for (sigma in c(1)){
    			for (rep in c(1:10)) {
  			      filename = paste0("plots/plot-", n, "-", obsnoise, "-", beta, "-", sigma, "-cmp-", rep, ".RData")
              env1 <- new.env()
  			      load(filename, env=env1)
  			      combined_fnm = paste0("plots/plot-", n, "-", obsnoise, "-", beta, "-", sigma, "-", rep)
              plot_policy_all(env1$adr_eval$Z,
                              env1$adr_eval$policy$month,
                              env1$q_opt_eval_res$Z,
                              env1$q_opt_eval_res$A_time,
                              H,
                              t_limit,
                              env1$oracle_pp,
                              env1$wadr_pp,
                              env1$wipw_pp,
                              combined_fnm)
          }
  	    }
  	}
  }
}
