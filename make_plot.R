rm(list=ls())
library(tidyverse)
library(RColorBrewer)
library(latex2exp)
library(grid)
library(tidyr)
library(dplyr)


pal = brewer.pal(n = 8, name = "Dark2")
plotsize = function(x,y) options(repr.plot.width=x, repr.plot.height=y)
plotsize(8,15)


make_regret_plot = function(data, fnm) {

  df = data%>%
      select(n, wadr_regret,  wipw_regret, fq_regret, wadr_regretupper, wipw_regretupper, fq_regretupper, wadr_regretlower, wipw_regretlower, fq_regretlower) %>%
        gather(v, value, wadr_regret:fq_regretlower) %>%
        separate( v,
                 into = c("method", "line"))%>%
        arrange(n) %>%
        spread(line, value)

      ggplot(df, aes(x=n,
                 y=log(regret),
                 colour=method)) +
       geom_point(size=3) +
       geom_line(size=2) +
       scale_x_log10()+
       theme_bw(base_size = 11) +
       geom_ribbon(aes(ymin=log(regretlower), ymax=log(regretupper), fill="grey"),
                       alpha=0.2) +
       guides(fill=FALSE) +
       scale_color_manual(name="method",
                          breaks=c("wadr", "wipw","fq"),
                          labels = c("ADR",
                                     "IPW",
                                     "Q-Opt"),
                          values =c(pal[3], pal[1], pal[2])) +
        scale_fill_manual(values = c(pal[3], pal[1], pal[2]))
    ggsave(fnm)
}
make_mse_plot = function(data, fnm) {
      df = data%>%
      select(n, wadr_mse,  wipw_mse, wadr_mseupper, wipw_mseupper,  wadr_mselower, wipw_mselower) %>%
        gather(v, value, wadr_mse:wipw_mselower) %>%
        separate( v,
                 into = c("method", "line"))%>%
        arrange(n) %>%
        spread(line, value)

      ggplot(df, aes(x=n,
                 y=log(mse),
                 colour=method)) +
       geom_point(size=3) +
       geom_line(size=2) +
       scale_x_log10()+
       theme_bw(base_size = 11) +
       geom_ribbon(aes(ymin=log(mselower), ymax=log(mseupper), fill=method),
                       alpha=0.2) +
       guides(fill=FALSE) +
       scale_color_manual(name="method",
                          breaks=c("wadr", "wipw"),
                          labels = c("ADR",
                                     "IPW"),
                          values =c(pal[1], pal[2]))+
        scale_fill_manual(values = c(pal[1], pal[2]))

    ggsave(fnm)

}

make_eval_mse_plot = function(data, fnm) {

      df = data%>%
      select(n, wadr_evalmse,  wipw_evalmse, waipw_evalmse, q_evalmse, wadr_evalmseupper, wipw_evalmseupper,  wadr_evalmselower, wipw_evalmselower, waipw_evalmseupper, waipw_evalmselower, q_evalmseupper, q_evalmselower) %>%
        gather(v, value, wadr_evalmse:q_evalmselower) %>%
        separate( v,
                 into = c("method", "line"))%>%
        arrange(n) %>%
        spread(line, value)

      ggplot(df, aes(x=n,
                 y=log(evalmse),
                 colour=method)) +
       geom_point(size=3) +
       geom_line(size=2) +
        ylab("log(mse)")+
       scale_x_log10()+
       theme_bw(base_size = 11) +
       geom_ribbon(aes(ymin=log(evalmselower), ymax=log(evalmseupper), fill=method),
                       alpha=0.2) +
       guides(fill=FALSE) +
       scale_color_manual(name="method",
                          breaks=c("wadr", "wipw", "waipw", "q"),
                          labels = c("ADR",
                                     "IPW",
                                     "AIPW",
                                     "Q-Eval"),
                          values =c(pal[3], pal[1], pal[4], pal[2]))+
        scale_fill_manual(values = c(pal[3], pal[1], pal[4], pal[2]))

    ggsave(fnm)
}

make_policy_eval_plot = function(data, fnm) {
  oracle_value = mean(data$ora_v)
      df = data%>%
      select(n, wadr_value, wipw_value, waipw_value, q_value, wadr_valueupper, wipw_valueupper, waipw_valueupper, q_valueupper, wadr_valuelower, wipw_valuelower, waipw_valuelower, q_valuelower) %>%
        gather(v, value, wadr_value:q_valuelower) %>%
        separate( v,
                 into = c("method", "line"))%>%
        arrange(n) %>%
        spread(line, value)

      ggplot(df, aes(x=n,
                 y=value,
                 colour=method)) +
       geom_point(size=3) +
       geom_hline(yintercept=oracle_value,
                     color = "black", size=2) +
       geom_line(size=2) +
       scale_x_log10()+
       theme_bw(base_size = 11) +
       geom_ribbon(aes(ymin=valuelower, ymax=valueupper, fill=method),
                       alpha=0.2) +
       guides(fill=FALSE) +
       scale_color_manual(name="method",
                          breaks=c("wadr", "wipw", "waipw", "q"),
                          labels = c("ADR",
                                     "IPW",
                                     "AIPW",
                                     "Q-Eval"),
                          values =c(pal[3], pal[1], pal[4], pal[2]))+
        scale_fill_manual(values = c(pal[3], pal[1], pal[4], pal[2]))

    ggsave(fnm)

}

alpha = 1
for (setup in c(1:2)) {
  out = read.csv(paste0("results-simu", setup, "-all.csv"))
  out$wadr_regretupper = out$wadr_regret + alpha*out$wadr_regret_se
  out$wadr_regretlower = out$wadr_regret - alpha*out$wadr_regret_se

  out$wipw_regretupper = out$wipw_regret + alpha*out$wipw_regret_se
  out$wipw_regretlower = out$wipw_regret - alpha*out$wipw_regret_se

  out$fq_regretupper = out$fq_regret + alpha*out$fq_regret_se
  out$fq_regretlower = out$fq_regret - alpha*out$fq_regret_se

  out$wadr_mseupper = out$wadr_mse + alpha*out$wadr_mse_se
  out$wadr_mselower = out$wadr_mse - alpha*out$wadr_mse_se

  out$wipw_mseupper = out$wipw_mse + alpha*out$wipw_mse_se
  out$wipw_mselower = out$wipw_mse - alpha*out$wipw_mse_se

  out$wadr_valueupper = out$wadr_value + alpha*out$wadr_value_se
  out$wadr_valuelower= out$wadr_value - alpha*out$wadr_value_se

  out$wipw_valueupper = out$wipw_value + alpha*out$wipw_value_se
  out$wipw_valuelower= out$wipw_value - alpha*out$wipw_value_se

  out$waipw_valueupper = out$waipw_value + alpha*out$waipw_value_se
  out$waipw_valuelower= out$waipw_value - alpha*out$waipw_value_se

  out$q_valueupper = out$q_value + alpha*out$q_value_se
  out$q_valuelower= out$q_value - alpha*out$q_value_se

  out$wadr_evalmseupper = out$wadr_evalmse + alpha*out$wadr_evalmse_se
  out$wadr_evalmselower= out$wadr_evalmse - alpha*out$wadr_evalmse_se

  out$wipw_evalmseupper = out$wipw_evalmse + alpha*out$wipw_evalmse_se
  out$wipw_evalmselower= out$wipw_evalmse - alpha*out$wipw_evalmse_se

  out$waipw_evalmseupper = out$waipw_evalmse + alpha*out$waipw_evalmse_se
  out$waipw_evalmselower= out$waipw_evalmse - alpha*out$waipw_evalmse_se

  out$q_evalmseupper = out$q_evalmse + alpha*out$q_evalmse_se
  out$q_evalmselower= out$q_evalmse - alpha*out$q_evalmse_se

   if (setup == 2) {
    for (noise in unique(out$noise)){
      for (beta in unique(out$beta)) {
        for (sigma in unique(out$sigma)) {

        data = out[out$beta==beta & out$sigma==sigma & out$noise ==noise,]

        fnm = paste("setup", setup, "noise", noise, "beta", beta, "sigma", sigma, sep="-")
        fnm = str_replace(fnm, "\\.", "_")
        make_regret_plot(data, paste0("plots/regret-", fnm, ".pdf"))
        make_mse_plot(data, paste0("plots/mse-", fnm, ".pdf"))
        make_policy_eval_plot(data, paste0("plots/eval-", fnm, ".pdf"))
        make_eval_mse_plot(data, paste0("plots/evalmse-", fnm, ".pdf"))
        }
      }
    }
  }
  if (setup == 1) {
      for (noise in unique(out$noise)) {

        data = out[out$noise==noise,]

        fnm = paste("setup", setup, "noise", noise, sep="-")
        fnm = str_replace(fnm, "\\.", "_")
        make_regret_plot(data, paste0("plots/regret-", fnm, ".pdf"))
        make_mse_plot(data, paste0("plots/mse-", fnm, ".pdf"))
        make_policy_eval_plot(data, paste0("plots/eval-", fnm, ".pdf"))
        make_eval_mse_plot(data, paste0("plots/evalmse-", fnm, ".pdf"))
      }
  }
}


