rm(list = ls())

library(xtable)
library(data.table)
library(scales)
library(dplyr)
library(readr)


results_dir = "results/simu1"

for (n in c(250, 500, 1000,5000,10000, 20000)) {
  for (obs_noise in c(0,0.5,1)) {
#      if (obs_noise == 0.5) {
#        filenames = list.files(results_dir, pattern=paste0(paste("output", n, "0\\.5", "-", sep="-"), "*"), full.names=TRUE)
#      } else {
        filenames = list.files(results_dir, pattern=paste0(paste("simu1", n, obs_noise,  "-", sep="-"), "*"), full.names=TRUE)
#      }
      df <- t(filenames  %>%
        lapply(read_csv) %>%
        bind_cols)
      rownames(df) = c()
      write.csv(df, file = paste0(paste(paste0(results_dir, "/simu1"), n, obs_noise, sep = "-"), "-combined.csv"))
  }
}

results_dir = "results/simu2"

for (n in c(250, 500, 1000,5000,10000, 20000,30000)) {
  for (obs_noise in c(0,0.5)) {
    for (beta in c(0.5,1)) {
      for (sigma in c(1,3)) {
        filenames = Sys.glob(paste0(results_dir, paste("/simu2", n, obs_noise, beta, sigma, sep="-"), "*"))
        df <- t(filenames  %>%
          lapply(read_csv) %>%
          bind_cols)
        rownames(df) = c()
        write.csv(df, file = paste0(paste(paste0(results_dir, "/simu2"), n, obs_noise, beta, sigma, sep = "-"), "-combined.csv"))
      }
    }
  }
}
