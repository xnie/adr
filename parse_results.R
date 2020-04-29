rm(list = ls())

library(xtable)
library(data.table)
library(scales)

setup=2 # change for setup = 1
beta = 0.5 # change for beta = 1 if setup = 2. Ignore this for setup = 1
filenames = list.files(paste0("results/simu", setup), pattern="*combined*", full.names=TRUE)

if (setup == 1) { # multiple
  param.names = c("n", "noise")
} else if (setup == 2) { # binary
  param.names = c("n", "noise", "beta", "sigma")
}

raw = data.frame(t(sapply(filenames, function(fnm) {

  output = as.matrix(read.table(fnm, header=TRUE, sep=","))
  output = output[, colnames(output)!="X"]
  output = cbind(output, output[,10]-output[,3], output[,10]-output[,6], output[,10]-output[,11], (output[,23]- output[,10])^2, (output[,26]-output[,10])^2, (output[,28] - output[,10])^2, (output[,29] - output[,10])^2)# corresponding to wadr, wipw, fq regret, and wadr, wipw, waipw, q
  if (setup == 1) {
    params = strsplit(fnm, "-")[[1]][2:3]
  } else if (setup == 2) {
    params = strsplit(fnm, "-")[[1]][2:5]
  }

  mean = colMeans(output)
  se = apply(output, 2, sd) / sqrt(dim(output)[1])

  c(params,
    mean,
    se)
})))

var_colnames = c("adr_v", "adr_reg_v", "wadr_v","wadr_reg_v","ipw_v","wipw_v","aipw_v","waipw_v","q_v", "ora_v","fq_v","adr_mse","adr_reg_mse","wadr_mse","wadr_reg_mse","ipw_mse","wipw_mse","aipw_mse","waipw_mse", "q_mse", "adr_value","adr_reg_value", "wadr_value", "wadr_reg_value", "ipw_value", "wipw_value", "aipw_value", "waipw_value", "q_value")
rownames(raw) = 1:nrow(raw)
colnames_original = c(var_colnames, "wadr_regret", "wipw_regret", "fq_regret", "wadr_eval_mse", "wipw_eval_mse", "waipw_eval_mse", "q_eval_mse")
names(raw) = c(param.names, colnames_original, paste(colnames_original, "_se", sep=""))

options(stringsAsFactors = FALSE)
if (setup == 1) {
  raw = raw[order(as.numeric(raw$noise)),]
  raw = raw[order(as.numeric(raw$n)),]
} else if (setup == 2) {
  raw = raw[order(as.numeric(raw$sigma)),]
  raw = raw[order(as.numeric(raw$beta)),]
  raw = raw[order(as.numeric(raw$noise)),]
  raw = raw[order(as.numeric(raw$n)),]
}

rownames(raw) = 1:nrow(raw)
mincol = length(param.names) + 1
maxcol = dim(raw)[2]

for (col in (mincol:maxcol)){
  raw[,col] = formatC(as.numeric(unlist(raw[,col])), format = "e", digits = 2)
}

# write raw csv output file
write.csv(raw, file=paste0("results-simu", setup, "-all.csv"))

if (setup == 2) {
  raw<-raw[raw$beta ==beta, ] # only look at beta at a certain value to make the table smaller
}
 if (setup == 1) {
   raw = raw[,names(raw) %in% c("n", "noise", "wadr_v", "ora_v", "fq_v", "wipw_v", "wadr_mse", "wipw_mse")]
   raw <- raw[c("n", "noise", "wadr_v", "wipw_v", "fq_v", "ora_v", "wadr_mse", "wipw_mse")]
 } else if (setup == 2) {
   raw = raw[,names(raw) %in% c("n", "noise", "beta", "sigma",  "wadr_v", "ora_v", "fq_v", "wipw_v", "wadr_mse", "wipw_mse")]
   raw <- raw[c("n", "noise", "beta", "sigma",  "wadr_v", "wipw_v", "fq_v", "ora_v", "wadr_mse", "wipw_mse")]
 }
 # write to latex tables
 tab.setup = cbind("", raw)
 if (setup == 1) {
   val.idx = c(4:6)
   mse.idx = c(8:9)
 } else if (setup == 2) {
   val.idx = c(6:8)
   mse.idx = c(10:11)
 }
 for(iter in 1:nrow(tab.setup)) {
   best.idx = val.idx[which(as.numeric(tab.setup[iter,val.idx]) == max(as.numeric(tab.setup[iter,val.idx])))]
   for (j in 1:length(best.idx)) {
     best.idx.j = best.idx[j]
     tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,best.idx.j])
   }
 }
 for(iter in 1:nrow(tab.setup)) {
   best.idx = mse.idx[which(as.numeric(tab.setup[iter,mse.idx]) == min(as.numeric(tab.setup[iter,mse.idx])))]
   for (j in 1:length(best.idx)) {
     best.idx.j = best.idx[j]
     tab.setup[iter,best.idx.j] = paste("\\bf", tab.setup[iter,best.idx.j])
   }
 }

 tab.setup = tab.setup[,-1]
 tab.setup
 if (setup == 1) {
   xtab.setup = xtable(tab.setup,  align="ccccccccc")
   names(xtab.setup) = c("n", "$\\sigma$", "ADR", "IPW", "Q-Opt", "Oracle", "MSE:ADR", "MSE:IPW")
   print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste0("tables/simulation_results_setup_", setup, ".tex"))
  } else if (setup == 2) {
   xtab.setup = xtable(tab.setup,  align="ccccccccccc")
   names(xtab.setup) = c("n", "$\\nu$", "$\\beta$", "$\\sigma$", "ADR",  "IPW", "Q-Opt","Oracle", "MSE:ADR", "MSE:IPW")
   print(xtab.setup, include.rownames = FALSE, include.colnames = TRUE, sanitize.text.function = identity, file = paste0("tables/simulation_results_setup_", setup, "_beta_", beta, ".tex"))
 }
