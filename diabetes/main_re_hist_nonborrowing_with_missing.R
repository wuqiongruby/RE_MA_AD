wk.dir <- "C:/Users/wuqio/Dropbox/master_thesis/"

# load libraries
library(rstan)
library(RBesT)
library(teigen)

# load data set
data.dir <- paste0(wk.dir, "reference/main paper/")
code.dir <- paste0(wk.dir, "code/")
# load data and use only complete case for main meta-analysis
data.paper <- read.csv(paste0(data.dir, "diabetes_data.csv"), sep = ";", na.strings = ".")
# convert planned trial duration from weeks to years
# assume 7 days in a week and 365.2425 days in a year
data.paper$tau <- data.paper$tau * 7 / 365.2425

# data set for main meta-analysis
data.ma.complete <- subset(data.paper, historical == 0)
# data.ma.complete <- data.ma.complete[complete.cases(data.ma.complete), ]
data.ma.complete$z[is.na(data.ma.complete$z)] <- 99999L
data.ma.complete$m[is.na(data.ma.complete$m)] <- 99999L
nstudy <- length(unique(data.ma.complete$trialn))
data.ma.complete$trialn <- rep(seq(nstudy), each = 2)

# load fitting of historial data for MAP
load(file = paste0(wk.dir, "results/with_missing_data/hist_vague_refined_long.RData"))

# extract posterior samples and fit approximate beta and mixture Gaussian
# extract samples for q and fit beta distribution

q.sample <- extract(fit.hist.vague, pars = 'q0')$q0
mix.beta <- mixfit(as.vector(q.sample), type = "beta", Nc = 2, constrain_gt1=FALSE)
beta.comp1 <- mix.beta[["comp1"]]
beta.comp2 <- mix.beta[["comp2"]]

mix.beta.weight <- c(beta.comp1[1, 1], beta.comp2[1, 1])
mix.beta.alpha <- c(beta.comp1[2, 1], beta.comp2[2, 1])
mix.beta.beta <- c(beta.comp1[3, 1], beta.comp2[3, 1])

# extract samples for loglambda and fit mixture t_3 with 3 components
loglambda <- as.vector(extract(fit.hist.vague, pars = 'lambda_sim')$lambda_sim)
fit.loglambda <- teigen(loglambda, Gs = 3, models = "univUC",
                        dfstart = 3, dfupdate = FALSE, scale = FALSE)

mix.prob.L <- fit.loglambda$parameters$pig
mix.mean.L <- fit.loglambda$parameters$mean[,1]
mix.sd.L <- sqrt(Reduce(c, fit.loglambda$parameters$sigma))

# extract samples for logmu and fit mixture t_3 with 3 components
logmu <- as.vector(extract(fit.hist.vague, pars = 'mu_sim')$mu_sim)
fit.logmu <- teigen(logmu, Gs = 3, models = "univUC",
                    dfstart = 3, dfupdate = FALSE, scale = FALSE)

mix.prob.M <- fit.logmu$parameters$pig
mix.mean.M <- fit.logmu$parameters$mean[,1]
mix.sd.M <- sqrt(Reduce(c, fit.logmu$parameters$sigma))

rm(fit.hist.vague, fit.loglambda, fit.logmu, q.sample, loglambda, logmu)
# set the rstan data and run stan file for model fitting
for (w.map in c(1, 0.8, 0.5)) {
  # data input for main meta-analysis
  madata <- list(I = nstudy,
                 NRec = length(data.ma.complete$trialn),
                 i = data.ma.complete$trialn,
                 j = data.ma.complete$test,
                 n = data.ma.complete$n,
                 y = data.ma.complete$y,
                 m = data.ma.complete$m,
                 z = data.ma.complete$z,
                 tau = data.ma.complete$tau,
                 mw = w.map,
                 wbeta = mix.beta.weight,
                 alphaQ = mix.beta.alpha,
                 betaQ = mix.beta.beta,
                 wL = mix.prob.L,
                 meansL = mix.mean.L,
                 sdL = mix.sd.L,
                 wM = mix.prob.M,
                 meansM = mix.mean.M,
                 sdM = mix.sd.M)
  
  # fixed effect model
  start <- Sys.time()
  fit.hist.nonborrow.037.re <- stan(file = paste0(code.dir, "re_hist_nonborrow_037_t3_map_with_missing.stan"), 
                                 data = madata,
                                 chains = 5,
                                 iter = 60000,
                                 warmup = 10000,
                                 thin = 5,
                                 algorithm = 'NUTS',
                                 control = list(adapt_delta = 0.99))
  run.time <- Sys.time() - start
  
  save(fit.hist.nonborrow.037.re, run.time,
       file = paste0(wk.dir, "results/re_hist_nonborrow_037_refined_", 100*w.map, "_v3.RData"))
  rm(fit.hist.nonborrow.037)
  
  # fixed effect model
  start <- Sys.time()
  fit.hist.nonborrow.250.re <- stan(file = paste0(code.dir, "re_hist_nonborrow_250_t3_map_with_missing.stan"), 
                              data = madata,
                              chains = 5,
                              iter = 60000,
                              warmup = 10000,
                              thin = 5,
                              algorithm = 'NUTS',
                              control = list(adapt_delta = 0.99))
  run.time <- Sys.time() - start
  
  save(fit.hist.nonborrow.250.re, run.time,
       file = paste0(wk.dir, "results/re_hist_nonborrow_250_refined_", 100*w.map, "_v3.RData"))
  rm(fit.hist.nonborrow.250)
}