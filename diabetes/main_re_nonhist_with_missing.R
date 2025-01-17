wk.dir <- "C:/Users/wuqio/Dropbox/master_thesis/"

# set up stan
library(rstan)

data.dir <- paste0(wk.dir, "reference/main paper/")
code.dir <- paste0(wk.dir, "code/")
# load data and use only complete case for main meta-analysis
data.paper <- read.csv(paste0(data.dir, "diabetes_data.csv"), sep = ";", na.strings = ".")

################################################################
## data pre-processing based on supplemented SAS codes #########
################################################################

# convert planned trial duration from weeks to years
# assume 7 days in a week and 365.2425 days in a year
data.paper$tau <- data.paper$tau * 7 / 365.2425
# data.paper$logtau <- log(data.paper$tau)

#################################################
## test run on ma data set with missing values ##
#################################################
data.ma.complete <- subset(data.paper, historical == 0)
data.ma.complete$z[is.na(data.ma.complete$z)] <- 99999L
data.ma.complete$m[is.na(data.ma.complete$m)] <- 99999L
nstudy <- length(unique(data.ma.complete$trialn))

# data list for rstan
madata <- list(I = nstudy,
               NRec = length(data.ma.complete$trialn),
               i = data.ma.complete$trialn,
               j = data.ma.complete$test,
               n = data.ma.complete$n,
               y = data.ma.complete$y,
               m = data.ma.complete$m,
               z = data.ma.complete$z,
               tau = data.ma.complete$tau)

# fixed effect model
start <- Sys.time()
fit.re.nonhist.037.vague <- stan(file = paste0(code.dir, "re_nonhist_037_vague_with_missing.stan"), 
                                 data = madata,
                                 chains = 5,
                                 iter = 60000,
                                 warmup = 10000,
                                 thin = 5,
                                 algorithm = 'NUTS',
                                 control = list(adapt_delta = 0.99))
run.time <- Sys.time() - start

save(fit.re.nonhist.037.vague, run.time, file = paste0(wk.dir, "results/re_nonhist_037_vague_with_missing.RData"))

# fixed effect model
start <- Sys.time()
fit.re.nonhist.250.vague <- stan(file = paste0(code.dir, "re_nonhist_250_vague_with_missing.stan"), 
                                  data = madata,
                                  chains = 5,
                                  iter = 60000,
                                  warmup = 10000,
                                  thin = 5,
                                  algorithm = 'NUTS',
                                  control = list(adapt_delta = 0.99))
run.time <- Sys.time() - start

save(fit.re.nonhist.250.vague, run.time, file = paste0(wk.dir, "results/re_nonhist_250_vague_with_missing.RData"))

# fixed effect model
start <- Sys.time()
fit.re.nonhist.037.weak <- stan(file = paste0(code.dir, "re_nonhist_037_weak_with_missing.stan"), 
                                  data = madata,
                                  chains = 5,
                                  iter = 60000,
                                  warmup = 10000,
                                  thin = 5,
                                  algorithm = 'NUTS',
                                  control = list(adapt_delta = 0.99))
run.time <- Sys.time() - start

save(fit.re.nonhist.037.weak, run.time, file = paste0(wk.dir, "results/re_nonhist_037_weak_with_missing.RData"))

# fixed effect model
start <- Sys.time()
fit.re.nonhist.250.weak <- stan(file = paste0(code.dir, "re_nonhist_250_weak_with_missing.stan"), 
                                  data = madata,
                                  chains = 5,
                                  iter = 60000,
                                  warmup = 10000,
                                  thin = 5,
                                  algorithm = 'NUTS',
                                  control = list(adapt_delta = 0.99))
run.time <- Sys.time() - start

save(fit.re.nonhist.250.weak, run.time, file = paste0(wk.dir, "results/re_nonhist_250_weak_with_missing.RData"))
