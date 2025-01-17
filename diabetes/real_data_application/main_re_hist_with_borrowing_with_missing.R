wk.dir <- "C:/Users/wuqio/Dropbox/master_thesis/"

# load libraries
library(rstan)
library(RBesT)
library(mclust)

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
data.ma.complete$z[is.na(data.ma.complete$z)] <- 99999L
data.ma.complete$m[is.na(data.ma.complete$m)] <- 99999L
nstudy <- length(unique(data.ma.complete$trialn))

# load fitting of historial data for MAP
load(file = paste0(wk.dir, "results/hist_vague_refined_long.RData"))

# extract posterior samples and fit approximate beta and mixture Gaussian
# extract samples for q and fit beta distribution

q.sample <- extract(fit.hist.vague, pars = 'q0')$q0
mix.beta <- mixfit(as.vector(q.sample), type = "beta", Nc = 2, constrain_gt1=FALSE)
beta.comp1 <- mix.beta[["comp1"]]
beta.comp2 <- mix.beta[["comp2"]]

mix.beta.weight <- c(beta.comp1[1, 1], beta.comp2[1, 1])
mix.beta.alpha <- c(beta.comp1[2, 1], beta.comp2[2, 1])
mix.beta.beta <- c(beta.comp1[3, 1], beta.comp2[3, 1])

# extract samples for muL and sigmaL and fit mixture Gaussian with 5 components
# sigmaL needs to transform to log-scale before approximation

musigL <- extract(fit.hist.vague, pars = c('muL', 'sigmaL'))
musigL <- Reduce(cbind, musigL)
musigL[, 2] <- log(musigL[, 2])

musigL.norm.approx <- Mclust(musigL, modelNames=c("VVV"), G = 5)

mix.prob.L <- musigL.norm.approx$parameters$pro
mix.mean.L <- t(musigL.norm.approx$parameters$mean)
mix.cov.L <- musigL.norm.approx$parameters$variance$sigma
mix.cov.L <- aperm(mix.cov.L, c(3,1,2))


# extract samples for muM0 and sigmaM0 and fit mixture Gaussian with 5 components
# sigmaM0 needs to transform to log-scale before approximation

musigM <- extract(fit.hist.vague, pars = c('muM0', 'sigmaM0'))
musigM <- Reduce(cbind, musigM)
musigM[, 2] <- log(musigM[, 2])

musigM.norm.approx <- Mclust(musigM, modelNames=c("VVV"), G = 5)

mix.prob.M <- musigM.norm.approx$parameters$pro
mix.mean.M <- t(musigM.norm.approx$parameters$mean)
mix.cov.M <- musigM.norm.approx$parameters$variance$sigma
mix.cov.M <- aperm(mix.cov.M, c(3,1,2))

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
                 covL = mix.cov.L,
                 wM = mix.prob.M,
                 meansM = mix.mean.M,
                 covM = mix.cov.M)
  
  # fixed effect model
  start <- Sys.time()
  fit.hist.borrow.037 <- stan(file = paste0(code.dir, "re_hist_borrow_037_refined.stan"), 
                              data = madata,
                              chains = 5,
                              iter = 60000,
                              warmup = 10000,
                              thin = 5,
                              algorithm = 'NUTS',
                              control = list(adapt_delta = 0.99))
  run.time <- Sys.time() - start
  
  save(fit.hist.borrow.037, run.time,
       file = paste0(wk.dir, "results/re_hist_borrow_037_refined_", 100*w.map, ".RData"))
  
  
  # fixed effect model
  start <- Sys.time()
  fit.hist.borrow.250 <- stan(file = paste0(code.dir, "re_hist_borrow_250_refined.stan"), 
                              data = madata,
                              chains = 5,
                              iter = 60000,
                              warmup = 10000,
                              thin = 5,
                              algorithm = 'NUTS',
                              control = list(adapt_delta = 0.99))
  run.time <- Sys.time() - start
  
  save(fit.hist.borrow.250, run.time,
       file = paste0(wk.dir, "results/re_hist_borrow_250_refined_", 100*w.map, ".RData"))
  
}