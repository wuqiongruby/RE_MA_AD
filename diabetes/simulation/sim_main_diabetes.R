# read-in array id ------------------------
# array id should within 0-7999
args <- commandArgs(trailingOnly = TRUE)
array.id <- as.numeric(args[1])

# set up directories --------------------------
wk.dir <- "~/"
code.dir <- paste0(wk.dir, "code/")
res.dir <- "/scratch1/projects/mthesis_qiong_wu/"

# load R packages ----------------------------------
library(rstan)
library(mclust)
library(RBesT)

# generate input variables for data generation ------------
num.sim <- 1000
scenario <- floor(array.id / num.sim) + 1
sim.id <- array.id %% num.sim
seed <- 12345 * sim.id + array.id

# number of trials ----------------------------------
num.hist.trials <- 12
num.ma.trials <- 6
num.trials.total <- num.hist.trials + num.ma.trials

# true parameters --------------------------------------
if(scenario %in% c(1,2,3,7)){
  trueLHR <- 0 # phi: true log-HR for an event of interest
}else{
  trueLHR <- 0.25 # phi: true log-HR for an event of interest
}

if(scenario %in% c(1, 4)){
  sdLHR <- 0.25*sqrt(2)/sqrt(pi) # mean of half-norm(0, 0.25)
}else if(scenario %in% c(2, 5)){
  sdLHR <- 0.5*sqrt(2)/sqrt(pi) # mean of half-norm(0, 0.5)
}else if(scenario %in% c(3, 6)){
  sdLHR <- sqrt(2)/sqrt(pi) # mean of half-norm(0, 1)
}else{
  sdLHR <- 0 # fixed effect model
}
q <- 0.35 # q_0 and q_1: prob of fatal event
evtrate <- 0.5 # lambda_0: event time Weibull rate parameter for control group
censrate <- 0.5 # mu_0: drop-out Weibull rate parameter for control group
censLHR <- 0 # true log-HR for drop-out
censshape <- 1 # drop-out Weibull shape parameter for control group; value 1 means exponential
evtshape <- 1 # event time Weibull shape parameter for control group; value 1 means exponential

# trial durations ---------------------------------------------
tau <- c(0.25, 0.25, 0.5, 0.5, 0.5, 1,
         0.5, 0.25, 0.25, 1, 1, 1,
         0.5, 0.25, 0.25, 1, 1, 1)
# number of patients in each arm of each trial -----------------
num.patients <- c(50, 300, 100, 100, 25, 50,
                  100, 200, 100, 200, 100, 100,
                  100, 100, 50, 300, 300, 300,
                  100, 100, 50, 300, 300, 300)
trial.no <- c(rep(1:num.ma.trials, each = 2), (num.ma.trials+1):num.trials.total)

# function to generate aggregated patient-level data for a single simulation
gen_patient_one_trial <- function(trialn){
    tau.trial <- tau[trialn]
    hist.trial <- (trialn > num.ma.trials)
    trueLHR.trial <- rnorm(1, mean = trueLHR, sd = sdLHR)
    
    patient.data <- NULL # simulated aggregate data for the trial
    for (test in 0:(1-hist.trial)) {  # test indicates arm
      n.patient <- num.patients[trial.no == trialn] # n_i
      n.patient <- n.patient[test+1] # n_ij
      y.total <- z.total <- m.total <- 0 # aggregate data
      for(pat in seq(n.patient)){
        t <- rweibull(1, shape = censshape, scale = 1/(censrate * exp(censLHR * test)))
        x <- rweibull(1, shape = evtshape, scale = 1/(evtrate * exp(trueLHR.trial * test)))
        y <- x<=min(t, tau.trial)
        # check fatal event
        if(y){
          m <- runif(1, min = 0, max = 1) < q
        }else{
          m <- FALSE
        }
        z <- m || (t<tau.trial)
        
        # update aggregate data
        y.total <- y.total + y
        z.total <- z.total + z
        m.total <- m.total + m
      }
      patient.data <- rbind(patient.data,
                            c(trialn, hist.trial, tau.trial, test, n.patient, y.total, z.total, m.total, trueLHR.trial))
    }
  colnames(patient.data) <- c("trialn", "historical", "tau", "test", "n", "y", "z", "m", "trueLHR")
  patient.data
}

#######################################################
######### simulation starts here ! ####################
#######################################################

# start timing
start <- Sys.time()

# set seed for reproducibility
set.seed(seed)

##### generate aggregated data for single simulation ####
data.sim <- sapply(seq(num.trials.total), gen_patient_one_trial)
data.sim <- as.data.frame(Reduce(rbind, data.sim))
  
################# fit historical data ###################
data.hist.complete <- subset(data.sim, historical == 1)

# data input for fitting on historical data to obtain MAP
histdata <- list(NRec = nrow(data.hist.complete),
                 n = data.hist.complete$n,
                 y = data.hist.complete$y,
                 m = data.hist.complete$m,
                 z = data.hist.complete$z,
                 tau = data.hist.complete$tau)
  
# fixed effect model on historical data to obtain MAP
fit.hist.vague <- stan(file = paste0(code.dir, "hist_vague_prior_with_missing.stan"), 
                       data = histdata,
                       chains = 4,
                       iter = 20000,
                       thin = 5,
                       algorithm = 'NUTS',
                       control = list(adapt_delta = 0.99))
  
########### approximation of MAP prior ################
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
  
rm(q.sample, musigL, musigM, musigL.norm.approx, musigM.norm.approx)

################################################
### Start to fit Bayesian hierachical models ###
################################################

############ Fix-effect models #################
#------------- Vague prior ------------
data.ma.complete <- subset(data.sim, historical == 0)
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
fit.nonhist.250.vague <- stan(file = paste0(code.dir, "replicate_nonhist_250_vague_with_missing.stan"), 
                              data = madata,
                              chains = 4,
                              iter = 20000,
                              warmup = 10000,
                              thin = 5,
                              algorithm = 'NUTS',
                              control = list(adapt_delta = 0.99))

### TBD: what summary of outputs should be saved ###
  
#------------- MAP with borrowing -------------
# data input for main meta-analysis
w.map <- 0.5
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
  
fit.hist.borrow.250 <- stan(file = paste0(code.dir, "replicate_hist_borrow_250_refined.stan"), 
                            data = madata,
                            chains = 4,
                            iter = 20000,
                            warmup = 10000,
                            thin = 5,
                            algorithm = 'NUTS',
                            control = list(adapt_delta = 0.99))

### TBD: what summary of outputs should be saved ###
  
############ Random-effect models ################
#------------- Vague prior ------------
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
  
# random effect model
fit.re.nonhist.250.vague <- stan(file = paste0(code.dir, "re_nonhist_250_vague_with_missing.stan"), 
                                 data = madata,
                                 chains = 4,
                                 iter = 20000,
                                 warmup = 10000,
                                 thin = 5,
                                 algorithm = 'NUTS',
                                 control = list(adapt_delta = 0.99))

### TBD: what summary of outputs should be saved ###
  
#------------- MAP with borrowing -------------
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
  
# random effect model
fit.re.hist.borrow.250 <- stan(file = paste0(code.dir, "re_hist_borrow_250_refined.stan"), 
                               data = madata,
                               chains = 4,
                               iter = 20000,
                               warmup = 10000,
                               thin = 5,
                               algorithm = 'NUTS',
                               control = list(adapt_delta = 0.99))

### TBD: what summary of outputs should be saved ###

# stop timing
run.time <- Sys.time() - start

################### outputs #####################
res <- list(scenario = scenario,
            seed = seed,
            sim_id = sim.id,
            data = data.sim,
            fit_hist = fit.hist.vague,
            fit_fe_nonhist = fit.nonhist.250.vague,
            fit_fe_hist = fit.hist.borrow.250,
            fit_re_nonhist = fit.re.nonhist.250.vague,
            fit_re_hist = fit.re.hist.borrow.250,
            run_time = run.time)
save(res, file = paste0(res.dir, "sim_res_sce_", scenario, "_sim_", sim.id, ".RData"))
