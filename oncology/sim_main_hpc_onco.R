# read-in array id ------------------------
# array id should within 0-1999
args <- commandArgs(trailingOnly = TRUE)
array.id <- as.numeric(args[1])

# set up directories --------------------------
wk.dir <- "~/"
code.dir <- paste0(wk.dir, "code/")
res.dir <- "/scratch1/projects/mthesis_qiong_wu/"

# load R packages ----------------------------------
library(rstan)

# generate input variables for data generation ------------
num.sim <- 1000
scenario <- floor(array.id / num.sim) + 1
sim.id <- array.id %% num.sim
seed <- 12345 * sim.id + array.id

# number of trials ----------------------------------
num.hist.trials <- 0
num.ma.trials <- 9
num.trials.total <- num.hist.trials + num.ma.trials

# true parameters --------------------------------------
if(scenario %in% c(1,2,3,7)){
  trueLHR <- 0 # phi: true log-HR for an event of interest
}else{
  trueLHR <- 0.5 # phi: true log-HR for an event of interest
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
q <- 0.01 # q_0 and q_1: prob of fatal event
evtrate <- 0.02 # lambda_1: event time Weibull rate parameter for control group
logevtrate_sd <- 1.2 # sd of loglambda_1
censrate <- 0.5 # mu_0: drop-out Weibull rate parameter for control group
censLHR <- 0 # true log-HR for drop-out
censshape <- 1 # drop-out Weibull shape parameter for control group; value 1 means exponential
evtshape <- 1 # event time Weibull shape parameter for control group; value 1 means exponential

# trial durations ---------------------------------------------
tau <- c(30, 20, 65, 65, 30, 30,
         15, 15, 35, 35, 25, 20,
         30, 25, 25, 25, 10, 10) / 12
# number of patients in each arm of each trial -----------------
num.patients <- c(400, 370, 680, 500, 250, 240,
                  190, 180, 350, 350, 440, 440,
                  190, 180, 330, 340, 350, 350)
# trial.no <- c(rep(1:num.ma.trials, each = 2), (num.ma.trials+1):num.trials.total)
trial.no <- rep(1:num.ma.trials, each = 2)

# function to generate aggregated patient-level data for a single simulation
gen_patient_one_trial <- function(trialn){
  hist.trial <- (trialn > num.ma.trials)
  trueLHR.trial <- rnorm(1, mean = trueLHR, sd = sdLHR)
  trueevtrate.trial <- exp(rnorm(1, mean = log(evtrate), sd = logevtrate_sd))
  
  patient.data <- NULL # simulated aggregate data for the trial
  for (test in (1-hist.trial):0) {  # test indicates arm
    tau.trial <- tau[trial.no == trialn]
    tau.trial <- tau.trial[2-test] # tau_ij
    n.patient <- num.patients[trial.no == trialn] # n_i
    n.patient <- n.patient[2-test] # n_ij
    y.total <- z.total <- m.total <- 0 # aggregate data
    for(pat in seq(n.patient)){
      t <- rweibull(1, shape = censshape, scale = 1/(censrate * exp(censLHR * test)))
      x <- rweibull(1, shape = evtshape, scale = 1/(trueevtrate.trial / exp(trueLHR.trial * (1-test))))
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
                          c(trialn, hist.trial, tau.trial, test, n.patient,
                            y.total, z.total, m.total, trueLHR.trial, trueevtrate.trial))
  }
  colnames(patient.data) <- c("trialn", "historical", "tau", "test", "n", "y", "z", "m", "trueLHR", "trueevtrate")
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
data.sim <- lapply(seq(num.trials.total), gen_patient_one_trial)
data.sim <- as.data.frame(Reduce(rbind, data.sim))

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
fit.nonhist.250.vague <- stan(file = paste0(code.dir, "fe_nonhist_250_vague.stan"), 
                              data = madata,
                              chains = 4,
                              iter = 20000,
                              warmup = 10000,
                              thin = 5,
                              algorithm = 'NUTS',
                              control = list(adapt_delta = 0.99))

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
fit.re.nonhist.250.vague <- stan(file = paste0(code.dir, "re_nonhist_250_vague.stan"), 
                                 data = madata,
                                 chains = 4,
                                 iter = 20000,
                                 warmup = 10000,
                                 thin = 5,
                                 algorithm = 'NUTS',
                                 control = list(adapt_delta = 0.99))

# stop timing
run.time <- Sys.time() - start

################### outputs #####################
res <- list(scenario = scenario,
            seed = seed,
            sim_id = sim.id,
            data = data.sim,
            fit_fe_nonhist = fit.nonhist.250.vague,
            fit_re_nonhist = fit.re.nonhist.250.vague,
            run_time = run.time)
save(res, file = paste0(res.dir, "sim_res_sce_", scenario, "_sim_", sim.id, ".RData"))
