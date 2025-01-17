wk.dir <- "C:/Users/wuqio/Dropbox/master_thesis/upload/code/oncology/"
code.dir <- paste0(wk.dir, "code/")

# load R packages
library(rstan)

# load data and use only complete case for main meta-analysis
data.anja <- readxl::read_xlsx(paste0(wk.dir, "oncology_data.xlsx"))
data.anja <- data.anja[, -4]
colnames(data.anja) <- c("trialn", "historical", "test", "n", "y", "z", "m", "tau")

# convert planned trial duration from months to years
# assume 12 month in a year
data.anja$tau <- data.anja$tau / 12

# each trial only retain two records
data.anja <- data.anja[-c(14, 17), ]   # can be further decided 

# reformat arm code
data.anja$test <- rep(c(1,0), max(data.anja$trialn))

# convert test, n, y, z, m into integer type
data.anja[, 2:7] <- type.convert(data.anja[, 2:7])

####################################################################
## model fitting
nstudy <- length(unique(data.anja$trialn))

# data list for rstan
madata <- list(I = nstudy,
               NRec = nrow(data.anja),
               i = data.anja$trialn,
               j = data.anja$test,
               n = data.anja$n,
               y = data.anja$y,
               m = data.anja$m,
               z = data.anja$z,
               tau = data.anja$tau)

# random effect model with vague prior
start <- Sys.time()
fit.re.nonhist.250.vague <- stan(file = paste0(code.dir, "re_nonhist_250_vague.stan"), 
                                 data = madata,
                                 chains = 5,
                                 iter = 20000,
                                 warmup = 10000,
                                 thin = 5,
                                 algorithm = 'NUTS',
                                 control = list(adapt_delta = 0.99))
run.time <- Sys.time() - start

save(fit.re.nonhist.250.vague, run.time, file = paste0(wk.dir, "results/re_nonhist_250_vague.RData"))
