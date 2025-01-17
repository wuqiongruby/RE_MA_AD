// Model with robust historical prior and borrowing of information
functions {
  real logL_mymodel(int n, // num of patients
                    int y, // num of events
                    int m, // num of fatal events
                    int z, // num of discontinuation
                    // int R1, // lower bound of w_3 event
                    // int R2, // upper bound of w_3 event
                    vector logpij){
           int ymin;
           int ymax;
           int mmin;
           int mmax;
           int zmin;
           int zmax;
           
           // range for y values
           if (y != 99999){
             ymin = y;
             ymax = y;
           }else{
             if (m != 99999){
               ymin = m;
             }else{
               ymin = 0;
             }
             ymax = n;
           }
           
           // range for m values
           if (m != 99999){
             mmin = m;
             mmax = m;
           }else{
             mmin = 0;
             mmax = ymax;
           }
           
           // range for z values
           if(z != 99999){
             zmin = z;
             zmax = z;
           }else{
             if(m != 99999){
               zmin = m;
             }else{
               zmin = 0;
             }
             zmax = n;
           }
           
           vector[(ymax-ymin+1)*(min(mmax, ymax)-mmin+1)*(zmax-zmin+1)] like;
           int counter;
           int R1;
           int R2;
           counter = 1;
           for(yX in ymin:ymax){
             for(mX in mmin:min(mmax, ymax)){
               for(zX in max(zmin, mX):zmax){
                 R1 = max(0, yX+zX-mX-n);
                 R2 = min(yX-mX, zX-mX);
                 
                 vector[R2-R1+1] temp;
                 for (r in R1:R2){
                  temp[r-R1+1] = lgamma(n+1) + mX*logpij[1] - lgamma(mX+1)
                          + (yX-mX-r)*logpij[2] - lgamma(yX-mX-r+1) + r*logpij[3] - lgamma(r+1)
                          + (n-yX-zX+mX+r)*logpij[4] - lgamma(n-yX-zX+mX+r+1) + (zX-mX-r)*logpij[5] - lgamma(zX-mX-r+1);
                 }
                 like[counter] = log_sum_exp(temp);
                 
                 counter = counter + 1;
               }
             }
           }
           return (log_sum_exp(like[1:(counter-1)]));
  }
}

data{
  int<lower=0> I; // number of studies
  int<lower=I> NRec; // number of records
  int<lower=1> i[NRec]; // study number
  int<lower=0> j[NRec]; // group (j=0 control, j=1 test)
  int<lower=0> n[NRec]; // patients in arm
  int<lower=0> y[NRec]; // patients with event
  int<lower=0> m[NRec]; // patients with fatal event
  int<lower=0> z[NRec]; // patients that did not complete
  real<lower=0> tau[NRec]; // planned duration
  
  real<lower=0,upper=1> mw; // weight given to MAP prior: same weight for all MAP priors
  real<lower=0,upper=1> wbeta[2]; // weights of mixture beta component
  real<lower=0> alphaQ[2]; // alphas of mixture beta component
  real<lower=0> betaQ[2]; // betas of mixture beta component
  real<lower=0, upper=1> wL[3]; // weights of mixture components
  vector[3] meansL; // mean vectors of mixture components
  vector[3] sdL; // covariance matrices of mixture components
  real<lower=0,upper=1> wM[3]; // weights of mixture components
  vector[3] meansM; // mean vectors of mixture components
  vector[3] sdM; // covariance matrices of mixture components
}

transformed data{
  // Derive possible number of patients that could have experienced a 
  // non-fatal event and then have been non-administratively censored thereafter
  // before completing the trial. The possible minimum number of such patients
  // is R1 and the maximum number R2
  // int<lower=0> R1[NRec];
  // int<lower=0> R2[NRec];
  int<lower=0> j1m[NRec];
  real<lower=0,upper=1> mw1m;
  real avgmeansL;
  real avgmeansM;
  // vector[3] zero_vector;
  // cov_matrix[3] identity_matrix;
  
  avgmeansL = 0;
  avgmeansM = 0;
  for (r in 1:3){
    avgmeansL = avgmeansL + wL[r]*meansL[r];
    avgmeansM = avgmeansM + wM[r]*meansM[r];
  }
  
  for (r in 1:NRec){
    // R1[r] = max(0, y[r]+z[r]-m[r]-n[r]);
    // R2[r] = min(y[r]-m[r], z[r]-m[r]);
    j1m[r] = 1-j[r];
  }
  // zero_vector = rep_vector(0,3);
  // identity_matrix = diag_matrix(rep_vector(1,3));
  mw1m = 1-mw; 
}

parameters{
  // Log-hazard ratio for event times
  real muPhi;
  real<lower=0> sigmaPhi;
  // Trivariate random trial effect
  // 1st component relates to control event log-hazard rate
  // 2nd component relates to control drop-out log-hazard rate
  // 3rd component relates to test drop-out log-hazard rate
  vector[I] eta;
  // P(event being fatal on control (q0) or test (q1) treatment)
  real<lower=0,upper=1> q0;
  real<lower=0,upper=1> q1;
  // mean and scale parameters for the random trial effect
  array[I] real loglambda0;
  array[I] real logmu0;
  array[I] real logmu1;
}

transformed parameters {
  // Log-probability of each of the 5 multinomial outcomes in each trial arm
  vector[5] logp[NRec];
  // Log-hazard rate for events and drop-out in each trial arm
  real loglambda[NRec];
  real logmu[NRec];
  vector[I] logHR; // log-HR for each trial
  // Probability of an event being fatal in each trial arm
  real q[NRec];
  
  // Derive the transformed parameters based on hyperparameters 
  // and random effect components
  for (r in 1:NRec){
    logHR[i[r]] = muPhi + eta[i[r]]*sigmaPhi;
    loglambda[r] = loglambda0[i[r]] + logHR[i[r]]*j[r];
    logmu[r] = logmu0[i[r]]*j1m[r] + logmu1[i[r]]*j[r];
    q[r] = q0*(1-j[r])+q1*j[r];
    logp[r,1] = log(q[r]) + loglambda[r]
                 - log_sum_exp(logmu[r], loglambda[r])
                 + log1m_exp(-tau[r]*exp(log_sum_exp(logmu[r],loglambda[r])));
    logp[r,2] = log1m(q[r]) + log_diff_exp(-exp(logmu[r])*tau[r],
                 -tau[r]*exp(log_sum_exp(logmu[r],loglambda[r])));
    logp[r,4] = -tau[r]*exp(log_sum_exp(logmu[r],loglambda[r]));
    logp[r,5] = logmu[r] - log_sum_exp(logmu[r],loglambda[r])
                 + log1m_exp(-tau[r]*exp(log_sum_exp(logmu[r],loglambda[r])));
    // Most convenient to calculate logp[r,3] as below
    logp[r,3] = log1m_exp(log_sum_exp(log_sum_exp(logp[r,1],logp[r,2]),
                 log_sum_exp(logp[r,4], logp[r,5])));
  }
}

model {
  muPhi ~ cauchy(0, 0.37); // scale = 0.37 or 2.5
  sigmaPhi ~ normal(0, 0.5);
  // eta ~ multi_normal(zero_vector, identity_matrix);
  
  target +=
    log_sum_exp(
      log(mw) + log_sum_exp(
            log(wbeta[1]) + beta_lupdf(q0 | alphaQ[1], betaQ[1]),
            log(wbeta[2]) + beta_lupdf(q0 | alphaQ[2], betaQ[2])
                ),
      log(mw1m) + beta_lupdf(q0 | 0.5, 0.5)
    );
  
  target +=
    log_sum_exp(
      log(mw) + log_sum_exp(
            log(wbeta[1]) + beta_lupdf(q1 | alphaQ[1], betaQ[1]),
            log(wbeta[2]) + beta_lupdf(q1 | alphaQ[2], betaQ[2])
                ),
      log(mw1m) + beta_lupdf(q1 | 0.5, 0.5)
    );
  
  for (r in 1:I){
    eta[r] ~ normal(0, 1);
    // robust hist. priors with weakly informative mixture components 
    target +=
      log_sum_exp(
        log(mw1m) + normal_lupdf(loglambda0[r] | avgmeansL, log(10)),
        log(mw) + log_sum_exp(
          log_sum_exp(
            log(wL[1]) + student_t_lupdf(loglambda0[r] | 3, meansL[1], sdL[1]),
            log(wL[2]) + student_t_lupdf(loglambda0[r] | 3, meansL[2], sdL[2])
          ),
          log(wL[3]) + student_t_lupdf(loglambda0[r] | 3, meansL[3], sdL[3])
        )
      );
  
    target += 
      log_sum_exp(
        log(mw1m) + normal_lupdf(logmu0[r] | avgmeansM, log(10)),
        log(mw) + log_sum_exp(
          log_sum_exp(
            log(wM[1]) + student_t_lupdf(logmu0[r] | 3, meansM[1], sdM[1]),
            log(wM[2]) + student_t_lupdf(logmu0[r] | 3, meansM[2], sdM[2])
          ),
          log(wM[3]) + student_t_lupdf(logmu0[r] | 3, meansM[3], sdM[3])
        )
      );
  
    target += 
      log_sum_exp(
        log(mw1m) + normal_lupdf(logmu1[r] | avgmeansM, log(10)),
        log(mw) + log_sum_exp(
          log_sum_exp(
            log(wM[1]) + student_t_lupdf(logmu1[r] | 3, meansM[1], sdM[1]),
            log(wM[2]) + student_t_lupdf(logmu1[r] | 3, meansM[2], sdM[2])
          ),
          log(wM[3]) + student_t_lupdf(logmu1[r] | 3, meansM[3], sdM[3])
        )
      );
  }
  
  
  // Increment the log-L by likelihood contribution of each trial arm 
  for (r in 1:NRec) {
    target += logL_mymodel(n[r], y[r], m[r], z[r], logp[r]);
  }
    
}
