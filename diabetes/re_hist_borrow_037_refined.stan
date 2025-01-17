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
  real<lower=0, upper=1> wL[5]; // weights of mixture components
  array[5] vector[2] meansL; // mean vectors of mixture components
  array[5] cov_matrix[2] covL; // covariance matrices of mixture components
  real<lower=0,upper=1> wM[5]; // weights of mixture components
  array[5] vector[2] meansM; // mean vectors of mixture components
  array[5] cov_matrix[2] covM; // covariance matrices of mixture components
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
  vector[4] zero_vector;
  cov_matrix[4] identity_matrix;
  
  avgmeansL = 0;
  avgmeansM = 0;
  for (r in 1:5){
    avgmeansL = avgmeansL + wL[r]*meansL[r,1];
    avgmeansM = avgmeansM + wM[r]*meansM[r,1];
  }
  
  for (r in 1:NRec){
    // R1[r] = max(0, y[r]+z[r]-m[r]-n[r]);
    // R2[r] = min(y[r]-m[r], z[r]-m[r]);
    j1m[r] = 1-j[r];
  }
  zero_vector = rep_vector(0,4);
  identity_matrix = diag_matrix(rep_vector(1,4));
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
  vector[4] eta[I];
  // P(event being fatal on control (q0) or test (q1) treatment)
  real<lower=0,upper=1> q0;
  real<lower=0,upper=1> q1;
  // Hierarchical mean and scale parameters for the random trial effect
  vector[2] musigL; // this is a vector of length two, 1st element is muL, 2nd element is log-sigmaL
  vector[2] musigM0;
  vector[2] musigM1;
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
    logHR[i[r]] = muPhi + eta[i[r], 4]*sigmaPhi;
    loglambda[r] = musigL[1] + eta[i[r],1]*exp(musigL[2]) + logHR[i[r]]*j[r];
    logmu[r] = (musigM0[1]+eta[i[r],2]*exp(musigM0[2]))*j1m[r]
                + (musigM1[1]+eta[i[r],3]*exp(musigM1[2]))*j[r];
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
  sigmaPhi ~ normal(0, 0.5); // test run with scale = 0.1, 0.5, log(10), 100
  eta ~ multi_normal(zero_vector, identity_matrix);
  // robust hist. priors with weakly informative mixture components 
  target +=
    log_sum_exp(
      log(mw) + log_sum_exp(
            log(wbeta[1]) + beta_lpdf(q0 | alphaQ[1], betaQ[1]),
            log(wbeta[2]) + beta_lpdf(q0 | alphaQ[2], betaQ[2])
                ),
      log(mw1m) + beta_lpdf(q0 | 0.5, 0.5)
    );
  
  target +=
    log_sum_exp(
      log(mw) + log_sum_exp(
            log(wbeta[1]) + beta_lpdf(q1 | alphaQ[1], betaQ[1]),
            log(wbeta[2]) + beta_lpdf(q1 | alphaQ[2], betaQ[2])
                ),
      log(mw1m) + beta_lpdf(q1 | 0.5, 0.5)
    );
  
  target +=
    log_sum_exp(
      log(mw1m) + normal_lpdf(musigL[1] | avgmeansL, log(10)) + log(2)
      + musigL[2] + normal_lpdf(exp(musigL[2]) | 0,log(10)),
      log(mw) + log_sum_exp(
        log_sum_exp(
          log_sum_exp(
            log(wL[1]) + multi_normal_lpdf(musigL | meansL[1], covL[1]),
            log(wL[2]) + multi_normal_lpdf(musigL | meansL[2], covL[2])
          ),
          log_sum_exp(
            log(wL[3]) + multi_normal_lpdf(musigL | meansL[3], covL[3]),
            log(wL[4]) + multi_normal_lpdf(musigL | meansL[4], covL[4])
          )
        ),
        log(wL[5]) + multi_normal_lpdf(musigL | meansL[5], covL[5])
      ));
  
  target += 
    log_sum_exp(
      log(mw1m) + normal_lpdf(musigM0[1] | avgmeansM,log(10)) + log(2)
       + musigM0[2] + normal_lpdf(exp(musigM0[2]) | 0,log(10)),
      log(mw) + log_sum_exp(
        log_sum_exp(
          log_sum_exp(
            log(wM[1]) + multi_normal_lpdf(musigM0 | meansM[1], covM[1]),
            log(wM[2]) + multi_normal_lpdf(musigM0 | meansM[2], covM[2])
          ),
          log_sum_exp(
            log(wM[3]) + multi_normal_lpdf(musigM0 | meansM[3], covM[3]),
            log(wM[4]) + multi_normal_lpdf(musigM0 | meansM[4], covM[4])
          )
        ),
        log(wM[5]) + multi_normal_lpdf(musigM0 | meansM[5], covM[5])
      )
    );
  
  target +=
    log_sum_exp(
      log(mw1m) + normal_lpdf(musigM1[1] | avgmeansM,log(10)) + log(2)
       + musigM1[2] + normal_lpdf(exp(musigM1[2]) | 0,log(10)),
      log(mw) + log_sum_exp(
        log_sum_exp(
          log_sum_exp(
            log(wM[1]) + multi_normal_lpdf(musigM1 | meansM[1], covM[1]),
            log(wM[2]) + multi_normal_lpdf(musigM1 | meansM[2], covM[2])
          ),
          log_sum_exp(
            log(wM[3]) + multi_normal_lpdf(musigM1 | meansM[3], covM[3]),
            log(wM[4]) + multi_normal_lpdf(musigM1 | meansM[4], covM[4])
          )
        ),
        log(wM[5]) + multi_normal_lpdf(musigM1 | meansM[5], covM[5])
      )
    );
  
  // Increment the log-L by likelihood contribution of each trial arm 
  for (r in 1:NRec) {
    target += logL_mymodel(n[r], y[r], m[r], z[r], logp[r]);
  }
    
}
