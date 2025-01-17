// Functions for calculating likelihood
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

data {
  int<lower=0> I; // num of studies
  int<lower=I> NRec; // num of records (usually one per trial and treatment group combination)
  array[NRec] int i; // study number
  array[NRec] int j; // group (j=0 control, j=1 treatment)
  array[NRec] int n; // patients in trial i=1,...,I
  array[NRec] int y; // patients with event in trial i
  array[NRec] int m; // patients with fatal event in trial i
  array[NRec] int z; // num of patients discontinued
  array[NRec] real tau; // planned duration of trial i
}

transformed data {
  // Derive possible number of patients that could have
  // experienced a non-fatal event and then have been
  // non-administratively censored thereafter before completing
  // the trial. The possible minimum number of such patients is
  // R1 and the maximum number R2 for each trial i=1,...,I
  // array[NRec] int R1;
  // array[NRec] int R2;

  // Indicator for being in the control group
  // array[NRec] int j1m;

  // Zero vector and identity matrix to be used as mean and covariance matrix
  // of the tri-variate normal random trial effect. We specify the random effect
  // using a non-centered parameterization to improve the efficiency
  // of the Hamiltonian Monte-Carlo sampling by Stan
  vector[4] zero_vector;
  cov_matrix[4] identity_matrix;
  // for(r in 1:NRec){
  //   R1[r] = max(0, y[r]+z[r]-m[r]-n[r]);
  //   R2[r] = min(y[r]-m[r], z[r]-m[r]);
  //   // j1m[r] = 1 - j[r];
  // }
  zero_vector = rep_vector(0, 4);
  identity_matrix = diag_matrix(rep_vector(1, 4));
}

parameters {
  // log-harzard ratio for event times
  real muPhi;
  real<lower=0> sigmaPhi;
  
  // tri-variate random trial effect
  // 1st component is control event log-harzard rate lambda_i0
  // 2nd component is control arm drop-out log-harzard rate mu_i0
  // 3rd component is treatment arm drop-out log-harzard rate mu_i1
  array[I] vector[4] eta;
  
  // probability of fatal event
  real<lower=0, upper=1> q0;
  real<lower=0, upper=1> q1;
  
  // hierarchical mean and scale parameters for the random trial effects
  real muL;
  real muM0;
  real muM1;
  real<lower=0> sigmaL;
  real<lower=0> sigmaM0;
  real<lower=0> sigmaM1;
}

transformed parameters {
  // log-probability of each of the 5 multinomial outcomes in each trial arm
  array[NRec] vector[5] logp;
  
  // log-harzard rate for events and drop-out in each trial arm
  vector[NRec] loglambda;
  vector[NRec] logmu;
  vector[I] logHR; // log-HR for each trial
  
  // probability of an event being fatal
  vector[NRec] q;
  
  // derive the transformed parameters based on hyper-parameters and random effect components
  for (r in 1:NRec) {
    logHR[i[r]] = muPhi + eta[i[r], 4]*sigmaPhi;
    loglambda[r] = muL + eta[i[r], 1]*sigmaL + logHR[i[r]]*j[r];
    logmu[r] = (muM0 + eta[i[r], 2]*sigmaM0)*(1-j[r]) + (muM1 + eta[i[r], 3]*sigmaM1)*j[r];
    q[r] = q0*(1-j[r]) + q1*j[r];
    // five log-probabilities
    logp[r, 1] = log(q[r]) + loglambda[r] - log_sum_exp(loglambda[r], logmu[r]) + log1m_exp(-tau[r] * exp(log_sum_exp(loglambda[r], logmu[r])));
    logp[r, 2] = log1m(q[r]) + log1m_exp(-tau[r]*exp(loglambda[r])) - exp(logmu[r])*tau[r];
    logp[r, 4] = -tau[r] * exp(log_sum_exp(loglambda[r], logmu[r]));
    logp[r, 5] = logmu[r] - log_sum_exp(loglambda[r], logmu[r]) + log1m_exp(-tau[r] * exp(log_sum_exp(loglambda[r], logmu[r])));
    logp[r, 3] = log1m_exp( log_sum_exp(log_sum_exp(logp[r, 1], logp[r, 2]), log_sum_exp(logp[r, 4], logp[r, 5])) );
  }
}


model {
  // priors
  muPhi ~ cauchy(0, 0.37); // scale = 0.37 or 2.5
  sigmaPhi ~ normal(0, 100);
  q0 ~ beta(0.5, 0.5);
  q1 ~ beta(0.5, 0.5);
  muL ~ normal(0, 1000); // normal(-4.27, log(10)) as weakly-informative prior
  muM0 ~ normal(0, 1000); // normal(log(0.22), log(10)) as weakly-informative prior
  muM1 ~ normal(0, 1000); // normal(log(0.22), log(10)) as weakly-informative prior
  sigmaL ~ normal(0, 100); // normal(0, log(10)) as weakly-informative prior
  sigmaM0 ~ normal(0, 100); // normal(0, log(10)) as weakly-informative prior
  sigmaM1 ~ normal(0, 100); // normal(0, log(10)) as weakly-informative prior
  
  // specify tri-variate random trial effects
  for (r in 1:I){
    eta[r] ~ multi_normal(zero_vector, identity_matrix);
  }
  // eta ~ multi_normal(zero_vector, identity_matrix);
  
  // logHR ~ normal(muPhi, sigmaPhi);
  
  // likelihood
  for(r in 1:NRec){
    target += logL_mymodel(n[r], y[r], m[r], z[r], logp[r]);
  }
}
