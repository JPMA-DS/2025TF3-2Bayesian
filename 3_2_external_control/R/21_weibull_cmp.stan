data {
  int <lower = 0> N;
  array[N] real time;
  array[N] int cnsr;
  array[N] int trt;
  array[N] int ext;
}

parameters {
  real beta_ext;
  real <lower = 0> tau;
  real beta_trt;
  real beta0;
  real <lower = 0> r;
}

transformed parameters {
  array[N] real lambda;
  for (n in 1:N){
    lambda[n] = exp(beta0*(1 - ext[n]) + beta_ext * ext[n] + beta_trt * trt[n]);
  }
}

model {
  beta_ext ~ normal(0, 1000);
  tau      ~ gamma(1, 0.001);
  beta0    ~ normal(beta_ext, (1/tau));
  beta_trt ~ normal(0, 1000);
  r        ~ exponential(1);

  for(n in 1:N){
    if(cnsr[n] == 1){
      target += weibull_lccdf(time[n]| r, (1/lambda[n])^(1/r));
    }else{
      target += weibull_lpdf(time[n]| r, (1/lambda[n])^(1/r));
    }
  }
}
