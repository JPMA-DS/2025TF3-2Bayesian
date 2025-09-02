/*----------------------------------------------------------------
#Title : Bayesian dynamic approach (CMP)
#Date  : 2024/10/30
#----------------------------------------------------------------*/

proc import out = data_CMP
		datafile = "****/data_CMP.csv" /* ****にはdata_CMP.csvを格納したフォルダのパスを指定する */
		dbms= csv replace;
run;

title 'Weibull Survival Model (CMP)';
proc mcmc data = data_CMP outpost = weisurvout nbi = 10000 nmc = 110000 thin = 10 seed = 1234
          monitor = (_parms_) stats = (summary intervals);

   parms beta_ext 0;
   parms tau 1;
   parms beta_trt 0;
   parms beta0 0;
   parms r 1; 
   
   hyperprior beta_ext: ~ normal(0, prec = 0.001);
   hyperprior tau:      ~ gamma(1, is = 0.001);
   prior beta_trt:      ~ normal(0, prec = 0.001);
   prior beta0:         ~ normal(beta_ext, prec = tau);
   prior r:             ~ expon(is = 1);

   lambda = beta0*(1-ext) + beta_ext*ext + beta_trt*trt;
   mu = (1 / exp(lambda) )**(1/r);

   llike = (1-cnsr)*logpdf('weibull', time, r, mu) + cnsr*logsdf('weibull', time, r, mu);
   model general(llike);
run;
