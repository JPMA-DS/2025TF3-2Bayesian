/*Bayesian MMRM with a mixture prior*/

/*
"Sim" macroの説明:
scenario:	1=Optimistic, 3=50% effect of optimistic scenario, 5=Tadalafil is no different than placebo (Null)
prior: 			1=diffuse prior, 2=mixture prior, 3=modified mixture prior
mixw: 			Weight for mixture prior(only prior=2)
suc: 			Threshold of success
nperg: 		Number of patient by group
*/

libname tf32 "格納先のフォルダパスを指定してください";

/*"Sim"マクロ内で自動実行するためマクロ変数の指定は不要*/
%macro MCMC(ess=);
%if &ess. = 0 %then %do;
	ods output PostSumInt=PostSumInt;
%end;
/*ESS用にimproper priorでmcmcを実行するための条件分岐*/
%if &ess. = 1 %then %do;
	%let prior_s=&prior.;
	%let prior=1;
	ods output PostSumInt=PostSumInt_ESS;
%end;

/*各priorの設定は適宜調整して使用してください。現在は"自己相関が高い"等のwarningが頻出しないように調整しております。*/
%if &prior.=1 %then %do;
proc mcmc data=H6DMCLVHV seed=4989 nbi=10000 nmc=50000 thin=10
	%if &ess. = 0 %then %do;
		outpost=POST
	%end;
	%if &ess. = 1 %then %do;
		outpost=POST_ESS
	%end;
monitor=(_parms_ lsm11-lsm16 lsm21-lsm26 lsmd1-lsmd6) diag=all plot(smooth)=all maxtune=50;
%end;
%if &prior.=2 %then %do;
proc mcmc data=H6DMCLVHV seed=4989 nbi=10000 nmc=200000 thin=40 outpost=POST monitor=(_parms_ lsm11-lsm16 lsm21-lsm26 lsmd1-lsmd6) diag=all plot(smooth)=all maxtune=50;
%end;
%if &prior.=3 %then %do;
proc mcmc data=H6DMCLVHV seed=4989 nbi=10000 nmc=400000 thin=80 outpost=POST monitor=(_parms_ lsm11-lsm16 lsm21-lsm26 lsmd1-lsmd6 w) diag=all plot(smooth)=all maxtune=50;
%end;

  array y [4] m1-m4;
  array mu[4] mu1-mu4;
  array R[4, 4];
  array S[4, 4];
  array beta_m [4] beta_m1-beta_m4;
  array beta_tm[4] beta_tm11-beta_tm14;
  array lsm1[4] lsm11-lsm14;
  array lsm2[4] lsm21-lsm24;
  array lsmd[4] lsmd1-lsmd4; 

  %if &prior. = 2 or &prior. = 3 %then %do;
  /*mixture prior*/
  /*adult component of the mixture prior*/
  /*論文中で示された事前分布が非正定値行列であったため（四捨五入とかの兼ね合いと思われる）4,4要素を82→83、5,5要素を89→90にして正定値行列に変更した*/
  array muA[4] muA1-muA4 (8 15 21 24);
  array SA[4, 4] SA1-SA16 (67 69 56 45,
						                		 69 78 71 63,
                                        56 71 78 79,
                                        45 63 79 90);

  /*skeptical component of the mixture prior*/
  array muS[4] muS1-muS4 (0 0 0 0);
  array SS[4, 4] SS1-SS16 (81 0 0 0 ,
						                		 0 81 0 0 ,
                                        0 0 81 0 ,
                                        0 0 0 81);
  /*for calc of mixture prior*/
  array DIFA[4];
  array tDIFA[1, 4];
  array invSA[4, 4];
  array XMA[1, 4];
  array XMXA[1];

  array DIFS[4];
  array tDIFS[1, 4];
  array invSS[4, 4];
  array XMS[1, 4];
  array XMXS[1];
  %end;
 
  begincnst;
    call identity(S);

	%if &prior. = 2 or &prior. = 3 %then %do;
	call det(SA, dSA);
	call inv(SA, invSA);

	call det(SS, dSS);
	call inv(SS, invSS);
	%end;
  endcnst;

  parms beta_m1 0 beta_m2 0 beta_m3 0 beta_m4 0;
  %if &prior.=1 %then %do;
  parms beta_tm11 0 beta_tm12 0 beta_tm13 0 beta_tm14 0;
  %end;

  %if &prior. = 2 or &prior. = 3 %then %do;
  /*自分で関数を設定した場合は初期値を設定する必要あり*/
  parms beta_tm {8 15 21 24};
  %end;
  %if &prior.=3 %then %do;
  parms w 0.5;
  %end;
  parms R;

  %if &prior. = 2 or &prior. = 3 %then %do;
  beginnodata;
  	n=4;
	%if &prior.=2 %then %do;
	w=&mixw.;
	%end;
	call subtractmatrix(beta_tm, muA, DIFA);
	call transpose(DIFA, tDIFA);
	call mult(tDIFA, invSA, XMA);
	call mult(XMA, DIFA, XMXA);

	call subtractmatrix(beta_tm, muS, DIFS);
	call transpose(DIFS, tDIFS);
	call mult(tDIFS, invSS, XMS);
	call mult(XMS, DIFS, XMXS);

	const_adlut=1/((CONSTANT('PI')**(n/2))*(dSA**(1/2)));
	adlut_pdf=const_adlut*exp(-0.5*XMXA[1]);

	const_skeptical=1/((CONSTANT('PI')**(n/2))*(dSS**(1/2)));
	skeptical_pdf=const_skeptical*exp(-0.5*XMXS[1]);

	lp=log(w*adlut_pdf + (1-w)*skeptical_pdf);
  endnodata;
  %end;

  prior beta_m: ~ general(0);
  /*Diffuse prior(本PRGではimproper priorを指定)*/
  %if &prior.=1 %then %do;
  prior beta_tm: ~ general(0);
  %end;
  %if &prior.=3 %then %do;
  /*Prior for the weight*/
  prior w ~ beta(0.5, 0.5);
  %end;
  %if &prior. = 2 or &prior. = 3 %then %do;
  /*Mixture prior*/
  prior beta_tm ~ general(lp);
  %end;
  prior R ~ iwish(4, S);

  beginnodata;
    do j=1 to 4;
      lsm1[j] = beta_m[j] + beta_tm[j];
      lsm2[j] = beta_m[j];
      lsmd[j] = beta_tm[j];
    end;
  endnodata;

 do j=1 to 4;
    if trt=1
      then mu[j] = beta_m[j] + beta_tm[j];
      else mu[j] = beta_m[j];
 end;

 model y ~ mvn(mu, R);
run;

/*再格納*/
%if &ess. = 1 %then %do;
	%let prior=&prior_s.;
%end;
%mend;

%macro Sim(nsim=, scenario=, prior=, mixw=, suc=, nperg=);

%let START_TIME = %sysfunc( datetime() );

%do sim=1 %to &nsim.;

ods noresults;
ods listing close;
ods html close;
options nosource nonotes;

/*Simulation data generation*/
PROC IML;
   %if &scenario.=1 %then %do;
   mean1={ 8 15 21 24};
   %end;
   %if &scenario.=3 %then %do;
   mean1={ 4 7.5 10.5 12};
   %end;
   %if &scenario.=5 %then %do;
   mean1={ 0 0 0 0};
   %end;
   mean2={ 0 0 0 0};

   cov={ 	2078 1463 1075 950,
 			1463 2266 1543 1450,
			1075 1543 2858 1950,
      		950 1450 1950 2858};

   rv1=RandNormal(&nperg., mean1, cov);
   rv2=RandNormal(&nperg., mean2, cov);

	trt1=J(&nperg.,1,1);
	trt2=J(&nperg.,1,2);
	rv1b=trt1||pat||rv1;
	rv2b=trt2||pat||rv2;

   CREATE H6DMCLVHV FROM rv1b[colname={'trt' 'm1' 'm2' 'm3' 'm4'}];
   APPEND FROM rv1b;
   APPEND FROM rv2b;
quit;

/*BMMRM*/
%MCMC(ess=0);

/*ESSの算出はmixture priroのみ実行*/
%if &prior. = 2 or &prior. = 3 %then %do;
%MCMC(ess=1);
%end;

/*Decision Criteria by posterior probability*/
/*Power and alpha-error*/
data post2;
	set post;
	if lsmd4 >0 then resp=1;
	else resp=0;
run;

proc freq data=post2 noprint;
	table resp / out=out;
run;

data out2;
	set out(where=(resp=1));
	if &suc.<=percent then suc=1;
	else suc=0;
	sim=&sim.;
run;

proc append base = base data = out2 force;
run;

/*Bias and MSE*/
proc append base = base_smry data = PostSumInt(where=(Parameter in("lsmd4", "w"))) force;
run;

/*ESS*/
%if &prior. = 2 or &prior. = 3 %then %do;
data ESS;
	merge 	PostSumInt(keep=Parameter StdDev where=(Parameter in("lsmd4"))) 
				PostSumInt_ESS(keep=Parameter StdDev rename=(StdDev=StdDev_noninfo) where=(Parameter in("lsmd4")));
	V = StdDev**2;
	V_noninfo = StdDev_noninfo**2;
	posterior_ESS = (&nperg. * 2) * (V_noninfo / V);
	prior_ESS = posterior_ESS - (&nperg. * 2); 
run;

proc append base = base_ESS data = ESS force;
run;
%end;

%end;

/*ESS(成果物では"prior ESSを出力")*/
%if &prior. = 2 or &prior. = 3 %then %do;
proc means data=base_ESS noprint;
	var prior_ESS posterior_ESS;
	output out=tf32.sim_p&prior._w%sysevalf(&mixw. * 10)_s&scenario._suc&suc._n&nperg._ESS mean= std= / autoname;
run;
%end;

/*Power and alpha error(suc=1の割合が検出力またはalpha error)*/
proc freq data=base;
	table suc / out=tf32.sim_p&prior._w%sysevalf(&mixw. * 10)_s&scenario._suc&suc._n&nperg.;
run;

/*Bias and MSE(prior=3の場合はmix weightの推定値も出力)*/
proc sort data=base_smry; by Parameter; run;
data base_smry2;
	set base_smry;
	%if &scenario.=1 %then %do;
	mu=24;
	%end;
	%if &scenario.=3 %then %do;
	mu=12;
	%end;
	%if &scenario.=5 %then %do;
	mu=0;
	%end;

	if Parameter = "lsmd4" then do;
		bias = mean -mu;
		mse = (mean -mu)**2;
	end;
run;
proc means data=base_smry2 noprint;
	var mean bias mse;
	output out=tf32.sim_p&prior._w%sysevalf(&mixw. * 10)_s&scenario._suc&suc._n&nperg._smry mean= std= / autoname;
	by Parameter;
run;

proc datasets library=work kill nolist;
quit;

ods results;
ods listing;
ods html;
options source notes;

%put 処理時間 = %sysevalf( %sysfunc( datetime() ) - &START_TIME. ); 

%mend;

/*macro実行例。mixwはprior=2の時のみ有効*/
%sim(nsim=3000, scenario=1, prior=1, mixw=1, suc=95, nperg=17);
%sim(nsim=3000, scenario=3, prior=2, mixw=0.8, suc=95, nperg=51);
%sim(nsim=3000, scenario=5, prior=3, mixw=1, suc=95, nperg=51);
