ods graphics off;
ods exclude all; 
ods results off; 
options nonotes nosource; 


/************************************************
シミュレーションプログラム
*************************************************/
/*マクロ変数の定義*/
/*nsim: シミュレーション回数*/
/*n_min: 最小例数*/
/*n_max: 最大例数*/
/*alpha_e, beta_e: 期待奏効割合の事前分布の超パラメータ*/
/*alpha_s, beta_s: 閾値奏効割合の事前分布の超パラメータ*/
/*d0: 試験治療に関する最低許容増分*/
/*p: 閾値奏功割合若しくは期待奏功割合の真値*/

%macro SIM(nsim, n_min, n_max, alpha_e, beta_e, alpha_s, beta_s, d0, p);
/******************** 
乱数データ生成 
*********************/
data D1;
    call streaminit(123456);
    retain sum 0 n 0;
    do study=1to &nsim.;
     do i = 1 to &n_max.;
                AVAL = rand('BERNOULLI',&p.);
                sum = sum + AVAL;
                n = n + 1;
                if i >= &n_min. then output;
       end;
      sum = 0; n=0; 
    end;
    drop i;
run;

/***********************
Thall and simonの積分計算
************************/
proc iml;
  aE=&alpha_e.; bE=&beta_e.; aS=&alpha_s.;  bS=&beta_s.; d0=&d0.;
  /*データセットを読み込む */ 
  use D1;
  read all var {study sum n} into input;/
  close D1;
  /* 同時密度関数fの定義 */
  start f(x) global(aS, bS, aE, bE, d0, n, y, delta);
        first=1-cdf('Beta', X+d0, aE+y, bE+n-y);
        second= pdf('Beta', X, aS, bS);
        fx = first*second;
        return(fx);
   finish;
   delta = 1-d0;
  /* IMLプロシジャ内にてDOループ */
  do i=1 to nrow(input);
    y = input[i,2];
    n = input[i,3];
    call quad(result, "f", 0 || delta);/*同時密度関数fを0からdeltaまで積分する*/
    t_results = t_results//result;
   end;
/* データセットに出力  */
  output = input || t_results;
  cname={"study" "sum" "n" "result"};
  create outdata from output [colname=cname];
  append from output;
  close outdata;
quit;

/******************************
データセットの加工
*******************************/
data outdata2 ;
   set outdata;
   if result>=0.95 then flg=1;else  flg=0;
run;

proc sort data=outdata2; by study;run;
proc means data=outdata2 noprint;
  var flg;
  by study;
  output out=OUT_flg  SUM=SUM;
 run;

data flg;
  set OUT_flg;
  if SUM >=1 then ef=1;else ef=0;
run;



/**************************
Power(alpha error)の計算
***************************/
proc freq data=flg;
  tables ef/out=Power;
run;

 /****************************
有効中止時点の平均被検者数の計算
******************************/
 data num1;
   set outdata2;
run;

proc sort data=num1; by study descending flg ;run;

data num2;
  set num1;
  by study descending flg;
  if first.flg=1;
  if flg=1 then NUM=n ;else NUM=15;
run;

proc sort data=num2 out=num3 nodupkey; by study; run; 

proc means data=num3;
  var num;
  output out=MeanN  Mean=Mean;
run;

%mend;


%SIM(nsim = 100000, n_min = 10, n_max = 15, alpha_e = 0.5, beta_e = 0.5, alpha_s = 34.4, beta_s = 137.6, d0 = 0, p = 0.5);/*対立仮説p=0.5を真値とした場合のPower*/
%SIM(nsim = 100000, n_min = 10, n_max = 15, alpha_e = 0.5, beta_e = 0.5, alpha_s = 34.4, beta_s = 137.6, d0 = 0, p = 0.2);/*帰無仮説p=0.2を真値とした場合のPower（＝αエラー）*/


ods graphics on;
ods exclude none;
ods results on;
options notes;
