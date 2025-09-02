ods graphics off;
ods exclude all; 
ods results off; 
options nonotes nosource; 


/************************************************
�V�~�����[�V�����v���O����
*************************************************/
/*�}�N���ϐ��̒�`*/
/*nsim: �V�~�����[�V������*/
/*n_min: �ŏ��ᐔ*/
/*n_max: �ő�ᐔ*/
/*alpha_e, beta_e: ���ґt�������̎��O���z�̒��p�����[�^*/
/*alpha_s, beta_s: 臒l�t�������̎��O���z�̒��p�����[�^*/
/*d0: �������ÂɊւ���Œዖ�e����*/
/*p: 臒l�t�������Ⴕ���͊��ґt�������̐^�l*/

%macro SIM(nsim, n_min, n_max, alpha_e, beta_e, alpha_s, beta_s, d0, p);
/******************** 
�����f�[�^���� 
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
Thall and simon�̐ϕ��v�Z
************************/
proc iml;
  aE=&alpha_e.; bE=&beta_e.; aS=&alpha_s.;  bS=&beta_s.; d0=&d0.;
  /*�f�[�^�Z�b�g��ǂݍ��� */ 
  use D1;
  read all var {study sum n} into input;/
  close D1;
  /* �������x�֐�f�̒�` */
  start f(x) global(aS, bS, aE, bE, d0, n, y, delta);
        first=1-cdf('Beta', X+d0, aE+y, bE+n-y);
        second= pdf('Beta', X, aS, bS);
        fx = first*second;
        return(fx);
   finish;
   delta = 1-d0;
  /* IML�v���V�W�����ɂ�DO���[�v */
  do i=1 to nrow(input);
    y = input[i,2];
    n = input[i,3];
    call quad(result, "f", 0 || delta);/*�������x�֐�f��0����delta�܂Őϕ�����*/
    t_results = t_results//result;
   end;
/* �f�[�^�Z�b�g�ɏo��  */
  output = input || t_results;
  cname={"study" "sum" "n" "result"};
  create outdata from output [colname=cname];
  append from output;
  close outdata;
quit;

/******************************
�f�[�^�Z�b�g�̉��H
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
Power(alpha error)�̌v�Z
***************************/
proc freq data=flg;
  tables ef/out=Power;
run;

 /****************************
�L�����~���_�̕��ϔ팟�Ґ��̌v�Z
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


%SIM(nsim = 100000, n_min = 10, n_max = 15, alpha_e = 0.5, beta_e = 0.5, alpha_s = 34.4, beta_s = 137.6, d0 = 0, p = 0.5);/*�Η�����p=0.5��^�l�Ƃ����ꍇ��Power*/
%SIM(nsim = 100000, n_min = 10, n_max = 15, alpha_e = 0.5, beta_e = 0.5, alpha_s = 34.4, beta_s = 137.6, d0 = 0, p = 0.2);/*�A������p=0.2��^�l�Ƃ����ꍇ��Power�i�����G���[�j*/


ods graphics on;
ods exclude none;
ods results on;
options notes;
