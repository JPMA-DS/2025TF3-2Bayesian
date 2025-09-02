# パッケージのロード
library(MASS)#多変量正規分布で乱数発生
library(tidyverse)#データ成形用
library(brms)
library(rstan)
library(bayesplot)

#仮想データの生成準備（scenarioごとにdata_mean1を変更する必要あり）とシミュレーション準備

#小児薬効★
#data_mean1 <- c(8, 15, 21,  24) #scenario1
#data_mean1 <- c(4, 7.5, 10.5,  12) #scenario3
data_mean1 <- c(0, 0, 0,  0) #scenario5
data_mean2 <- c(0, 0,  0, 0)

cov <- matrix(c(2078, 1463, 1075,  950,
                1463, 2266,  1543,  1450,
                1075, 1543,  2858,  1950,
                950, 1450,  1950,  2858), nrow = 4, byrow = TRUE)


#シミュレーション用の値
w<-0.8               #★mixw
#w<-1

n_simulations <- 1 

#事前分布形成用
mean1 <- c(8, 15, 21,  24)
mean2 <- c(0, 0,  0, 0)
var1 <- matrix(c(67, 69,  56,  45,
                 69, 78,  71,  63,
                 56, 71,  78,  79,
                 45, 63,  79,  90), nrow = 4, byrow = TRUE)

var2 <- matrix(c(81, 0, 0, 0, 
                 0, 81, 0, 0, 
                 0, 0, 81, 0, 
                 0, 0, 0, 81 ), nrow = 4, byrow = TRUE)

#仮想データの生成準備（scenarioごとにdata_mean1を変更する必要あり）とシミュレーション準備

#小児薬効★
#data_mean1 <- c(8, 15, 21,  24) #scenario1
#data_mean1 <- c(4, 7.5, 10.5,  12) #scenario3
data_mean1 <- c(0, 0, 0,  0) #scenario5
data_mean2 <- c(0, 0,  0, 0)

cov <- matrix(c(2078, 1463, 1075,  950,
                1463, 2266,  1543,  1450,
                1075, 1543,  2858,  1950,
                950, 1450,  1950,  2858), nrow = 4, byrow = TRUE)


#シミュレーション用の値
w<-0.8               #★mixw
#w<-1

n_simulations <- 3000

#事前分布形成用
mean1 <- c(8, 15, 21,  24)
mean2 <- c(0, 0,  0, 0)
var1 <- matrix(c(67, 69,  56,  45,
                 69, 78,  71,  63,
                 56, 71,  78,  79,
                 45, 63,  79,  90), nrow = 4, byrow = TRUE)

var2 <- matrix(c(81, 0, 0, 0, 
                 0, 81, 0, 0, 
                 0, 0, 81, 0, 
                 0, 0, 0, 81 ), nrow = 4, byrow = TRUE)

# 乱数生成(17例/群)
rv1 <- MASS::mvrnorm(17, data_mean1, cov)
rv2 <- MASS::mvrnorm(17, data_mean2, cov)

# 治療群の準備：1:active, 0:placebo
trt1 <- rep(1, 17)
trt2 <- rep(0, 17)

# データフレームの作成
rv1b <- data.frame(trt = trt1, rv1)
rv2b <- data.frame(trt = trt2, rv2)


#モデルの設定
formula <- bf(value ~ 1 + time * trt, autocor = ~unstr(time=time, gr=subject_id))

#仮の事前分布の設定
priors <- c(
  set_prior("normal(0,1e10)", class = "b",coef = "time2"),
  set_prior("normal(0,1e10)", class = "b",coef = "time3"),
  set_prior("normal(0,1e10)", class = "b",coef = "time4"),
  set_prior("normal(9,81)", class = "b", coef = "trt"),
  set_prior("normal(9,81)", class = "b", coef = "time2:trt"),
  set_prior("normal(9,81)", class = "b", coef = "time3:trt"),
  set_prior("normal(9,81)", class = "b", coef = "time4:trt"),
  set_prior("lkj(1)", class = "Lcortime")
)


results1 <- data.frame(simulation = integer(n_simulations),
                       mean_n7 = logical(n_simulations), sd_b7 = logical(n_simulations), '2.5%' = logical(n_simulations), '5%' = logical(n_simulations), '10%' = logical(n_simulations), '50%' = logical(n_simulations)) 
results2 <- data.frame(simulation = integer(n_simulations), lower_bound_zero = logical(n_simulations)) 
for (i in 1:n_simulations){
  
  
  # データの結合
  H6DMCLVHV <- rbind(rv1b, rv2b)
  
  # データの整形
  H6DMCLVHV <- H6DMCLVHV %>%
    mutate(subject_id = rep(1:34,1))
  
  data1 <- H6DMCLVHV %>%
    pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "value") %>%
    mutate(time = as.character(str_remove(time, "X")))
  
  #stan_code:仮のモデルのstanコードを格納＝＞コードを更新してシミュレーションを実施．
  stan_code <- stancode(
    formula,
    data = data1,
    family = gaussian,
    prior = priors
  )
  
  #データをstanで使えるようにする
  stan_data<-standata(
    formula,
    data = data1,
    family = gaussian
  )
  #データ以外で使いたい情報をstan_dataに組み込む*stanコードへの追記と併せて必要
  stan_data$mean1 <- mean1
  stan_data$mean2 <- mean2
  stan_data$var1 <- var1
  stan_data$var2 <- var2
  stan_data$w <- w 
  stan_data$prior_only <- 0 #1：prior onlyでsampling
  
  #stancode2:stan_codeを元にモデルを改造．内容はstancode2.Rを参照
  stanmodel2 <- stan_model(model_code = stan_code2)
  fit_custom <- rstan::sampling(stanmodel2, data= stan_data,
                                iter=200000, warmup=10000, thin=50, chain =1, cores=1)
  
  # パラメータ b[7] trt:time4 を抽出
  b7_samples <- extract(fit_custom, pars = "b[7]")$`b[7]`
  # 統計量を計算
  mean_b7 <- mean(b7_samples)
  sd_b7 <- sd(b7_samples)
  
  quantiles_b7 <- quantile(b7_samples, probs = c(0.025,0.05,0.10, 0.5))
  # 結果を表示
  summary_bayes <- cbind(mean_b7, sd_b7, t(quantiles_b7))
  
  lower_bound_zero <- quantiles_b7[2] > 0
  
  results1[i,] <-c(i,summary_bayes)
  results2[i,] <-c(i,lower_bound_zero )
}

power00<-results2[,2]
power <- sum(power00)/n_simulations

#effect size確認
actual <- subset(H6DMCLVHV,H6DMCLVHV[,1]==1)
placebo <- subset(H6DMCLVHV,H6DMCLVHV[,1]==0)
trt_diff <- mean(actual[,5]) - mean(placebo[,5])

var_a<-var(actual[,5])
n_a<-length(actual[,5])
var_p<-var(placebo[,5])
n_p<-length(placebo[,5])

Sc <-sqrt((n_a*var_a+n_p*var_p)/(n_a+n_p))
ES <-trt_diff/Sc
c(trt_diff, Sc, ES)