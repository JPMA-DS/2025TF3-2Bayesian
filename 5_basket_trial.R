library(mvtnorm); library(dplyr); library(rstan); library(trialr); library(DoseFinding)

#### 設定 ####
#中間解析前及び, 中間解析後のデータ数
N_total <- 100
N_beforeIA <-100/2
N_afterIA <-N_total - N_beforeIA

#がん種ごとの奏効率
p = c(0.35, 0.35, 0.35, 0.35, 0.35)

#がん種ごとの発現割合
Ctype_prop <- c(0.3, 0.3, 0.3, 0.05, 0.05)

#帰無仮説の奏効割合
null_prop <- 0.15

#階層ベイズモデルでTargetとする奏効割合
target_prop <- 0.25

#階層ベイズモデルの設定
mu <- -2.19 　#muの超事前分布の平均-
imu <- 4.00　 #muの超事前分布の標準偏差
alpha <- 2.00　 #tauの超事前分布のパラメータ
beta <- 2.00   #tauの超事前分布のパラメータ

#中間解析で無益性中止する基準: Targetの奏効割合以上の奏効割合である事後確率
target_posterior_IA <- 0.4

#最終解析での有効性判定基準: Targetの奏効割合以上の奏効割合である事後確率
target_posterior_FA <- 0.8

#シミュレーションの繰り返し回数
reprep <- 10000

#各がん種・各症例での奏効の有無を発生させる関数の定義
idata <- function(i, N, p) { 
  #i:シミュレーションでの試行回数 
  #N:各がん種での例数 
  #p:各がん種での真の奏効割合
  
  set.seed(0827 + i)
  N_max <- max(N)
  y <- lapply(1:length(N), function(j) {
    c(rbinom(n = N[j], size = 1, prob = p[j]), rep(NA, N_max - N[j]))
  })
  d <- as.data.frame(y)
  names(d) <- paste0("Y", 1:length(N))
  return(d)
}

### 階層ベイズの推定
anhBay <- function(i, data, mu, imu, alpha, beta, final) {
  fith <<- stan_hierarchical_response_thall(
    
    #各がん種での奏効数
    group_responses <- c(sum(data$Y1, na.rm=TRUE),
                         sum(data$Y2, na.rm=TRUE),
                         sum(data$Y3, na.rm=TRUE),
                         sum(data$Y4, na.rm=TRUE),
                         sum(data$Y5, na.rm=TRUE)),
    
    #各がん種での例数
    group_sizes <- c(sum(!is.na(data$Y1)),
                     sum(!is.na(data$Y2)),
                     sum(!is.na(data$Y3)),
                     sum(!is.na(data$Y4)),
                     sum(!is.na(data$Y5))),
    mu_mean = mu,
    mu_sd = sqrt(1 / imu),
    tau_alpha = alpha,
    tau_beta = beta)
  
  fith2 <<- summary(fith)
  res_fith2 <- as.data.frame(fith2$summary)
  res_fith2 <- res_fith2 %>% mutate(paramn = row_number())
  post25 <- colMeans(abs(as.data.frame(fith, pars = 'prob_response')) > target_prop) #各がん種での奏効率が25%以上である事後確率を算出
  res_post2 <- data.frame(post25 = post25)
  res_post2$itr <- i
  res_post2$flag <- NA
  res_post2$flag[res_post2$post25 >= target_posterior_IA] <- 1　#奏効率が25%以上である事後確率が40%以上・未満に基づく無益性判断
  res_post2$flag[res_post2$post25 < target_posterior_IA] <- 0　 #奏効率が25%以上である事後確率が40%以上・未満に基づく無益性判断
  res_post2$final <- NA
  res_post2$final <- final
  if (i == 1 & final == 0) {res_fith <<- res_fith2} else {res_fith <<- rbind(res_fith, res_fith2)}
  if (i == 1 & final == 0) {res_post <<- res_post2} else {res_post <<- rbind(res_post, res_post2)}
  return(list(res_post, res_post2))
}

#中間から最終までの各がん種各症例でのレスポンスデータを発生させる関数
fdata <- function(i, N, p, flg) {
  set.seed(0827 + i)
  N_max <- max(N)
  y <- lapply(1:length(N), function(j) {
    if (flg[j] == 1) {
      c(rbinom(n = N[j], size = 1, prob = p[j]), rep(NA, N_max - N[j]))
    } else {
      rep(NA, N_max)
    }
  })
  d <- as.data.frame(y)
  names(d) <- paste0("Y", 1:length(N))
  return(d)
}



#がん種ごとの発現例数
N=ceiling(N_beforeIA*Ctype_prop)
data_final_list <- list()



#無益性中止に引っ掛からなかったデータのみを対象に中間から最終までに各がん種で追加される例数を各がん種の発現割合に基づき算出する関数。
f <- function(res_post2, N_add) {
  # flagが1のがん種だけcytpe_propをそのままにして
  # flagが0のがん種のcytpe_propは0にする
  d_remain <- res_post2 |> 
    # 発現割合の列を追加
    mutate(Ctype_prop = Ctype_prop) |> 
    mutate(Ctype_prop = if_else(flag == 1, Ctype_prop, 0))
  
  # すべて落ちている場合
  if (sum(d_remain$Ctype_prop) == 0) {
    return(data.frame(aN1 = 0, aN2 = 0, aN3 = 0, aN4 = 0, aN5 = 0))
  }
  # N_addをd_remain$pに応じて割り振る
  N_adds <- DoseFinding::rndDesign(d_remain$Ctype_prop, N_add)
  return(data.frame(aN1 = N_adds[1], aN2 = N_adds[2], aN3 = N_adds[3], aN4 = N_adds[4], aN5 = N_adds[5]))
}


for (i in 1:reprep) {
  data1 <- idata(i, N, p)　#各がん種・各症例での奏効の有無を発生させる
  data2 <- anhBay(i,data1, mu, imu, alpha, beta, 0)　#階層ベイズモデルによる推定（中間解析）
  data3 <- f(data2[[2]], N_afterIA) #中間から最終までに各がん種で追加される例数
  add_data <- fdata(i, N= c(data3$aN1, data3$aN2, data3$aN3, data3$aN4, data3$aN5), p, data2[[2]]$flag) #中間から最終までの各がん種各症例でのレスポンスデータ
  data_final <- rbind(data1, add_data)　#最終解析でのデータ
  data_final$itr <- i
  data_final_list[[i]] <- data_final
  bayes_all_result <- anhBay(i,data_final, mu, imu, alpha, beta, 1)　#階層ベイズモデルによる推定（最終解析）
}

#全データ
data_final_combined <- do.call(rbind, data_final_list) 

#中間解析と最終解析の結果を横並びにしたデータの作成
bayes_interim <- filter(bayes_all_result[[1]], final==0)
bayes_final <- filter(bayes_all_result[[1]], final==1)
bayes_interim <- data.frame(bayes_interim)
bayes_interim <-rename(bayes_interim, i_itr=itr, i_flag=flag, i_post25=post25) #Interim analysisの結果変数名にi_を付加
bayes_final <- data.frame(bayes_final)
bayes_final$flag <- NULL

#各がん種での成功判定基準: 今回の成果物ではPr(ORR > 25%) > 80%を追加
bayes_final$flag <- ifelse(bayes_final$post25 >= target_posterior_FA, 1, 0)
analysis<-cbind(bayes_interim, bayes_final)

#Proportion of Success
#1つ以上のがん種で成功と判定される確率
analysis$success_flag <- ifelse(analysis$i_flag == 1 & analysis$flag == 1, 1, 0)
grouped1 <- split(analysis$success_flag, analysis$itr)
group_success1 <- lapply(grouped1, function(x) any(x == 1))
prop_success <- sum(unlist(group_success1)) / length(group_success1)

##### Family-wise Error Rate (FWER) ####
#活性がないがん種で一つ以上偽陽性で有効性を主張してしまう確率

if (sum(p <= null_prop ) != 0){　#活性がないがん種がある場合
  #活性がないがん種名の取得
  rows_to_keep1 <- p <= null_prop 
  non_active    <- analysis[rows_to_keep1, ]

  #活性がないがん種での有効性の主張確率
  grouped2 <- split(non_active$success_flag, non_active$itr)
  group_success2 <- lapply(grouped2, function(x) any(x == 1))
  FWRE <- sum(unlist(group_success2)) / length(group_success2)
  
} else {
  #活性がないがん種がない場合はFWREはNA
  FWRE <- NA
}
  
#### Selection Power ####
#活性があるがん種のうち80%異常が最後の解析に含まれる確率

if (sum(p > null_prop ) != 0){
  #活性があるがん種名の取得
  rows_to_keep2 <- p > null_prop 
  active　<- analysis[rows_to_keep2, ]
  num_rows_itr1 <- nrow(active[active$itr == 1, ])
  
  count <-ceiling(num_rows_itr1*0.8)
  flag_counts <- table(active$itr[active$i_flag == 1])
  extracted_rows <- flag_counts[flag_counts >= count]
  selection_power<-nrow(extracted_rows)/reprep
} else {
  selection_power <- NA
}



#### Power ####
#主要解析の併合解析で有効性がいえる確率
new_data_final_list <- data_final_list
cols <- c("Y1", "Y2", "Y3", "Y4", "Y5")

#無益性中止となったがん種を主要解析から除外
for (i in seq_along(new_data_final_list)) {
  df <- new_data_final_list[[i]]
  max_row <- max(N) + 1
  #max_row <- max(N1, N2, N3, N4, N5) + 1
  df[, cols] <- lapply(df[, cols], function(x) {
    if (is.na(x[max_row])) NA else x
  })
  new_data_final_list[[i]] <- df
}
na.omit_counts <- lapply(new_data_final_list, function(df) colSums(!is.na(df[, 1:5])))
one_counts <- lapply(new_data_final_list, function(df) colSums(df[, 1:5] == 1, na.rm = TRUE))
sum_counts <- function(counts) {
  sum(counts[!is.na(counts)])
}
all_sum<-lapply(na.omit_counts, sum_counts)
all_sum2 <- unlist(all_sum)
all_sum2 <- ifelse(all_sum2 == 0, 1, all_sum2)
all_sum2 <- split(all_sum2, rep(1:length(all_sum), lengths(all_sum)))
one_sum<-lapply(one_counts, sum_counts)

run_binom_test <- function(x, n, p = null_prop, alternative = "two.sided") {
  binom.test(x, n, p, alternative = alternative)
}

binom_results <- lapply(seq_along(all_sum), function(i) {run_binom_test(one_sum[[i]], all_sum2[[i]])})
extract_pvalue <- function(binom_result) {binom_result$p.value}
p_values <- lapply(binom_results, extract_pvalue)
add_flag <- function(p_value) {
  flag <- ifelse(p_value < 0.05, 1, 0)
  return(flag)
}
flags <- lapply(p_values, add_flag)
power <- sum(unlist(flags)) / length(unlist(flags))

#### Mixing Nonactive ####
#最終解析で一つ以上活性がないがん種が含まれる中で、主要解析の併合解析で有効性がいえる確率
if (sum(p <= null_prop ) != 0){　#活性がないがん種がある場合
  levels <- factor(c(0, 1))
  flag_counts_non_active <- table(non_active$itr, factor(non_active$i_flag, levels = levels))
  flag_counts_non_active <- flag_counts_non_active[, 2] >= 1
  Mixing_Nonactive       <- sum(flag_counts_non_active * unlist(flags))/reprep
} else {
  Mixing_Nonactive <- NA
  }

sink("output_3_1B.txt")
cat("Summary of Design Setting:\n",
    "Ntotal:", N_total, "\n",
    "N at IA:", N_beforeIA, "\n",
    "True Response Rate for Each Cancer:", round(p, 2), "\n",
    "Prevalence for Each Cancer:", round(Ctype_prop, 2), "\n",
    "Null Response Rate:", round(null_prop, 2), "\n", 
    "BHM parameter: mu:", mu, " imu:", imu, " alpha:", alpha, " beta:", beta, "\n",
    "Target Response Rate for BHM:", round(target_prop, 2), "\n",
    "Target Posterior Probability for BHM at IA:", round(target_posterior_IA, 2), "\n",
    "Target Posterior Probability for BHM at FA:", round(target_posterior_FA, 2), "\n",
    "\n",
    "Summary of SimulationResult:\n",
    "Selection power:", round(selection_power, 2), "\n",
    "Proportion of Success:", round(prop_success, 2), "\n",
    "FWER:", round(FWRE, 2), "\n",
    "Power:", round(power, 2), "\n",
    "Mixing Nonactive:", round(Mixing_Nonactive , 2), "\n"
    )
sink() 
