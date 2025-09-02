library(rstan)
library(rstudioapi)
library(ggplot2)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
setwd("C:\\Rwork")

# Simulation setting 

seedvalue <- 1001   # Seed value for random numbers
iter <- 10          # Number of iterations
numia <- 1          # Flag for interim analysis (IA) (0: without IA, 1: with IA)
samplesize_fa <- 35 # Number of subjects per arm in final analysis (FA)
timing_ia <- 0.5    # Timing of IA (percentage of subjects)
stddev <- 2.3       # Standard deviation of outcome data
scenario <- "2a"    # Dose response scenario

# Pre-processing

set.seed(seedvalue)

if (scenario == "1" ){truemean <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00)}
if (scenario == "2a"){truemean <- c(0.00, 0.12, 0.24, 0.36, 0.48, 0.73, 1.09, 1.45)}
if (scenario == "2b"){truemean <- c(0.00, 0.95, 1.17, 1.27, 1.32, 1.38, 1.43, 1.45)}
if (scenario == "2c"){truemean <- c(0.00, 0.92, 1.21, 1.36, 1.43, 1.42, 1.12, 0.57)}
if (scenario == "2d"){truemean <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 1.02, 1.45)}
if (scenario == "3a"){truemean <- c(0.00, 0.18, 0.37, 0.55, 0.73, 1.10, 1.65, 2.20)}
if (scenario == "3b"){truemean <- c(0.00, 1.44, 1.78, 1.93, 2.01, 2.10, 2.17, 2.20)}
if (scenario == "3c"){truemean <- c(0.00, 1.40, 1.83, 2.07, 2.18, 2.15, 1.69, 0.87)}
if (scenario == "3d"){truemean <- c(0.00, 0.00, 0.00, 0.00, 0.00, 0.08, 1.54, 2.20)}
if (scenario == "4a"){truemean <- c(0.00, 0.24, 0.48, 0.73, 0.97, 1.45, 2.18, 2.90)}
if (scenario == "4b"){truemean <- c(0.00, 1.90, 2.34, 2.54, 2.65, 2.77, 2.85, 2.90)}
if (scenario == "4c"){truemean <- c(0.00, 1.85, 2.42, 2.72, 2.87, 2.84, 2.23, 1.15)}
if (scenario == "4d"){truemean <- c(0.00, 0.00, 0.00, 0.00, 0.01, 0.10, 2.04, 2.90)}

# Stan code for Bayesian hierarchical model

stanmodelcode <- "
data {
  int<lower=0> N[8];
  vector[N[1]] y0;
  vector[N[2]] y1;
  vector[N[3]] y2;
  vector[N[4]] y3;
  vector[N[5]] y4;
  vector[N[6]] y5;
  vector[N[7]] y6;
  vector[N[8]] y7;
}

parameters {
  real theta0;
  real mu;
  vector[7] theta;
  real<lower=0> sigma;
  real<lower=0> tau;
}

model {
  theta ~ normal(mu, tau);
  y0 ~ normal(theta0, sigma);
  y1 ~ normal(theta[1], sigma);
  y2 ~ normal(theta[2], sigma);
  y3 ~ normal(theta[3], sigma);
  y4 ~ normal(theta[4], sigma);
  y5 ~ normal(theta[5], sigma);
  y6 ~ normal(theta[6], sigma);
  y7 ~ normal(theta[7], sigma);
}
"

stanmodel1 <- stan_model(model_code = stanmodelcode)

# Simulation iterations

result_estimate <- NULL
result_probsuccess <- NULL
result_probmaxeffect <- NULL
result_probfutility <- NULL
result_ss <- NULL

for (i in 1:iter){

  estimate <- numeric(8)
  probsuccess <- numeric(7)
  probmaxeffect <- numeric(7)
  probfutility <- numeric(7)
  
  # Set number of subjects in FA
  n_fa <- rep(samplesize_fa,8)
  
  # Data generation
  y0 <- rnorm(n_fa[1], truemean[1], stddev)
  y1 <- rnorm(n_fa[2], truemean[2], stddev)
  y2 <- rnorm(n_fa[3], truemean[3], stddev)
  y3 <- rnorm(n_fa[4], truemean[4], stddev)
  y4 <- rnorm(n_fa[5], truemean[5], stddev)
  y5 <- rnorm(n_fa[6], truemean[6], stddev)
  y6 <- rnorm(n_fa[7], truemean[7], stddev)
  y7 <- rnorm(n_fa[8], truemean[8], stddev)

  # Interim analysis
  if (numia == 1){
    # Set number of subjects in IA 
    n_ia <-rep(floor(samplesize_fa * timing_ia), 8)
    n_ia[sample(c(1:8), 4, replace=FALSE)] <- ceiling(samplesize_fa * timing_ia)
    
    # MCMC sampling 
    data_list_ia <- list(N = n_ia, y0 = y0[1:n_ia[1]], y1 = y1[1:n_ia[2]], y2 = y2[1:n_ia[3]], y3 = y3[1:n_ia[4]], y4 = y4[1:n_ia[5]], y5 = y5[1:n_ia[6]], y6 = y6[1:n_ia[7]], y7 = y7[1:n_ia[8]])
    fit0 <- sampling(stanmodel1, data = data_list_ia, seed = seedvalue, chains = 2, iter = 3000, warmup = 1000, thin =1)
  
    mcmcsample <- rstan::extract(fit0, permuted = FALSE)
  
    theta_ia <- matrix(c(
      as.vector(mcmcsample[, , "theta0"]),
      as.vector(mcmcsample[, , "theta[1]"]),
      as.vector(mcmcsample[, , "theta[2]"]),
      as.vector(mcmcsample[, , "theta[3]"]),
      as.vector(mcmcsample[, , "theta[4]"]),
      as.vector(mcmcsample[, , "theta[5]"]),
      as.vector(mcmcsample[, , "theta[6]"]),
      as.vector(mcmcsample[, , "theta[7]"])
    ), nrow = 4000, ncol = 8)
  
    # Calculate probability of futility and update number of subjects in FA
    for (j  in 1:7){
      probfutility[j] <- mean(ifelse(theta_ia[, j + 1] - theta_ia[, 1] < 1.5, 1, 0))
      if (probfutility[j] > 0.8){ n_fa[j + 1] <- n_ia[j + 1] }
    }
    if (max(n_fa[2:8]) < samplesize_fa){ n_fa[1] <- n_ia[1] }
  }

  # Final analysis

  # MCMC sampling
  data_list_fa <- list(N = n_fa, y0 = y0[1:n_fa[1]], y1 = y1[1:n_fa[2]], y2 = y2[1:n_fa[3]], y3 = y3[1:n_fa[4]], y4 = y4[1:n_fa[5]], y5 = y5[1:n_fa[6]], y6 = y6[1:n_fa[7]], y7 = y7[1:n_fa[8]])
  fit1 <- sampling(stanmodel1, data = data_list_fa, seed = seedvalue, chains = 2, iter = 3000, warmup = 1000, thin =1)

  mcmcsample <- rstan::extract(fit1, permuted = FALSE)

  theta <- matrix(c(
    as.vector(mcmcsample[, , "theta0"]),
    as.vector(mcmcsample[, , "theta[1]"]),
    as.vector(mcmcsample[, , "theta[2]"]),
    as.vector(mcmcsample[, , "theta[3]"]),
    as.vector(mcmcsample[, , "theta[4]"]),
    as.vector(mcmcsample[, , "theta[5]"]),
    as.vector(mcmcsample[, , "theta[6]"]),
    as.vector(mcmcsample[, , "theta[7]"])
  ), nrow = 4000, ncol = 8)

  # Calculate posterior mean, probability of efficacy and probability of maximum effect
  maxdiff <- apply(theta[,2:8], 1, max) - theta[, 1]

  for (j  in 1:8){
    estimate[j] <- mean(theta[,j])
    if (j!=8){
      probsuccess[j] <- mean(ifelse(theta[, j + 1] - theta[, 1] > 1.5, 1, 0))
      if (probfutility[j] > 0.8) { probsuccess[j] <- 0 }
      probmaxeffect[j] <- mean(ifelse(theta[, j + 1] - theta[, 1] == maxdiff, 1, 0))
    }
  }

  result_estimate <- rbind(result_estimate, t(estimate))
  result_probsuccess <- rbind(result_probsuccess, t(probsuccess))
  result_probmaxeffect <- rbind(result_probmaxeffect, t(probmaxeffect))
  result_probfutility <- rbind(result_probfutility, t(probfutility))
  result_ss <-rbind(result_ss, t(n_fa))

}

# Calculate bias and squared error   
result_bias <- (result_estimate[,2:8] - result_estimate[,1]) - rep(1,iter) %*% t(truemean[2:8])
result_mse <- result_bias^2

# Summarize simulation results

meanestimate <- apply(result_estimate, 2, mean) # Mean estimate (posterior mean)
quarestimate <- apply(result_estimate, 2, quantile, c(0.25,0.5,0.75)) # Quantile estimate (posterior mean)
meanbias <- apply(result_bias, 2, mean) # Mean bias
rmse <- sqrt(apply(result_mse, 2, mean)) # Root mean squared error (RMSE)
power <- apply(ifelse(result_probsuccess > 0.8, 1, 0),2,mean) # Power for each dose
poweroverall <- mean(ifelse(apply(result_probsuccess, 1, max) > 0.8, 1, 0))  # Power for overall (at least one dose show efficacy)
selectmax <- mean(ifelse(apply(result_estimate[,2:8], 1, which.max) == which.max(truemean[2:8]), 1, 0)) # Probability of maximum effect 
meanss <- apply(result_ss, 2, mean) # Mean sample size

#Output

filenam1 <- sprintf("BHM_simulation_est_%s_%d.csv", scenario, numia)
filenam2 <- sprintf("BHM_simulation_psuccess_%s_%d.csv", scenario, numia)
filenam3 <- sprintf("BHM_simulation_pmaxeff_%s_%d.csv", scenario, numia)
filenam4 <- sprintf("BHM_simulation_ss_%s_%d.csv", scenario, numia)
filenam5 <- sprintf("BHM_simulation_pfutility_%s_%d.csv", scenario, numia)
filenam6 <- sprintf("BHM_simulation_out_%s_%d.txt", scenario, numia)

write.csv(result_estimate, filenam1)
write.csv(result_probsuccess, filenam2)
write.csv(result_probmaxeffect, filenam3)
write.csv(result_ss, filenam4)
write.csv(result_probfutility, filenam5)

cat("scenario:\n", scenario, "\n", file = filenam6, append = FALSE) 
cat("ture mean:\n", noquote(paste(round(truemean, 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("estimate mean:\n", noquote(paste(round(meanestimate, 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("estimate 25%:\n", noquote(paste(round(quarestimate[1,], 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("estimate 50%:\n", noquote(paste(round(quarestimate[2,], 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("estimate 75%:\n", noquote(paste(round(quarestimate[3,], 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("bias mean:\n", noquote(paste(round(meanbias, 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("rmse:\n", noquote(paste(round(rmse, 3), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("power:\n", noquote(paste(round(power, 2), collapse = ', ')), "\n", file = filenam6, append = TRUE)
cat("overall power:\n", round(poweroverall, 2), "\n", file = filenam6, append = TRUE)
cat("selectmax:\n", round(selectmax,2), "\n", file = filenam6, append = TRUE)
cat("expected sample size:\n", noquote(paste(round(meanss, 1), collapse = ', ')), "\n", file = filenam6, append = TRUE)

#Graph

dose <- c(0,50,100,150,200,300,450,600)
pctss <-meanss/sum(meanss)
DAT1  <- data.frame(dose, truemean, meanestimate, p25 = quarestimate[1,], p50 = quarestimate[2,], p75 = quarestimate[3,], pctss)
DAT2 <- data.frame(dose[2:8], meanbias, rmse, power)

ggplot(data = DAT1) +
  geom_point(aes(x = dose, y = truemean), shape = 5, size = 3) + 
  geom_line(aes(x = dose, y = p25), linetype = 2) +
  geom_line(aes(x = dose, y = p50), linetype = 1, linewidth = 1) +
  geom_line(aes(x = dose, y = p75), linetype = 2) +
  scale_x_continuous(breaks = seq(0,600,50)) +
  labs(x = "Dose", y = "Change from Baseline")
ggsave(sprintf("BHM_simulation_fig1_%s_%d.png", scenario, numia), width = 16, height = 9)

ggplot(data = DAT2) +
  geom_bar(aes(x = factor(dose.2.8.), y = meanbias), stat = "identity") + 
  labs(x = "Dose", y = "Mean Bias")
ggsave(sprintf("BHM_simulation_fig2_%s_%d.png", scenario, numia), width = 16, height = 9)

ggplot(data = DAT2) +
  geom_bar(aes(x = factor(dose.2.8.), y = rmse), stat = "identity") + 
  labs(x = "Dose", y = "Root Mean Squared Error")
ggsave(sprintf("BHM_simulation_fig3_%s_%d.png", scenario, numia), width = 16, height = 9)

ggplot(data = DAT2) +
  geom_bar(aes(x = factor(dose.2.8.), y = power), stat = "identity") + 
  labs(x = "Dose", y = "Power")
ggsave(sprintf("BHM_simulation_fig4_%s_%d.png", scenario, numia), width = 16, height = 9)

ggplot(data = DAT1) +
  geom_bar(aes(x = factor(dose), y = pctss), stat = "identity") + 
  labs(x = "Dose", y = "Allocation Proportion")
ggsave(sprintf("BHM_simulation_fig5_%s_%d.png", scenario, numia), width = 16, height = 9)