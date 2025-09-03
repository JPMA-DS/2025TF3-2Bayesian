#----------------------------------------------------------------#
#Title : TS design operation characteristics
#Date  : 2024/11/2
#----------------------------------------------------------------#

boundary <- boundary # from "02_TSdesign_Boundary.R"
nsim <- 100000

TS_sim <- function(sim, p, boundary, n_min, n_max){
  success_FL <- 0
  endflg <- 0
  true_p <- p
  responder <- rbinom(1, n_min, true_p)
  act_num <- n_min

  if (responder >= boundary[n_min, ]$bound){
    success_FL <- 1
    endflg <-  1
  }else{
    while(all(endflg == 0, act_num < n_max)){
      x <- rbinom(1, 1, true_p)
      responder <- responder + x
      act_num <- act_num + 1
      if (responder >= boundary[act_num, ]$bound){
        success_FL <- 1
        endflg <- 1
        break()
      }else {
        endflg <- 0
      }
    }
  }
  tidyr::tibble(sim, true_p, act_num, responder, success_FL)
}

TS_OC <- function(nsim, theta, n_min, n_max, alpha_e, beta_e, alpha_s, beta_s){
  p_H0 <- alpha_s / (alpha_s + beta_s)
  p_H1 <- alpha_e / (alpha_e + beta_e)
  result_H0 <- purrr::map2(c(1:nsim), p_H0, TS_sim, boundary = boundary, n_min = n_min, n_max = n_max) %>% do.call(rbind,.)
  result_H1 <- purrr::map2(c(1:nsim), p_H1, TS_sim, boundary = boundary, n_min = n_min, n_max = n_max) %>% do.call(rbind,.)

  #H0
  ASS_H0    <- mean(result_H0$act_num)        #Average sample size
  alpha_err <- sum(result_H0$success_FL)/nsim #Type I error

  #H1
  ASS_H1 <- mean(result_H1$act_num)           #Average sample size
  power  <- sum(result_H1$success_FL)/nsim    #Statistical power

  tidyr::tibble(theta, n_min, n_max, alpha_e, beta_e, alpha_s, beta_s, ASS_H0, alpha_err, ASS_H1, power)
}

result <- TS_OC(nsim = nsim, theta = theta, n_min = n_min, n_max = n_max, alpha_e = alpha_e, beta_e = beta_e, alpha_s = alpha_s, beta_s = beta_s)
result %<>% knitr::kable()

#Save results
capture.output(result, file = "./output/03_TSdesign_OC.txt")
