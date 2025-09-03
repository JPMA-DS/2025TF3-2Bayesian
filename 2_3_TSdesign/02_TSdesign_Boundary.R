#----------------------------------------------------------------#
#Title : TS design calculating boundary for case example
#Date  : 2024/11/2
#----------------------------------------------------------------#

n_min <- 10
n_max <- 15
alpha_e <- 0.5
beta_e <- 0.5
theta <- 0.95
alpha_s <- dat$alpha # from "01_TSdesign_PriorElicitation.R"
beta_s <- dat$beta   # see above

boundary_ <- ph2bayes::stopbound_post(theta     = theta,
                                      type      = "superiority",
                                      nmax      = n_max,
                                      alpha_e   = alpha_e,
                                      beta_e    = beta_e,
                                      alpha_s   = alpha_s,
                                      beta_s    = beta_s,
                                      delta     = 0)

boundary <- dplyr::full_join(data.frame(n = c(1:n_max)), boundary_, by = "n") %>%
  tidyr::fill(bound)

capture.output(boundary, file = "./output/02_TSdesign_Boundary.txt")

