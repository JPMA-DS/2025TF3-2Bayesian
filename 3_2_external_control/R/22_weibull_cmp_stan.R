#----------------------------------------------------------------#
# Title : Bayesian dynamic approach (CMP approach) for stan
#----------------------------------------------------------------#

data_CMP <- read.csv("data_CMP.csv") %>% as.data.frame()

data <- list(
    N = nrow(data_CMP),
    time = data_CMP$time,
    cnsr = data_CMP$cnsr,
    trt = data_CMP$trt,
    ext = data_CMP$ext
)

fit1 <- stan(
  file    = "./R/21_weibull_cmp.stan",
  data    = data,   # A named list of data
  chains  = 4,      # The number of Markov chains
  warmup  = 10000,  # The number of warmup iterations per chain
  iter    = 100000, # Total number of iterations per chain
  thin    = 10,     # The period for saving samples
  cores   = 1,      # The number of cores
  refresh = 0,      # No progress shown
  seed    = 2024    # The seed for random number generation
)

capture.output(fit1, file = "./output/21_model_run_cmp_stan.txt")
