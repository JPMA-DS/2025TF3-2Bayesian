#----------------------------------------------------------------#
# Title : Bayesian dynamic approach (CMP approach) for Jags
#----------------------------------------------------------------#

# 1.Pre-processing
data <- read.csv("data_CMP.csv") %>%
  as.data.frame() %>%
  arrange(cnsr)

# 2.Definition of the model (Jags code)
model_code <- "
  model{
    for (i in 1:n.obs){
      time[i] ~ dweib(r, mu.obs[i])
      mu.obs[i] <- exp(beta0 * (1 - ext[i]) + beta.ext * ext[i] + beta.trt * trt[i])
    }

    for (i in (n.obs + 1):n){
      cnsr[i] ~ dbern(S[i])
      S[i] <- 1 - pweib(time[i], r, mu.cens[i])
      mu.cens[i] <- exp(beta0 * (1 - ext[i]) + beta.ext * ext[i] + beta.trt * trt[i])
    }

    #Priors
    beta.ext ~ dnorm(0.0, 0.001)
    beta0    ~ dnorm(beta.ext, tau)
    beta.trt ~ dnorm(0.0, 0.001)
    tau      ~ dgamma(1, 0.001)
    r        ~ dexp(1)
  }
"

# 3.Set up data for MCMC
model_data <- list(
  n.obs = sum(data[, "cnsr"] == 0), n = nrow(data), time = data[, "time"], trt = data[, "trt"], cnsr = data[, "cnsr"],
  ext = data[, "ext"]
)

# 4.Choose parameter to monitor in MCMC
model_parameters <- c("beta0",
                      "beta.trt",
                      "beta.ext",
                      "tau",
                      "r")

# 5.Run the model
set.seed(2024)
model_run_cmp <- jags(
  data               = model_data,
  parameters.to.save = model_parameters,
  model.file         = textConnection(model_code),
  n.chains           = 4,
  n.iter             = 100000,
  n.burnin           = 10000,
  n.thin             = 10
)

# 6.Diagnosis of convergence
png("./output/11_traceplot_cmp.png", width = 1500, height = 800)
R2jags::traceplot(model_run_cmp, varname = c("beta.trt"), ask = FALSE)
dev.off()

# 7.Summary of results
post <- print(model_run_cmp)
PostBetaTrt <- data.frame(betatrt = post$sims.matrix[, "beta.trt"])
betatrt_plot <- ggplot(PostBetaTrt, aes(x = betatrt)) + geom_density(adjust = 1)

ggsave("./output/11_betatrt_posterior_cmp.png", plot = betatrt_plot)
capture.output(model_run_cmp, file = "./output/11_model_run_cmp.txt")

