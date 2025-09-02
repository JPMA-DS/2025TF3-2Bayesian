#----------------------------------------------------------------#
# Title : Effective historical sample size for Jags
#----------------------------------------------------------------#

# 1.Pre-processing
data.nb <- read.csv("data_CMP.csv") %>%
  as.data.frame() %>%
  filter(ext == 0) %>%
  arrange(cnsr)

# 2.Definition of the model (Jags code)
model_code.nb <- "
  model{
    for (i in 1:n.obs){
      time[i] ~ dweib(r, mu.obs[i])
      mu.obs[i] <- exp(beta0 + beta.trt * trt[i])
    }

    for (i in (n.obs + 1):n){
      cnsr[i] ~ dbern(S[i])
      S[i] <- 1 - pweib(time[i], r, mu.cens[i])
      mu.cens[i] <- exp(beta0 + beta.trt * trt[i])
    }

    #Priors
    beta0    ~ dnorm(0.0, 0.001)
    beta.trt ~ dnorm(0.0, 0.001)
    r        ~ dexp(1)
  }
"

# 3.Set up data for MCMC
model_data.nb <- list(
  n.obs = sum(data.nb[, "cnsr"] == 0), n = nrow(data.nb), time = data.nb[, "time"], trt = data.nb[, "trt"], cnsr = data.nb[, "cnsr"]
)

# 4.Choose parameter to monitor in MCMC
model_parameters.nb <- c("beta0",
                      "beta.trt",
                      "r")

# 5.Run the model
set.seed(2024)
model_run_ehss <- jags(
  data               = model_data.nb,
  parameters.to.save = model_parameters.nb,
  model.file         = textConnection(model_code.nb),
  n.chains           = 4,
  n.iter             = 100000,
  n.burnin           = 10000,
  n.thin             = 10
)

# 6.Diagnosis of convergence
png("./output/12_traceplot.png", width = 1500, height = 800)
R2jags::traceplot(model_run_ehss, varname = c("beta.trt"), ask = FALSE)
dev.off()

# 7.Summary of results
post <- print(model_run_ehss)
PostBetaTrt <- data.frame(betatrt = post$sims.matrix[, "beta.trt"])
betatrt_plot <- ggplot(PostBetaTrt, aes(x = betatrt)) + geom_density(adjust = 1)

# 8.Calculating EHSS
sd_ct <- model_run_ehss$BUGSoutput$summary %>% as.data.frame()

sd_ct %<>%
  filter(rownames(sd_ct) == "beta.trt") %>%
  select(sd) %>%
  as.numeric()

sd_all <- model_run_cmp$BUGSoutput$summary %>% as.data.frame()

sd_all %<>%
  filter(rownames(sd_all) == "beta.trt") %>%
  select(sd) %>%
  as.numeric()

ehss <- num_CT * ( (1/sd_all^2) / (1/sd_ct^2) - 1)

result.ehss <- list(
  model_run_ehss,
  ehss
)
names(result.ehss) <- c("MCMC summary", "EHSS")

ggsave("./output/12_betatrt_posterior_ehss.png", plot = betatrt_plot)
capture.output(result.ehss, file = "./output/12_model_run_ehss.txt")
