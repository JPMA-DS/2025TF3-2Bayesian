#----------------------------------------------------------------#
#Title : Design prior and Analysis prior (simple example under H0)
#Date  : 2024/11/2
#----------------------------------------------------------------#

# Design parameters
N <- 60             #Total sample size
a <- b <- 0.5       #Hyper parameter (Jeffery's prior)
SOC <- 0.2          #ORR of SOC setting
futility <- 0.95    #Futility boundary

# Cumulative distribution for posterior distribution
post_prob <- function(crit, alpha, beta, num, resp){
  prob <- 1 - pbeta(crit, alpha + resp, beta + num - resp)
  return(prob)
}

# Calculating threshold
prob_thresh <- function(x){
  responder <- sprintf("%d", x)
  orr <- sprintf("%.2f", 100 * (x / N))
  prob_SOC <- (100 * post_prob(crit = SOC, alpha = a, beta = b, num = N, resp = x)) %>% sprintf("%.2f", .)
  result <- data.frame(responder, orr, prob_SOC)
  return(result)
}

thresh <- lapply(0:N, prob_thresh) %>% do.call(rbind,.) %>%
  dplyr::filter(prob_SOC > (100 * futility)) %>%
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::select(responder) %>%
  as.numeric()

# Errorness probability
p <- SOC
x <- 0:N

## Ddirac distribution of design prior
prob <- dbinom(x, N, p)
names(prob) <- x

df_dirac <- data.frame(x, prob) %>%
  dplyr::mutate(Decision = ifelse(x >= thresh, "Go", "No-go")) %>%
  dplyr::mutate(Decision = factor(Decision , levels = c("No-go", "Go")))

plot_dirac <- ggplot(df_dirac, aes(x = x, y = prob, fill = Decision)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  coord_cartesian(xlim = c(0, 55), ylim = c(0, 0.13)) +
  xlab("Number of responders") +
  ylab("Probability") +
  ggtitle("Dirac design prior")

prob_go <- df_dirac %>%
  dplyr::filter(Decision == "Go")

prob_go <- sum(prob_go$prob)
prob_nogo <- 1 - prob_go

result_dirac <- data.frame(prob_go, prob_nogo) %>% dplyr::mutate(design = "Dirac")

## Informative Beta distribution of design prior
prob_infbeta <- extraDistr::dbbinom(x, N, alpha = 20, beta = 80, log = FALSE)

df_infbeta <- data.frame(x, prob_infbeta) %>%
  dplyr::mutate(Decision = ifelse(x >= thresh, "Go", "No-go")) %>%
  dplyr::mutate(Decision = factor(Decision , levels = c("No-go", "Go")))

plot_infbeta <- ggplot(df_infbeta, aes(x = x, y =  prob_infbeta, fill = Decision)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  coord_cartesian(xlim = c(0, 55), ylim = c(0, 0.13)) +
  xlab("Number of responders") +
  ylab("Probability") +
  ggtitle("Informative beta design prior")

prob_go <- df_infbeta %>%
  dplyr::filter(Decision == "Go")

prob_go <- sum(prob_go$prob_infbeta)
prob_nogo <- 1 - prob_go

result_infbeta <- data.frame(prob_go, prob_nogo) %>% dplyr::mutate(design = "Informative Beta")

## Vague Beta distribution of design prior
prob_vagbeta <- extraDistr::dbbinom(x, N, alpha = 2, beta = 8, log = FALSE)

df_vagbeta <- data.frame(x, prob_vagbeta) %>%
  dplyr::mutate(Decision = ifelse(x >= thresh, "Go", "No-go")) %>%
  dplyr::mutate(Decision = factor(Decision , levels = c("No-go", "Go")))

plot_vagbeta <- ggplot(df_vagbeta, aes(x = x, y =  prob_vagbeta, fill = Decision)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  coord_cartesian(xlim = c(0, 55), ylim = c(0, 0.13)) +
  xlab("Number of responders") +
  ylab("Probability") +
  ggtitle("Vague beta design prior")

prob_go <- df_vagbeta %>%
  dplyr::filter(Decision == "Go")

prob_go <- sum(prob_go$prob_vagbeta)
prob_nogo <- 1 - prob_go

result_vagbeta <- data.frame(prob_go, prob_nogo) %>% dplyr::mutate(design = "Vague Beta")

result <- rbind(result_dirac, result_infbeta, result_vagbeta)

# Save results
ggsave("./output/04_plot_dirac_H0.png", plot = plot_dirac)
ggsave("./output/05_plot_infbeta_H0.png", plot = plot_infbeta)
ggsave("./output/06_plot_vagbeta_H0.png", plot = plot_vagbeta)

capture.output(result, file = "./output/07_probability_H0.txt")
