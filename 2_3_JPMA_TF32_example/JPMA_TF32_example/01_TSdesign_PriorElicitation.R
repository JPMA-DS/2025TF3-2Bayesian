#----------------------------------------------------------------#
#Title : TS design prior elicitation
#Date  : 2024/11/2
#----------------------------------------------------------------#

alpha <- seq(1, 100, 0.1)
beta <- alpha * 4
q5 <- q95 <- interval <- diff <- numeric(length(alpha))
for (i in 1:length(alpha)){
  q5[i]  <- qbeta(p = 0.05, shape1 = alpha[i], shape2 = beta[i], lower.tail = TRUE, log.p = FALSE)
  q95[i] <- qbeta(p = 0.95, shape1 = alpha[i], shape2 = beta[i], lower.tail = TRUE, log.p = FALSE)
  interval[i] <- q95[i] - q5[i]
}
diff <- abs(interval - 0.10)
dat <- data.frame(alpha, beta, interval, diff)
dat %<>% dplyr::filter(dat$diff == min(abs(dat$diff)))

#Save results
capture.output(dat, file = "./output/01_TSdesign_Prior.txt")
