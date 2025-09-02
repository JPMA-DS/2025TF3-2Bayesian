#----------------------------------------------------------------#
# Title : Kaplan-Meier and PH estimation (HC vs. CC)
#----------------------------------------------------------------#

# Kaplan-Meier estimate and plot
est_km.C <- survfit(Surv(time, 1 - cnsr) ~ ext , data = data.C)
kmplot.C <- ggsurvplot(data       = data.C,
                       fit        = est_km.C,
                       size       = 1,
                       conf.int   = TRUE,
                       censor     = TRUE,
                       risk.table = TRUE,
                       ggtheme    = theme_bw()
)

# AFT model
res.aft.C <- survreg(Surv(time, 1 - cnsr) ~ ext, dist = "weibull", data = data.C)
res.aft.C_sum <- summary(res.aft.C)

# Hazard ratio scale
res.ph.C <- exp(- res.aft.C$coefficients / res.aft.C$scale)[2]

# Save results
result.C <- list(
  est_km.C,
  res.aft.C_sum,
  res.ph.C
)
names(result.C) <- c("Kaplan-Meier estimate (CC vs. HC)",
                   "AFT model (CC vs. HC)",
                   "Hazard ratio scale (CC vs. HC)")

capture.output(result.C, file = "./output/01_result_c.txt")
ggsave("./output/01_kmplot_c.png", plot = survminer:::.build_ggsurvplot(kmplot.C))
