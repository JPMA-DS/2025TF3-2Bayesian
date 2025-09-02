#----------------------------------------------------------------#
# Title : Kaplan-Meier and PH estimation (No borrowing)
#----------------------------------------------------------------#

# Kaplan-Meier estimate
est_km.nb <- survfit(Surv(time, 1 - cnsr) ~ trt , data = data.nb)
kmplot.nb <- ggsurvplot(data       = data.nb,
                        fit        = est_km.nb,
                        size       = 1,
                        conf.int   = TRUE,
                        censor     = TRUE,
                        risk.table = TRUE,
                        ggtheme    = theme_bw()
)

# AFT model
res.aft.nb <- survreg(Surv(time, 1 - cnsr) ~ trt, dist = "weibull", data = data.nb)
res.aft.nb_sum <- summary(res.aft.nb)

# Hazard ratio scale
res.ph.nb <- exp(- res.aft.nb$coefficients / res.aft.nb$scale)[2]

# Save results
result.nb <- list(
  est_km.nb,
  res.aft.nb_sum,
  res.ph.nb
)
names(result.nb) <- c("Kaplan-Meier estimate (No borrowing)",
                     "AFT model (No borrowing)",
                     "Hazard ratio scale (No borrowing)")

capture.output(result.nb, file = "./output/02_result_nb.txt")
ggsave("./output/02_kmplot_nb.png", plot = survminer:::.build_ggsurvplot(kmplot.nb))
