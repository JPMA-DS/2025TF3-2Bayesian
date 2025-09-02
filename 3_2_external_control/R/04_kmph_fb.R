#----------------------------------------------------------------#
# Title : Kaplan-Meier and PH estimation (Full borrowing)
#----------------------------------------------------------------#

# Kaplan-Meier estimate and plot
est_km.fb <- survfit(Surv(time, 1 - cnsr) ~ trt , data = data)
kmplot.fb <- ggsurvplot(data       = data,
                        fit        = est_km.fb,
                        size       = 1,
                        conf.int   = TRUE,
                        censor     = TRUE,
                        risk.table = TRUE,
                        ggtheme    = theme_bw()
)

# AFT model
res.aft.fb <- survreg(Surv(time, 1 - cnsr) ~ trt, dist = "weibull", data = data)
res.aft.fb_sum <- summary(res.aft.fb)

# Hazard ratio scale
res.ph.fb <- exp(- res.aft.fb$coefficients / res.aft.fb$scale)[2]

# Save results
result.fb <- list(
  est_km.fb,
  res.aft.fb_sum,
  res.ph.fb
)
names(result.fb) <- c("Kaplan-Meier estimate (Full borrowing)",
                     "AFT model (Full borrowing)",
                     "Hazard ratio scale (Full borrowing)")

capture.output(result.fb, file = "./output/03_result_fb.txt")
ggsave("./output/03_kmplot_fb.png", plot = survminer:::.build_ggsurvplot(kmplot.fb))
