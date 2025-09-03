library(OncoBayes2)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doRNG) # for seeds used OC 

set.seed(12345)

## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(OncoBayes2.MC.warmup=1000, OncoBayes2.MC.iter=5000, OncoBayes2.MC.chains=4,
                            OncoBayes2.MC.save_warmup=FALSE)


###################################################################
########### study design
###################################################################
group_ids = c("QD","BID") # list of all groups of historical and new studies
stratum_ids = c("single_stratum") # factors to separately model logistic model e.g. one model for each drug of combo therapy 
group_id  = factor("BID", group_ids)
stratum_id= factor("single_stratum",stratum_ids)


num_comp = 1 # one investigational drug
num_inter = 0 # no drug-drug interactions need to be modeled
num_groups = nlevels(group_id)
num_strata = nlevels(stratum_id)

hist_data = tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "QD",	10 ,3,0,0,
  "single_stratum", "QD",	20 ,3,0,0,
  "single_stratum", "QD",	40 ,3,0,0,
  "single_stratum", "QD",	80 ,6,1,0,
  
  # end of list without ","
)
hist_data = mutate(hist_data,
                   stratum_id = factor(stratum_id, levels=stratum_ids),
                   group_id = factor(group_id,levels=group_ids)
)

drug_info = tibble(
  drug_name = "drug_A", 
  dose_ref  = 80, # set reference dose d*
  dose_unit = "mg/day"  # dosing units
)

dose_info = tibble(
  stratum_id = factor("single_stratum", stratum_ids),
  group_id = factor(c("BID"), group_ids),
  drug_A = c(10, 20, 40, 80, 120, 160), # be sure to replace the variable name "drug_A" with the name of your study drug defined in `drug_info`
  arm_id = group_id
)

starting_doses = tibble(
  starting_dose = TRUE, # leave this as TRUE
  stratum_id = factor("single_stratum", stratum_ids),
  group_id = factor("BID", group_ids),
  drug_A = 20, # replace with your starting dose
  arm_id = group_id
)

dose_info = suppressMessages(replace_na(left_join(dose_info, starting_doses),list(starting_dose = FALSE)))


########## Parameter Setting for prior distributions #####
prior_mu_mean = c(logit(0.3), log(1)) #c(mean of intercept on logit scale, mean of log-slope on logit scale)
prior_mu_sd = c(2, 1) # c(sd of intercept, sd of log-slope)
prior_tau_mean = c(log(0.25), log(0.125))
prior_tau_sd =  c(log(2)/1.96, log(2)/1.96)
prior_prob = 0.8 
#prior_prob=1
#prior_prob=0

#for EXNEX ###
prior_glmu_mean = c(logit(0.3), log(1))
prior_glmu_sd = c(2,1)

###################################################################
########### BLRM setup
###################################################################
blrm_setup <- blrm_trial(
  data = hist_data,
  dose_info = dose_info,
  drug_info = drug_info,
  simplified_prior = FALSE,
  formula_generator = blrm_formula_saturating # functional form of interaction model (alternate is blrm_formula_linear); has no impact for single-agent models.
)


blrmfit <- update(
  blrm_setup,
  # formula=cbind(num_toxicities, num_patients - num_toxicities) ~
  #   1 + log(drug_A / drug_info$dose_ref) |
  #   0 |
  #   group_id,
  # data = hist_data
  # ,
  prior_EX_mu_mean_comp = matrix(
    prior_mu_mean,          # mean of intercept on logit scale
    nrow = num_comp,        # mean of log-slope on logit scale
    ncol = 2
  ),
  prior_EX_mu_sd_comp = matrix(
    prior_mu_sd,            # sd of intercept
    nrow = num_comp,        # sd of log-slope
    ncol = 2
  )
  ,
  
  prior_EX_tau_mean_comp = matrix(
    prior_tau_mean,
    nrow = num_comp,
    ncol = 2
  ),
  prior_EX_tau_sd_comp = matrix(
    prior_tau_sd,
    nrow = num_comp,
    ncol = 2
  ),
  prior_EX_corr_eta_comp = array(
    1,
    dim = c(num_strata, num_comp, 1)
  ),
  prior_tau_dist = 1, # log-normal priors for tau's

 prior_EX_prob_comp = matrix(
    prior_prob, 
    nrow = num_groups, 
    ncol = num_comp
    ),  # for exnex
  prior_is_EXNEX_comp = rep(TRUE,num_comp), # for exnex
  prior_NEX_mu_mean_comp = matrix( # for nex
      prior_glmu_mean, # expected value of global-mean intercept on logit scale, expected value of global-mean log-slope
      nrow = num_comp,
      ncol = 2
    ),
  prior_NEX_mu_sd_comp = matrix( # for nex
    prior_glmu_sd,              # sd of global-mean intercept on logit scale,  sd of global-mean log-slope
    nrow = num_comp,
    ncol = 2
  )
  # prior_PD = FALSE
)

###################################################################
########### check model and prior settings 
###################################################################

# # check model and prior parameter setting
# prior_summary(blrmfit)
# 
# ###################################################################
# ## Posterior based on historical data will be used as a prior for new study
# ###################################################################
# # check prior tox curve based on historical data 
# plot_toxicity_curve(blrmfit,x = "drug_A",group = ~ group_id)
# # plot_toxicity_intervals_stacked(blrmfit,x = "drug_A",group = ~ group_id,predictive = TRUE,
# #                                 interval_prob = c(-1, 0, 1, 6),
# #                                 num_patients = cohort_size,
# #                                 num_toxicities = 0
# # )
# # check posterior based on historical data
#print(blrmfit)
# # check all posterior distribution of all parameters
# print(blrmfit$blrmfit$stanfit)
# # summary of posterior for DLT rate by dose for new set of dose levels
# summary(blrmfit, "dose_prediction")[,-(5:6)]
# # another way for customised result


 # pred_range <- tibble(expand.grid(
 #   stratum_id = "single_stratum", group_id = c("BID"),
 #   drug_A = dose_info$drug_A,
 #   stringsAsFactors = FALSE
 # ))
 # summ_pred <- summary(blrmfit, newdata = pred_range, interval_prob = c(0, 0.16, 0.33, 1))
 # print(cbind(pred_range, summ_pred))

