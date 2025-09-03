library(here) # to get path where R project is located
### echo = NULL if plot for prior curve is not needed.
set.seed(12345)
source(here::here("BHLRM_exnex/bhlrm_setup.R"), echo = TRUE)  


demo_scenario.fun=function(demo_study_data){
set.seed(12345)
curr_dose=demo_study_data$drug_A[nrow(demo_study_data)]  
  
## update BLRM
blrmfit_demo <- update(blrmfit, data = demo_study_data)
## prediction 
pred_range <- tibble(expand.grid(
  stratum_id = "single_stratum", group_id = c("BID"),
  drug_A = dose_info$drug_A, # available dose levels for next
  stringsAsFactors = FALSE
))
summ_pred <- summary(blrmfit_demo, newdata = pred_range, interval_prob = c(0, 0.16, 0.33, 1))


#plot_toxicity_curve(blrmfit_demo,x = "drug_A",group = ~ group_id)
# p = plot_toxicity_curve(blrmfit_demo,x = "drug_A",group = ~ group_id)
# p + geom_line(aes(x=drug_A, y=prob1), data=senarios)

# check all posterior distribution of all parameters
# summary(blrmfit_demo, "dose_prediction")

#print(cbind(pred_range, summ_pred))

if(any(summ_pred$ewoc_ok)){ND=summ_pred[max(which(summ_pred$ewoc_ok)),c(1,2,3,16,17)]}
else{ND=NULL}

return(
  list(Data=demo_study_data,
       CD=summ_pred[which(summ_pred$drug_A==curr_dose),c(1,2,3,16,17)],
       ND=ND)
)
}


### new data ###

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,0,1,
)
S1=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,1,1,
)
S2=demo_scenario.fun(demo_study_data)

##############
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,2,1,
)
S3=demo_scenario.fun(demo_study_data)

###############

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,0,1,
)
S4=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,1,1,
)
S5=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,5,1,1,
)
S6=demo_scenario.fun(demo_study_data)

##############
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,6,1,1,
)
S7=demo_scenario.fun(demo_study_data)

##############
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,6,2,1,
)
S8=demo_scenario.fun(demo_study_data)



sink("Out_demo1.txt")
print(S1)
cat("\n")
print(S2)
cat("\n")
print(S3)
cat("\n")
print(S4)
cat("\n")
print(S5)
cat("\n")
print(S6)
cat("\n")
print(S7)
cat("\n")
print(S8)
sink()

####################################


### new data ###

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,0,1,
  "single_stratum", "BID",	40 ,3,0,2,
)
S9=demo_scenario.fun(demo_study_data)

############

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,0,1,
  "single_stratum", "BID",	40 ,3,1,2,
)
S10=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,1,1,
  "single_stratum", "BID",	10 ,3,0,2,
)
S11=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,1,1,
  "single_stratum", "BID",	10 ,3,1,2,
)
S12=demo_scenario.fun(demo_study_data)

###############

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,0,1,
  "single_stratum", "BID",	40 ,3,0,2,
)
S13=demo_scenario.fun(demo_study_data)

###############

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,0,1,
  "single_stratum", "BID",	40 ,3,1,2,
)
S14=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,1,1,
  "single_stratum", "BID",	20 ,3,0,2,
)
S15=demo_scenario.fun(demo_study_data)

#########
demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,4,1,1,
  "single_stratum", "BID",	20 ,3,1,2,
)
S16=demo_scenario.fun(demo_study_data)



sink("Out_demo2.txt")
print(S9)
cat("\n")
print(S10)
cat("\n")
print(S11)
cat("\n")
print(S12)
cat("\n")
print(S13)
cat("\n")
print(S14)
cat("\n")
print(S15)
cat("\n")
print(S16)
sink()

####################################

set.seed(12345)
source(here::here("BHLRM_exnex/bhlrm_setup.R"), echo = TRUE) 

demo_study_data <- tribble(
  ~stratum_id,~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time,
  "single_stratum", "BID",	20 ,3,1,1,
  "single_stratum", "BID",	10 ,3,0,2,
  "single_stratum", "BID",	20 ,3,1,3,
)
S_=demo_scenario.fun(demo_study_data)