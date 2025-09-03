library(here) # to get path where R project is located
### echo = FALSE if plot for prior curve is not needed.
source(here::here("BHLRM_exnex/bhlrm_setup.R"), echo = FALSE)


########## Cohort Setting #####

max_percent_inc = 2 # percentage of maximum increment for next dose 1 = +100% increment is ok

max_patients= 36
min_patients = 9

# seeds=sample(1:10000000, 1)
# print(seeds)

### example 1 fixed cohort size
# cohort_sizes=3
# cohort_sizes_probs=rep(1/length(cohort_sizes),length(cohort_sizes)) ## user-defined probability can be defined here

### example 2 cohort size = 3 to 6 with equal probability
cohort_sizes=3:6 # cohort_sizes=3 if cohort size is fixed at 3
cohort_sizes_probs=rep(1/length(cohort_sizes),length(cohort_sizes)) ## user-defined probability can be defined here


### example 3
## cohort_size =3 if specific cohort size is prefered for a dose level which is never treated
## cohort_size =2 otherwise
# cohort_size_first=3
# cohort_size_add=2

### example 4
## safe start with fixed size
# cohort_size_initial=1
# cohort_size_others=3


patients_on_dose_for_MTD=6
prob_target_cond=0.5

  

########## Scenario Setting #####
trueDLTs1=c(
  0.05,
  0.1,
  0.2,
  0.3,
  0.4,
  0.6
)

trueDLTs2=c(
  0.33,
  0.4,
  0.5,
  0.6,
  0.7,
  0.8
)

trueDLTs3=c(
  0.025,
  0.05,
  0.1,
  0.15,
  0.2,
  0.3
)

trueDLTs=trueDLTs2

####################################


if (length(dose_info$drug_A)!=length(trueDLTs)) print("error")
scenarios=tibble(
  drug_A=dose_info$drug_A,
  prob1=trueDLTs
)


###################################################################
########### begin trial
########### simulation start
###################################################################
bhlrm_oc <- function(scenarios,seeds){
  ## initialization of iteration
  decision=NULL
  cohort_time=0
  MTD=NA
  current_dose=starting_doses$drug_A
  cumulative_n = 0
  available_patients=max_patients
  
  ## argmument blrmfit and simulated_study_data 
  simulated_study_data <- tribble(
    ~stratum_id,         ~group_id, ~drug_A, ~num_patients, ~num_toxicities, ~cohort_time, ~seed
  )
  
  ################## iteration start ###########
  while(is.null(decision)){
    next_dose=NULL
    cohort_time = cohort_time + 1
    
    # "L'Ecuyer-CMRG" generator needs a set of 6 random numbers
    # "seeds" is a list of sets of 6 random numbers.
    current_seed = seeds[[cohort_time]]
    .Random.seed=current_seed
    
    ################# START: cohort size (update the following if needed) ############
    #### for example 1 and 2
    if (available_patients<max(cohort_sizes)) {
      cohort_size = available_patients
    } else if (length(cohort_sizes)==1) {
      cohort_size = cohort_sizes
    } else if (max(cohort_sizes)<=available_patients) {
      cohort_size = sample(cohort_sizes, 1, replace = TRUE, prob = cohort_sizes_probs)
    }
    
    ### for example 3
    ## if the dose level was not treated for patients 
    # if (nrow(filter(simulated_study_data, drug_A == current_dose , num_patients > 0)) == 0) {
    #   if (available_patients<=cohort_size_first) cohort_size = available_patients 
    #   else cohort_size = cohort_size_first
    # } else {
    #   if (available_patients<=cohort_size_add) cohort_size = available_patients
    #   else cohort_size = cohort_size_add
    # } 
    
    ### for example 4
    ## safe start with starting dose  
    # if (cohort_time <=2) {
    #   cohort_size = cohort_size_initial
    #   current_dose=starting_doses$drug_A # overwrite current_dose by starting dose
    # } else {
    #   cohort_size = cohort_size_others
    # } 
    
    ################# END: cohort size (update the following if needed) ############
    
    
    available_patients=available_patients-cohort_size
    cumulative_n = cumulative_n + cohort_size
    
    
    simulated_study_data <- simulated_study_data %>% 
      add_row(
        stratum_id = stratum_id, 
        group_id = group_id, 
        drug_A = current_dose,
        num_patients = cohort_size,
        num_toxicities = rbinom(1,cohort_size,filter(scenarios,drug_A == current_dose)$prob1),
        cohort_time = cohort_time,
        seed = paste(current_seed, collapse = ",")
      )
    # print(simulated_study_data)
    
    ## update BLRM
    blrmfit_new <- update(blrmfit, data = simulated_study_data)
    ## prediction 
    pred_range <- tibble(expand.grid(
      stratum_id = stratum_id, group_id = c(group_id),
      drug_A = dose_info$drug_A, # available dose levels for next
      stringsAsFactors = FALSE
    ))
    
    ################# START: MTD definition (update the following if needed) ############
    pred_range = pred_range %>%
      mutate(increment_ok = drug_A <= (1 + max_percent_inc) * current_dose)
    
    summ_pred = summary(blrmfit_new, newdata = pred_range, interval_prob = c(0, 0.16, 0.33, 1))
    
    # choose next dose
    candidate_doses=filter(summ_pred, ewoc_ok==TRUE &  increment_ok)
    next_dose=max(candidate_doses$drug_A)
    
    # number of patients treated with next dose in the past data of the new study
    past_next_dose=filter(simulated_study_data, drug_A == next_dose)
    n_nextdose=sum(past_next_dose$num_patients)
    n_MTD_ok = (n_nextdose>=patients_on_dose_for_MTD)
    
    #check if target prob > 50%
    pred_next_dose=filter(summ_pred, drug_A == next_dose)
    prob_targ_ok = (pred_next_dose$prob_target>prob_target_cond)
    
    
    ## check if early termination criteria met or not first
    ## if not, propose the next dose satisfies EOWC criteria and maximum increments
    ## additional condition can be added here if needed
    # check max patients
    if (available_patients > 0) {
      # check if MTD met or not
      if (nrow(candidate_doses) == 0){ ## add another criterion for trial stop if needed
        decision="No dose selected"
        MTD = NA
      } else if ( n_MTD_ok & (prob_targ_ok | cumulative_n >= min_patients)){
        decision ="MTD achieved"
        MTD=next_dose
      }
    } else if (available_patients ==0) {
      if(n_MTD_ok){
        decision ="MTD achieved"
        MTD = next_dose
      } else {
        decision="Reached max num of pts"
        MTD = NA
      }
    } else decision =NULL
    
    current_dose=next_dose
    
  }
  
  
  
  OC_results=list(
    blrmfit_new=blrmfit_new,
    summ_pred=summ_pred,
    simulated_study_data=simulated_study_data,
    samplesize=cumulative_n,
    MTD=MTD,
    decision=decision
  )
  
  return(OC_results)
}


# OC_results=bhlrm_oc(scenarios)
# #final results at the end of the trial
# print(OC_results$summ_pred,width=Inf)
# OC_results$simulated_study_data
# OC_results$decision
# OC_results$samplesize
# OC_results$MTD

## to visualise DLT curve
# plot_toxicity_curve(OC_results$blrmfit_new,x = "drug_A",group = ~ group_id)
# p = plot_toxicity_curve(blrmfit_new,x = "drug_A",group = ~ group_id)
# p + geom_line(aes(x=drug_A, y=prob1), data=scenarios)

# summary(blrmfit_new, "dose_prediction")
# check all posterior distribution of all parameters
# print(cbind(pred_range, summ_pred))


## Recover user set sampling defaults
options(.user_mc_options)
