
rm(list = ls())
gc()

library(here) # to get path where R project is located
### echo = FALSE if plot for prior curve is not needed.
source(here::here("BHLRM_exnex/bhlrm_oc.R"), echo = FALSE)
source(here::here("BHLRM_exnex/bhlrm_summarize.R"), echo = FALSE)

# rm(OC_results)


## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
# .user_mc_options <- options(OncoBayes2.MC.warmup=10, OncoBayes2.MC.iter=50, OncoBayes2.MC.chains=4,
#                             OncoBayes2.MC.save_warmup=FALSE)

###################################################

library(doParallel)
library(progress)

t<-proc.time()

cores <- getOption("mc.cores", detectCores())/2
# cores = 8
cl <- makeCluster(cores)
registerDoParallel(cl)

nTrials=1000
n_rounded_trials = floor(nTrials / cores) * cores
trials_dropped = nTrials - n_rounded_trials
n_itr = n_rounded_trials/cores

if (trials_dropped != 0){
  cat("Original nTrials:", nTrials, "\n")
  cat("Rounded nTrials:", n_rounded_trials, "\n")
  cat("Number of dropped trials due to rounding:", trials_dropped, "\n")
  
}

## The number of decision making <= max number of patients for each trial
## Create seeds for each decision making in all trials using RNGseq with an initial seed 123
## Note: 6 integers are needed to create a seed with RNG
RNGkind("L'Ecuyer-CMRG")
all_seeds = RNGseq(nTrials * max_patients * 6, seed = 123) 
seed_list1 = split(all_seeds[1:(n_rounded_trials*max_patients*6)], rep(1:n_rounded_trials, each = max_patients * 6))



### uncomment the following if one needs to visualise progress of simulation
pb <- progress_bar$new(
  total = n_rounded_trials,
  format = "  Processing [:bar] :percent :elapsedfull ETA: :eta"
)

OC_results = list()
OC_itr_results <- foreach(i=1:n_itr, 
                          # .packages=c("OncoBayes2", "tidyr", "dplyr"), 
                          .packages=c("OncoBayes2", "tidyr", "dplyr","doRNG") ) %dopar% { 
                            
                            for (j in 1:cores) {
                              trial_idx = (i - 1) * cores + j
                              OC_result = bhlrm_oc(scenarios, seeds = seed_list1[[trial_idx]])
                              OC_results = append(OC_results, list(OC_result))
                              
                              ### uncomment the following if one needs to visualise progress of simulation
                              pb$tick()  # Update the progress bar
                            }
                            
                            OC_summary=summarize.OC(OC_results, cores)
                            OC_results = list()
                            return(OC_summary)
                          }

stopCluster(cl)

if(nTrials != n_rounded_trials){
  seed_list2 = split(all_seeds[(n_rounded_trials*max_patients*6+1):(nTrials * max_patients * 6) ], 
                     rep((n_rounded_trials+1):nTrials, each = max_patients * 6))
  for (j in 1:(nTrials-n_rounded_trials)) {
    trial_idx =  j
    OC_result = bhlrm_oc(scenarios, seeds = seed_list2[[trial_idx]])
    OC_results = append(OC_results, list(OC_result))
  }
  OC_summary_sub=summarize.OC(OC_results, (nTrials-n_rounded_trials) )
}

OC_itr_results = append(OC_itr_results, list(OC_summary_sub))

sim_time=proc.time()-t




gc()
###########################################

interval_prob = OC_itr_results[[1]]$interval_prob

MTD.prop = 
  bind_rows(lapply(OC_itr_results, \(x) x$MTD.dist)) %>%
  group_by(Interval) %>%
  summarize(prop.MTD = sum(n)/nTrials)

prop.noselect = sum(sapply(OC_itr_results, \(x) x$tot.noselect))/nTrials

prop.maxpts = sum(sapply(OC_itr_results, \(x) x$tot.maxpts))/nTrials

n_pts.prop = 
  bind_rows(lapply(OC_itr_results, \(x) x$n_pts.prop)) %>%
  group_by(Interval) %>%
  summarize(average_Nprop = sum(sum_Nprop)/nTrials)

total.n_ave = mean(sapply(OC_itr_results, \(x) x$total.n_ave))

MTD.dist =   
  bind_rows(lapply(OC_itr_results, \(x) x$MTD.dist)) %>%
  group_by(drug_A) %>%
  summarize(MTD_selected = sum(n)) %>%
  mutate(MTD.prop = round(MTD_selected/nTrials, 4))

n_pts.count = 
  bind_rows(lapply(OC_itr_results, \(x) x$n_pts.count)) %>%
  group_by(drug_A) %>%
  summarize(average_N = sum(N/nTrials)) 

n_pts.count_int = 
  bind_rows(lapply(OC_itr_results, \(x) x$n_pts.count)) %>%
  group_by(Interval) %>%
  summarize(average_N = sum(N/nTrials)) 

n_tox.count = 
  bind_rows(lapply(OC_itr_results, \(x) x$n_tox.count)) %>%
  group_by(drug_A) %>%
  summarize(average_N = sum(N/nTrials)) 

n_tox.count_int = 
  bind_rows(lapply(OC_itr_results, \(x) x$n_tox.count)) %>%
  group_by(Interval) %>%
  summarize(average_N = sum(N/nTrials))

total.n_tox_ave = bind_rows(lapply(OC_itr_results, \(x) x$n_tox.count)) %>%
  summarize(average_N = sum(N/nTrials)) 


#######################


sink("Out_summary_OC_s1.txt")

cat(paste("Proportion of trials with MTD within under dose region (< ",
          interval_prob[2]*100, "%) = ",
          round(MTD.prop %>% filter(Interval == "UD") %>% select(prop.MTD),3), sep=""))
cat("\n")
cat(paste("Proportion of trials with MTD within target dose region (>= ", 
          interval_prob[2]*100,"%", "-", interval_prob[3]*100, "%) = ",
          round(MTD.prop %>% filter(Interval == "TD") %>% select(prop.MTD),3), sep=""))
cat("\n")
cat(paste("Proportion of trials with MTD within over dose region (>= ",  
          interval_prob[3]*100, "%) = ",
          round(MTD.prop %>% filter(Interval == "OD") %>% select(prop.MTD),3), sep=""))
cat("\n")
cat(paste("Proportion of trials with no dose selected : ", 
          round(prop.noselect,3), sep=""))
cat("\n")
cat(paste("Proportion of trials with reached max number of patients : ", 
          round(prop.maxpts,3), sep=""))
cat("\n")
cat("\n")
#cat(paste("Average proportion of patients in under dose (<", 
#          interval_prob[2]*100, "%) = ",
#          round(n_pts.prop %>% filter(Interval == "UD") %>% select(average_Nprop),3), sep=""))
#cat("\n")
#cat(paste("Average proportion of patients in target dose (>=", 
#          interval_prob[2]*100,"%", "-", interval_prob[3]*100, "%) = ",
#          round(n_pts.prop %>% filter(Interval == "TD") %>% select(average_Nprop),3),sep=""))
#cat("\n")
#cat(paste("Average proportion of patients in over dose (>=",  
#          interval_prob[3]*100, "%) = ",
#          round(n_pts.prop %>% filter(Interval == "OD") %>% select(average_Nprop),3),sep=""))
#cat("\n")
#cat("\n")
cat(paste("Average Number of patients in under dose (<", 
          interval_prob[2]*100, "%) = ",
          round(n_pts.count_int %>% filter(Interval == "UD") %>% select(average_N),3), sep=""))
cat("\n")
cat(paste("Average Number of patients in target dose (>=", 
          interval_prob[2]*100,"%", "-", interval_prob[3]*100, "%) = ",
          round(n_pts.count_int %>% filter(Interval == "TD") %>% select(average_N),3),sep=""))
cat("\n")
cat(paste("Average Number of patients in over dose (>=",  
          interval_prob[3]*100, "%) = ",
          round(n_pts.count_int %>% filter(Interval == "OD") %>% select(average_N),3),sep=""))
cat("\n")
cat("\n")
cat(paste("Average Number of DLT in under dose (<", 
          interval_prob[2]*100, "%) = ",
          round(n_tox.count_int %>% filter(Interval == "UD") %>% select(average_N),3), sep=""))
cat("\n")
cat(paste("Average Number of DLT in target dose (>=", 
          interval_prob[2]*100,"%", "-", interval_prob[3]*100, "%) = ",
          round(n_tox.count_int %>% filter(Interval == "TD") %>% select(average_N),3),sep=""))
cat("\n")
cat(paste("Average Number of DLT in over dose (>=",  
          interval_prob[3]*100, "%) = ",
          round(n_tox.count_int %>% filter(Interval == "OD") %>% select(average_N),3),sep=""))

cat("\n")
cat("\n")
cat(paste("Average sample size : ", round(total.n_ave,3), sep=""))
cat("\n")
cat("\n")
cat(paste("Average Number of DLT : ", round(total.n_tox_ave,3), sep=""))
cat("\n")
cat("\n")
cat("Distribution of MTD")
cat("\n")
print(as.data.frame(MTD.dist))
cat("\n")
cat("Average Number of patients in each dose: ")
cat("\n")
print(round(as.data.frame(n_pts.count),3))
cat("\n")
cat("Average Number of DLT in each dose:")
cat("\n")
print(round(as.data.frame(n_tox.count),3))

cat("\n")
cat("\n")
cat(paste("(# of trials: ", nTrials, ")",sep=""))
cat("\n")
cat(paste("(Time [sec]:  ", sim_time[3], ")", sep=""))
sink()

## Recover user set sampling defaults
options(.user_mc_options)


