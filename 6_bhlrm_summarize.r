summarize.OC=function(OC_itr_results, cores){
  
  interval_prob = summary(OC_itr_results[[1]]$blrmfit_new, summarize="interval_prob")
  dose_lv = scenarios$drug_A
  
  UD = scenarios[scenarios$prob1 < interval_prob[2],]$drug_A
  TD = scenarios[scenarios$prob1 >= interval_prob[2] & scenarios$prob1 < interval_prob[3],]$drug_A
  OD = scenarios[scenarios$prob1 >= interval_prob[3],]$drug_A
  
  #Distribution of MTD in trials
  MTD.dist = tibble(drug_A=sapply(OC_itr_results, \(x) x$MTD))  %>% 
    count(drug_A) %>%
    mutate(Interval = case_when(
      drug_A %in% UD ~ "UD",
      drug_A %in% TD ~ "TD",
      drug_A %in% OD ~ "OD",
      TRUE ~ NA
    )) %>%
    arrange(drug_A)
  
  
  #number of subjects in each dose in trials
  n_pts.count = bind_rows(lapply(OC_itr_results, 
                                 \(x){x$simulated_study_data %>% 
                                     group_by(drug_A) %>% 
                                     summarize(N=sum(num_patients)) 
                                 }))%>%
    group_by(drug_A) %>%
    mutate(Interval = case_when(
      drug_A %in% UD ~ "UD",
      drug_A %in% TD ~ "TD",
      drug_A %in% OD ~ "OD",
      TRUE ~ NA
    ))
  
  #number of DLT in each dose in trials
  n_tox.count = bind_rows(lapply(OC_itr_results, 
                                 \(x){x$simulated_study_data %>% 
                                     group_by(drug_A) %>% 
                                     summarize(N=sum(num_toxicities))
                                 })) %>%
    group_by(drug_A) %>%
    mutate(Interval = case_when(
      drug_A %in% UD ~ "UD",
      drug_A %in% TD ~ "TD",
      drug_A %in% OD ~ "OD",
      TRUE ~ NA
    ))
  
  #Average proportion of patients
  n_pts.prop =  bind_rows(lapply(OC_itr_results, 
                                 \(x){x$simulated_study_data %>% 
                                     group_by(drug_A) %>% 
                                     summarize(N=sum(num_patients)) %>%
                                     mutate(Nprop = N/x$samplesize)
                                 }))%>%
    mutate(Interval = case_when(
      drug_A %in% UD ~ "UD",
      drug_A %in% TD ~ "TD",
      drug_A %in% OD ~ "OD",
      TRUE ~ NA
    )) %>%
    group_by(Interval) %>%
    summarize(sum_Nprop = sum(Nprop))
  
  #Average sample size
  total.n_ave = mean(sapply(OC_itr_results, \(x) x$samplesize))   
  
  #Proportion of no dose selected (stopped by too toxic)
  tot.noselect = sum(sapply(OC_itr_results, \(x) x$decision == "No dose selected")) 
  
  #Proportion of reached max number of pts
  tot.maxpts = sum(sapply(OC_itr_results, \(x) x$decision == "Reached max num of pts")) 
  
  
  summary = list( interval_prob = interval_prob,
                  n_pts.prop = n_pts.prop, 
                  tot.noselect = tot.noselect, tot.maxpts = tot.maxpts, 
                  total.n_ave = total.n_ave,
                  MTD.dist=MTD.dist, n_pts.count = n_pts.count, n_tox.count = n_tox.count
  )
  return(summary)
}

