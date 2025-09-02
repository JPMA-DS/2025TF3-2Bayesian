#----------------------------------------------------------------#
# Title : Sample data preparation
#----------------------------------------------------------------#

# Parameters
set.seed(2024)             #Seed
acc_CT <- 24               #Recruitment period for clinical trial
fol_CT <- 36               #Follow-up period for clinical trial
acc_HC <- 24               #Recruitment period for historical control
fol_HC <- 36               #Follow-up period for historical control
HR <- 0.65                 #Hazard ratio for clinical trial
driftHR <- 1.2             #Hazard ratio for external control arm
median_CC <- 24            #MST of concurrent control arm
median_E <- 24 / HR        #MST of experimental arm
median_HC <- 24 / driftHR  #MST of historical control arm
num_CT <- 300              #Number of patients for clinical trial
num_HC <- 100              #Number of patients for historical control

# Arm
rnd <- runif(num_CT)
trt <- ifelse(rnd < 0.5, 0, 1)

# Accrual time
dummy_CT <- runif(num_CT)

# Time to event
time_CT <- rexp(num_CT)
time_CT <- ifelse(trt == 0, time_CT / (-log(0.5) / median_CC), time_CT /(-log(0.5) / median_E) )

# Time to event (Calendar time) and censored data
c_t_CT <- dummy_CT * acc_CT + time_CT
cnsr_CT <- numeric(num_CT)

monitor_CT <- tibble(trt, time_CT, cnsr_CT) %>%
  mutate(ext = 0) %>%
  mutate(cnsr = ifelse(c_t_CT > acc_CT + fol_CT, 1, 0)) %>%  #cnsr: 0(event), 1(censor)
  mutate(time = ifelse(c_t_CT > acc_CT + fol_CT, acc_CT + fol_CT - acc_CT * dummy_CT, time_CT)) %>%
  mutate(trt = factor(trt, levels = c(0, 1))) %>%  #Reference: control arm (0)
  select(trt, ext, time, cnsr)

# Generation dataset for histrical control arm
time_HC <- rexp(num_HC) / (-log(0.5) / median_HC)
cnsr_HC <- numeric(num_HC)

# Accrual time
dummy_HC <- runif(num_HC)

# Time to event (Calendar time) and censored data
c_t_HC <- dummy_HC * acc_HC + time_HC
cnsr_HC <- numeric(num_HC)

monitor_HC <- tibble(time_HC, cnsr_HC) %>%
  mutate(ext = 1) %>%
  mutate(trt = 0) %>%
  mutate(cnsr = ifelse(c_t_HC > acc_HC + fol_HC, 1, 0)) %>%
  mutate(time = ifelse(c_t_HC > acc_HC + fol_HC, acc_HC + fol_HC - acc_HC * dummy_HC, time_HC)) %>%
  select(trt, ext, time, cnsr)

# Integrating dataset
data <- rbind(monitor_CT, monitor_HC) %>% as.data.frame()
patient <- c(1:(num_CT + num_HC))
data <- data.frame(patient, data)

# Export data
write.csv(data, "data_CMP.csv", row.names = FALSE)
