#----------------------------------------------------------------#
# Title : Execution file
#----------------------------------------------------------------#

# Package install
#renv::install(c("magrittr", "dplyr", "survival", "survminer", "R2jags", "ggplot2", "rstan", "readr", "logr", "sessioninfo"), type = "win.binary")

# Package loading
library(magrittr)
library(dplyr)
library(survival)
library(R2jags)
library(ggplot2)
library(survminer)
library(readr)
library(sessioninfo)
library(rstan)
library(logr)

log_open("TF32_ext")

# Data preparation
data <- readr::read_csv(file = "data_CMP.csv",
                        col_types = cols(trt = col_factor(levels = c("0", "1")),
                                         ext = col_factor(levels = c("0", "1")))
)

data.nb <- data %>% filter(ext == 0) # No borrowing dataset
data.C  <- data %>% filter(trt == 0) # Control arm dataset

# Execution
source("./R/01_data_preparation.R")
source("./R/02_kmph_c.R")
source("./R/03_kmph_nb.R")
source("./R/04_kmph_fb.R")
source("./R/11_weibull_cmp.R")
source("./R/12_weibull_ehss.R")
source("./R/22_weibull_cmp_stan.R")

#Session
session <- session_info()
capture.output(session, file = "./output/00_session.txt")

log_close()
