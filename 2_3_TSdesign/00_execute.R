#----------------------------------------------------------------#
#Title : TS design prior elicitation
#Date  : 2024/11/2
#----------------------------------------------------------------#

#Setting
#renv::install(c("dplyr", "magrittr", "purrr", "tidyr", "knitr", "extraDistr", "ph2bayes", "ggplot2" ,"sessioninfo", "logr"), type = "win.binary")

logr::log_open("TF32_simple_example")

library(magrittr)
library(ggplot2)

#Parameters for TS design
n_min <- 10
n_max <- 15
alpha_e <- 0.5
beta_e <- 0.5

set.seed(12345)
source("01_TSdesign_PriorElicitation.R")
source("02_TSdesign_Boundary.R")
source("03_TSdesign_OC.R")
source("04_DesignPrior_AnalysisPrior_example_H0.R")
source("05_DesignPrior_AnalysisPrior_example_H1.R")

#Session
session <- sessioninfo::session_info()
capture.output(session, file = "./output/00_session.txt")

logr::log_close()
