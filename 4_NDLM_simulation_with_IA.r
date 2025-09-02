

library(rstan)
library(tidyverse)

#set work space
base_dir <- "//"
setwd(base_dir)


# iteration number for simulation 
simu_num=

# study with IA (1) or without IA (0)
IAYN<-

#dose response scenario
DR1<- c(0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00,	0.00)
DR2a<-c(0.00,	0.12,	0.24,	0.36,	0.48,	0.73,	1.09,	1.45)
DR2b<-c(0.00,	0.95,	1.17,	1.27,	1.32,	1.38,	1.43,	1.45)
DR2c<-c(0.00,	0.92,	1.21,	1.36,	1.43,	1.42,	1.12,	0.57)
DR2d<-c(0.00,	0.00,	0.00,	0.00,	0.00,	0.05,	1.02,	1.45)
DR3a<-c(0.00,	0.18,	0.37,	0.55,	0.73,	1.10,	1.65,	2.20)
DR3b<-c(0.00,	1.44,	1.78,	1.93,	2.01,	2.10,	2.17,	2.20)
DR3c<-c(0.00,	1.40,	1.83,	2.07,	2.18,	2.15,	1.69,	0.87)
DR3d<-c(0.00,	0.00,	0.00,	0.00,	0.00,	0.08,	1.54,	2.20)
DR4a<-c(0.00,	0.24,	0.48,	0.73,	0.97,	1.45,	2.18,	2.90)
DR4b<-c(0.00,	1.90,	2.34,	2.54,	2.65,	2.77,	2.85,	2.90)
DR4c<-c(0.00,	1.85,	2.42,	2.72,	2.87,	2.84,	2.23,	1.15)
DR4d<-c(0.00,	0.00,	0.00,	0.00,	0.01,	0.10,	2.04,	2.90)

DR_list<-list(DR1,DR2a,DR2b,DR2c,DR2d,DR3a,DR3b,DR3c,DR3d,DR4a,DR4b,DR4c,DR4d)
scenario_name<-c("1","2a","2b","2c","2d","3a","3b","3c","3d","4a","4b","4c","4d")


#Specify NDLM
stanmodelcode <- "
data {
  int<lower=0> N; 
  int<lower=0> K; 
  real y[N]; 
  int<lower=0> arm[N]; 
}
parameters {
  real theta[K];
  real delta[K];
  real<lower=0> sigma2;
  real<lower=0> W;
}

model {
  for(k in 2:K){
      theta[k] ~ normal(theta[k-1]+delta[k-1], sqrt(W*sigma2));
      delta[k] ~ normal(delta[k-1], sqrt(W*sigma2));
    }
    theta[1] ~ normal(0, 1000);
    delta[1] ~ normal(0, 1000);
    for(j in 1:N){
      y[j] ~ normal(theta[arm[j]], sqrt(sigma2));
    }
  W~uniform(0.001,100);
  sigma2~inv_gamma(0.001, 0.001);
  }
"
# Model information is transferred to STAN
stanmodel1 <- stan_model(model_code = stanmodelcode)


for(scenario in 1 :13){ #repetitions for scenarios
#set the number of seed
set.seed(scenario) 
#select dose response scenario
DR<-unlist(DR_list[scenario])

data <- as.data.frame(matrix(ncol=2, nrow=280))
colnames(data) <- c("y", "arm")
result <- as.data.frame(matrix(ncol=48, nrow=simu_num))
colnames(result) <- c("theta_1", "theta_2", "theta_3", "theta_4", "theta_5", "theta_6", "theta_7", "theta_8", "sigma2", "W",
                      "max_2","max_3","max_4","max_5","max_6","max_7","max_8",
                      "power_2","power_3","power_4","power_5","power_6","power_7","power_8", "power",
                      "selected_2","selected_3","selected_4","selected_5","selected_6","selected_7","selected_8", 
                      "N1","N2","N3","N4","N5","N6","N7","N8","Ntotal",
                      "futile_2","futile_3","futile_4","futile_5","futile_6","futile_7","futile_8")


  for (i in 1:simu_num){ #simulation repetition
    # data generation
    
    data<-data %>%
      mutate(arm = rep(1:8, 35)) %>%
      mutate(y=rnorm(280, mean=DR[arm], sd=2.3)) %>%
      mutate(ran=runif(280, min = 0, max = 1)) %>%
      arrange(by=ran) %>%
      mutate(rep=row_number()) %>%
      mutate(IA_F=if_else(rep<=140,1,0)) #randomly select 50% of patients for IA data set
  
    if(IAYN==1){  
      # IA data set
      data_IA <-data %>%
        filter(IA_F==1)
      data_list_IA<-list(K=8, N=140, y=data_IA$y, arm=data_IA$arm)

      #IA
      # MCMC by STAN
      fit1 <- sampling(stanmodel1, data = data_list_IA, chains = 1, thin=10, iter = 30000, warmup = 3000)
        rstan::extract(fit1,pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]","theta[7]","theta[8]","sigma2", "W")) %>% 
        data.frame() ->MC_sample_IA 
      colnames(MC_sample_IA) <- c("theta_1", "theta_2", "theta_3", "theta_4", "theta_5", "theta_6", "theta_7", "theta_8", "sigma2", "W")
  
      #Drop arm
      MC_sample_IA <-MC_sample_IA %>% 
        mutate(MCdiff_2=theta_2-theta_1, #difference between active arm 2,3,4,5,6,7 vs PLB (arm=0) for each MCMC sample
               MCdiff_3=theta_3-theta_1,
               MCdiff_4=theta_4-theta_1,
               MCdiff_5=theta_5-theta_1,
               MCdiff_6=theta_6-theta_1,
               MCdiff_7=theta_7-theta_1,
               MCdiff_8=theta_8-theta_1,
               MCID=row_number()) %>%
        mutate(futile_2f=if_else(MCdiff_2<1.5, 1, 0), # inferiority beyond margin(=1.5)          
               futile_3f=if_else(MCdiff_3<1.5, 1, 0), 
               futile_4f=if_else(MCdiff_4<1.5, 1, 0),
               futile_5f=if_else(MCdiff_5<1.5, 1, 0),
               futile_6f=if_else(MCdiff_6<1.5, 1, 0),
               futile_7f=if_else(MCdiff_7<1.5, 1, 0),
               futile_8f=if_else(MCdiff_8<1.5, 1, 0))
  
        Post_IA<-bind_rows(colMeans(MC_sample_IA[,])) %>% #calculate futility prob P(diff between active arms vs PLB < margin=1.5|data) by average out MCMC samples
          mutate(futile_2=if_else(futile_2f>0.8,1,0), 
                 futile_3=if_else(futile_3f>0.8,1,0),
                 futile_4=if_else(futile_4f>0.8,1,0),
                 futile_5=if_else(futile_5f>0.8,1,0),
                 futile_6=if_else(futile_6f>0.8,1,0),
                 futile_7=if_else(futile_7f>0.8,1,0),
                 futile_8=if_else(futile_8f>0.8,1,0))
        #Increment of data by FA from IA(selection by flagging by IA_F=0)
        data_FA <-data %>% 
          filter(IA_F==0) %>% 
          filter(  arm==1                        | #PLB is not dropped (always kept)
                  (Post_IA$futile_2==0 & arm==2) | #drop active doses whose futility probabilities are >0.8
                  (Post_IA$futile_3==0 & arm==3) |
                  (Post_IA$futile_4==0 & arm==4) |              
                  (Post_IA$futile_5==0 & arm==5) |
                  (Post_IA$futile_6==0 & arm==6) |
                  (Post_IA$futile_7==0 & arm==7) |
                  (Post_IA$futile_8==0 & arm==8) )
    
          #If study is stopped at IA (all active arms are dropped) data for NDLM estimation will be IA data set
          #otherwise, the data will be combined data set IA data set and incremented data by FA
          if(Post_IA$futile_2==1 &        
             Post_IA$futile_3==1 & 
             Post_IA$futile_4==1 & 
             Post_IA$futile_5==1 &
             Post_IA$futile_6==1 &
             Post_IA$futile_7==1 &
             Post_IA$futile_8==1){
             data_all<-data_IA
          }else{
               data_all<-rbind(data_IA, data_FA)  # in the case of study with IA, the data used for dose response estimate "data_all" is dat_IA plus data_FA      
          }
    }else{
      data_all<-data # in the case of study without IA, the data used for dose response estimate "data_all" is data
    }
        data_list<-list(K=8, N=nrow(data_all), y=data_all$y, arm=data_all$arm)
    
        # MCMC by STAN
        fit2 <- sampling(stanmodel1, data = data_list, chains =  1, thin=10, iter = 30000, warmup = 3000)
        rstan::extract(fit2, pars=c("theta[1]","theta[2]","theta[3]","theta[4]","theta[5]","theta[6]","theta[7]","theta[8]","sigma2", "W")) %>% 
          data.frame() ->MC_sample 
        colnames(MC_sample) <- c("theta_1", "theta_2", "theta_3", "theta_4", "theta_5", "theta_6", "theta_7", "theta_8", "sigma2", "W")
      
        #test
        MC_sample <-MC_sample %>% 
          mutate(MCdiff_2=theta_2-theta_1,
                 MCdiff_3=theta_3-theta_1,
                 MCdiff_4=theta_4-theta_1,
                 MCdiff_5=theta_5-theta_1,
                 MCdiff_6=theta_6-theta_1,
                 MCdiff_7=theta_7-theta_1,
                 MCdiff_8=theta_8-theta_1,
                 MCID=row_number()) %>%
          mutate(test_2f=if_else(MCdiff_2>1.5, 1, 0), # superiority beyond margin(=1.5)                
                 test_3f=if_else(MCdiff_3>1.5, 1, 0),
                 test_4f=if_else(MCdiff_4>1.5, 1, 0),
                 test_5f=if_else(MCdiff_5>1.5, 1, 0),
                 test_6f=if_else(MCdiff_6>1.5, 1, 0),
                 test_7f=if_else(MCdiff_7>1.5, 1, 0),
                 test_8f=if_else(MCdiff_8>1.5, 1, 0)) %>%
          group_by(MCID) %>%
            mutate(max_2=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_2,1,0), # flagging max effect dose 
                 max_3=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_3,1,0),
                 max_4=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_4,1,0),
                 max_5=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_5,1,0),
                 max_6=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_6,1,0),
                 max_7=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_7,1,0),
                 max_8=if_else(max(MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8)==MCdiff_8,1,0)) %>%
          ungroup()  %>%
          select(-c(MCID,MCdiff_2,MCdiff_3,MCdiff_4,MCdiff_5,MCdiff_6,MCdiff_7,MCdiff_8))
      
        Post<-bind_rows(colMeans(MC_sample[,]))  #calculate superiority prob P(diff between active arms vs PLB > margin=1.5|data) by average out MCMC samples
  
          if(IAYN==1){
            Post <-Post %>%                          #prob of with max effect prob by average out MCMC samples  
              mutate(power_2=if_else(Post_IA$futile_2==0 & test_2f>0.8, 1,0),         
                     power_3=if_else(Post_IA$futile_3==0 & test_3f>0.8, 1,0),
                     power_4=if_else(Post_IA$futile_4==0 & test_4f>0.8, 1,0),
                     power_5=if_else(Post_IA$futile_5==0 & test_5f>0.8, 1,0),
                     power_6=if_else(Post_IA$futile_6==0 & test_6f>0.8, 1,0),
                     power_7=if_else(Post_IA$futile_7==0 & test_7f>0.8, 1,0),
                     power_8=if_else(Post_IA$futile_8==0 & test_8f>0.8, 1,0))
        }else{
            Post <-Post %>% 
              mutate(power_2=if_else(test_2f>0.8, 1,0),         
                    power_3=if_else(test_3f>0.8, 1,0),
                    power_4=if_else(test_4f>0.8, 1,0),
                    power_5=if_else(test_5f>0.8, 1,0),
                    power_6=if_else(test_6f>0.8, 1,0),
                    power_7=if_else(test_7f>0.8, 1,0),
                    power_8=if_else(test_8f>0.8, 1,0))
          }
          Post <-Post %>%
          mutate(power=if_else(max(power_2,power_3,power_4,power_5,power_6,power_7,power_8)>0,1,0))
        
          if(IAYN==1){
            Post <-Post %>%
              mutate(selected_2=if_else(Post_IA$futile_2==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_2,1,0), #select the dose with max effect as the dose which have the largest prob of max effect
                   selected_3=if_else(Post_IA$futile_3==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_3,1,0),
                   selected_4=if_else(Post_IA$futile_4==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_4,1,0),
                   selected_5=if_else(Post_IA$futile_5==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_5,1,0),
                   selected_6=if_else(Post_IA$futile_6==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_6,1,0),
                   selected_7=if_else(Post_IA$futile_7==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_7,1,0),
                   selected_8=if_else(Post_IA$futile_8==0 & max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_8,1,0)) %>%
              select(-c(test_2f,test_3f,test_4f,test_5f,test_6f,test_7f,test_8f))
          }else{
            Post <-Post %>%
              mutate(selected_2=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_2,1,0), #select the dose with max effect as the dose which have the largest prob of max effect
                   selected_3=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_3,1,0),
                   selected_4=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_4,1,0),
                   selected_5=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_5,1,0),
                   selected_6=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_6,1,0),
                   selected_7=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_7,1,0),
                   selected_8=if_else(max(max_2,max_3,max_4,max_5,max_6,max_7,max_8)==max_8,1,0)) %>%
              select(-c(test_2f,test_3f,test_4f,test_5f,test_6f,test_7f,test_8f))
          }
        
        Post <- Post %>%
         # number of enrolled patients for each dose, and total number
          mutate(
            N1=sum(data_all$arm==1),
            N2=sum(data_all$arm==2),
            N3=sum(data_all$arm==3),
            N4=sum(data_all$arm==4),
            N5=sum(data_all$arm==5),
            N6=sum(data_all$arm==6),
            N7=sum(data_all$arm==7),
            N8=sum(data_all$arm==8),
            Ntotal=nrow(data_all) )

        # Save futility results
         if(IAYN==1){      
           Post <- Post %>%
            mutate(
              futile_2=Post_IA$futile_2,
              futile_3=Post_IA$futile_3,
              futile_4=Post_IA$futile_4,
              futile_5=Post_IA$futile_5,
              futile_6=Post_IA$futile_6,
              futile_7=Post_IA$futile_7,
              futile_8=Post_IA$futile_8)
        }else{
           Post <- Post %>%
             mutate(
               futile_2=NA,
               futile_3=NA,
               futile_4=NA,
               futile_5=NA,
               futile_6=NA,
               futile_7=NA,
               futile_8=NA)
         }
         result[i,] <-Post
  }

  #Bias, MSE, 25%, 50%, 75% tiles
  result<-result %>%
    mutate(diff_2=theta_2-theta_1,
           diff_3=theta_3-theta_1,
           diff_4=theta_4-theta_1,
           diff_5=theta_5-theta_1,
           diff_6=theta_6-theta_1,
           diff_7=theta_7-theta_1,
           diff_8=theta_8-theta_1) %>%
    mutate(bias_2=diff_2-(DR[2]-DR[1]),
           bias_3=diff_3-(DR[3]-DR[1]),
           bias_4=diff_4-(DR[4]-DR[1]),
           bias_5=diff_5-(DR[5]-DR[1]),
           bias_6=diff_6-(DR[6]-DR[1]),
           bias_7=diff_7-(DR[7]-DR[1]),
           bias_8=diff_8-(DR[8]-DR[1])) %>%
  mutate(MSE_2= (diff_2-(DR[2]-DR[1]))^2,
         MSE_3= (diff_3-(DR[3]-DR[1]))^2,
         MSE_4= (diff_4-(DR[4]-DR[1]))^2,
         MSE_5= (diff_5-(DR[5]-DR[1]))^2,
         MSE_6= (diff_6-(DR[6]-DR[1]))^2,
         MSE_7= (diff_7-(DR[7]-DR[1]))^2,
         MSE_8= (diff_8-(DR[8]-DR[1]))^2) 
  # calculate RMSE 
  summary_result<-bind_rows(colMeans(result[,])) %>%
    mutate(RMSE_2= MSE_2^0.5,
           RMSE_3= MSE_3^0.5,
           RMSE_4= MSE_4^0.5,
           RMSE_5= MSE_5^0.5,
           RMSE_6= MSE_6^0.5,
          RMSE_7= MSE_7^0.5,
          RMSE_8= MSE_8^0.5) %>%
  #calculate percentiles of posterior mean of patient outcome in each arm over simulated studies for graphical display 
    mutate(theta_1_P25=quantile(result$theta_1, probs=0.25),
           theta_1_P50=quantile(result$theta_1, probs=0.50),
           theta_1_P75=quantile(result$theta_1, probs=0.75),
         
           theta_2_P25=quantile(result$theta_2, probs=0.25),
           theta_2_P50=quantile(result$theta_2, probs=0.50),
           theta_2_P75=quantile(result$theta_2, probs=0.75),
         
           theta_3_P25=quantile(result$theta_3, probs=0.25),
           theta_3_P50=quantile(result$theta_3, probs=0.50),
           theta_3_P75=quantile(result$theta_3, probs=0.75),
         
           theta_4_P25=quantile(result$theta_4, probs=0.25),
           theta_4_P50=quantile(result$theta_4, probs=0.50),
           theta_4_P75=quantile(result$theta_4, probs=0.75),
         
           theta_5_P25=quantile(result$theta_5, probs=0.25),
           theta_5_P50=quantile(result$theta_5, probs=0.50),
           theta_5_P75=quantile(result$theta_5, probs=0.75),
         
           theta_6_P25=quantile(result$theta_6, probs=0.25),
           theta_6_P50=quantile(result$theta_6, probs=0.50),
           theta_6_P75=quantile(result$theta_6, probs=0.75),
         
           theta_7_P25=quantile(result$theta_7, probs=0.25),
           theta_7_P50=quantile(result$theta_7, probs=0.50),
           theta_7_P75=quantile(result$theta_7, probs=0.75),
         
           theta_8_P25=quantile(result$theta_8, probs=0.25),
           theta_8_P50=quantile(result$theta_8, probs=0.50),
           theta_8_P75=quantile(result$theta_8, probs=0.75)) %>%
    mutate(scenario_name=scenario_name[scenario])

  #specify file name for simulation result csv file
  file_name <- paste0(base_dir,"/result/","simulation_result", "IA_", IAYN, "_Scenario_", scenario_name[scenario],".csv")
  write_csv(summary_result, file_name) 
  

  #save the graph of dose response
  summary_mean1<-summary_result %>%
    select(c(theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, theta_7)) %>% 
    gather(key=variable, value= PosteriorMean) %>%
    mutate(dose=row_number()) %>% 
    select(-c(variable))
  
  summary_mean2<-summary_result %>%
    select(c(theta_1_P25, theta_2_P25, theta_3_P25, theta_4_P25, theta_5_P25, theta_6_P25, theta_7_P25)) %>% 
    gather(key=variable, value= P25) %>%
    mutate(dose=row_number()) %>% 
    select(-c(variable))
  summary_mean3<-summary_result %>%
    select(c(theta_1_P75, theta_2_P75, theta_3_P75, theta_4_P75, theta_5_P75, theta_6_P75, theta_7_P75)) %>% 
    gather(key=variable, value= P75) %>%
    mutate(dose=row_number()) %>% 
    select(-c(variable))

  dose_N <- data.frame(dose_N=c(0,50,100,150,200,300,450,600))
  dose_N <-dose_N %>% 
    mutate(dose=row_number())

  summary_mean <- left_join(summary_mean1, summary_mean2, by="dose") %>%
    left_join(summary_mean3, by="dose") %>%
    mutate(DR=DR[dose]) %>% 
    left_join(dose_N, by="dose")

  p<-ggplot(data = summary_mean) +
    geom_point(aes(x = dose_N, y = DR), shape = 5, size = 3) + 
    geom_line(aes(x = dose_N,  y = P25), linetype = 2) +
    geom_line(aes(x = dose_N,  y = PosteriorMean), linetype = 1, linewidth = 1) +
    geom_line(aes(x = dose_N,  y = P75), linetype = 2) +
    scale_x_continuous(breaks = c(0,50,100,150,200,300,450)) +
    labs(x = "Dose", y = "Change from Baseline")

  #specify file name for simulation result png file
  pic_name <- paste0(base_dir,"/result/","simulation_graph", "IA_", IAYN, "_Scenario_", scenario_name[scenario],".png")
  ggsave(file = pic_name, plot = p)
}
