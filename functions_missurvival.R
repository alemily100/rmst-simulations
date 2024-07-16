library(survival)
library(tidyverse)
library(survminer)
#library(Rlab)
library(calculus)
library(dplyr)
library(eha)
#library(rgl)
library(survRM2)
library(parallel)

###########################generate_ps##########################################
#Description: generate time to observed failure for a sample of patients 

#Input:
#n - trial sample size
#t1 - knot point 
#beta 1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta12 - interaction  coefficient  between covariate 1 and 2   
#beta3 - covariate 3  coefficient 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)
#lambda_censor - censoring parameter
#recruit_end - time until recruitment ends (months)
#end_time - end of trial (months)

#Output: data frame containing patient ID, covariate indicators, time until 
#failure time, censoring indicator and arrival time within trial.

generate_ps <- function(n,
                        t1,
                        beta1a,
                        beta1b,
                        beta2a,
                        beta2b,
                        beta12,
                        beta3,
                        h0a,
                        h0b,
                        lambda_censor,
                        recruit_end,
                        end_time){
  
  
  
  #fixed covariates
  Z <- cbind(cov_1<-rbinom(n, 1, prob = 0.5),
             cov_2<-rbinom(n, 1, prob = 0.5), cov12<- ifelse(cov_1==1 & cov_2==1, 1, 0), cov_3<-rbinom(n, 1, prob = 0.5))
  
  #Generate survival times for each patient
  T_val <- sapply(1:n, function(i){
    
    #Cumulative hazard function
    H <- function(t){
      #Before knot-point, t1
      if(t<=t1){
        H_val <- t*exp(h0a+beta1a*Z[i,1]+beta2a*Z[i,2]+beta12*Z[i,3]+beta3*Z[i,4])
      }
      #After knot-point, t1
      if(t>t1){
        H_val <- t1*exp(h0a+beta1a*Z[i,1]+beta2a*Z[i,2]+beta12*Z[i,3]+beta3*Z[i,4]) + (t-t1)*exp(h0b+beta1b*Z[i,1]+beta2b*Z[i,2]+beta12*Z[i,3]+beta3*Z[i,4])
      }
      return(H_val)
    }
    
    
    #### inversion method e.g T = H^-1(u) ##########
    #Set a while loop to avoid having Ts=0 when generated uniform RV is tiny.
    Ts<-0
    while(Ts==0){
      #generate uniform random variable for use in the inversion method
      u <- runif(1)
      Ts <- uniroot(function(t) H(t) + log(1-u), lower=0,
                    upper = 100,
                    extendInt = "upX")$root  
    }
    return(Ts)})
  
  ####Censoring times ######
  C_val <- rexp(n, lambda_censor)
  status <-  as.numeric(as.numeric(T_val) < as.numeric(C_val))
  T_val[status==0] <- C_val[status==0] #replace event times with censoring times for
  T_val <- as.numeric(T_val)           #censored observations
  
  #Generate arrival times
  arrival_times <- runif(n, 0, recruit_end)
  
  #All survival data
  surv.data <- data.frame(patient=c(1:n),
                          treatment=Z[,1],
                          inherit = Z[,2],
                          treatment.inherit = Z[,3],
                          sex = Z[,4],
                          surv_times = T_val,
                          status = status,
                          arrival_times = arrival_times)
  
  ###End-of-study censoring
  #Remove patients who haven't arrived yet (relevant for group sequential tests only)
  surv.data <- surv.data[surv.data$arrival_times <= end_time,]
  #Censor unobserved observations and change observation time
  un_obs <- surv.data$surv_times+surv.data$arrival_times > end_time
  surv.data$surv_times[un_obs] <- end_time-surv.data$arrival_times[un_obs]
  surv.data$status[un_obs] <- 0
  
  return(surv.data)
  return(C_val)
}



###############################nparam_rmst###################################### 
#Description: Find non-parametric area below Kaplan-Meier plot 

#Input: 
#kaplan_meier - a kaplan-meier fitted using rmst2
#t_star - upper bound for calculating RMST (t*)

#Output: non-parametric RMST estimate (numeric)
nparam_rmst<- function(kaplan_meier, t_star){
  survival:::survmean(kaplan_meier, rmean=t_star)[[1]][2,5] - 
    survival:::survmean(kaplan_meier, rmean=t_star)[[1]][1,5]
}

#######################nparam_rmst_var##########################################
#Description: Find variance for non-parametric RMST estimate 

#Input: 
#kaplan_meier - a kaplan-meier fitted using rmst2
#t_star - upper bound for calculating RMST (t*)

#Output: variance of non-parametric RMST estimate (numeric)
nparam_rmst_var<- function(kaplan_meier, t_star){
  survival:::survmean(kaplan_meier, rmean=t_star)[[1]][1,6]^2 + survival:::survmean(kaplan_meier, rmean=t_star)[[1]][2,6]^2
}

##########################################area#################################
#Description: Find parametric area below cox-exponential survival curve 

#Input: 
#t_star - upper bound for calculating RMST (t*)
#base_hazard - baseline hazard (\lambda)
#vec_cov_coef - vector of covariate coefficients 
#vec_cov_param - vector of covariate indicators 

#Ouput - area below cox-exponential survival curve (numeric)
area<- function(t_star, base_hazard, vec_cov_coef, vec_cov_param){
  x_t_b<- vec_cov_coef%*%vec_cov_param
  area<-(1/(base_hazard*exp(x_t_b)))*(1-exp(-base_hazard*t_star*exp(x_t_b)))
  return(area)
}

###############################h_scale_t########################################
#Description: Find fully specified parametric RMST estimate 

#Input:
#log_scale - log transformed base hazard (log(\lambda))
#b_1 - covariate 1 coefficient 
#b_2 - covariate 2 coefficient 
#b_12 - interaction  coefficient  between covariate 1 and 2  
#b_3 - covariate 3 coefficient 
#t_star - upper bound for calculating RMST (t*)

#Output: Fully specified parametric RMST estimate (numeric)

h_scale_t<- function(log_scale, b_1, b_2, b_12, b_3, t_star){
  0.25*((exp(log_scale)/exp(b_1+b_2+b_12+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_2+b_12+b_3))) + 
          (exp(log_scale)/exp(b_1+b_2+b_12))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_2+b_12)))+
          (exp(log_scale)/exp(b_1+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_3)))+
          (exp(log_scale)/exp(b_1))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1)))) - 
    0.25*((exp(log_scale)/exp(b_2+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_2+b_3))) + 
            (exp(log_scale)/exp(b_2))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_2)))+ 
            (exp(log_scale)/exp(b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_3)))+
            (exp(log_scale))*(1-exp(-(1/exp(log_scale))*t_star)))}


##############################param_rmst_var####################################
#Description: Find variance for fully speicified parametric RMST estimate 

#Input: 
#log_scale - log transformed base hazard (log(\lambda))
#beta1 - covariate 1 coefficient 
#beta2 - covariate 2 coefficient 
#beta12 - interaction  coefficient  between covariate 1 and 2  
#beta3 - covariate 3 coefficient 
#t - upper bound for calculating RMST (t*)
#simulated_dataset - dataframe of observation times (generated using generate_ps)

#Output: variance of fully specified parametric RMST estimate (numeric)

param_rmst_var<- function(log_of_scale, beta1, beta2, beta12, beta3, t, simulated_dataset){
  h_scale<- function(log_scale, b_1, b_2, b_12, b_3){
    0.25*((exp(log_scale)/exp(b_1+b_2+b_12+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_2+b_12+b_3))) + 
            (exp(log_scale)/exp(b_1+b_2+b_12))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_2+b_12)))+
            (exp(log_scale)/exp(b_1+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1+b_3)))+
            (exp(log_scale)/exp(b_1))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_1)))) - 
      0.25*((exp(log_scale)/exp(b_2+b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_2+b_3))) + 
              (exp(log_scale)/exp(b_2))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_2))) + 
              (exp(log_scale)/exp(b_3))*(1-exp(-(1/exp(log_scale))*t_star*exp(b_3)))+
              (exp(log_scale))*(1-exp(-(1/exp(log_scale))*t_star)))
  }
  t_star<-t
  #Calculate Jacobian for this expression 
  grad_h_beta<-calculus::jacobian(h_scale, c(log_scale=log_of_scale, b_1=beta1, b_2=beta2, b_12=beta12, b_3=beta3))
  #Rearrange vector to gradient(beta_1, beta_2, log_scale)
  grad_h_beta<- c(grad_h_beta[2], grad_h_beta[3],grad_h_beta[5], grad_h_beta[4],grad_h_beta[1])
  #Calculate covariance matrix using phreg function 
  new<- phreg(Surv(surv_times,status)~treatment+inherit+sex+treatment*inherit, data=simulated_dataset, shape=1)
  covariance_phreg<-new$var
  #This covariance matrix has entries (1,1)=Var(beta_1), (2,2)=Var(beta_2), (3,3)=Var(log_scale)
  #Implement Delta method
  var<-t(grad_h_beta)%*%covariance_phreg%*%grad_h_beta
  return(var)
}

#####################################miss_param_rmst_var########################
#Description: Find variance for mis-speicified parametric RMST estimate 

#Input: 
#log_scale - log transformed base hazard (log(\lambda))
#beta1 - covariate 1 coefficient 
#beta2 - covariate 2 coefficient 
#beta12 - interaction  coefficient  between covariate 1 and 2  
#t - upper bound for calculating RMST (t*)
#simulated_dataset - dataframe of observation times (generated using generate_ps)

#Output: variance of mis-specified parametric RMST estimate (numeric)


miss_param_rmst_var<- function(log_of_scale, beta1, beta2,beta_12, t_star, simulated_dataset){
  area<- function(log_scale, b_1, b_2, b_12){
    area<- 0.5*((exp(log_scale)/(exp(b_1+b_2+b_12)))*(1-exp(-exp(-log_scale)*t_star*exp(b_1+b_2+b_12)))+
                  (exp(log_scale)/(exp(b_1)))*(1-exp(-exp(-log_scale)*t_star*exp(b_1))))  -
      0.5*((exp(log_scale)/(exp(b_2)))*(1-exp(-exp(-log_scale)*t_star*exp(b_2)))+
             exp(log_scale)*(1-exp(-exp(-log_scale)*t_star)))
  }
  #Calculate Jacobian for this expression 
  grad_h_beta<-calculus::jacobian(area, c(log_scale=log_of_scale, b_1=beta1, b_2=beta2, b_12=beta_12))
  #Rearrange vector to gradient(beta_1, log_scale)
  grad_h_beta<- c(grad_h_beta[2], grad_h_beta[3], grad_h_beta[4], grad_h_beta[1])
  #Calculate covariance matrix using phreg function 
  new<- phreg(Surv(surv_times,status)~treatment + inherit + treatment*inherit, data=simulated_dataset, shape=1)
  covariance_phreg<-new$var
  #This covariance matrix has entries (1,1)=Var(beta_1), (2,2)=Var(log_scale)
  #Implement Delta method
  var<-t(grad_h_beta)%*%covariance_phreg%*%grad_h_beta
  return(var)
}

######################################generate_z_values#########################
#Description: Generate z-values for the non-parametric, fully specified parametric
#and mis-specified parametric RMST estimate for one generated trial

#Input:
#sample_size - trial sample size
#base_hazard - exp(\lambda)
#vec_cov_coef - vector of covariate coefficients c(\beta_1, \beta_2, \beta_12, \beta_3)
#censoring rate - censoring parameter
#arrival_interval - time until recruitment ends (months)
#end_time - end of trial (months)
#t_star - time of RMST estimate (t*)

#Output: vector of z-values and their breakdown
generate_z_values<- function(sample_size, base_hazard, vec_cov_coef, censoring_rate, arrival_interval, end_time, t_star){
  data.set<-generate_ps(sample_size,1,vec_cov_coef[1],vec_cov_coef[1],vec_cov_coef[2],vec_cov_coef[2],vec_cov_coef[3], vec_cov_coef[4], log(base_hazard),log(base_hazard),censoring_rate,arrival_interval,end_time)
  data_columns <- apply(as.matrix(data.set[, c(2, 3, 5, 7)]), 2, as.numeric)
  while(nrow(unique(data_columns[data_columns[,4]==1,]))<8){
    data.set<-generate_ps(sample_size,1,vec_cov_coef[1],vec_cov_coef[1],vec_cov_coef[2],vec_cov_coef[2],vec_cov_coef[3], vec_cov_coef[4], log(base_hazard),log(base_hazard),censoring_rate,arrival_interval,end_time)
    data_columns <- apply(as.matrix(data.set[, c(2, 3, 5, 7)]), 2, as.numeric)
    }
  #Z-value for non-parametric RMST estimate 
  if(t_star<=min(max((data.set%>%filter(treatment==1))$surv_times), max((data.set%>%filter(treatment==0))$surv_times))){
    model.2 <- rmst2(time=data.set$surv_times,
                     status=data.set$status,
                     arm=data.set$treatment,
                     tau=t_star)
    RMST.val <- model.2$RMST.arm1$rmst[1]-model.2$RMST.arm0$rmst[1]
    RMST.var <- model.2$RMST.arm1$rmst[2]^2+model.2$RMST.arm0$rmst[2]^2
    z.nonparam <- RMST.val/sqrt(RMST.var)
  #Z-value for fully specified parametric RMST estimate 
  cox_model<- phreg(Surv(surv_times,status)~treatment+inherit+sex+treatment*inherit, data=data.set, shape=1)
  logscale<-as.numeric(cox_model$coefficients[5])
  b1<-as.numeric(cox_model$coefficients[1])
  b2<- as.numeric(cox_model$coefficients[2])
  b3<- as.numeric(cox_model$coefficients[3])
  b12<- as.numeric(cox_model$coefficients[4])
  full_te<-h_scale_t(logscale, b1, b2, b12, b3, t_star)
  full_se<-(param_rmst_var(logscale, b1,b2,b12, b3, t_star, data.set))^(1/2)
  z.param<- full_te/full_se
  #Z-value for miss specified parametric RMST estimate  
  cox_model.miss<- phreg(Surv(surv_times,status)~treatment+inherit+treatment*inherit, data=data.set, shape=1)
  logscale.miss<-as.numeric(cox_model.miss$coefficients[4])
  b1.miss<-as.numeric(cox_model.miss$coefficients[1])
  b2.miss<-as.numeric(cox_model.miss$coefficients[2])
  b12.miss<-as.numeric(cox_model.miss$coefficients[3])
  se.miss<-(miss_param_rmst_var(logscale.miss, b1.miss,b2.miss, b12.miss, t_star, data.set))^(1/2)
  #Find treatment effect using only cov_1, cov_2 and cov_3 parameter value
  treatment_effect<- 0.5*(area(t_star, exp(-(logscale.miss)), c(b1.miss, b2.miss, b12.miss),c(1,1,1))+area(t_star, exp(-(logscale.miss)), c(b1.miss, b2.miss, b12.miss),c(1,0,0)))-
    0.5*(area(t_star, exp((-logscale.miss)), c(b1.miss, b2.miss, b12.miss),c(0,1,0)) + area(t_star, exp(-(logscale.miss)), c(b1.miss, b2.miss, b12.miss),c(0,0,0)))
  z.miss.param<- treatment_effect/se.miss
  #cox ph model 
  cox<- phreg(Surv(surv_times, status)~treatment + inherit + sex+ treatment*inherit, data=data.set, shape=1)
  z.cox<- cox$coefficients[1]/sqrt(cox$var[1,1])
  return(c(z.nonparam, RMST.val, sqrt(RMST.var), z.param, full_te, full_se, z.miss.param, treatment_effect, se.miss, z.cox))
  }else{
  return(rep(NA, times=10))
  }
}





