###########################generate_ps##########################################
#Description: generate time to observed failure for a sample of patients 

#Input:
#n - trial sample size
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)
#lambda_censor - censoring parameter
#recruit_end - time until recruitment ends (months)
#end_time - end of trial (months)

#Output: data frame containing patient ID, covariate indicators, time until 
#failure time, censoring indicator and arrival time within trial.

generate_ps <- function(n = 1e3,
                        t1 = 40,
                        beta1a = -1.116749,
                        beta1b = -1.116749,
                        beta2a = 0.09437255,
                        beta2b = 0.09437255,
                        beta3a = -0.4021883,
                        beta3b = -0.4021883,
                        beta12a = 0.4754451,
                        beta12b = 0.4754451,
                        h0a = log(0.01577724),
                        h0b = log(0.01577724),
                        lambda_censor=0.001,
                        recruit_end = 20,
                        end_time = 120){
  
  
  #fixed covariates
  Z <- cbind(rbinom(n, 1, prob = 0.5),
             rbinom(n, 1, prob = 0.5),
             rbinom(n, 1, prob = 0.5))
  
  #Generate survival times for each patient
  T_val <- sapply(1:n, function(i){
    
    #Cumulative hazard function
    H <- function(t){
      #Before knot-point, t1
      if(t<=t1){
        H_val <- t*exp(h0a+beta1a*Z[i,1]+beta2a*Z[i,2]+beta3a*Z[i,3]+beta12a*Z[i,1]*Z[i,2])
      }
      #After knot-point, t1
      if(t>t1){
        H_val <- t1*exp(h0a+beta1a*Z[i,1]+beta2a*Z[i,2]+beta3a*Z[i,3]+beta12a*Z[i,1]*Z[i,2]) +
          (t-t1)*exp(h0a+beta1b*Z[i,1]+beta2b*Z[i,2]+beta3b*Z[i,3]+beta12b*Z[i,1]*Z[i,2])
      }
      return(H_val)
    }
    
    
    #### inversion method e.g T = H^-1(u) ##########
    #Set a while loop to avoid having Ts=0 when generated uniform RV is tiny.
    Ts <- 0
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
                          gender = Z[,3],
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
}


###############################h.i###################################### 
#Description: Find hazard rate for patient i

#Input: 
#t - time
#Z - treatment indicator
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: hazard rate at time t (numeric)
h.i <- function(t, Z, beta1a, beta1b, beta2a, beta2b, 
                beta3a, beta3b, beta12a, beta12b,
                h0a, h0b, t1){
  #Before knot-point, t1
  if(t<=t1){
    h_val <- exp(h0a+beta1a*Z[1]+beta2a*Z[2]+beta3a*Z[3]+beta12a*Z[1]*Z[2])
  }
  #After knot-point, t1
  if(t>t1){
    h_val <- exp(h0a+beta1b*Z[1]+beta2b*Z[2]+beta3b*Z[3]+beta12b*Z[1]*Z[2])
  }
  return(h_val)
}


###############################H.i###################################### 
#Description: Find cumulative hazard for patient i

#Input: 
#t - time
#Z - treatment indicator
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: cumulative hazard at time t (numeric)
H.i <- function(t, Z, beta1a, beta1b, beta2a, beta2b, 
                beta3a, beta3b, beta12a, beta12b,
                h0a, h0b, t1){
  #Before knot-point, t1
  if(t<=t1){
    H_val <- t*exp(h0a+beta1a*Z[1]+beta2a*Z[2]+beta3a*Z[3]+beta12a*Z[1]*Z[2])
  }
  #After knot-point, t1
  if(t>t1){
    H_val <- t1*exp(h0a+beta1a*Z[1]+beta2a*Z[2]+beta3a*Z[3]+beta12a*Z[1]*Z[2]) +
      (t-t1)*exp(h0a+beta1b*Z[1]+beta2b*Z[2]+beta3b*Z[3]+beta12b*Z[1]*Z[2])
  }
  return(H_val)
}

###############################S.i###################################### 
#Description: Find survival probability for patient i

#Input: 
#t - time
#Z - treatment indicator
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: survival probability at time t (numeric)
S.i <- function(t, Z, beta1a, beta1b, beta2a, beta2b, 
                beta3a, beta3b, beta12a, beta12b,
                h0a, h0b, t1){
  exp(-H.i(t, Z, beta1a, beta1b, beta2a, beta2b, 
           beta3a, beta3b, beta12a, beta12b,
           h0a, h0b, t1))
}
Si.Vec <- Vectorize(S.i, vectorize.args = "t")

###############################f.i###################################### 
#Description: Find density function for patient i

#Input: 
#t - time
#Z - treatment indicator
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: density at time t (numeric)
f.i <- function(t, Z, beta1a, beta1b, beta2a, beta2b, 
                beta3a, beta3b, beta12a, beta12b,
                h0a, h0b, t1){
  h.i(t, Z, beta1a, beta1b, beta2a, beta2b, 
      beta3a, beta3b, beta12a, beta12b,
      h0a, h0b, t1)*
    S.i(t, Z, beta1a, beta1b, beta2a, beta2b, 
        beta3a, beta3b, beta12a, beta12b,
        h0a, h0b, t1)
}


###############################RMST.func###################################### 
#Description: Difference (between treatment and control arms) in parametric 
#restricted mean survival at time t*

#Input: 
#t_star - truncation time
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: RMST differencee at time t* (numeric)
RMST.func <- function(t_star, beta1a, beta1b, beta2a, beta2b, 
                      beta3a, beta3b, beta12a, beta12b,
                      h0a, h0b, t1){
  integrate(function(x){
    0.25*Si.Vec(x,c(1,0,0), beta1a, beta1b, beta2a, beta2b, 
                beta3a, beta3b, beta12a, beta12b,
                h0a, h0b, t1)+
      0.25*Si.Vec(x,c(1,0,1), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)+
      0.25*Si.Vec(x,c(1,1,0), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)+
      0.25*Si.Vec(x,c(1,1,1), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)-
      0.25*Si.Vec(x,c(0,0,0), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)-
      0.25*Si.Vec(x,c(0,0,1), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)-
      0.25*Si.Vec(x,c(0,1,0), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)-
      0.25*Si.Vec(x,c(0,1,1), beta1a, beta1b, beta2a, beta2b, 
                  beta3a, beta3b, beta12a, beta12b,
                  h0a, h0b, t1)
  }, lower = 0, upper = t_star)$value
}

###############################d.RMST###################################### 
#Description: Differentiate the RMST function with respect to vector of parameters

#Input: 
#t_star - truncation time
#beta1a - covariate 1 coefficient before knot point  
#beta1b - covariate 1 coefficient after knot point  
#beta2a - covariate 2 coefficient before knot point  
#beta2b - covariate 2 coefficient after knot point  
#beta3a - covariate 3  coefficient before knot point  
#beta3b - covariate 3  coefficient after knot point  
#beta12a - interaction  coefficient  between covariate 1 and 2 before knot point
#beta12b - interaction  coefficient  between covariate 1 and 2  after knot point
#t1 - knot point 
#h0a - baseline hazard before knot point (\lambda)
#h0b - baseline hazard after knot point  (\lambda)

#Output: gradient of RMST.func at time t* with respect to beta vector,
#vector of length 6
d.RMST <- function(t_star, beta1a, beta1b, beta2a, beta2b, 
                   beta3a, beta3b, beta12a, beta12b,
                   h0a, h0b, t1){
  grad(function(x){RMST.func(t_star, x[2], x[7], x[3], x[8], x[4],
                             x[9], x[5], x[10], x[1], x[6], t1)},
       x = c(h0a, beta1a, beta2a, beta3a, beta12a,
             h0b, beta1b, beta2b, beta3b, beta12b))
}

