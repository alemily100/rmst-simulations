library(pch)
library(numDeriv)
library(survRM2)
source("functions_crossing.R")


#Parameters to fix
N <- 1e5
n <- 200
t1 <- 40
beta2a <- 0.09437255
beta2b <- 0.09437255
beta3a <- -0.4021883
beta3b <- -0.4021883
h0a <- log(0.01577724)
h0b <- log(0.01577724)
lambda_censor <- 0.001
recruit_end <- 20
stop_time <- 90

#Varying parameters
t_star.grid <- seq(20,150,length = 30)
t1_miss.grid <- seq(15,60,length = 30)

#Parallelise
print("Make parallel")
library(parallel)
cl <- makeCluster(56)  
invisible(clusterEvalQ(cl,{
  library(pch)
  library(numDeriv)
  library(survRM2)
  source("functions_crossing.R")
}))
clusterExport(cl, c("n", "t1", "beta3a","beta3b", "beta2a", "beta2b",
                    "h0a", "h0b","lambda_censor", "recruit_end", "stop_time"))

#### Under HA, varying t*
beta1a <- -1.116749
beta1b <- 0.75
beta12a <- 0.4754451
beta12b <- 0.4754451
t1_miss <- 50
clusterExport(cl, c("beta1a", "beta1b", "beta12a", "beta12b", "t1_miss"))
result.data <- data.frame(t_star=NULL, t1 = NULL,
                          power.fullpara.dir1=NULL, power.misspara.dir1=NULL,
                          power.nonpara.dir1=NULL, power.fullpara.dir2=NULL,
                          power.misspara.dir2=NULL, power.nonpara.dir2=NULL)
print("Simulations under HA varying t*")
for(i in 1:length(t_star.grid)){
  
  #Set parameter values for this replicate
  t_star <- t_star.grid[i]
  clusterExport(cl, c("t_star"))
  
  #Simulation study at this set of parameters
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
    
    #Simulate patients with crossing survival curves
    surv.data <- generate_ps(n,t1,beta1a,beta1b,beta2a,beta2b,
                             beta3a,beta3b,beta12a,beta12b,
                             h0a,h0b,lambda_censor,recruit_end,stop_time)
    
    
    #Fully specified parametric model
    model.1 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.1$beta
    Sigma <- model.1$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.fullpara <- RMST.val/sqrt(RMST.var)
    
    #Miss-specified parametric model
    model.2 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1_miss, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.2$beta
    Sigma <- model.2$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1_miss)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1_miss)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.misspara <- RMST.val/sqrt(RMST.var)
    
    #Non-parametric RMST
    if(max(surv.data[surv.data$status==1&surv.data$treatment==1,
                     "surv_times"]) < t_star|
       max(surv.data[surv.data$status==1&surv.data$treatment==0,
                     "surv_times"]) < t_star){
      Z.nonpara <- NA
    }else{
      model.2 <- rmst2(time=surv.data$surv_times,
                       status=surv.data$status,
                       arm=surv.data$treatment,
                       tau=t_star)
      
      RMST.val <- model.2$RMST.arm1$rmst[1]-model.2$RMST.arm0$rmst[1]
      RMST.var <- model.2$RMST.arm1$rmst[2]^2+model.2$RMST.arm0$rmst[2]^2
      Z.nonpara <- RMST.val/sqrt(RMST.var)
    }
    
    return(c(Z.fullpara, Z.misspara, Z.nonpara))
    
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  
  #Gather results
  fullpara_power_dir1 <- mean(vals[1,]>qnorm(1-0.025))
  misspara_power_dir1 <- mean(vals[2,]>qnorm(1-0.025))
  non_power_dir1 <- mean(vals[3,]>qnorm(1-0.025))
  fullpara_power_dir2 <- mean(-vals[1,]>qnorm(1-0.025))
  misspara_power_dir2 <- mean(-vals[2,]>qnorm(1-0.025))
  non_power_dir2 <- mean(-vals[3,]>qnorm(1-0.025))
  result.data <- rbind(result.data, data.frame(
    t_star=t_star, t1 = t1_miss,
    power.fullpara.dir1=fullpara_power_dir1,
    power.misspara.dir1=misspara_power_dir1,
    power.nonpara.dir1=non_power_dir1,
    power.fullpara.dir2=fullpara_power_dir2,
    power.misspara.dir2=misspara_power_dir2,
    power.nonpara.dir2=non_power_dir2
  ))
}

#### Under HA, varying t1_miss
t_star <- 40
clusterExport(cl, "t_star")
print("Simulations under HA varying t1_miss")
for(i in 1:length(t1_miss.grid)){
  
  #Set parameter values for this replicate
  t1_miss <- t1_miss.grid[i]
  clusterExport(cl, c("t1_miss"))
  
  #Simulation study at this set of parameters
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
    
    #Simulate patients with crossing survival curves
    surv.data <- generate_ps(n,t1,beta1a,beta1b,beta2a,beta2b,
                             beta3a,beta3b,beta12a,beta12b,
                             h0a,h0b,lambda_censor,recruit_end,stop_time)
    
    
    #Fully specified parametric model
    model.1 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.1$beta
    Sigma <- model.1$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.fullpara <- RMST.val/sqrt(RMST.var)
    
    #Miss-specified parametric model
    model.2 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1_miss, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.2$beta
    Sigma <- model.2$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1_miss)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1_miss)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.misspara <- RMST.val/sqrt(RMST.var)
    
    #Non-parametric RMST
    if(max(surv.data[surv.data$status==1&surv.data$treatment==1,
                     "surv_times"]) < t_star|
       max(surv.data[surv.data$status==1&surv.data$treatment==0,
                     "surv_times"]) < t_star){
      Z.nonpara <- NA
    }else{
      model.2 <- rmst2(time=surv.data$surv_times,
                       status=surv.data$status,
                       arm=surv.data$treatment,
                       tau=t_star)
      
      RMST.val <- model.2$RMST.arm1$rmst[1]-model.2$RMST.arm0$rmst[1]
      RMST.var <- model.2$RMST.arm1$rmst[2]^2+model.2$RMST.arm0$rmst[2]^2
      Z.nonpara <- RMST.val/sqrt(RMST.var)
    }
    
    
    return(c(Z.fullpara, Z.misspara, Z.nonpara))
    
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  
  #Gather results
  fullpara_power_dir1 <- mean(vals[1,]>qnorm(1-0.025))
  misspara_power_dir1 <- mean(vals[2,]>qnorm(1-0.025))
  non_power_dir1 <- mean(vals[3,]>qnorm(1-0.025))
  fullpara_power_dir2 <- mean(-vals[1,]>qnorm(1-0.025))
  misspara_power_dir2 <- mean(-vals[2,]>qnorm(1-0.025))
  non_power_dir2 <- mean(-vals[3,]>qnorm(1-0.025))
  result.data <- rbind(result.data, data.frame(
    t_star=t_star, t1 = t1_miss,
    power.fullpara.dir1=fullpara_power_dir1,
    power.misspara.dir1=misspara_power_dir1,
    power.nonpara.dir1=non_power_dir1,
    power.fullpara.dir2=fullpara_power_dir2,
    power.misspara.dir2=misspara_power_dir2,
    power.nonpara.dir2=non_power_dir2
  ))
}
print("Save dataset")
write.csv(result.data, "cross_sims_HA.csv")


#### Under H0, varying t*
beta1a <- 0
beta1b <- 0
beta12a <- 0
beta12b <- 0
t1_miss <- 50
clusterExport(cl, c("beta1a", "beta1b", "beta12a", "beta12b", "t1_miss"))
result.data <- data.frame(t_star=NULL, t1 = NULL,
                          power.fullpara.dir1=NULL, power.misspara.dir1=NULL,
                          power.nonpara.dir1=NULL, power.fullpara.dir2=NULL,
                          power.misspara.dir2=NULL, power.nonpara.dir2=NULL)
print("Simulations under H0 varying t*")
for(i in 1:length(t_star.grid)){
  
  #Set parameter values for this replicate
  t_star <- t_star.grid[i]
  clusterExport(cl, c("t_star"))
  
  #Simulation study at this set of parameters
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
    
    #Simulate patients with crossing survival curves
    surv.data <- generate_ps(n,t1,beta1a,beta1b,beta2a,beta2b,
                             beta3a,beta3b,beta12a,beta12b,
                             h0a,h0b,lambda_censor,recruit_end,stop_time)
    
    
    #Fully specified parametric model
    model.1 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.1$beta
    Sigma <- model.1$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.fullpara <- RMST.val/sqrt(RMST.var)
    
    #Miss-specified parametric model
    model.2 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1_miss, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.2$beta
    Sigma <- model.2$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1_miss)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1_miss)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.misspara <- RMST.val/sqrt(RMST.var)
    
    #Non-parametric RMST
    if(max(surv.data[surv.data$status==1&surv.data$treatment==1,
                     "surv_times"]) < t_star|
       max(surv.data[surv.data$status==1&surv.data$treatment==0,
                     "surv_times"]) < t_star){
      Z.nonpara <- NA
    }else{
      model.2 <- rmst2(time=surv.data$surv_times,
                       status=surv.data$status,
                       arm=surv.data$treatment,
                       tau=t_star)
      
      RMST.val <- model.2$RMST.arm1$rmst[1]-model.2$RMST.arm0$rmst[1]
      RMST.var <- model.2$RMST.arm1$rmst[2]^2+model.2$RMST.arm0$rmst[2]^2
      Z.nonpara <- RMST.val/sqrt(RMST.var)
    }
    
    
    return(c(Z.fullpara, Z.misspara, Z.nonpara))
    
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  
  #Gather results
  fullpara_power_dir1 <- mean(vals[1,]>qnorm(1-0.025))
  misspara_power_dir1 <- mean(vals[2,]>qnorm(1-0.025))
  non_power_dir1 <- mean(vals[3,]>qnorm(1-0.025))
  fullpara_power_dir2 <- mean(-vals[1,]>qnorm(1-0.025))
  misspara_power_dir2 <- mean(-vals[2,]>qnorm(1-0.025))
  non_power_dir2 <- mean(-vals[3,]>qnorm(1-0.025))
  result.data <- rbind(result.data, data.frame(
    t_star=t_star, t1 = t1_miss,
    power.fullpara.dir1=fullpara_power_dir1,
    power.misspara.dir1=misspara_power_dir1,
    power.nonpara.dir1=non_power_dir1,
    power.fullpara.dir2=fullpara_power_dir2,
    power.misspara.dir2=misspara_power_dir2,
    power.nonpara.dir2=non_power_dir2
  ))
  
  write.csv(result.data, "cross_sims_H0.csv")
}

#### Under H0, varying t1_miss
t_star <- 40
clusterExport(cl, "t_star")
print("Simulations under HA varying t1_miss")
for(i in 1:length(t1_miss.grid)){
  
  #Set parameter values for this replicate
  t1_miss <- t1_miss.grid[i]
  clusterExport(cl, c("t1_miss"))
  
  #Simulation study at this set of parameters
  start_time <- Sys.time()
  vals <- parSapply(cl, 1:N, function(sim.rep){
    
    #Simulate patients with crossing survival curves
    surv.data <- generate_ps(n,t1,beta1a,beta1b,beta2a,beta2b,
                             beta3a,beta3b,beta12a,beta12b,
                             h0a,h0b,lambda_censor,recruit_end,stop_time)
    
    
    #Fully specified parametric model
    model.1 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.1$beta
    Sigma <- model.1$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.fullpara <- RMST.val/sqrt(RMST.var)
    
    #Miss-specified parametric model
    model.2 <- pch:::pchreg(Surv(surv_times, status)~treatment + inherit + gender + treatment*inherit, 
                            breaks=c(0, t1_miss, Inf), surv.data)
    #MLE and covariance matrix
    MLE.val <- model.2$beta
    Sigma <- model.2$covar
    if(dim(MLE.val)[2]==1){
      MLE.val <- cbind(MLE.val, MLE.val)
      Sigma <- rbind(cbind(Sigma, matrix(0, dim(Sigma)[1], dim(Sigma)[2])),
                     cbind(matrix(0, dim(Sigma)[1], dim(Sigma)[2]), Sigma))
    }
    #RMST
    RMST.val <- RMST.func(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                          MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                          MLE.val[1,1], MLE.val[1,2], t1_miss)
    d.RMST.val <- d.RMST(t_star, MLE.val[2,1], MLE.val[2,2], MLE.val[3,1], MLE.val[3,2], 
                         MLE.val[4,1], MLE.val[4,2], MLE.val[5,1], MLE.val[5,2],
                         MLE.val[1,1], MLE.val[1,2], t1_miss)
    RMST.var <- t(as.matrix(d.RMST.val, ncol = 1))%*%Sigma%*%as.matrix(d.RMST.val, ncol = 1)
    Z.misspara <- RMST.val/sqrt(RMST.var)
    
    #Non-parametric RMST
    if(max(surv.data[surv.data$status==1&surv.data$treatment==1,
                     "surv_times"]) < t_star|
       max(surv.data[surv.data$status==1&surv.data$treatment==0,
                     "surv_times"]) < t_star){
      Z.nonpara <- NA
    }else{
      model.2 <- rmst2(time=surv.data$surv_times,
                       status=surv.data$status,
                       arm=surv.data$treatment,
                       tau=t_star)
      
      RMST.val <- model.2$RMST.arm1$rmst[1]-model.2$RMST.arm0$rmst[1]
      RMST.var <- model.2$RMST.arm1$rmst[2]^2+model.2$RMST.arm0$rmst[2]^2
      Z.nonpara <- RMST.val/sqrt(RMST.var)
    }
    
    
    return(c(Z.fullpara, Z.misspara, Z.nonpara))
    
  })
  end_time <- Sys.time()
  print(end_time-start_time)
  
  #Gather results
  fullpara_power_dir1 <- mean(vals[1,]>qnorm(1-0.025))
  misspara_power_dir1 <- mean(vals[2,]>qnorm(1-0.025))
  non_power_dir1 <- mean(vals[3,]>qnorm(1-0.025))
  fullpara_power_dir2 <- mean(-vals[1,]>qnorm(1-0.025))
  misspara_power_dir2 <- mean(-vals[2,]>qnorm(1-0.025))
  non_power_dir2 <- mean(-vals[3,]>qnorm(1-0.025))
  result.data <- rbind(result.data, data.frame(
    t_star=t_star, t1 = t1_miss,
    power.fullpara.dir1=fullpara_power_dir1,
    power.misspara.dir1=misspara_power_dir1,
    power.nonpara.dir1=non_power_dir1,
    power.fullpara.dir2=fullpara_power_dir2,
    power.misspara.dir2=misspara_power_dir2,
    power.nonpara.dir2=non_power_dir2
  ))
  
  write.csv(result.data, "cross_sims_H0.csv")
}
print("Save dataset")
write.csv(result.data, "cross_sims_H0.csv")

