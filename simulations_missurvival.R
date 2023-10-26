library(survival)
library(tidyverse)
library(survminer)
library(Rlab)
library(calculus)
library(dplyr)
library(eha)
library(rgl)
library(survRM2)
library(parallel
        )
source(functions_interaction.R)

#######################RUN BEFORE EVERY CALCULATION#############################
#Before each power/type I error calculation, re-run this code chunk to ensure 
#all parameters are reset to original values   

#MOTIVATING EXAMPLE - CGD
data(cgd)
cgd<-cgd%>%filter(tstart==0)%>%select(id,tstop, status,treat, sex, inherit)

cgd$treat<-as.factor(cgd$treat)
cgd$sex<-as.factor(cgd$sex)
cgd$inherit<-as.factor(cgd$inherit)
cgd$time<- cgd$tstop/7
fit3 <- phreg(Surv(time, status) ~ treat
              + inherit +treat*inherit+sex, data=cgd, shape=1)

treat.coef<- fit3$coefficients[1]
inherit.coef<- fit3$coefficients[2]
treat.inherit.coef<- fit3$coefficients[4]
sex.coef<- fit3$coefficients[3]
base.hazard1<- exp(-fit3$coefficients[5])
cens<- 0.001
arrival<-20
end<-120

N<-100000
sample<-100
base_haz<- base.hazard1
vec_cov<-c(treat.coef, inherit.coef, treat.inherit.coef, sex.coef)


############################POWER CALCULATIONS##################################
################################################################################

################################VARYING t*######################################
cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl,{
  library(survival)
  library(tidyverse)
  library(survminer)
  library(Rlab)
  library(calculus)
  library(dplyr)
  library(eha)
  library(rgl)
  library(survRM2)
  library(parallel)
  setwd("~/MRC/FINAL")
  source("functions.r")
}))
clusterExport(cl, c("sample", "base_haz", "vec_cov", "cens",
                    "arrival", "end"))
result.data <- data.frame(t_star=NULL,
                          nonparam.power=NULL, missparam.power=NULL, param.power=NULL)
stat.data <- data.frame(t_star=NULL,
                          nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
t<-seq(from=1, to=150, length.out=15)
error<-c()
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(t)){
  t_star<-t[i]
  clusterExport(cl, c("t_star"))
  val<-parSapply(cl, rep(t_star, times=N),  function(k) generate_z_values(sample, base_haz, vec_cov, cens, arrival, end, k))
  result.data <- rbind(result.data, data.frame(t_star=t_star, nonparam.power=sum(val[1,]>qnorm(0.975),na.rm=TRUE)/(N-sum(is.na(val[1,]))),param.power=sum(val[4,]>qnorm(0.975))/N, missparam.power=sum(val[7,]>qnorm(0.975))/N))
  val_subset<-val[,val[6,]<1000]
  stat.data<- rbind(stat.data, data.frame(t_star=t_star, nonparam.z=mean(val_subset[1,], na.rm=TRUE),nonparam.se=mean(val_subset[3,],na.rm = TRUE), nonparam.te=mean(val_subset[2,],na.rm=TRUE),
                                          fullparam.z=mean(val_subset[4,]), fullparam.se=mean(val_subset[6,]), fullparam.te=mean(val_subset[5,]),
                                          missparam.z=mean(val_subset[7,]), missparam.se=mean(val_subset[9,]), missparam.te=mean(val_subset[8,])))
  no_z[i]<-sum(is.na(val[1,]))
  large_var[i]<-sum(val[6,]>=1000)
}
stopCluster(cl)
end_time<-Sys.time()
total_time<-end_time-start_time

write.csv(result.data, "power_varyt_star.csv")
write.csv(stat.data, "power_varyt_starbreakdown.csv")
write.csv(large_var, "power_varyt_starlargevariance.csv")
write.csv(no_z, "power_varyt_star_nonparamnotpossible.csv")

################################VARYING COV_3###################################
t_star<-100
cl <- makeCluster(detectCores()-1)
invisible(clusterEvalQ(cl,{
  library(survival)
  library(tidyverse)
  library(survminer)
  library(Rlab)
  library(calculus)
  library(dplyr)
  library(eha)
  library(rgl)
  library(survRM2)
  library(parallel)
  setwd("~/MRC/FINAL")
  source("functions.r")
}))
clusterExport(cl, c("sample", "treat.coef", "inherit.coef", "treat.inherit.coef", "base_haz", "cens",
                    "arrival", "end", "t_star"))
result.data <- data.frame(beta3=NULL,
                          nonparam.power=NULL, missparam.power=NULL)
stat.data <- data.frame(t_star=NULL,
                        nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
b<-seq(from=-2, to=2, length.out=10)
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(b)){
  beta3<-b[i]
  clusterExport(cl, c("beta3"))
  start_time <- Sys.time()
  val<-parSapply(cl, rep(beta3, times=N),  function(k) generate_z_values(sample, base_haz, c(treat.coef, inherit.coef, treat.inherit.coef, k), cens, arrival, end, t_star))
  result.data <- rbind(result.data, data.frame(beta3=beta3, nonparam.power=sum(val[1,]>qnorm(0.975),na.rm=TRUE)/(N-sum(is.na(val[1,]))),param.power=sum(val[4,]>qnorm(0.975),na.rm=TRUE)/N, missparam.power=sum(val[7,]>qnorm(0.975),na.rm=TRUE)/N))
  val_subset<-val[,val[6,]<1000]
  stat.data<- rbind(stat.data, data.frame(beta3=beta3, nonparam.z=mean(val_subset[1,], na.rm=TRUE),nonparam.se=mean(val_subset[3,],na.rm=TRUE), nonparam.te=mean(val_subset[2,],na.rm=TRUE),
                                          fullparam.z=mean(val_subset[4,]), fullparam.se=mean(val_subset[6,]), fullparam.te=mean(val_subset[5,]),
                                          missparam.z=mean(val_subset[7,]), missparam.se=mean(val_subset[9,]), missparam.te=mean(val_subset[8,])))
  large_var[i]<-sum(val[6,]>=1000)
  no_z[i]<-sum(is.na(val[1,]))
  }
end_time<-Sys.time()
total_time<-end_time-start_time

write.csv(result.data, "power_varycov_3.csv")
write.csv(stat.data, "power_varycov_3breakdown.csv")
write.csv(large_var, "power_varycov_3largevariance.csv")
write.csv(no_z, "power_varycov3_nonparamnotpossible.csv")

###########################TYPE I ERROR CALCULATIONS############################
################################################################################

##############################Varying t^*#######################################
vec_cov<-c(0, inherit.coef, 0, sex.coef)

cl <- makeCluster(detectCores())
invisible(clusterEvalQ(cl,{
  library(survival)
  library(tidyverse)
  library(survminer)
  library(Rlab)
  library(calculus)
  library(dplyr)
  library(eha)
  library(rgl)
  library(survRM2)
  library(parallel)
  setwd("~/MRC/FINAL")
  source("functions.r")
}))
clusterExport(cl, c("sample", "base_haz", "vec_cov", "cens",
                    "arrival", "end"))
result.data <- data.frame(t_star=NULL,
                          nonparam.power=NULL, missparam.power=NULL, param.power=NULL)
stat.data <- data.frame(t_star=NULL,
                        nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
t<-seq(from=0, to=150, length.out=15)
error<-c()
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(t)){
  t_star<-t[i]
  clusterExport(cl, c("t_star"))
  val<-parSapply(cl, rep(t_star, times=N),  function(k) generate_z_values(sample, base_haz, vec_cov, cens, arrival, end, k))
  #result.data <- rbind(result.data, data.frame(t_star=t_star, nonparam.power=sum(val[1,]>qnorm(0.975))/N,param.power=sum(val[2,]>qnorm(0.975))/N, missparam.power=sum(val[3,]>qnorm(0.975))/N))
  result.data <- rbind(result.data, data.frame(t_star=t_star, nonparam.power=sum(val[1,]>qnorm(0.975), na.rm = TRUE)/(N-sum(is.na(val[1,]))),param.power=sum(val[4,]>qnorm(0.975))/N, missparam.power=sum(val[7,]>qnorm(0.975))/N))
  val_subset<-val[,val[6,]<1000]
  stat.data<- rbind(stat.data, data.frame(t_star=t_star, nonparam.z=mean(val_subset[1,], na.rm=TRUE),nonparam.se=mean(val_subset[3,], na.rm=TRUE), nonparam.te=mean(val_subset[2,], na.rm=TRUE),
                                          fullparam.z=mean(val_subset[4,]), fullparam.se=mean(val_subset[6,]), fullparam.te=mean(val_subset[5,]),
                                          missparam.z=mean(val_subset[7,]), missparam.se=mean(val_subset[9,]), missparam.te=mean(val_subset[8,])))
  no_z[i]<-sum(is.na(val[1,]))
  large_var[i]<-sum(val[6,]>=1000)
}
stopCluster(cl)
end_time<-Sys.time()
total_time<-end_time-start_time

write.csv(result.data, "type1_varyt_star.csv")
write.csv(stat.data, "type1_varyt_starbreakdown.csv")
write.csv(large_var, "type1_varyt_starlargevariance.csv")
write.csv(no_z, "type1_vary_t_star_nonparamnotpossible.csv")


################################Varying beta3###################################
treat.coef<- 0
inherit.coef<- fit3$coefficients[2]
treat.inherit.coef<- 0
t_star<-100

cl <- makeCluster(detectCores()-1)
invisible(clusterEvalQ(cl,{
  library(survival)
  library(tidyverse)
  library(survminer)
  library(Rlab)
  library(calculus)
  library(dplyr)
  library(eha)
  library(rgl)
  library(survRM2)
  library(parallel)
  setwd("~/MRC/FINAL")
  source("functions.r")
}))
clusterExport(cl, c("sample", "treat.coef", "inherit.coef", "treat.inherit.coef", "base_haz", "cens",
                    "arrival", "end", "t_star"))
result.data <- data.frame(beta3=NULL,
                          nonparam.power=NULL, missparam.power=NULL)
stat.data <- data.frame(t_star=NULL,
                        nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
b<-seq(from=-2, to=2, length.out=10)
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(b)){
  beta3<-b[i]
  clusterExport(cl, c("beta3"))
  start_time <- Sys.time()
  val<-parSapply(cl, rep(beta3, times=N),  function(k) generate_z_values(sample, base_haz, c(treat.coef, inherit.coef, treat.inherit.coef, k), cens, arrival, end, t_star))
  result.data <- rbind(result.data, data.frame(beta3=beta3, nonparam.power=sum(val[1,]>qnorm(0.975), na.rm=TRUE)/(N-sum(is.na(val[1,]))),param.power=sum(val[4,]>qnorm(0.975))/N, missparam.power=sum(val[7,]>qnorm(0.975))/N))
  val_subset<-val[,val[6,]<1000]
  stat.data<- rbind(stat.data, data.frame(beta3=beta3, nonparam.z=mean(val_subset[1,], na.rm=TRUE),nonparam.se=mean(val_subset[3,], na.rm=TRUE), nonparam.te=mean(val_subset[2,], na.rm=TRUE),
                                          fullparam.z=mean(val_subset[4,]), fullparam.se=mean(val_subset[6,]), fullparam.te=mean(val_subset[5,]),
                                          missparam.z=mean(val_subset[7,]), missparam.se=mean(val_subset[9,]), missparam.te=mean(val_subset[8,])))
  large_var[i]<-sum(val[6,]>=1000)
  no_z[i]<-sum(is.na(val[1,]))
}
end_time<-Sys.time()
total_time<-end_time-start_time

write.csv(result.data, "type1_varycov_3.csv")
write.csv(stat.data, "type1_varycov_3breakdown.csv")
write.csv(large_var, "type1_varycov_3largevariance.csv")
write.csv(no_z, "type1_varycov_3_nonparamnotpossible.csv")