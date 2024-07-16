setwd("~/Library/CloudStorage/OneDrive-TheAlanTuringInstitute/RMST_paper")
source("functions_missurvival.R")
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
#sex.coef<- fit3$coefficients[3]
sex.coef<- -1
base.hazard1<- exp(-fit3$coefficients[5])
cens<- 0.001
arrival<-20
end<-120
base_haz<- base.hazard1
vec_cov<-c(treat.coef, inherit.coef, treat.inherit.coef, sex.coef)

N<- 10000
sample<-100

############################POWER CALCULATIONS##################################
################################################################################


################################VARYING t*######################################
cl <- makeCluster(detectCores()-1)
clusterSetRNGStream(cl,1999)
invisible(clusterEvalQ(cl,{
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
  setwd("~/Library/CloudStorage/OneDrive-TheAlanTuringInstitute/RMST_paper")
  source("functions_missurvival.R")
}))
clusterExport(cl, c("sample", "base_haz", "vec_cov", "cens",
                    "arrival", "end"))
result.data <- data.frame(t_star=NULL,
                          nonparam.power=NULL, param.power=NULL, missparam.power=NULL, cox.power=NULL)
stat.data <- data.frame(t_star=NULL,
                          nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
t<-seq(from=1, to=120, length.out=20)
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(t)){
  #t_star<-t[i]
  t_star<-t[i]
  clusterExport(cl, c("t_star"))
  val<-parSapply(cl, rep(t_star, times=N),  function(k) generate_z_values(sample, base_haz, vec_cov, cens, arrival, end, k))
  n<- N-sum(is.na(val[4,]))
  te<- data.frame(val[1,], val[2,], val[3,], val[4,], val[5,], val[6,], val[7,], val[8,], val[9,])
  te[,10]<- rep(t_star, times =nrow(te))
  colnames(te)<-c("npz","npt","nps", "fpz","fpt", "fps", "mpz", "mpt", "mps","t")
  
  result.data <- rbind(result.data, data.frame(t_star=t_star, nonparam.power=sum(val[1,]>qnorm(0.975),na.rm=TRUE)/n,param.power=sum(val[4,]>qnorm(0.975),na.rm=TRUE)/n, missparam.power=sum(val[7,]>qnorm(0.975),na.rm=TRUE)/n, 
                                               cox.power=sum(val[10,]<qnorm(0.025),na.rm=TRUE)/n))
  no_z[i]<-N-n
  write.csv(te, paste0(i, "tpbreak.csv"))
}
stopCluster(cl)
write.csv(result.data, "tp.csv")
write.csv(no_z, "tpimpossible.csv")

################################VARYING COV_3###################################
t_star<-100
cl <- makeCluster(detectCores()-1)
clusterSetRNGStream(cl,1999)
invisible(clusterEvalQ(cl,{
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
  setwd("~/Library/CloudStorage/OneDrive-TheAlanTuringInstitute/RMST_paper")
  source("functions_missurvival.R")
}))
clusterExport(cl, c("sample", "treat.coef", "inherit.coef", "treat.inherit.coef", "base_haz", "cens",
                    "arrival", "end", "t_star"))
result.data <- data.frame(beta3=NULL,
                          nonparam.power=NULL, missparam.power=NULL)
stat.data <- data.frame(t_star=NULL,
                        nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
b<-seq(from=-2, to=2, length.out=20)
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(b)){
  beta3<-b[i]
  clusterExport(cl, c("beta3"))
  start_time <- Sys.time()
  val<-parSapply(cl, rep(beta3, times=N),  function(k) generate_z_values(sample, base_haz, c(treat.coef, inherit.coef, treat.inherit.coef, k), cens, arrival, end, t_star))
  te<- data.frame(val[1,], val[2,], val[3,], val[4,], val[5,], val[6,], val[7,], val[8,], val[9,])
  te[,10]<- rep(t_star, times =nrow(te))
  colnames(te)<-c("npz","npt","nps", "fpz","fpt", "fps", "mpz", "mpt", "mps","t")
  n<- N-sum(is.na(val[4,]))
  result.data <- rbind(result.data, data.frame(beta3=beta3, nonparam.power=sum(val[1,]>qnorm(0.975),na.rm=TRUE)/n,param.power=sum(val[4,]>qnorm(0.975),na.rm=TRUE)/n, missparam.power=sum(val[7,]>qnorm(0.975),na.rm=TRUE)/n,
                                               cox.power=sum(val[10,]<qnorm(0.025), na.rm=TRUE)/n))
  no_z[i]<-N-n
  write.csv(te, paste0(i, "bpbreak.csv"))
  }

write.csv(result.data, "bp.csv")
write.csv(no_z, "bpimpossible.csv")

###########################TYPE I ERROR CALCULATIONS############################
################################################################################

##############################Varying t^*#######################################
vec_cov<-c(0, inherit.coef, 0, sex.coef)

cl <- makeCluster(detectCores()-1)
clusterSetRNGStream(cl,1999)
invisible(clusterEvalQ(cl,{
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
  setwd("~/Library/CloudStorage/OneDrive-TheAlanTuringInstitute/RMST_paper")
  source("functions_missurvival.R")
}))
clusterExport(cl, c("sample", "base_haz", "vec_cov", "cens",
                    "arrival", "end"))
result.data <- data.frame(t_star=NULL,
                          nonparam.power=NULL, param.power=NULL, missparam.power=NULL, cox.power=NULL)
stat.data <- data.frame(t_star=NULL,
                        nonparam.z=NULL,nonparam.se=NULL, nonparam.te=NULL,
                        fullparam.z=NULL, fullparam.se=NULL, fullparam.te=NULL,
                        missparam.z=NULL, missparam.se=NULL, missparam.te=NULL)
t<-seq(from=1, to=120, length.out=20)
error<-c()
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(t)){
  t_star<-t[i]
  clusterExport(cl, c("t_star"))
  val<-parSapply(cl, rep(t_star, times=N),  function(k) generate_z_values(sample, base_haz, vec_cov, cens, arrival, end, k))
  te<- data.frame(val[1,], val[2,], val[3,], val[4,], val[5,], val[6,], val[7,], val[8,], val[9,])
  te[,10]<- rep(t_star, times =nrow(te))
  colnames(te)<-c("npz","npt","nps", "fpz","fpt", "fps", "mpz", "mpt", "mps","t")
  n<-N-sum(is.na(val[4,]))
  result.data <- rbind(result.data, data.frame(t_star=t_star, nonparam.power=sum(val[1,]>qnorm(0.975), na.rm = TRUE)/n, param.power=sum(val[4,]>qnorm(0.975), na.rm=TRUE)/n, missparam.power=sum(val[7,]>qnorm(0.975), na.rm=TRUE)/n, 
                                               cox.power=sum(val[10,]<qnorm(0.025), na.rm=TRUE)/n))
  no_z[i]<-N-n
  write.csv(te, paste0(i, "tt1break.csv"))
}
stopCluster(cl)
end_time<-Sys.time()
total_time<-end_time-start_time

write.csv(result.data, "tt1.csv")
write.csv(no_z, "tt1impossible.csv")

################################Varying beta3###################################
treat.coef<- 0
inherit.coef<- fit3$coefficients[2]
treat.inherit.coef<- 0
t_star<-100

cl <- makeCluster(detectCores()-1)
clusterSetRNGStream(cl,1999)
invisible(clusterEvalQ(cl,{
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
  setwd("~/Library/CloudStorage/OneDrive-TheAlanTuringInstitute/RMST_paper")
  source("functions_missurvival.R")
}))
clusterExport(cl, c("sample", "treat.coef", "inherit.coef", "treat.inherit.coef", "base_haz", "cens",
                    "arrival", "end", "t_star"))
result.data <- data.frame(beta3=NULL,
                          nonparam.power=NULL, param.power=NULL, missparam.power=NULL)
b<-seq(from=-2, to=2, length.out=20)
large_var<-c()
no_z<-c()
start_time<-Sys.time()
for(i in 1:length(b)){
  beta3<-b[i]
  clusterExport(cl, c("beta3"))
  start_time <- Sys.time()
  val<-parSapply(cl, rep(beta3, times=N),  function(k) generate_z_values(sample, base_haz, c(treat.coef, inherit.coef, treat.inherit.coef, k), cens, arrival, end, t_star))
  te<- data.frame(val[1,], val[2,], val[3,], val[4,], val[5,], val[6,], val[7,], val[8,], val[9,])
  te[,10]<- rep(t_star, times =nrow(te))
  colnames(te)<-c("npz","npt","nps", "fpz","fpt", "fps", "mpz", "mpt", "mps","t")
  n<- N-sum(is.na(val[4,]))
  result.data <- rbind(result.data, data.frame(beta3=beta3, nonparam.power=sum(val[1,]>qnorm(0.975), na.rm=TRUE)/n,param.power=sum(val[4,]>qnorm(0.975), na.rm=TRUE)/n, missparam.power=sum(val[7,]>qnorm(0.975), na.rm=TRUE)/n))
  no_z[i]<-N-n
  write.csv(te, paste0(i, "bt1break.csv"))
}

write.csv(result.data, "bt1.csv")
write.csv(no_z, "bt1impossible.csv")
