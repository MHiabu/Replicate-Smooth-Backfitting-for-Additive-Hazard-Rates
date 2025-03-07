library(mets)
library(parallel)
library(tikzDevice)
source("my-functions-backfit.r")
source("SBF_reg_Additive.R")
source("SBF_MH_LL_Additive.R")
source("SBF_MH_CLL.R")


source("predict_hazard.add.R")
source("predict_hazard.R")

### calc surv fucntion for fixed covariate points
source("surv_pred_add.R")
source("surv_pred_mult.R")


##################################################
### Running the trace data 
##################################################

### setting up the data  
data(TRACE)
dsort(TRACE) <- ~id
set.seed(100)
## make sure we do not have ties 
#TRACE$time <- TRACE$time+runif(nrow(TRACE))*0.01
TRACE$age0 <- TRACE$age
TRACE$age <- TRACE$age+TRACE$time
TRACE$time0 <- 0
## start age and start time=0
TRACE$age.null <- TRACE$age0
TRACE$time.null <- TRACE$time0
head(TRACE)
head(TRACE)
dcut(TRACE) <- ~age
TRACE$vff <- factor(TRACE$vf)




## decide where to look at model here box 
t0 <- 0
mt <- 5
a0 <- 40
ma <- 85




dtable(TRACE,~vf+chf+diabetes+sex,level=2)

## ct is data that is only using relevant box where t \in [0,5] a \in [40,85]
## model \beta(a) + \alpha(t) 
ct <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct$statusd <- ct$status!=0



ct$status[ct$status!=0] <- 1

ctsex1 <- ct[ct$vf==1,]
ctsex0 <- ct[ct$vf==0,]

ctsex11 <- ct[ct$vf==1&ct$chf==1,]
ctsex01 <- ct[ct$vf==0&ct$chf==1,]
ctsex10 <- ct[ct$vf==1&ct$chf==0,]
ctsex00 <- ct[ct$vf==0&ct$chf==0,]


xgrid <- list()
xgrid[[1]]<- sort(ct$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct$wmi) )
xgrid[[4]]<- unique(sort(ct$vf) )
xgrid[[5]]<- unique(sort(ct$extra) )

# 
# ####
# #### ################################################################
# #### first on full data
# ############################################################################

 alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
 
 
tikz( "full.tex" ,width = 4, height = 3)
#quartz(width = 11.3, height = 5.7)
par(mfrow = c(1, 1))
plot(alpha_backfit.LL$x.grid[[1]], alpha_backfit.LL$alpha_backfit[[1]], type="l", lwd=2, xlab="duration (years)", ylab="alpha1")
#quartz.save("full.pdf", type = "pdf")
dev.off()
 

# ####
# #### ################################################################
# #### Now considering first three month and later three months separately
# ############################################################################

# ############################################################################
# late part: get crps scores from 200 simulations with 80/20 training test split
# ############################################################################

t0 <- 1/4
mt <-5
a0 <- 40
ma <- 85
ct_late <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_late$statusd <- ct_late$status!=0
ct_late$status[ct_late$status!=0] <- 1
ct_latesex1 <- ct_late[ct_late$vf==1,]
ct_latesex0 <- ct_late[ct_late$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_late$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_late$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_late$wmi) )
xgrid[[4]]<- unique(sort(ct_late$vf) )
xgrid[[5]]<- unique(sort(ct_late$extra) )



sim_res_late <- mclapply(1:200,
                    function(s){
                      
                      
           
                      
                      train0 <- sample(1:nrow(ct_latesex0), 0.8*nrow(ct_latesex0))
                      train1 <- sample(1:nrow(ct_latesex1), 0.8*nrow(ct_latesex1))
                      ct_latesex0_train <- ct_latesex0[train0,] 
                      ct_latesex0_test <- ct_latesex0[-train0,] 
                      ct_latesex1_train <- ct_latesex1[train1,] 
                      ct_latesex1_test <- ct_latesex1[-train1,] 
                      
                      test0 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex0_train,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
                      test1 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex1_train,c(1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
                      test0.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex0_train,c(1.5,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
                      test1.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex1_train,c(1.5,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
                      
                      pred0_add <- predict.hazard.add(test0 ,cbind(ct_latesex0_test$time, ct_latesex0_test$age.null,ct_latesex0_test$wmi))
                      pred0_mult <- predict.hazard(test0.mult,cbind(ct_latesex0_test$time, ct_latesex0_test$age.null,ct_latesex0_test$wmi))
                      
                      pred1_add <- predict.hazard.add(test1 ,cbind(ct_latesex1_test$time, ct_latesex1_test$age.null,ct_latesex1_test$wmi))
                      pred1_mult <- predict.hazard(test1.mult,cbind(ct_latesex1_test$time, ct_latesex1_test$age.null,ct_latesex1_test$wmi))
                      
                      print(s)
                      return(c(mean(pred0_add$crps_1 +ct_latesex0_test$status* pred0_add$crps_2), mean(pred0_add$crps_1_pos +ct_latesex0_test$status* pred0_add$crps_2_pos)
                               , mean(pred0_mult$crps_1 +ct_latesex0_test$status* pred0_mult$crps_2),mean(pred1_add$crps_1 +ct_latesex1_test$status* pred1_add$crps_2), mean(pred1_add$crps_1_pos +ct_latesex1_test$status* pred1_add$crps_2_pos)
                               , mean(pred1_mult$crps_1 +ct_latesex1_test$status* pred1_mult$crps_2)
                      ))
                    },
                    mc.cores=1)


# ############################################################################
# early  part: get crps scores from 200 simulations with 80/20 training test split
# ############################################################################

t0 <-0
mt <-1/4
a0 <- 40
ma <- 85
ct_late <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_late$statusd <- ct_late$status!=0
ct_late$status[ct_late$status!=0] <- 1
ct_latesex1 <- ct_late[ct_late$vf==1,]
ct_latesex0 <- ct_late[ct_late$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_late$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_late$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_late$wmi) )
xgrid[[4]]<- unique(sort(ct_late$vf) )
xgrid[[5]]<- unique(sort(ct_late$extra) )



sim_res_early <- mclapply(1:200,
                    function(s){
                      
                      
                      
                      
                      train0 <- sample(1:nrow(ct_latesex0), 0.8*nrow(ct_latesex0))
                      train1 <- sample(1:nrow(ct_latesex1), 0.8*nrow(ct_latesex1))
                      ct_latesex0_train <- ct_latesex0[train0,] 
                      ct_latesex0_test <- ct_latesex0[-train0,] 
                      ct_latesex1_train <- ct_latesex1[train1,] 
                      ct_latesex1_test <- ct_latesex1[-train1,] 
                      
                      test0 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex0_train,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
                      test1 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex1_train,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
                      test0.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex0_train,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
                      test1.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex1_train,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
                      
                      pred0_add <-  predict.hazard.add(test0 ,cbind(ct_latesex0_test$time, ct_latesex0_test$age.null,ct_latesex0_test$wmi))
                      pred0_mult <- predict.hazard(test0.mult,cbind(ct_latesex0_test$time, ct_latesex0_test$age.null,ct_latesex0_test$wmi))
                      
                      pred1_add <- predict.hazard.add(test1 ,cbind(ct_latesex1_test$time, ct_latesex1_test$age.null,ct_latesex1_test$wmi))
                      pred1_mult <- predict.hazard(test1.mult,cbind(ct_latesex1_test$time, ct_latesex1_test$age.null,ct_latesex1_test$wmi))
                      
                      print(s)
                      return(c(mean(pred0_add$crps_1 +ct_latesex0_test$status* pred0_add$crps_2), mean(pred0_add$crps_1_pos +ct_latesex0_test$status* pred0_add$crps_2_pos)
                               , mean(pred0_mult$crps_1 +ct_latesex0_test$status* pred0_mult$crps_2),mean(pred1_add$crps_1 +ct_latesex1_test$status* pred1_add$crps_2), mean(pred1_add$crps_1_pos +ct_latesex1_test$status* pred1_add$crps_2_pos)
                               , mean(pred1_mult$crps_1 +ct_latesex1_test$status* pred1_mult$crps_2)
                      ))
                    },
                    mc.cores=10)

# ############################################################################
# Box plot for crps results
# ############################################################################

tikz( "boxplots.tex" ,
      width = 6.3, height = 3.15)
par(mfrow = c(1, 4))
par(mar = c(8,4,4,2) + 0.1)

bp <- boxplot(data.frame("additive"=as.numeric(lapply(sim_res_early, `[[`, 1)),"additive.adj"=as.numeric(lapply(sim_res_early, `[[`, 2)),"multiplicative"=as.numeric(lapply(sim_res_early, `[[`, 3))), xaxt = "n",main="vf=0, first three months", ylim=c(0.009,0.024), ylab="CRPS score")
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.004, bp$names, srt = 75, xpd = TRUE)



bp <- boxplot(data.frame("additive"=as.numeric(lapply(sim_res_early, `[[`, 4)),"additive.adj"=as.numeric(lapply(sim_res_early, `[[`, 5)),"multiplicative"=as.numeric(lapply(sim_res_early, `[[`, 6))), xaxt = "n",main="vf=1, first three months")
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.01, bp$names, srt = 75, xpd = TRUE)


bp <- boxplot(data.frame("additive"=as.numeric(lapply(sim_res_late, `[[`, 1)),"additive.adj"=as.numeric(lapply(sim_res_late, `[[`, 2)),"multiplicative"=as.numeric(lapply(sim_res_late, `[[`,3))), xaxt = "n", main="vf=0, after three months")
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.05, bp$names, srt = 75, xpd = TRUE)

bp <- boxplot(data.frame("additive"=as.numeric(lapply(sim_res_late, `[[`, 4)),"additive.adj"=as.numeric(lapply(sim_res_late, `[[`, 5)),"multiplicative"=as.numeric(lapply(sim_res_late, `[[`, 6))), xaxt = "n", main="vf=1, after three months")
tick <- seq_along(bp$names)
axis(1, at = tick, labels = FALSE)
text(tick, par("usr")[3] - 0.23, bp$names, srt = 75, xpd = TRUE)
dev.off()


######################################################################
################## splitting data late/early without train/test
####################################################################################################################################################################################

################## first early 

t0 <- 0
mt <-1/4
a0 <- 40
ma <- 85

ct_early <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_early$statusd <- ct_early$status!=0




ct_early$status[ct_early$status!=0] <- 1

ct_earlysex1 <- ct_early[ct_early$vf==1,]
ct_earlysex0 <- ct_early[ct_early$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_early$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_early$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_early$wmi) )
xgrid[[4]]<- unique(sort(ct_early$vf) )
xgrid[[5]]<- unique(sort(ct_early$extra) )




## bandwidth for later part: 1: c(1,20,0.8), 0: c(1,15,0.8), mult: 1: c(1.5,20,0.8), 0: c(1.5,15,0.8)
## bandwidth for earlier part: same but duration is 0.1 for all

test1_early <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_earlysex1,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
test0_early <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_earlysex0,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)


test1_early.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_earlysex1,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
test0_early.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_earlysex0,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

##### make plots


tikz( "early_add.tex" ,
      width = 6.3, height = 4)
par(mfrow=c(2,3))

plot(test0_early$x.grid[[1]],test0_early$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)",ylim=c(min(-3*sqrt(test0_early$err[[1]])+test0_early$alpha_backfit[[1]],na.rm=T), max(3*sqrt(test0_early$err[[1]])+test0_early$alpha_backfit[[1]],na.rm=T)))
lines(test0_early$x.grid[[1]],1.96*sqrt(test0_early$err[[1]])+test0_early$alpha_backfit[[1]], lty=2)
lines(test0_early$x.grid[[1]],-1.96*sqrt(test0_early$err[[1]])+test0_early$alpha_backfit[[1]], lty=2)


plot(test0_early$x.grid[[2]],test0_early$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)",ylim=c(min(-3*sqrt(test0_early$err[[2]])+test0_early$alpha_backfit[[2]],na.rm=T), max(3*sqrt(test0_early$err[[2]])+test0_early$alpha_backfit[[2]],na.rm=T)), main="vf=0")
lines(test0_early$x.grid[[2]],1.96*sqrt(test0_early$err[[2]])+test0_early$alpha_backfit[[2]], lty=2)
lines(test0_early$x.grid[[2]],-1.96*sqrt(test0_early$err[[2]])+test0_early$alpha_backfit[[2]], lty=2)

plot(test0_early$x.grid[[3]],test0_early$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi",ylim=c(min(-3*sqrt(test0_early$err[[3]])+test0_early$alpha_backfit[[3]],na.rm=T), max(3*sqrt(test0_early$err[[3]])+test0_early$alpha_backfit[[3]],na.rm=T)))
lines(test0_early$x.grid[[3]],1.96*sqrt(test0_early$err[[3]])+test0_early$alpha_backfit[[3]], lty=2)
lines(test0_early$x.grid[[3]],-1.96*sqrt(test0_early$err[[3]])+test0_early$alpha_backfit[[3]], lty=2)


plot(test1_early$x.grid[[1]],test1_early$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)",ylim=c(min(-3*sqrt(test1_early$err[[1]])+test1_early$alpha_backfit[[1]],na.rm=T), max(3*sqrt(test1_early$err[[1]])+test1_early$alpha_backfit[[1]],na.rm=T)))
lines(test1_early$x.grid[[1]],1.96*sqrt(test1_early$err[[1]])+test1_early$alpha_backfit[[1]], lty=2)
lines(test1_early$x.grid[[1]],-1.96*sqrt(test1_early$err[[1]])+test1_early$alpha_backfit[[1]], lty=2)

plot(test1_early$x.grid[[2]],test1_early$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)",ylim=c(min(-3*sqrt(test1_early$err[[2]])+test1_early$alpha_backfit[[2]],na.rm=T), max(3*sqrt(test1_early$err[[2]])+test1_early$alpha_backfit[[2]],na.rm=T)), main="vf=1")
lines(test1_early$x.grid[[2]],1.96*sqrt(test1_early$err[[2]])+test1_early$alpha_backfit[[2]], lty=2)
lines(test1_early$x.grid[[2]],-1.96*sqrt(test1_early$err[[2]])+test1_early$alpha_backfit[[2]], lty=2)

plot(test1_early$x.grid[[3]],test1_early$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi",ylim=c(min(-3*sqrt(test1_early$err[[3]])+test1_early$alpha_backfit[[3]],na.rm=T), max(3*sqrt(test1_early$err[[3]])+test1_early$alpha_backfit[[3]],na.rm=T)))
lines(test1_early$x.grid[[3]],1.96*sqrt(test1_early$err[[3]])+test1_early$alpha_backfit[[3]], lty=2)
lines(test1_early$x.grid[[3]],-1.96*sqrt(test1_early$err[[3]])+test1_early$alpha_backfit[[3]], lty=2)
#mtext("Hazard function estimates conditional on surviving the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title
mtext("Hazard function estimates for the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title
dev.off()





tikz( "early_mult.tex" ,
      width = 6.3, height = 3.15)
par(mfrow=c(2,3))

plot(test0_early.mult$x.grid[[1]],test0_early.mult$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)")


plot(test0_early.mult$x.grid[[2]],test0_early.mult$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)", main="vf=0")

plot(test0_early.mult$x.grid[[3]],test0_early.mult$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi")


plot(test1_early.mult$x.grid[[1]],test1_early.mult$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)")

plot(test1_early.mult$x.grid[[2]],test1_early.mult$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)", main="vf=1")

plot(test1_early.mult$x.grid[[3]],test1_early.mult$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi")
mtext("Hazard function estimates for the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title


dev.off()


########### ########### ########### ########### ########### ########### ########### 
########### now late
########### ########### ########### ########### ########### 

t0 <- 1/4
mt <-5
a0 <- 40
ma <- 85

ct_late <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_late$statusd <- ct_late$status!=0




##### now sbf


ct_late$status[ct_late$status!=0] <- 1

ct_latesex1 <- ct_late[ct_late$vf==1,]
ct_latesex0 <- ct_late[ct_late$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_late$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_late$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_late$wmi) )
xgrid[[4]]<- unique(sort(ct_late$vf) )
xgrid[[5]]<- unique(sort(ct_late$extra) )




####on full data
## bandwidth for later part: 1: c(1,20,0.8), 0: c(1,15,0.8), mult: 1: c(1.5,20,0.8), 0: c(1.5,15,0.8)
## bandwidth for earlier part: same but duration is 0.1 for all

test1_late <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex1,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
test0_late <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex0,c(1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=FALSE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)


test1_late.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex1,c(1.5,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
test0_late.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex0,c(1.5,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

##### make plots


tikz( "late_add.tex" ,
      width = 6.3, height = 4)
par(mfrow=c(2,3))

plot(test0_late$x.grid[[1]],test0_late$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)",ylim=c(min(-3*sqrt(test0_late$err[[1]])+test0_late$alpha_backfit[[1]],na.rm=T), max(3*sqrt(test0_late$err[[1]])+test0_late$alpha_backfit[[1]],na.rm=T)))
lines(test0_late$x.grid[[1]],1.96*sqrt(test0_late$err[[1]])+test0_late$alpha_backfit[[1]], lty=2)
lines(test0_late$x.grid[[1]],-1.96*sqrt(test0_late$err[[1]])+test0_late$alpha_backfit[[1]], lty=2)


plot(test0_late$x.grid[[2]],test0_late$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)",ylim=c(min(-3*sqrt(test0_late$err[[2]])+test0_late$alpha_backfit[[2]],na.rm=T), max(3*sqrt(test0_late$err[[2]])+test0_late$alpha_backfit[[2]],na.rm=T)), main="vf=0")
lines(test0_late$x.grid[[2]],1.96*sqrt(test0_late$err[[2]])+test0_late$alpha_backfit[[2]], lty=2)
lines(test0_late$x.grid[[2]],-1.96*sqrt(test0_late$err[[2]])+test0_late$alpha_backfit[[2]], lty=2)

plot(test0_late$x.grid[[3]],test0_late$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi",ylim=c(min(-3*sqrt(test0_late$err[[3]])+test0_late$alpha_backfit[[3]],na.rm=T), max(3*sqrt(test0_late$err[[3]])+test0_late$alpha_backfit[[3]],na.rm=T)))
lines(test0_late$x.grid[[3]],1.96*sqrt(test0_late$err[[3]])+test0_late$alpha_backfit[[3]], lty=2)
lines(test0_late$x.grid[[3]],-1.96*sqrt(test0_late$err[[3]])+test0_late$alpha_backfit[[3]], lty=2)


plot(test1_late$x.grid[[1]],test1_late$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)",ylim=c(min(-3*sqrt(test1_late$err[[1]])+test1_late$alpha_backfit[[1]],na.rm=T), max(3*sqrt(test1_late$err[[1]])+test1_late$alpha_backfit[[1]],na.rm=T)))
lines(test1_late$x.grid[[1]],1.96*sqrt(test1_late$err[[1]])+test1_late$alpha_backfit[[1]], lty=2)
lines(test1_late$x.grid[[1]],-1.96*sqrt(test1_late$err[[1]])+test1_late$alpha_backfit[[1]], lty=2)

plot(test1_late$x.grid[[2]],test1_late$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)",ylim=c(min(-3*sqrt(test1_late$err[[2]])+test1_late$alpha_backfit[[2]],na.rm=T), max(3*sqrt(test1_late$err[[2]])+test1_late$alpha_backfit[[2]],na.rm=T)), main="vf=1")
lines(test1_late$x.grid[[2]],1.96*sqrt(test1_late$err[[2]])+test1_late$alpha_backfit[[2]], lty=2)
lines(test1_late$x.grid[[2]],-1.96*sqrt(test1_late$err[[2]])+test1_late$alpha_backfit[[2]], lty=2)

plot(test1_late$x.grid[[3]],test1_late$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi",ylim=c(min(-3*sqrt(test1_late$err[[3]])+test1_late$alpha_backfit[[3]],na.rm=T), max(3*sqrt(test1_late$err[[3]])+test1_late$alpha_backfit[[3]],na.rm=T)))
lines(test1_late$x.grid[[3]],1.96*sqrt(test1_late$err[[3]])+test1_late$alpha_backfit[[3]], lty=2)
lines(test1_late$x.grid[[3]],-1.96*sqrt(test1_late$err[[3]])+test1_late$alpha_backfit[[3]], lty=2)
#mtext("Hazard function estimates conditional on surviving the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title
mtext("Hazard function estimates conditional on surviving the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title
dev.off()





tikz( "late_mult.tex" ,
      width = 6.3, height = 3.15)
par(mfrow=c(2,3))

plot(test0_late.mult$x.grid[[1]],test0_late.mult$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)")


plot(test0_late.mult$x.grid[[2]],test0_late.mult$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)", main="vf=0")

plot(test0_late.mult$x.grid[[3]],test0_late.mult$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi")


plot(test1_late.mult$x.grid[[1]],test1_late.mult$alpha_backfit[[1]],type="l",lwd=2, col="grey",ylab="$alpha_1$",xlab="duration (years)")

plot(test1_late.mult$x.grid[[2]],test1_late.mult$alpha_backfit[[2]],type="l",lwd=2, col="grey",ylab="$alpha_2$", xlab="age (years)", main="vf=1")

plot(test1_late.mult$x.grid[[3]],test1_late.mult$alpha_backfit[[3]],type="l",lwd=2, col="grey",ylab="$alpha_3$", xlab="wmi")
mtext("Hazard function estimates conditional on surviving the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) # added mtext to show how I would use it to create a title

dev.off()
########### till here split part



######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### simulation from multiplicative fit
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### 

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ct_earlysex0

t0 <- 0
mt <-1/4
a0 <- 40
ma <- 85

ct_early <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_early$statusd <- ct_early$status!=0
ct_early$status[ct_early$status!=0] <- 1

ct_earlysex1 <- ct_early[ct_early$vf==1,]
ct_earlysex0 <- ct_early[ct_early$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_early$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_early$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_early$wmi) )
xgrid[[4]]<- unique(sort(ct_early$vf) )
xgrid[[5]]<- unique(sort(ct_early$extra) )




test0.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_earlysex0,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

test0 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_earlysex0,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)



xxgrid<-c(test0.mult$x.grid[[1]],6)
dx<-lapply(1:length(test0.mult$x.grid),function(k){  ddx<-diff(test0.mult$x.grid[[k]])
c(  (test0.mult$x.grid[[k]][1]-x.min[k])  ,ddx)
}
)

###### start sim

res<- mclapply(1:200, function(s)
{
x_index <- sample(nrow(ct_earlysex0),10000, replace=TRUE)
x <- cbind(ct_earlysex0$age.null[x_index],ct_earlysex0$wmi[x_index])

t <-numeric(10000)
for(i in 1:10000){


pred0_surv_mult80 <- predict.surv.mult(test0.mult,cbind(ct_earlysex0$time, ct_earlysex0$age.null,ct_earlysex0$wmi),test0.mult$x.grid[[1]],x[i,])

haz<- predict.hazard(test0.mult, cbind(test0.mult$x.grid[[1]][1],x[i,1], x[i,2]))
haz <- prod(haz$hazard[,2:3])
haz <- test0.mult$alpha_backfit[[1]]*haz

#probs <- -diff(pred0_surv_mult80 )
probs <- pred0_surv_mult80[-1]*haz[-1]*dx[[1]][-1]
probs <- probs*((1-tail(pred0_surv_mult80,1))/(sum(probs)))
probs <- c(probs, 1-sum(probs))
t_index <- sample(length(test0.mult$x.grid[[1]]),1,prob=probs)
t[i] <- test0.mult$x.grid[[1]][t_index]+ runif(1,0,xxgrid[t_index+1]-xxgrid[t_index])
}


mydata <- data.frame(time=t,age.null=x[,1], wmi=x[,2])
mydata$status <- (mydata$time<1/4)*1
mydata$time <- pmin(mydata$time,1/4)

mydata1<- mydata
mydata <- mydata[1:1655,]

test0.mult_A <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

test0_A <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)

test0_AB <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=TRUE)


mydata1$time <- sapply(mydata1$time, function(l) runif(1,0,l))

myhazard <- predict.hazard(test0.mult, cbind(mydata1$time,mydata1$age.null,mydata1$wmi))
mydata1$haz <- apply(myhazard$hazard[,1:3],  1, prod)  


res_reg <- SBF.reg.LL(haz~time+age.null+wmi,data=mydata1,bandwidth=c(0.1,15,0.8),x.grid=xgrid,n.grid.additional=200, integral.approx='midd',it=15,kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=FALSE)

return(list(res_reg=res_reg,test0_A=test0_A,test0.mult=test0.mult_A,test0_AB=test0_AB ))
},
mc.cores=8)






######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ct_earlysex1


test1.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_earlysex1,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

test1 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_earlysex1,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)



xxgrid<-c(test1.mult$x.grid[[1]],6)
dx<-lapply(1:length(test1.mult$x.grid),function(k){  ddx<-diff(test1.mult$x.grid[[k]])
c(  (test1.mult$x.grid[[k]][1]-x.min[k])  ,ddx)
}
)

###### start sim

res1<- mclapply(1:200, function(s)
{
  x_index <- sample(nrow(ct_earlysex1),10000, replace=TRUE)
  x <- cbind(ct_earlysex1$age.null[x_index],ct_earlysex1$wmi[x_index])
  
  t <-numeric(10000)
  for(i in 1:10000){
    
    
    pred1_surv_mult80 <- predict.surv.mult(test1.mult,cbind(ct_earlysex1$time, ct_earlysex1$age.null,ct_earlysex1$wmi),test1.mult$x.grid[[1]],x[i,])
    
    haz<- predict.hazard(test1.mult, cbind(test1.mult$x.grid[[1]][1],x[i,1], x[i,2]))
    haz <- prod(haz$hazard[,2:3])
    haz <- test1.mult$alpha_backfit[[1]]*haz
    
    #probs <- -diff(pred0_surv_mult80 )
    probs <- pred1_surv_mult80[-1]*haz[-1]*dx[[1]][-1]
    probs <- probs*((1-tail(pred1_surv_mult80,1))/(sum(probs)))
    probs <- c(probs, 1-sum(probs))
    t_index <- sample(length(test1.mult$x.grid[[1]]),1,prob=probs)
    t[i] <- test1.mult$x.grid[[1]][t_index]+ runif(1,0,xxgrid[t_index+1]-xxgrid[t_index])
  }
  
  
  mydata <- data.frame(time=t,age.null=x[,1], wmi=x[,2])
  mydata$status <- (mydata$time<1/4)*1
  mydata$time <- pmin(mydata$time,1/4)
  
  mydata1<- mydata
  mydata <- mydata[1:132,]
  
  test1.mult_A <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
  
  test1_A <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
  
  test1_AB <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(0.1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=TRUE)
  
  
  
  mydata1$time <- sapply(mydata1$time, function(l) runif(1,0,l))
  
  myhazard <- predict.hazard(test1.mult, cbind(mydata1$time,mydata1$age.null,mydata1$wmi))
  mydata1$haz <- apply(myhazard$hazard[,1:3],  1, prod)  
  #  mydata2 <- mydata1[mydata1$status==1,]

  
  
  res_reg <- SBF.reg.LL(haz~time+age.null+wmi,data=mydata1,bandwidth=c(0.05,15,0.8),x.grid=xgrid,n.grid.additional=200, integral.approx='midd',it=15,kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=FALSE)
  
  return(list(res_reg=res_reg,test1_A=test1_A,test1.mult=test1.mult_A,test1_AB=test1_AB ))
},
mc.cores=8)


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ct_latesex0

t0 <- 1/4
mt <-5
a0 <- 40
ma <- 85

ct_late <- timescale.box(t0,mt,a0,ma,TRACE,gplot=1)
ct_late$statusd <- ct_late$status!=0
ct_late$status[ct_late$status!=0] <- 1

ct_latesex1 <- ct_late[ct_late$vf==1,]
ct_latesex0 <- ct_late[ct_late$vf==0,]

xgrid <- list()
xgrid[[1]]<- sort(ct_late$time); xgrid[[1]]<- unique(xgrid[[1]][xgrid[[1]]<mt&xgrid[[1]]>t0]) #
xgrid[[2]]<- sort(ct_late$age.null); xgrid[[2]]<- unique(xgrid[[2]][xgrid[[2]]<ma&xgrid[[2]]>a0]) #
xgrid[[3]]<- unique(sort(ct_late$wmi) )
xgrid[[4]]<- unique(sort(ct_late$vf) )
xgrid[[5]]<- unique(sort(ct_late$extra) )




test0.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex0,c(1.5,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

test0 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex0,c(1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)



xxgrid<-c(test0.mult$x.grid[[1]],6)
dx<-lapply(1:length(test0.mult$x.grid),function(k){  ddx<-diff(test0.mult$x.grid[[k]])
c(  (test0.mult$x.grid[[k]][1]-x.min[k])  ,ddx)
}
)

###### start sim

res2<- mclapply(1:200, function(s)
{
  x_index <- sample(nrow(ct_latesex0),10000, replace=TRUE)
  x <- cbind(ct_latesex0$age.null[x_index],ct_latesex0$wmi[x_index])
  
  t <-numeric(10000)
  for(i in 1:10000){
    
    
    pred0_surv_mult80 <- predict.surv.mult(test0.mult,cbind(ct_latesex0$time, ct_latesex0$age.null,ct_latesex0$wmi),test0.mult$x.grid[[1]],x[i,])
    
    haz<- predict.hazard(test0.mult, cbind(test0.mult$x.grid[[1]][1],x[i,1], x[i,2]))
    haz <- prod(haz$hazard[,2:3])
    haz <- test0.mult$alpha_backfit[[1]]*haz
    
    #probs <- -diff(pred0_surv_mult80 )
    probs <- pred0_surv_mult80[-1]*haz[-1]*dx[[1]][-1]
    probs <- probs*((1-tail(pred0_surv_mult80,1))/(sum(probs)))
    probs <- c(probs, 1-sum(probs))
    t_index <- sample(length(test0.mult$x.grid[[1]]),1,prob=probs)
    t[i] <- test0.mult$x.grid[[1]][t_index]+ runif(1,0,xxgrid[t_index+1]-xxgrid[t_index])
  }
  
  
  mydata <- data.frame(time=t,age.null=x[,1], wmi=x[,2])
  mydata$status <- (mydata$time<5)*1
  mydata$time <- pmin(mydata$time,5)
  
  mydata1<- mydata
  mydata <- mydata[1:1482,]
  
  test0.mult_A <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,mydata,c(1.5,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
  
  test0_A <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
  
  test0_AB <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(1,15,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=TRUE)
  
  
  
  mydata1$time <- sapply(mydata1$time, function(l) runif(1,1/4,l))
  
  myhazard <- predict.hazard(test0.mult, cbind(mydata1$time,mydata1$age.null,mydata1$wmi))
  mydata1$haz <- apply(myhazard$hazard[,1:3],  1, prod)  
  
  
  res_reg <- SBF.reg.LL(haz~time+age.null+wmi,data=mydata1,bandwidth=c(1,15,0.8),x.grid=xgrid,n.grid.additional=200, integral.approx='midd',it=15,kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=FALSE)
  
  return(list(res_reg=res_reg,test0_A=test0_A,test0.mult=test0.mult_A,test0_AB=test0_AB ))
},
mc.cores=200)



######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ct_latesex1


test1.mult <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,ct_latesex1,c(1.5,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)

test1 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex1,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)

test1 <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,ct_latesex1,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=TRUE)



xxgrid<-c(test1.mult$x.grid[[1]],6)
dx<-lapply(1:length(test1.mult$x.grid),function(k){  ddx<-diff(test1.mult$x.grid[[k]])
c(  (test1.mult$x.grid[[k]][1]-x.min[k])  ,ddx)
}
)


res3<- mclapply(1:200, function(s)
{
  x_index <- sample(nrow(ct_latesex1),10000, replace=TRUE)
  x <- cbind(ct_latesex1$age.null[x_index],ct_latesex1$wmi[x_index])
  
  t <-numeric(10000)
  for(i in 1:10000){
    
    
    pred1_surv_mult80 <- predict.surv.mult(test1.mult,cbind(ct_latesex1$time, ct_latesex1$age.null,ct_latesex1$wmi),test1.mult$x.grid[[1]],x[i,])
    
    haz<- predict.hazard(test1.mult, cbind(test1.mult$x.grid[[1]][1],x[i,1], x[i,2]))
    haz <- prod(haz$hazard[,2:3])
    haz <- test1.mult$alpha_backfit[[1]]*haz
    
    #probs <- -diff(pred0_surv_mult80 )
    probs <- pred1_surv_mult80[-1]*haz[-1]*dx[[1]][-1]
    probs <- probs*((1-tail(pred1_surv_mult80,1))/(sum(probs)))
    probs <- c(probs, 1-sum(probs))
    t_index <- sample(length(test1.mult$x.grid[[1]]),1,prob=probs)
    t[i] <- test1.mult$x.grid[[1]][t_index]+ runif(1,0,xxgrid[t_index+1]-xxgrid[t_index])
  }
  
  
  mydata <- data.frame(time=t,age.null=x[,1], wmi=x[,2])
  mydata$status <- (mydata$time<5)*1
  mydata$time <- pmin(mydata$time,5)
  
  mydata1<- mydata
  mydata <- mydata[1:75,]
  
  test1.mult_A <- alpha_backfit.LL<-SBF.MH.CLL(Surv(time, status) ~ age.null+wmi,mydata,c(1.5,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=xgrid,integral.approx='midd',kcorr=TRUE,LC=TRUE)
  
  test1_A <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
  
  test1_AB <- alpha_backfit.LL<-SBF.MH.CLL.add(Surv(time, status) ~ age.null+wmi,mydata,c(1,20,0.8),weight='sw',it=15,n.grid.additional=200,x.grid=NULL,integral.approx='right',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=TRUE)
  
  
  
  mydata1$time <- sapply(mydata1$time, function(l) runif(1,1/4,l))
  
  myhazard <- predict.hazard(test1.mult, cbind(mydata1$time,mydata1$age.null,mydata1$wmi))
  mydata1$haz <- apply(myhazard$hazard[,1:3],  1, prod)  
  #  mydata2 <- mydata1[mydata1$status==1,]
  
  
  
  res_reg <- SBF.reg.LL(haz~time+age.null+wmi,data=mydata1,bandwidth=c(1,20,0.8),x.grid=xgrid,n.grid.additional=200, integral.approx='midd',it=15,kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=FALSE)
  
  return(list(res_reg=res_reg,test1_A=test1_A,test1.mult=test1.mult_A,test1_AB=test1_AB ))
},
mc.cores=1)

######### generating simplots used in paper
######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### 



quartz(width = 11.3, height = 5.7)
par(mfrow=c(2,3))

# 
# plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[1]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), xlab="duration (years)", ylab="alpha_1")
plot(res[[1]][[2]]$x.grid[[1]], res[[1]][[2]]$alpha_backfit[[1]],type='l', ylim=c(0,2), xlab="duration (years)", ylab="alpha_0", col=adjustcolor("grey", alpha = 0.3),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 1:200){
  lines(res[[s]][[2]]$x.grid[[1]],  res[[s]][[2]]$alpha_backfit[[1]], col=adjustcolor("grey", alpha = 0.3))
}
#lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[1]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), col="green")
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[1]],  res[[s]][[1]]$f_backfit[[1]], col=adjustcolor("yellow", alpha = 0.3))
}
#for(s in 1:200){
#  lines(res[[s]][[4]]$x.grid[[1]],  res[[s]][[4]]$alpha_backfit[[1]], col=adjustcolor("blue", alpha = 0.3))
#}



#plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[2]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[2]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res[[1]][[2]]$x.grid[[2]],res[[1]][[2]]$alpha_backfit[[2]],type='l', ylim=c(-0.4,1.3),col=adjustcolor("grey", alpha = 0.3),,ylab="alpha_1", xlab="age (years)", main="vf=0",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res[[s]][[2]]$x.grid[[2]],  res[[s]][[2]]$alpha_backfit[[2]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[2]],  res[[s]][[1]]$f_backfit[[2]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[2]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[2]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res[[s]][[4]]$x.grid[[2]],  res[[s]][[4]]$alpha_backfit[[2]], col=adjustcolor("blue", alpha = 0.3))
# }


#plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[3]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[3]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res[[1]][[2]]$x.grid[[3]], res[[1]][[2]]$alpha_backfit[[3]],type='l',col=adjustcolor("grey", alpha = 0.3),  ylim=c(-0.4,1.3), ylab="alpha_2", xlab="wmi",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res[[s]][[2]]$x.grid[[3]],  res[[s]][[2]]$alpha_backfit[[3]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[3]],  res[[s]][[1]]$f_backfit[[3]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[3]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[3]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res[[s]][[4]]$x.grid[[3]],  res[[s]][[4]]$alpha_backfit[[3]], col=adjustcolor("blue", alpha = 0.3))
# }
mtext("Hazard function estimates for the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) #





# 
# plot(smooth.spline(unlist(lapply(res1, function(l) l$test0_A$x.grid[[1]])),unlist(lapply(res1, function(l) l$test0_A$alpha_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), xlab="duration (years)", ylab="alpha_1")
plot(res1[[1]][[2]]$x.grid[[1]], res1[[1]][[2]]$alpha_backfit[[1]],type='l', ylim=c(0,15), xlab="duration (years)", ylab="alpha_0", col=adjustcolor("grey", alpha = 0.3),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 1:200){
  lines(res1[[s]][[2]]$x.grid[[1]],  res1[[s]][[2]]$alpha_backfit[[1]], col=adjustcolor("grey", alpha = 0.3))
}
#lines(smooth.spline(unlist(lapply(res1, function(l) l$res1_reg$x.grid[[1]])),unlist(lapply(res1, function(l) l$res1_reg$f_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), col="green")
for(s in 1:200){
  lines(res1[[s]][[1]]$x.grid[[1]],  res1[[s]][[1]]$f_backfit[[1]], col=adjustcolor("yellow", alpha = 0.3))
}
#for(s in 1:200){
#  lines(res1[[s]][[4]]$x.grid[[1]],  res1[[s]][[4]]$alpha_backfit[[1]], col=adjustcolor("blue", alpha = 0.3))
#}



#plot(smooth.spline(unlist(lapply(res1, function(l) l$test0_A$x.grid[[2]])),unlist(lapply(res1, function(l) l$test0_A$alpha_backfit[[2]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res1[[1]][[2]]$x.grid[[2]],res1[[1]][[2]]$alpha_backfit[[2]],type='l', ylim=c(-4,5.3),col=adjustcolor("grey", alpha = 0.3),,ylab="alpha_1", xlab="age (years)", main="vf=1",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res1[[s]][[2]]$x.grid[[2]],  res1[[s]][[2]]$alpha_backfit[[2]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res1[[s]][[1]]$x.grid[[2]],  res1[[s]][[1]]$f_backfit[[2]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res1, function(l) l$res1_reg$x.grid[[2]])),unlist(lapply(res1, function(l) l$res1_reg$f_backfit[[2]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res1[[s]][[4]]$x.grid[[2]],  res1[[s]][[4]]$alpha_backfit[[2]], col=adjustcolor("blue", alpha = 0.3))
# }


#plot(smooth.spline(unlist(lapply(res1, function(l) l$test0_A$x.grid[[3]])),unlist(lapply(res1, function(l) l$test0_A$alpha_backfit[[3]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res1[[1]][[2]]$x.grid[[3]], res1[[1]][[2]]$alpha_backfit[[3]],type='l',col=adjustcolor("grey", alpha = 0.3),  ylim=c(-5.4,7.3), ylab="alpha_2", xlab="wmi",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res1[[s]][[2]]$x.grid[[3]],  res1[[s]][[2]]$alpha_backfit[[3]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res1[[s]][[1]]$x.grid[[3]],  res1[[s]][[1]]$f_backfit[[3]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res1, function(l) l$res1_reg$x.grid[[3]])),unlist(lapply(res1, function(l) l$res1_reg$f_backfit[[3]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res1[[s]][[4]]$x.grid[[3]],  res1[[s]][[4]]$alpha_backfit[[3]], col=adjustcolor("blue", alpha = 0.3))
# }
mtext("Hazard function estimates for the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) #
quartz.save("sim1.pdf", type = "pdf")
#dev.off()


quartz(width = 11.3, height = 5.7)
par(mfrow=c(2,3))

# 
# plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[1]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), xlab="duration (years)", ylab="alpha_1")
plot(res[[1]][[2]]$x.grid[[1]], res[[1]][[2]]$alpha_backfit[[1]],type='l', ylim=c(0,2.3), xlab="duration (years)", ylab="alpha_0", col=adjustcolor("grey", alpha = 0.3),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 1:200){
  lines(res[[s]][[2]]$x.grid[[1]],  res[[s]][[2]]$alpha_backfit[[1]], col=adjustcolor("grey", alpha = 0.3))
}
#lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[1]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), col="green")
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[1]],  res[[s]][[1]]$f_backfit[[1]], col=adjustcolor("yellow", alpha = 0.3))
}
#for(s in 1:200){
#  lines(res[[s]][[4]]$x.grid[[1]],  res[[s]][[4]]$alpha_backfit[[1]], col=adjustcolor("blue", alpha = 0.3))
#}



#plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[2]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[2]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res[[1]][[2]]$x.grid[[2]],res[[1]][[2]]$alpha_backfit[[2]],type='l', ylim=c(-0.4,0.8),col=adjustcolor("grey", alpha = 0.3),,ylab="alpha_1", xlab="age (years)", main="vf=0",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res[[s]][[2]]$x.grid[[2]],  res[[s]][[2]]$alpha_backfit[[2]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[2]],  res[[s]][[1]]$f_backfit[[2]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[2]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[2]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res[[s]][[4]]$x.grid[[2]],  res[[s]][[4]]$alpha_backfit[[2]], col=adjustcolor("blue", alpha = 0.3))
# }


#plot(smooth.spline(unlist(lapply(res, function(l) l$test0_A$x.grid[[3]])),unlist(lapply(res, function(l) l$test0_A$alpha_backfit[[3]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res[[1]][[2]]$x.grid[[3]], res[[1]][[2]]$alpha_backfit[[3]],type='l',col=adjustcolor("grey", alpha = 0.3),  ylim=c(-1,1.7), ylab="alpha_2", xlab="wmi",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res[[s]][[2]]$x.grid[[3]],  res[[s]][[2]]$alpha_backfit[[3]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res[[s]][[1]]$x.grid[[3]],  res[[s]][[1]]$f_backfit[[3]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res, function(l) l$res_reg$x.grid[[3]])),unlist(lapply(res, function(l) l$res_reg$f_backfit[[3]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res[[s]][[4]]$x.grid[[3]],  res[[s]][[4]]$alpha_backfit[[3]], col=adjustcolor("blue", alpha = 0.3))
# }
mtext("Hazard function estimates for the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) #





# 
# plot(smooth.spline(unlist(lapply(res3, function(l) l$test0_A$x.grid[[1]])),unlist(lapply(res3, function(l) l$test0_A$alpha_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), xlab="duration (years)", ylab="alpha_1")
plot(res3[[1]][[2]]$x.grid[[1]], res3[[1]][[2]]$alpha_backfit[[1]],type='l', ylim=c(-0.05,0.6), xlab="duration (years)", ylab="alpha_0", col=adjustcolor("grey", alpha = 0.3),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 1:200){
  lines(res3[[s]][[2]]$x.grid[[1]],  res3[[s]][[2]]$alpha_backfit[[1]], col=adjustcolor("grey", alpha = 0.3))
}
#lines(smooth.spline(unlist(lapply(res3, function(l) l$res3_reg$x.grid[[1]])),unlist(lapply(res3, function(l) l$res3_reg$f_backfit[[1]])), tol=0.01),type='l', ylim=c(0,2), col="green")
for(s in 1:200){
  lines(res3[[s]][[1]]$x.grid[[1]],  res3[[s]][[1]]$f_backfit[[1]], col=adjustcolor("yellow", alpha = 0.3))
}
#for(s in 1:200){
#  lines(res3[[s]][[4]]$x.grid[[1]],  res3[[s]][[4]]$alpha_backfit[[1]], col=adjustcolor("blue", alpha = 0.3))
#}



#plot(smooth.spline(unlist(lapply(res3, function(l) l$test0_A$x.grid[[2]])),unlist(lapply(res3, function(l) l$test0_A$alpha_backfit[[2]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res3[[1]][[2]]$x.grid[[2]],res3[[1]][[2]]$alpha_backfit[[2]],type='l', ylim=c(-0.2,.35),col=adjustcolor("grey", alpha = 0.3),,ylab="alpha_1", xlab="age (years)", main="vf=1",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res3[[s]][[2]]$x.grid[[2]],  res3[[s]][[2]]$alpha_backfit[[2]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res3[[s]][[1]]$x.grid[[2]],  res3[[s]][[1]]$f_backfit[[2]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res3, function(l) l$res3_reg$x.grid[[2]])),unlist(lapply(res3, function(l) l$res3_reg$f_backfit[[2]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res3[[s]][[4]]$x.grid[[2]],  res3[[s]][[4]]$alpha_backfit[[2]], col=adjustcolor("blue", alpha = 0.3))
# }


#plot(smooth.spline(unlist(lapply(res3, function(l) l$test0_A$x.grid[[3]])),unlist(lapply(res3, function(l) l$test0_A$alpha_backfit[[3]])), tol=0.01),type='l', ylim=c(-0.4,1.3))
plot(res3[[1]][[2]]$x.grid[[3]], res3[[1]][[2]]$alpha_backfit[[3]],type='l',col=adjustcolor("grey", alpha = 0.3),  ylim=c(-0.2,0.5), ylab="alpha_2", xlab="wmi",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
for(s in 2:200){
  lines(res3[[s]][[2]]$x.grid[[3]],  res3[[s]][[2]]$alpha_backfit[[3]], col=adjustcolor("grey", alpha = 0.3))
}
for(s in 1:200){
  lines(res3[[s]][[1]]$x.grid[[3]],  res3[[s]][[1]]$f_backfit[[3]], col=adjustcolor("yellow", alpha = 0.3))
}
# lines(smooth.spline(unlist(lapply(res3, function(l) l$res3_reg$x.grid[[3]])),unlist(lapply(res3, function(l) l$res3_reg$f_backfit[[3]])), tol=0.01),type='l', ylim=c(0,2), col="green")
# for(s in 1:200){
#   lines(res3[[s]][[4]]$x.grid[[3]],  res3[[s]][[4]]$alpha_backfit[[3]], col=adjustcolor("blue", alpha = 0.3))
# }
mtext("Hazard function conditional on surviving the first three months", side = 3, line = -1.3, outer = TRUE, col = 1) #
quartz.save("sim2.pdf", type = "pdf")
#dev.off()