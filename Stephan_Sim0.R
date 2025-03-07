library(parallel)
library(MASS)
library(survival)
#library(VGAM)
library(eha)
library(timereg)
# load some functions 
source("mhdata_Additive.R")
source("SBF_MH_LL_Additive.R")

plot.pdf <- TRUE
transp.level<-75   # transperency level of curves in figure 100=full transperancy. 

n.sim <- 70
n.cores <- 70  ## runs over different simulations (n.sim)
n.cores.main <-1 ## runs over different jobs (parameter settings)


################# n.range and d.range and bandwidths need to be manually changed for every desired setting



n.range<-c(5000)
d.range<-c(30)

##### optimal parameters
##### d=4 n=500
# b1.range<-c(0.35)
# bt.range<-c(0.2)
# 
# b1LC.range<-c(0.15) #(also showcase b=0.35)
# btLC.range<-c(0.2)
# 
# b1BF.range<-c(0.6)
# btBF.range<-c(0.2)
# 
# b1LCBF.range<-c(0.5)
# btLCBF.range<-c(0.2)

##### d=4 n=5000
# b1.range<-c(0.2) #or with worse mse (0.25)
# bt.range<-c(0.1)

#b1LC.range<-c(0.08)
#btLC.range<-c(0.1)

#b1LLBF.range<-c(0.25)
#btLLBF.range<-c(0.1)

#b1LCBF.range<-c(0.18)
#btLCBF.range<-c(0.1)

##### d=10 n=500
# b1.range<-c(0.42+)
# bt.range<-c(0.2)

#b1LC.range<-c(0.22+)
#btLC.range<-c(0.2)


##### d=10 n=5000
# b1.range<-c(0.25)
# bt.range<-c(0.2)

#b1LC.range<-c(0.10)
#btLC.range<-c(0.2)

#b1LCBF.range<-c(0.17)
#btLCBF.range<-c(0.2)

# b1BF.range<-c(0.25)
# btBF.range<-c(0.2)

##### d=30 n=500
# b1.range<-c(0.5)
# bt.range<-c(0.2)

#b1LC.range<-c(0.10)
#btLC.range<-c(0.2)

##### d=30 n=5000





b1.range<-c(0.3)
bt.range<-c(0.2)

b1LC.range<-c(0.18)
btLC.range<-c(0.2)

b1BF.range<-c(0.4)
btBF.range<-c(0.2)

b1LCBF.range<-c(0.25)
btLCBF.range<-c(0.2)

sparse<-c(FALSE)
rho.range<-c(0.5)
shift <- c(0)
trunc <- c(0)
n.grid <- c(200) 


job<-expand.grid(list(n=n.range,d=d.range,b1=b1.range, bt=bt.range,b1LC=b1LC.range, btLC=btLC.range,b1BF=b1BF.range,btBF=btBF.range,b1LCBF=b1LCBF.range, btLCBF=btLCBF.range, rho=rho.range,shift=shift,trunc=trunc,n.grid=n.grid,sparse=sparse ))


  
my_simulation <- function(i,job){

myfun<-function(ss,n,d,rho,b1,bt,b1LC,btLC,b1BF,btBF,b1LCBF,btLCBF,shift,trunc,n.grid,sparse ){
  
#seed1<-156
seed1<-sample(1:999999,1) 
set <- list(n=n,d=d,rho=rho,model=2,violate.cox=TRUE,seed=seed1,shift=shift,trunc=trunc,sparse)
train.data <- do.call("mhdata",set)


phi <-train.data$phi
train.data <- train.data$data

pred <- paste0("V", 1:(set$d-1))
frmla<-reformulate(pred,"Surv(time, status)")
# plot(prodlim(Hist(time,status)~1,data=large.data),xlim=c(0,100))

pred2 <- paste0("pspline(V", 1:(set$d-1),")")
frmla2<-reformulate(pred2,"Surv(time, status)")
# ranger
#forest1 <-ranger(frmla, data=train.data,importance = "none",num.trees=1000,mtry=set$d,verbose = TRUE)
# randomForestSRC
#forest2 <-rfsrc(frmla,data=train.data,importance = "none",num.trees=1000,mtry=set$d)
# cox
#cox.fit<-coxph(frmla , data=train.data)
# munir's stuff
it=100



xgrid <-NULL
#xgrid<-lapply(1:d, function(z) seq(-1,1,length=2))
#xgrid[[1]]<-seq(0,1,length=2)


b.LL.grid<-b.LC.grid<-b.LLBF.grid <- b.LCBF.grid <-numeric(set$d)
b.LL.grid[1]<-bt
b.LC.grid[1]<-btLC
b.LLBF.grid[1]<-btBF
b.LCBF.grid[1]<-btLCBF

b.LL.grid[2:set$d]<-rep(b1,set$d-1)
b.LC.grid[2:set$d]<-rep(b1LC,set$d-1)
b.LLBF.grid[2:set$d]<-rep(b1BF,set$d-1)
b.LCBF.grid[2:set$d]<-rep(b1LCBF,set$d-1)



n.size<-n.grid 

if (n.size==0) xgrid<-NULL

alpha_backfit.LL<-SBF.MH.CLL(frmla,train.data,b.LL.grid,weight='sw',it=it,x.grid=xgrid,n.grid.additional=n.size,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=FALSE)
alpha_backfit.LC<-SBF.MH.CLL(frmla,train.data,b.LC.grid,weight='sw',it=it,x.grid=xgrid,n.grid.additional=n.size,integral.approx='midd',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=FALSE)
 alpha_backfit.BF.LL<-SBF.MH.CLL(frmla,train.data,b.LLBF.grid,weight='sw',it=it,x.grid=xgrid,n.grid.additional=n.size,integral.approx='midd',kcorr=TRUE,LC=FALSE,wrong=FALSE,classic.backfit=TRUE)
 alpha_backfit.BF.LC<-SBF.MH.CLL(frmla,train.data,b.LCBF.grid,weight='sw',it=it,x.grid=xgrid,n.grid.additional=n.size,integral.approx='midd',kcorr=TRUE,LC=TRUE,wrong=FALSE,classic.backfit=TRUE)

#<-alpha_backfit.BF.LC<-NULL
# splines <- aalen(frmla,train.data,resample.iid=1)
# termplot(splines , term=2, se=TRUE, col.term=1, col.se=1)#predict(splines,X=rbind(c(1,0,0,0,0),rep(1,5)))                 
                 
return(list(seed1=seed1,alpha_backfit.LL=alpha_backfit.LL,alpha_backfit.LC=alpha_backfit.LC,alpha_backfit.BF.LC=alpha_backfit.BF.LC,alpha_backfit.BF.LL=alpha_backfit.BF.LL,phi=phi))



}

jobi<- as.list(job[i,])

res <- do.call( mclapply, c( list(1:n.sim,myfun), jobi,list(mc.cores=n.cores))   )

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

blueT<-t_col("blue", perc = transp.level, name = "blueT")
redT<-t_col("red", perc = transp.level, name = "redT")
greenT<-t_col("green", perc = transp.level, name = "greenT")
yellowT<-t_col("yellow", perc = transp.level, name = "yellowT")
greyT<-t_col("grey", perc = transp.level, name = "greyT")

### phi

# phi <- vector(length=jobi$d,mode="list")
# # Cox violated
#   for(k in 1:jobi$d){
#     if ((k%%2)==1)  {phi[[k]]<- function(z) 2/sqrt(jobi$d-1)*sin(pi*z)} else {phi[[k]]<-function(z) -2/sqrt(jobi$d-1)*sin(pi*z)
#     }}

#### plot result
phi <- res[[1]]$phi




if (plot.pdf==TRUE) pdf(file = paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid,".pdf", sep=""),  width = 10, height = 7) 

for(k in 2:jobi$d){
  
par(mfrow=c(2,2),oma = c(0,0,0,0), mar=c(2.1, 4.1, 2.1, 1.1))
  
plot(res[[1]]$alpha_backfit.LL$x.grid[[k]],phi[[k-1]](res[[1]]$alpha_backfit.LL$x.grid[[k]]),type="l",lwd=3, ylim=4/sqrt(jobi$d-1)*c(-1.2,1.2), xlab=expression("x"[k]), ylab=expression(alpha[k]), main=c(paste("Local Linear Smooth Backfitting",'b_t=',jobi$bt,'b_k=',jobi$b1,sep=" ")))
  
#points( res[[1]]$alpha_backfit.LL$X[,k], rep(0, length(res[[1]]$alpha_backfit.LL$X[,k])), pch=4)
  
lapply(1:n.sim, function(ss){
  lines( res[[ss]]$alpha_backfit.LL$x.grid[[k]], res[[ss]]$alpha_backfit.LL$alpha_backfit[[k]], col=greyT)})


plot(res[[1]]$alpha_backfit.LL$x.grid[[k]],phi[[k-1]](res[[1]]$alpha_backfit.LL$x.grid[[k]]),type="l",lwd=3, ylim=4/sqrt(jobi$d-1)*c(-1.2,1.2), xlab=expression("x"[k]), ylab=expression(alpha[k]), main=c(paste("Local Constant Smooth Backfitting",'b_t=',jobi$btLC,'b_k=',jobi$b1LC,sep=" ")))

lapply(1:n.sim, function(ss){
  lines( res[[ss]]$alpha_backfit.LC$x.grid[[k]], res[[ss]]$alpha_backfit.LC$alpha_backfit[[k]], col=greyT)})


 plot(res[[1]]$alpha_backfit.LL$x.grid[[k]],phi[[k-1]](res[[1]]$alpha_backfit.LL$x.grid[[k]]),type="l",lwd=3, ylim=4/sqrt(jobi$d-1)*c(-1.2,1.2), xlab=expression("x"[k]), ylab=expression(alpha[k]), main=c(paste("Local Linear Classical Backfitting",'b_t=',jobi$btBF,'b_k=',jobi$b1BF,sep=" ")))
 
 lapply(1:n.sim, function(ss){
   lines( res[[ss]]$alpha_backfit.BF.LL$x.grid[[k]], res[[ss]]$alpha_backfit.BF.LL$alpha_backfit[[k]], col=greyT)
#   
 })

# 
 plot(res[[1]]$alpha_backfit.LL$x.grid[[k]],phi[[k-1]](res[[1]]$alpha_backfit.LL$x.grid[[k]]),type="l",lwd=3, ylim=4/sqrt(jobi$d-1)*c(-1.2,1.2), xlab=expression("x"[k]), ylab=expression(alpha[k]), main=c(paste("Local Constant Classical Backfitting",'b_t=',jobi$btLCBF,'b_k=',jobi$b1LCBF,sep=" ")))
# 
 lapply(1:n.sim, function(ss){
   lines( res[[ss]]$alpha_backfit.BF.LC$x.grid[[k]], res[[ss]]$alpha_backfit.BF.LC$alpha_backfit[[k]], col=greyT)})
# 
# 


# legend("topleft", legend=c("True",'LL SBF',"LC SBF", "LC BF","LL BF"), col = c("black","blue","red","green","yellow"), lty=1, lwd=2)

}
if (plot.pdf==TRUE) dev.off()


for(k in 2:jobi$d){
  
MSE_LL <- mean(sapply(1:n.sim, function(ss) { mean((phi[[k-1]](res[[ss]]$alpha_backfit.LL$x.grid[[k]])-res[[ss]]$alpha_backfit.LL$alpha_backfit[[k]])^2) }))
MSE_LC <- mean(sapply(1:n.sim, function(ss) { mean((phi[[k-1]](res[[ss]]$alpha_backfit.LC$x.grid[[k]])-res[[ss]]$alpha_backfit.LC$alpha_backfit[[k]])^2) }))
 MSE_BF_LC <- mean(sapply(1:n.sim, function(ss) { mean((phi[[k-1]](res[[ss]]$alpha_backfit.LL$x.grid[[k]])-res[[ss]]$alpha_backfit.BF.LC$alpha_backfit[[k]])^2) }))
 MSE_BF_LL <- mean(sapply(1:n.sim, function(ss) { mean((phi[[k-1]](res[[ss]]$alpha_backfit.LL$x.grid[[k]])-res[[ss]]$alpha_backfit.BF.LL$alpha_backfit[[k]])^2) }))


x.grid.all <- as.numeric(unlist(sapply(1:n.sim,  function(ss) {as.numeric(round(res[[ss]]$alpha_backfit.LL$X[,k],15))})))


x.grid.all <- sample(x.grid.all, min(5000, length(x.grid.all)), replace=FALSE)

alpha1.all<-matrix(nrow=length(x.grid.all), ncol=n.sim)
for(ss in 1:n.sim){
for(i in 1:length(x.grid.all)){
    index <- which.min(abs(x.grid.all[i]-res[[ss]]$alpha_backfit.LL$x.grid[[k]]))
    error <-  res[[ss]]$alpha_backfit.LL$x.grid[[k]][index] - x.grid.all[i]
    if (error>0&index>=2) {index2<-index-1} else index2<- index
    if (error<0) {index2<-index+1} 
    error2 <-   res[[ss]]$alpha_backfit.LL$x.grid[[k]][index2] -   x.grid.all[i]
    
    if (error!=0&index2>=1&index2<=length(res[[ss]]$alpha_backfit.LL$x.grid[[k]]))  {
      alpha1.all[i,ss] <- (abs(error)*res[[ss]]$alpha_backfit.LL$alpha_backfit[[k]][index]+abs(error2)*res[[ss]]$alpha_backfit.LL$alpha_backfit[[k]][index2])/(abs(error)+abs(error2))} else{
        alpha1.all[i,ss] <- res[[ss]]$alpha_backfit.LL$alpha_backfit[[k]][index]
      }
}
}
alpha1.all.all <- rowMeans(alpha1.all)

MSE_LL2<- mean((alpha1.all - phi[[k-1]](x.grid.all))^2) 

Bias_LL <- mean( (phi[[k-1]](x.grid.all)- alpha1.all.all)^2)             

Var_LL <-   mean( (alpha1.all- alpha1.all.all)^2)  



# write((phi[[k-1]](x.grid.all)- alpha1.all.all), file= paste("BIAS_LL"," ",paste("n=",jobi$n," ","d=",jobi$d," ","rho=",jobi$rho,sep="")," ","b1=",jobi$b1," ", "bt=",jobi$bt,'shift=',jobi$shift," ","trunc=",jobi$trunc, ".txt", sep=""), append=TRUE)
# 
# 
# write(x.grid.all, file= paste("XGRID"," ",paste("n=",jobi$n," ","d=",jobi$d," ","rho=",jobi$rho,sep="")," ","b1=",jobi$b1," ", "bt=",jobi$bt,'shift=',jobi$shift," ","trunc=",jobi$trunc, ".txt", sep=""), append=TRUE)
# 
# 
# write(( alpha1.all.all), file= paste("LL"," ",paste("n=",jobi$n," ","d=",jobi$d," ","rho=",jobi$rho,sep="")," ","b1=",jobi$b1," ", "bt=",jobi$bt,'shift=',jobi$shift," ","trunc=",jobi$trunc, ".txt", sep=""), append=TRUE)
# 
# 
# write((phi[[k-1]](x.grid.all)), file= paste("PHI"," ",paste("n=",jobi$n," ","d=",jobi$d," ","rho=",jobi$rho,sep="")," ","b1=",jobi$b1," ", "bt=",jobi$bt,'shift=',jobi$shift," ","trunc=",jobi$trunc, ".txt", sep=""), append=TRUE)




alpha1.all<-matrix(nrow=length(x.grid.all), ncol=n.sim)
for(ss in 1:n.sim){
  for(i in 1:length(x.grid.all)){
    index <- which.min(abs(x.grid.all[i]-res[[ss]]$alpha_backfit.LC$x.grid[[k]]))
    error <-  res[[ss]]$alpha_backfit.LC$x.grid[[k]][index] - x.grid.all[i]
    if (error>0&index>=2) {index2<-index-1} else index2<- index
    if (error<0) {index2<-index+1} 
    error2 <-   res[[ss]]$alpha_backfit.LC$x.grid[[k]][index2] -   x.grid.all[i]
    
    if (error!=0&index2>=1&index2<=length(res[[ss]]$alpha_backfit.LC$x.grid[[k]]))  {
      alpha1.all[i,ss] <- (abs(error)*res[[ss]]$alpha_backfit.LC$alpha_backfit[[k]][index]+abs(error2)*res[[ss]]$alpha_backfit.LC$alpha_backfit[[k]][index2])/(abs(error)+abs(error2))} else{
        alpha1.all[i,ss] <- res[[ss]]$alpha_backfit.LC$alpha_backfit[[k]][index]
      }
  }
}
alpha1.all.all <- rowMeans(alpha1.all)

MSE_LC2<- mean((alpha1.all - phi[[k-1]](x.grid.all))^2) 

Bias_LC <- mean( (phi[[k-1]](x.grid.all)- alpha1.all.all)^2)             

Var_LC <-   mean( (alpha1.all- alpha1.all.all)^2)  

alpha1.all<-matrix(nrow=length(x.grid.all), ncol=n.sim)
for(ss in 1:n.sim){
  for(i in 1:length(x.grid.all)){
    index <- which.min(abs(x.grid.all[i]-res[[ss]]$alpha_backfit.BF.LC$x.grid[[k]]))
    error <-  res[[ss]]$alpha_backfit.BF.LC$x.grid[[k]][index] - x.grid.all[i]
    if (error>0&index>=2) {index2<-index-1} else index2<- index
    if (error<0) {index2<-index+1}
    error2 <-   res[[ss]]$alpha_backfit.BF.LC$x.grid[[k]][index2] -   x.grid.all[i]

    if (error!=0&index2>=1&index2<=length(res[[ss]]$alpha_backfit.BF.LC$x.grid[[k]]))  {
      alpha1.all[i,ss] <- (abs(error)*res[[ss]]$alpha_backfit.BF.LC$alpha_backfit[[k]][index]+abs(error2)*res[[ss]]$alpha_backfit.BF.LC$alpha_backfit[[k]][index2])/(abs(error)+abs(error2))} else{
        alpha1.all[i,ss] <- res[[ss]]$alpha_backfit.BF.LC$alpha_backfit[[k]][index]
      }
  }
}
alpha1.all.all <- rowMeans(alpha1.all)



# ,
MSE_BF_LC2<- mean((alpha1.all - phi[[k-1]](x.grid.all))^2) 
# 
 Bias_BF_LC <- mean( (phi[[k-1]](x.grid.all)- alpha1.all.all)^2)             
# 
 Var_BF_LC <-   mean( (alpha1.all- alpha1.all.all)^2)  
# # 
 # write((phi[[k-1]](x.grid.all)- alpha1.all.all), file= paste("BIAS_BF_LC"," ",paste("n=",jobi$n," ","d=",jobi$d," ","rho=",jobi$rho,sep="")," ","b1=",jobi$b1," ", "bt=",jobi$bt,'shift=',jobi$shift," ","trunc=",jobi$trunc, ".txt", sep=""), append=TRUE)
# 
alpha1.all<-matrix(nrow=length(x.grid.all), ncol=n.sim)
 for(ss in 1:n.sim){
   for(i in 1:length(x.grid.all)){
     index <- which.min(abs(x.grid.all[i]-res[[ss]]$alpha_backfit.BF.LL$x.grid[[k]]))
     error <-  res[[ss]]$alpha_backfit.BF.LL$x.grid[[k]][index] - x.grid.all[i]
     if (error>0&index>=2) {index2<-index-1} else index2<- index
     if (error<0) {index2<-index+1} 
     error2 <-   res[[ss]]$alpha_backfit.BF.LL$x.grid[[k]][index2] -   x.grid.all[i]
     
     if (error!=0&index2>=1&index2<=length(res[[ss]]$alpha_backfit.BF.LL$x.grid[[k]]))  {
       alpha1.all[i,ss] <- (abs(error)*res[[ss]]$alpha_backfit.BF.LL$alpha_backfit[[k]][index]+abs(error2)*res[[ss]]$alpha_backfit.BF.LL$alpha_backfit[[k]][index2])/(abs(error)+abs(error2))} else{
         alpha1.all[i,ss] <- res[[ss]]$alpha_backfit.BF.LL$alpha_backfit[[k]][index]
       }
   }
 }
 alpha1.all.all <- rowMeans(alpha1.all)
# 
# 
# 
 MSE_BF_LL2<- mean((alpha1.all - phi[[k-1]](x.grid.all))^2) 
# 
 Bias_BF_LL <- mean( (phi[[k-1]](x.grid.all)- alpha1.all.all)^2)             
# 
 Var_BF_LL <-   mean( (alpha1.all- alpha1.all.all)^2)  
# 

 write(c("Method","MSE","MSE2", "Bias^2", "Var"), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid, ".txt", sep=""),ncolumns=5, append=TRUE)

 write(c(paste("k=",k),"","","",""), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid, ".txt", sep=""),ncolumns=5, append=TRUE)

 write(c("LL",MSE_LL,MSE_LL2,Bias_LL,Var_LL), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid, ".txt", sep=""),ncolumns=5, append=TRUE)
 
 write(c("LC",MSE_LC,MSE_LC2,Bias_LC,Var_LC), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid, ".txt", sep=""),ncolumns=5, append=TRUE)
 
  write(c("BF-LL",MSE_BF_LL,MSE_BF_LL2,Bias_BF_LL,Var_BF_LL), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid,".txt", sep=""),ncolumns=5, append=TRUE)
 # 
  write(c("BF-LC",MSE_BF_LC,MSE_BF_LC2,Bias_BF_LC,Var_BF_LC), file= paste("model=2","Stephan_Sim"," ",paste("n=",jobi$n," ","d=",jobi$d," ", "sparse=",jobi$sparse,"rho=",jobi$rho,sep="")," ","b.LL=",jobi$b1,'b.LC=',jobi$b1LC,'b.LL.BF=',jobi$b1BF,'b.LC.BF=',jobi$b1LCBF," ","bt.LL=",jobi$bt,'bt.LC=',jobi$btLC,'bt.LL.BF=',jobi$btBF,'bt.LC.BF=',jobi$btLCBF,"n.grid=",jobi$n.grid,".txt", sep=""),ncolumns=5, append=TRUE)
 
}
}

do.call( mclapply, c( list(1:nrow(job), my_simulation), list(job),list(mc.cores=n.cores.main))   )







# if(ss==1) 
# pdf(file = "/Users/ndphillips/Desktop/My Plot.pdf",   # The directory you want to save the file in
#     width = 4, # The width of the plot in inches
#     height = 4) # The height of the plot in inches
# 
# 


# plot(res[[1]]$alpha_backfit.LL$x.grid[[2]],1*sin(pi*(res[[1]]$alpha_backfit.LL$x.grid[[2]])))
# lines(alpha_backfit.LL$x.grid[[2]],alpha_backfit.LL$alpha_backfit[[2]], col='blue')
# lines(alpha_backfit.LC$x.grid[[2]],alpha_backfit.LC$alpha_backfit[[2]], col='red')
# lines(alpha_backfit.BF$x.grid[[2]],alpha_backfit.BF$alpha_backfit[[2]], col='green')
# lines(alpha_backfit2$x.grid[[2]],alpha_backfit2$O1[[2]]/alpha_backfit2$D[[2]]-2, col='red')
# 
# 
# plot(alpha_backfit2$x.grid[[3]],1*-sin(pi*(alpha_backfit2$x.grid[[3]])))
# lines(alpha_backfit2$x.grid[[3]],alpha_backfit2$alpha_backfit[[3]], col='green')
# lines(alpha_backfit2$x.grid[[3]],alpha_backfit2$O1[[3]]/alpha_backfit2$D[[3]]-2, col='red')
# 
# #Model2
# ##plot(alpha_backfit2$x.grid[[1]],alpha*exp(beta*alpha_backfit2$x.grid[[1]])+2, col='blue')
# #Model1
# plot(alpha_backfit2$x.grid[[1]],rep(4, length(alpha_backfit2$x.grid[[1]])), col='blue')
# lines(alpha_backfit2$x.grid[[1]],alpha_backfit2$alpha_backfit[[1]], col='red')
# plot(alpha_backfit2$x.grid[[1]],alpha_backfit2$O1[[1]]/alpha_backfit2$D[[1]], col='green')
# boxplot(train.data[,1])
