###########################################################################
# {{{
getAalen<-function(data,t0,mt,a0,ma,m,
                   time="time",time0="time0",age="age",age0="age0",status="status")
{# {{{
  #xt <- seq(t0,mt,length=m)
#xx<-seq(a0,ma,length=m)
  Atx<-function(x)
  {if (length(which(data$time<=x&data$status==1))>0)
  {
    At<-sum(sapply(data$time[which(data$time<=x&data$status==1)], 
                   function(y){ if (sum(data$time>=y&y>=data$time0)==0) return(0)
                     else return(1/sum(data$time>=y&y>=data$time0))}
    ) 
    )
  }
    else At<-0
    return(as.numeric(At)) }
  
  Aax<-function(x)
  {if (length(which(data$age<=x&data$status==1))>0)
  {
    Aa<-sum(sapply(data$age[which(data$age<=x&data$status==1)], function(y){
      if (sum(((data$age)>=y)&(data$age0<=y))==0) return(0) 
      else return(1/(sum((data$age>=y)&(data$age0<=y))))}))
  }
    else Aa<-0
    return(as.numeric(Aa)) }
  
  At<-Aa<-Nt<-numeric(m)
  #Nt<-sapply(xt, function(x){sum(data$time<=x&data$status==1)})
  #Yt<-sapply(xt, function(x){sum(data$time>x)})
  At<-sapply(xt, Atx)
  Aa<-sapply(xx, Aax)
  return(list(Aa=Aa,At=At))
}# }}}

getEaEt2.old<-function(data,t0,mt,a0,ma,
                       time="time",time0="time0",age="age",age0="age0",
                       status="status")# {{{
{
  xx <- seq(a0,ma,length=m)
  xt <- seq(t0,mt,length=m)
  #Yt<-sapply(xt, function(x){sum(data$time>x)})
  #Ya<-sapply(xx, function(x){sum((data$age)>x&(data$age0<=x))})
  E1<-tta<-matrix(nrow=length(xt),ncol=length(xx))
  E2<-matrix(nrow=length(xx),ncol=length(xt))
  
  for(i in 1:length(xx)){
    for(j in 1:length(xt)){
      Yt<-data$time>xt[j]
      Ya<-(data$age>xx[i])#&(data$age0<=xx[i])
      It<-0-data$age0<xt[j]&xx[i]-data$age0>=xt[j]
      Ia<-data$age0<xx[i]&data$age0+xt[j]>=xx[i]
      YIt<-Yt*It
      YIa<-Ya*Ia
      
      # Y.ai<-matrix(nrow=n,ncol=n)
      # for (l in 1:n)
      # {
      # Y.ai[l,]<-(data$time>xx[i]-data$age0[l])*((xx[i]-data$age0[l])>=data$time0)
      # }
      # Y.aii<-rowSums(Y.ai)
      Y.aii<-sapply(1:n,function (l) sum((data$time>xx[i]-data$age0[l])*((xx[i]-data$age0[l])>=0)))
      Y.tii<-sapply(1:n,function (l) sum((data$age>xt[j]+data$age0[l])*((xt[j]+data$age0[l])>=data$age0)))
      #   aaa<-xx[i]-0
      #   yaai <- Cpred(cbind(xt,1/cumat),rep(aaa,100))[,2]
      #   ww <- (yaai>0)
      # tta[j,i] <- sum(YIa[ww]/yaai[ww])
      
      
      # if (sum(Yt)==0) E2[i,j]<-0 else E2[i,j]<-sum(YIt)/sum(Yt)
      if (sum(Yt)==0) E2[i,j]<-0 else E2[i,j]<-sum(max(0,YIt/Y.tii,na.rm=TRUE))
      if (sum(Ya)==0) E1[j,i]<-0 else E1[j,i]<-sum(max(0,YIa/Y.aii,na.rm=TRUE)) #E1[j,i]<-sum(YIa)/sum(Ya) 
    }}
  # E1<-(xt/(diff(xx)[1]*rowSums(E1)))*E1
  #E2<-(xx/(diff(xt)[1]*rowSums(E2)))*E2
  #E1[is.nan(E1)]<-0
  #E2[is.nan(E2)]<-0
  return(list(E1b=E1,E1=tta,E2=E2))
}# }}}

getEaEt2.new<-function(data,t0,mt,a0,ma,xt,xx,
                       time="time",time0="time0",age="age",age0="age0",
                       status="status")
{# {{{
  #xx <- seq(a0,ma,length=m)
  #xt <- seq(t0,mt,length=m)
  
  Yit<-lapply(xt, function (j) data$time>j&j>=data$time0)
  Yia<-lapply(xx, function (i) {data$age>i&(data$age0<=i)})
  #Y.ia<-lapply(xx, function (i) {lapply(i-data$age0, function (j) {sum(data$time>j&j+ grid.length.t/2>=data$time0)} )})
  #Y.it<-lapply(xt, function (j) {lapply(j+data$age0, function(i){ sum((data$age>i)&(data$age0<=i+ grid.length.a/2))  })} )  
  Y.ia<-lapply(xx, function (i) {lapply(i-data$age0, 
				function (j) {sum(data$time>j&j>=data$time0)} )})
  Y.it<-lapply(xt, function (j) {lapply(j+data$age0, 
			function(i){ sum((data$age>i)&(data$age0<=i) ) })} )  
  
###  print(summary(unlist(Y.ia)))
###  print(summary(unlist(Y.it)))

  E1<-matrix(nrow=length(xt),ncol=length(xx))
  E2<-matrix(nrow=length(xx),ncol=length(xt))
  
  for(i in 1:length(xx)){
    #print(i)
    for(j in 1:length(xt)){
      It<-((a0-data$age0)<xt[j])&(xx[i]-data$age0>=xt[j])
      Ia<-t0+data$age0<xx[i]&data$age0+xt[j]>=xx[i]
      
      E2[i,j]<-sum(Yit[[j]]/as.numeric(Y.it[[j]])*It,na.rm=TRUE)
      E1[j,i]<-sum(Yia[[i]]/as.numeric(Y.ia[[i]])*Ia,na.rm=TRUE)
    }}
   #  E1<-(xt/(diff(xx)[1]*rowSums(E1)))*E1
   # E2<-(xx/(diff(xt)[1]*rowSums(E2)))*E2
   #  E1[is.nan(E1)]<-0
   # E2[is.nan(E2)]<-0
  return(list(E1=E1,E2=E2))
}# }}}


backfitEaEtMM <- function(E1,E2,cumaa,cumat,t0,mt,a0,ma,it=500,after=0,
                        time="time",time0="time0",age="age",age0="age0")
{# {{{
  At<-At2<-cumat
  Ba<-Ba2<-cumaa
  for(i in 1:it) {
    AAt<-At
    BBa<-Ba
###    At<-cumat - c(E1[,-1]%*%diff(c(Ba)))
###    Ba<-cumaa- c(E2[,-1]%*%diff(c(At))) + (tail(At,1)/tail(xt,1))*xx
    At<-cumat - c(E1%*%diff(c(0,Ba))) ## - (tail(At,1)/tail(xt,1))*xt
    Ba<-cumaa- c(E2%*%diff(c(0,At)))  + (tail(At,1)/tail(xt,1))*xx
###    At2<-cumat - c(E1[,-m]%*%diff(c(Ba2)))
###    Ba2<-cumaa- c(E2[,-m]%*%diff(c(At2)))+ (tail(At2,1)/tail(xt,1))*xx
    AAAt<-At-AAt
    BBBa<-Ba-BBa
    if ( (max(abs(AAAt)) <0.00001) & (max(abs(BBBa))<0.00001))  i <- it
  }
  
  if (after==1) {
	  gamma <- tail(At,1)/tail(xt,1)
	  At <- At- xt*gamma
	  Ba<-Ba+(xx-a0)*gamma
  }
  
###  gamma <-tail(At2,1)/tail(xt,1)
###  At2 <- At2- xt*gamma
###  Ba2<-Ba2+xx*gamma
###  At<-(At+At2)/2
###  Ba<-(Ba+Ba2)/2
  
###  AAAt <- BBBa <- 0

  return(list(At=At,Ba=Ba,AAAt=AAAt,BBBa=BBBa))
}# }}}

backfitEaEt <- function(E1,E2,cumaa,cumat,xt,xx,it=10)
{# {{{
  At<-At2<-cumat
  Ba<-Ba2<-cumaa
  m1 <- length(cumaa)
  m2 <- length(cumat)
  for(i in 1:it){
    AAt<-At
    BBa<-Ba
    
    ###right side version
  # At<-cumat - c(E1[,-1]%*%diff(c(Ba)))
  # Ba<-cumaa- c(E2[,-1]%*%diff(c(At))) + (tail(At,1)/tail(xt,1))*xx
    
   ###left-side version
   At<-cumat - c(E1[,-m1]%*%diff(c(Ba)))
   Ba<-cumaa- c(E2[,-m2]%*%diff(c(At)))+ (tail(At,1)/tail(xt,1))*xx
    
   AAAt<-At-AAt
   BBBa<-Ba-BBa
    if ( (max(abs(AAAt)) <0.00001) & (max(abs(BBBa))<0.00001))  i <- it
  }
  
  return(list(At=At,Ba=Ba,AAAt=AAAt,BBBa=BBBa))
}# }}}

getE<-function(E1,E2,m)
{# {{{
  
  dM<-matrix(ncol=m,nrow=m,0)
  dM[col(dM)==row(dM)]<--1
  dM[col(dM)==row(dM)+1]<-1
  dM[m,]<-0
  
###
#If one wishes to approximate the integral with the right values 
  # dM<-matrix(ncol=m,nrow=m,0)
   # dM[col(dM)==row(dM)-1]<--1
   # dM[col(dM)==row(dM)]<-1
   # dM[1,]<-0
### 
  E2<-E2%*%dM
  E1<-E1%*%dM
  
  
  E<-matrix(ncol=2*m,nrow=2*m)
  zeromatrix<-matrix(ncol=m,nrow=m,0)
  E[1:m,]<- -cbind(zeromatrix,E1)
  E[(m+1):(2*m),]<- -cbind(E2,zeromatrix)
  
  lambda<-numeric(2)
  lambda[1]<-max(abs(eigen(E)$value))  #max eigenvalue of  E before correction
  
  
  ##### Add Identification constraint for E2 and then into E 
  
  correction<-zeromatrix
  correction[,m]<--rep(1/xt[m],m)*xx
  E2<-E2+correction
  E[(m+1):(2*m),]<- -cbind(E2,zeromatrix)
  lambda[2]<-max(abs(eigen(E)$value))
  

  ##############
  print(c('Eigenvalues before and after correction:',lambda))
  return(list(E=E,lambda=lambda))
} # }}}

getAB<-function(E,At,Aa,summands=10)
{# {{{
  
  sol=0
  for (i in 0:summands)
  {
    sol<-sol+E%^%i%*%c(At,Aa)
  }
  return(sol)
}# }}}


getAalen.boot<-function(data,t0,mt,a0,ma,m,n.sim=1000,
                   time="time",time0="time0",age="age",age0="age0",status="status")
{# {{{
  #xt <- seq(t0,mt,length=m)
  #xx<-seq(a0,ma,length=m)
  At<-list()
  Aa<-list()
  for(i in 1:n.sim){
    print(i)
  g<-rnorm(nrow(data))
  
  Atx<-function(x)
  {if (length(which(data$time<=x&data$status==1))>0)
  {
    At<-sum(sapply(1:nrow(data), 
                   function(y){ if (sum(data$time>data$time[y]&data$time[y]>=data$time0)==0 ) return(0)
                     else return(g[y]*(data$time[y]<=x&data$status[y])/sum(data$time>data$time[y]&data$time[y]>=data$time0))})
             ) 
  }
    else At<-0
    return(as.numeric(At)) }
  
  Aax<-function(x)
  {if (length(which(data$age<=x&data$status==1))>0)
  {
    Aa<-sum(sapply(1:nrow(data),function(y){
      if (sum(((data$age)>data$age[y])*(data$age0<=data$age[y]))==0) return(0) 
      else return(g[y]*(data$age[y]<=x&data$status[y])/(sum(((data$age)>data$age[y])*(data$age0<=data$age[y]))))}))
  }
    else Aa<-0
    return(as.numeric(Aa)) }
    #Nt<-sapply(xt, function(x){sum(data$time<=x&data$status==1)})
  #Yt<-sapply(xt, function(x){sum(data$time>x)})
  At[[i]]<-sapply(xt, Atx)
  Aa[[i]]<-sapply(xx, Aax)
  }
  return(list(Aa=Aa,At=At))
}# }}}
# }}}
#######################################################################

timescale.box <- function(t0,mt,a0,ma,data,
			  time="time",time.start="time0",
			  age="age",age.start="age0",
		          age.null="age.null",time.null="time.null",
			  status="status",gplot=0)
{# {{{

# if (gplot==1) {
#    TT <- fast.reshape(data,list(c(age.start,age),c(time.start,time)))
#    par(mfrow=c(1,2))
#    ff <- toformula(age.start,time.start)
#    spaghetti(ff,id="id",data=TT)#,col=(TT$dead)+1,lwd=0.1)
#    abline(h=c(a0,ma),col=2)
#    abline(v=c(t0,mt),col=2)
# }

###data=TRACE; t0 <- 0; mt <- 5; a0 <- 55 ; ma <- 72; 
ccnames <- names(data)

data <- event.split(data,time=time,status=status,cuts=t0,name.start=time.start)
keep <- data[,time.start]>=t0
data <- data[keep,] 
## update start age
data[,age.start] <-  data[,age.null]+(data[,time.start]-data[,time.null])
keep <- data[,age.start]>0
data <- data[keep,] 
data <- data[,ccnames]
keep <- data[,age.start]<data[,age]
data <- data[keep,] 

data <- event.split(data,time=time,status=status,cuts=mt,name.start=time.start)
keep <- (data[,time] <= mt)
data <-  data[keep,]
###data
## update age
data[,age] <-  data[,age.null]+(data[,time]-data[,time.null])
data <- data[,ccnames]
###data

###print(with(data,lifetable(Surv(time0,time,status!=0)~+1)))
###print(with(data,lifetable(Surv(age0,age,status!=0)~+1)))

data <- event.split(data,time=age,status=status,cuts=a0,name.start=age.start)
keep <- data[,age.start]>=a0
data <- data[keep,] 
###data
## update start time
data[,time.start] <-  data[,time.null]+(data[,age.start]-data[,age.null])
rrr <- data[,time.start]<0
data[rrr,time.start] <- 0
data <- data[,ccnames]
keep <- data[,time.start]<data[,time]
data <- data[keep,] 

###print(with(data,lifetable(Surv(time0,time,status!=0)~+1)))
###print(with(data,lifetable(Surv(age0,age,status!=0)~+1)))

data <- event.split(data,time=age,status=status,cuts=ma,name.start=age.start)
keep <- (data[,age] <= ma)
data <-  data[keep,]
## update time 
data[,time] <-  data[,time.null]+(data[,age]-data[,age.null])
data <- data[,ccnames]

###print(with(data,lifetable(Surv(time0,time,status!=0)~+1)))
###print(with(data,lifetable(Surv(age0,age,status!=0)~+1)))

# if (gplot==1) {
# 	TT <- fast.reshape(data,list(c(age.start,age),c(time.start,time)))
#         spaghetti(ff,id="id",data=TT)
# 	abline(h=c(a0,ma),col=2)
# 	abline(v=c(t0,mt),col=2)
# }

return(data)
}# }}}



