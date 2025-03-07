


mhcovariate <- function(n,d=5,rho=0,seed){
    if (!missing(seed)) set.seed(seed)
    stddev<-rep(1,d-1)
    corMat<-matrix(rho, nrow=d-1,ncol=d-1)
    corMat[col(corMat)==row(corMat)]<-1
    covMat<-stddev %*% t(stddev) * corMat
    Z<-mvrnorm(n=n,  mu=rep(0,d-1), Sigma = covMat, empirical = FALSE)
    Z<-2.5*atan(Z)/pi
    #Z<-matrix(runif(n*(d-1),-1,1), nrow=n, ncol=d-1)
    Z
   
}



# data generating mechanism
mhrate <- function(Z,model=1,violate.cox=TRUE,sparse){
    Z<-as.matrix(Z)
    d <- NCOL(Z)
    phi <- vector(length=d,mode="list")
    if (violate.cox==TRUE){
        # Cox violated
        for(k in 1:d){
            if ((k%%2)==1)  {phi[[k]]<- function(z) (4/sqrt(d))*sin(pi*z)} else {phi[[k]]<-function(z) -(4/sqrt(d))*sin(pi*z)}
         if(sparse==TRUE){ if (k>=4) {phi[[k]]<-function(z) 0*z} }## sparse setting
        }
    }else{
        # Cox satisfied
        for(k in 1:d){
            if ((k%%2)==1)  {phi[[k]]<- function(z) -2*z} else {phi[[k]]<-function(z) 2*z}
        }
    }
    
    
    
    
    top <- rep(0,NROW(Z))
    for (k in 1:d)
    {
        top<-  top+phi[[k]](Z[,k])
    }
    
    alpha_prod <- rep(1,NROW(Z))
     for (k in 1:d)
     {
       alpha_prod <-  alpha_prod * (phi[[k]](Z[,k])+ (4/sqrt(d))+0.001)
     }
     alpha_prod <- alpha_prod/10
    
    
    phi_i <- list()
    for (k in 1:d)
    {
      phi_i[[k]] <- phi[[k]](Z[,k])
    }
  
    
    if (model==1)  surv.function <- function(t){ pexp(t,rate=top,lower.tail = FALSE)}
    if (model==2)  surv.function <- function(t) {pmakeham(t,scale=1, shape=c(1,1),epsilon=top,lower.tail =  FALSE)}
    if (model==3)  surv.function <- function(t) {pmakeham(t,scale=1, shape=c(1,1),epsilon=top,lower.tail =  FALSE)}     
    if (model==4)  surv.function <- function(t) {pmakeham(t,scale=1, shape=c(1,1),epsilon=top,lower.tail =  FALSE)}     
    
    # true_par<-0
    # for (k in 1:d){
    #     true_par<-  true_par+phi[[k]](Z[[k]])
    # }
    # true_par<-exp(true_par)

  return(list(top=top, surv.function=surv.function, phi=phi, alpha_prod=alpha_prod, phi_i=phi_i))
}



mhdata <- function(n=200,d=5,rho=0,model=1,violate.cox=TRUE,seed,shift,trunc,sparse){
    if (!missing(seed)) set.seed(seed)
    Z <- mhcovariate(n=10*n,d=d,rho=rho)
    # regression coefficients
    temp <- mhrate(Z,model,violate.cox=violate.cox,sparse)
    
    top <- temp$top
    phi <- temp$phi
    phi_i <- temp$phi_i
    alpha_prod <- temp$alpha_prod
    
    top <- top + shift
    Z<-Z[top>trunc,]
    top<-top[top>trunc]
    phi_i <- lapply(1:length(phi_i), function(j) phi_i[[j]][top>trunc]    )
    alpha_prod <- alpha_prod[top>trunc]
    
    phi_i <- phi_i[1:n]
    alpha_prod <- alpha_prod[1:n]
    top<-top[1:n]
    Z<-Z[1:n,]
    
    if (model==1){
    Time<-rexp(length(top),top)
    C<-rexp(length(top),top/1.75)
 
    }

  if (model==2)
  #survivaldistr=='makeham' 
  { 
     beta<-0.01
     alpha<- 1
    #Time<-rmakeham(n,beta, alpha, epsilon=top)
   # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
    Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha,top[i]),1/beta))
    C<-sapply(1:length(top), function(i) rmakeham(1,c(  alpha/1.75,top[i]),1/ beta))
    haz <- sapply(1:n, function(i) {alpha*exp(beta*Time[i])+top[i]})
    
  }
    
    
    if (model==3)
      #survivaldistr=='makeham' 
    { 
      beta<-0.01
      alpha<- 1 
      #Time<-rmakeham(n,beta, alpha, epsilon=top)
      # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
      Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha+pmin(3,1/alpha_prod)[i],top[i]),1/beta))
      C<-sapply(1:length(top), function(i) rmakeham(1,c(  (alpha+pmin(3,1/alpha_prod)[i])/1.75,top[i]),1/ beta))
      haz <- sapply(1:length(top), function(i) {(alpha+pmin(3,1/alpha_prod)[i])*exp(beta*Time[i])+top[i]})
      }
    
    if (model==4)
      #survivaldistr=='makeham' 
    { 
      beta<-0.01
      alpha<- 1 
      #Time<-rmakeham(n,beta, alpha, epsilon=top)
      # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
      Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha+pmin(3,1/top)[i],top[i]),1/beta))
      C<-sapply(1:length(top), function(i) rmakeham(1,c(  (alpha+pmin(3,1/top)[i])/1.75,top[i]),1/ beta))
      haz <- sapply(1:length(top), function(i) {(alpha+pmin(3,1/top)[i])*exp(beta*Time[i])+top[i]})
    }
    
  TT<-Time*(Time<=C)+C*(Time>C)
  status<-(Time<=C)*1  ## censoring indicator
  data<-data.frame(time=TT,status=status,as.data.frame(Z))

return(list(data=data,phi=phi, haz=haz))
}


##### function below outputse more than n instances in contrast to mhdata
mhdata2 <- function(n=200,d=5,rho=0,model=1,violate.cox=TRUE,seed,shift,trunc,sparse){
  if (!missing(seed)) set.seed(seed)
  Z <- mhcovariate(n=10*n,d=d,rho=rho)
  # regression coefficients
  temp <- mhrate(Z,model,violate.cox=violate.cox,sparse)
  
  top <- temp$top
  phi <- temp$phi
  phi_i <- temp$phi_i
  alpha_prod <- temp$alpha_prod
  
  top <- top + shift
  Z<-Z[top>trunc,]
  top<-top[top>trunc]
  phi_i <- lapply(1:length(phi_i), function(j) phi_i[[j]][top>trunc]    )
  alpha_prod <- alpha_prod[top>trunc]
  
  
  if (model==1){
    Time<-rexp(length(top),top)
    C<-rexp(length(top),top/1.75)
    
  }
  
  if (model==2)
    #survivaldistr=='makeham' 
  { 
    beta<-0.01
    alpha<- 1
    #Time<-rmakeham(n,beta, alpha, epsilon=top)
    # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
    Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha,top[i]),1/beta))
    C<-sapply(1:length(top), function(i) rmakeham(1,c(  alpha/1.75,top[i]),1/ beta))
    haz <- sapply(1:length(top), function(i) {alpha*exp(beta*Time[i])+top[i]})
    
  }
  
  
  if (model==3)
    #survivaldistr=='makeham' 
  { 
    beta<-0.01
    alpha<- 1 
    #Time<-rmakeham(n,beta, alpha, epsilon=top)
    # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
    Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha+pmin(3,1/alpha_prod)[i],top[i]),1/beta))
    C<-sapply(1:length(top), function(i) rmakeham(1,c(  (alpha+pmin(3,1/alpha_prod)[i])/1.75,top[i]),1/ beta))
    haz <- sapply(1:length(top), function(i) {(alpha+pmin(3,1/alpha_prod)[i])*exp(beta*Time[i])+top[i]})
  }
  
  
  if (model==4)
    #survivaldistr=='makeham' 
  { 
    beta<-0.01
    alpha<- 1 
    #Time<-rmakeham(n,beta, alpha, epsilon=top)
    # C<-rmakeham(n,beta, alpha/1.75, epsilon=top)
    Time<-sapply(1:length(top), function(i) rmakeham(1,c(alpha+pmin(3,1/top)[i],top[i]),1/beta))
    C<-sapply(1:length(top), function(i) rmakeham(1,c(  (alpha+pmin(3,1/top)[i])/1.75,top[i]),1/ beta))
    haz <- sapply(1:length(top), function(i) {(alpha+pmin(3,1/top)[i])*exp(beta*Time[i])+top[i]})
  }
  
  
  TT<-Time*(Time<=C)+C*(Time>C)
  status<-(Time<=C)*1  ## censoring indicator
  data<-data.frame(time=TT,status=status,as.data.frame(Z))
  
  return(list(data=data,phi=phi, haz=haz))
}

# 
# top=0.000000000000001
# alpha=1
# beta=0.08
# x=seq(0,5,length.out=100)
# dx<-c(0,diff(x))
# 
# plot(seq(0,5,length.out=100),dmakeham(seq(0,5,length.out=100), shape = c(0.000000000000001, 1), scale = 1/0.08, log = FALSE))
# 
# lines(seq(0,5,length.out=100),(top*exp(beta*x)+alpha)*exp(-alpha*x- (top/beta)*(exp(x*beta))), col="red")


######################################################################
### mhdata.R ends here
