predict.hazard<-function(result, data){
  data<-as.matrix(data)
  if (ncol(data)==1) data<-t(data)
  x.grid<-result$x.grid
  dx<-result$dx
  
  alpha_sbf<-result$alpha_backfit
  alpha.1_sbf<-result$alpha.1_backfit
  
 # mine <- c(30/365,70,2)
  
  
  hazard <- matrix(nrow=nrow(data), ncol=(length( alpha_sbf)+1))
  hazard.1 <- matrix(nrow=nrow(data), ncol=(length( alpha.1_sbf)+1))
  crps_1 <- crps_2 <- crps_2_pos <- crps_1_pos <- numeric(nrow(data))
  myhazard <- matrix(nrow=nrow(data), ncol=length(dx[[1]]))
  
  for(i in 1:nrow(data)){
    for (j in 1:length(alpha_sbf)){
      index<- which.min(abs(as.numeric(data[i,j])-x.grid[[j]]))
      error<-  x.grid[[j]][index] - as.numeric(as.numeric(data[i,j]))
      if (index==1|index==length(x.grid[[j]])|error==0) {index2 <- index}
      if (error>0&index>=2) {index2<-index-1} 
      if (error<0&index<length(x.grid[[j]])) {index2<-index+1} 
      error2<-  result$x.grid[[j]][index2] -   as.numeric(data[i,j])
      
      if (error!=0&index2>=1&index2<=length(x.grid[[j]]))  {
        hazard[i,j] <- (abs(error)*alpha_sbf[[j]][index]+abs(error2)*alpha_sbf[[j]][index2])/(abs(error)+abs(error2))
        hazard.1[i,j] <- (abs(error)*alpha.1_sbf[[j]][index]+abs(error2)*alpha.1_sbf[[j]][index2])/(abs(error)+abs(error2))
        if (j==1) {
          index_1 <- index
          index_2 <- index2
          temp <- matrix(NA,nrow=length(dx[[1]]),ncol=ncol(data))
        }
        } else{
          hazard[i,j] <- alpha_sbf[[j]][index]
          hazard.1[i,j] <- alpha.1_sbf[[j]][index]
          if (j==1) {
            index_1 <- index
            temp <- matrix(NA,nrow=length(dx[[1]]),ncol=ncol(data))
          }
        }
    }
    
    
    for (j in 1:length(alpha_sbf)){
      if (j==1){temp[,j] <-cumsum((dx[[1]]*alpha_sbf[[j]]))
      myhazard[i,] <- alpha_sbf[[j]]
      } else{    temp[,j] <- hazard[i,j]
      myhazard[i,] <- myhazard[i,] * hazard[i,j]
      }
    }
    myhazard[i,]  <-   pmax(0, myhazard[i,])
    myhazard[i,] <- cumsum(dx[[1]]*myhazard[i,])
    

    crps_1[i] <- sum(dx[[1]][1:index_1]*(1-exp(-apply(as.matrix(temp[1:index_1,]),1,prod)))^2)
    crps_2[i] <- sum(dx[[1]][index_1:length(dx[[1]])]*(exp(-apply(as.matrix(temp[index_1:length(dx[[1]]),]),1,prod)))^2)
    
    crps_1_pos[i] <- sum(dx[[1]][1:index_1]*(1-exp(-((myhazard[i,1:index_1]))))^2)
    crps_2_pos[i] <- sum(dx[[1]][index_1:length(dx[[1]])]*(exp(-((myhazard[i,index_1:length(dx[[1]])]))))^2)
    
    hazard[i,length( alpha_sbf)+1]<-prod(hazard[i,-c(1,(length( alpha_sbf)+1))])
    hazard.1[i,length( alpha.1_sbf)+1]<-prod(hazard.1[i,-(length( alpha.1_sbf)+1)])
    
  }
  
  # ##### one specific conditional survival pred
  # 
  # for (j in 1:length(alpha_sbf)){
  #   
  #   
  #   index<- which.min(abs(as.numeric(mine[j])-x.grid[[j]]))
  #   error<-  x.grid[[j]][index] - as.numeric(as.numeric(mine[j]))
  #   if (index==1|index==length(x.grid[[j]])) {index2 <- index}
  #   if (error>0&index>=2) {index2<-index-1} 
  #   if (error<0) {index2<-index+1} 
  #   error2<-  result$x.grid[[j]][index2] -   as.numeric(mine[j])
  #   
  #   if (error!=0&index2>=1&index2<=length(x.grid[[j]]))  {
  #     myhazard[j] <- (abs(error)*alpha_sbf[[j]][index]+abs(error2)*alpha_sbf[[j]][index2])/(abs(error)+abs(error2))
  #   } else{
  #     myhazard[j] <- alpha_sbf[[j]][index]
  #   }#}
  #   
  #   if (j==1) {
  #     index_1 <- index
  #   } 
  #   if(j==1) {surv_pred  <-sum((dx[[1]]*alpha_sbf[[j]])[1:index_1])} else{
  #     surv_pred <- surv_pred *myhazard[j] 
  #   }
  #   
  # } 
  # surv_pred <- exp(-surv_pred)
  return(list(hazard=hazard,hazard.1=hazard.1,crps_1=crps_1,crps_2=crps_2,crps_1_pos=crps_1_pos,crps_2_pos=crps_2_pos))
}
