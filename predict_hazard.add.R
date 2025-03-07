predict.hazard.add<-function(result, data){
  data<-as.matrix(data)
  if (ncol(data)==1) data<-t(data)
  x.grid<-result$x.grid
  dx<-result$dx
  
  dx<-lapply(1:ncol(data),function(k){  ddx<-diff(x.grid[[k]])
  c( 0  ,ddx)
  }
  )
  


  alpha_sbf<-result$alpha_backfit
  alpha.1_sbf<-result$alpha.1_backfit
  
  hazard <- matrix(nrow=nrow(data), ncol=(length( alpha_sbf)+1))
  hazard.1 <- matrix(nrow=nrow(data), ncol=(length( alpha.1_sbf)+1))
  cumhazard <- numeric(nrow(data))
  crps_1 <- crps_2 <- crps_2_pos <- crps_1_pos <- numeric(nrow(data))
  myhazard <- matrix(nrow=nrow(data), ncol=length(dx[[1]]))
  
  for(i in 1:nrow(data)){
    temp <- matrix(NA,nrow=length(dx[[1]]),ncol=ncol(data))
    for (j in 1:length(alpha_sbf)){
     # if (j==1) l= 8
     # if (j==2) l= 12
     # if (j==3) l= 2
      
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
         cumhazard[i] <- (abs(error)*sum((dx[[j]]*alpha_sbf[[j]])[1:index])+abs(error2)*sum((dx[[j]]*alpha_sbf[[j]])[1:index2]))/(abs(error)+abs(error2))
         } #else{
         #cumhazard[i] <-  (abs(error)*sum((dx[[1]][1:index_1]*alpha_sbf[[j]][index]))+abs(error2)*sum((dx[[1]][1:index_2]*alpha_sbf[[j]][index2])))/(abs(error)+abs(error2))
      # }
        } else{
          hazard[i,j] <- alpha_sbf[[j]][index]
          hazard.1[i,j] <- alpha.1_sbf[[j]][index]
          if (j==1) {
            index_1 <- index
            cumhazard[i] <- sum((dx[[1]]*alpha_sbf[[j]])[1:index])} #else {
            #cumhazard[i] <- sum((dx[[1]][1:index_1]*alpha_sbf[[j]][index]))
          
    } 
    }
    
    for (j in 1:length(alpha_sbf)){
      if (j==1){temp[,j] <-cumsum((dx[[1]]*alpha_sbf[[j]]))
      myhazard[i,] <- alpha_sbf[[j]]
      } else{    temp[,j] <-cumsum(dx[[1]]*hazard[i,j])
      myhazard[i,] <- myhazard[i,] + hazard[i,j]
      }
    }
     myhazard[i,]  <-   pmax(0, myhazard[i,])
      myhazard[i,] <- cumsum(dx[[1]]*myhazard[i,])
    
  

    crps_1[i] <- sum(dx[[1]][1:index_1]*(1-exp(-rowSums(as.matrix(temp[1:index_1,]))))^2)
    crps_2[i] <- sum(dx[[1]][index_1:length(dx[[1]])]*(exp(-rowSums(as.matrix(temp[index_1:length(dx[[1]]),]))))^2)
    
    crps_1_pos[i] <- sum(dx[[1]][1:index_1]*(1-exp(-((myhazard[i,1:index_1]))))^2)
    crps_2_pos[i] <- sum(dx[[1]][index_1:length(dx[[1]])]*(exp(-((myhazard[i,index_1:length(dx[[1]])]))))^2)
    
    hazard[i,length( alpha_sbf)+1]<-sum(hazard[i,-c(1,(length( alpha_sbf)+1))])
    hazard.1[i,length( alpha.1_sbf)+1]<-sum(hazard.1[i,-c(1,(length( alpha.1_sbf)+1))])
    
  }
  
 
  return(list(hazard=hazard,hazard.1=hazard.1,cumhazard=cumhazard,crps_1=crps_1,crps_2=crps_2,crps_1_pos=crps_1_pos,crps_2_pos=crps_2_pos))
}
