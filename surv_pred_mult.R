predict.surv.mult<-function(result, data, t_grid, other_points){
  data<-as.matrix(data)
  if (ncol(data)==1) data<-t(data)
  x.grid<-result$x.grid
 
  
  dx<-lapply(1:ncol(data),function(k){  ddx<-diff(x.grid[[k]])
  c(  (x.grid[[k]][1]-x.min[k])  ,ddx)
  }
  )
 #dx<-result$dx
  
  #dx<-lapply(1:ncol(data),function(k){  ddx<-diff(x.grid[[k]])
 # c( 0  ,ddx)
  #}
 # )
  
  
  
  alpha_sbf<-result$alpha_backfit

  cumhazard <- numeric(nrow(data))
  crps_1 <- crps_2 <- numeric(nrow(data))
  myhazard <- numeric(nrow(data))
  surv_pred <- numeric(length(t_grid))
  myhazard_f <- numeric(length(dx[[1]]))
  
  ##### one specific conditional survival pred
  
    mine <- c(1,other_points)
    
  for (j in 1:length(alpha_sbf)){
    
    
    index<- which.min(abs(as.numeric(mine[j])-x.grid[[j]]))
    error<-  x.grid[[j]][index] - as.numeric(as.numeric(mine[j]))
    if (index==1|index==length(x.grid[[j]])) {index2 <- index}
    if (error>0&index>=2) {index2<-index-1} 
    if (error<0) {index2<-index+1} 
    error2<-  result$x.grid[[j]][index2] -   as.numeric(mine[j])
    
    if (error!=0&index2>=1&index2<=length(x.grid[[j]]))  {
      myhazard[j] <- (abs(error)*alpha_sbf[[j]][index]+abs(error2)*alpha_sbf[[j]][index2])/(abs(error)+abs(error2))
    } else{
      myhazard[j] <- alpha_sbf[[j]][index]
    }#}
  }
    
    for (j in 1:length(alpha_sbf)){
      if (j==1){
        myhazard_f <- alpha_sbf[[j]]
      } else{    
        myhazard_f <- myhazard_f * myhazard[j]
      }
    }
    
    myhazard_f  <-   pmax(0, myhazard_f)

  for(t in 1:length(t_grid)){
 
    index<- which.min(abs(as.numeric(t_grid[t])-x.grid[[1]]))
  
   surv_pred[t]  <-sum((dx[[1]]*myhazard_f)[1:index])
  surv_pred[t] <- exp(-surv_pred[t])
  }
    
    
  return( surv_pred)
}
