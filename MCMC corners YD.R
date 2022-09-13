#some settings                          #integer of partitions
z<-0.8    # 0.9        #P(n-1)/P(n)                     #the value of z
N<-10000  #1000   #400000             #iterations of MCMC

#starting point  
lambda.start<-c(5,4,3,1,1,1)


#some functions 
# in this function how do I set 
add_y <- function(y, i) {
  sigma <- length(unique(y))
  if (i > sigma+1) {
    return(list(valid=FALSE))
  }
  if (i > sigma){ 
    y=c(y,1)
    return(list(valid=TRUE, y=y))
  } else
    pos <- which(y==unique(y)[i])[1]
  y[pos]<-y[pos]+1
  y<-y
  return(list(valid=TRUE, y=y))
} 


remove_y <- function(y, i) {
  sigma <- length(unique(y))
  pos <- max(which(y==unique(y)[i]))
  if (i >= sigma){
    return(list(valid=FALSE))
  } else if (y[1]==0){
    return(list(valid=FALSE))}
  if (y[1]!=0){
    y[pos]<-y[pos]-1
    y <- y}
    if (y[1]==0){ 
      y<-y
    } else {y= y[-which(y==0)] }
    return(list(valid=TRUE, y=y))
}

#MCMC method 
lambda<-vector("list",N)
lambda[[1]]<-lambda.start
lambda.abs<-matrix(0,nrow =N , ncol=1 )
lambda.abs[1,]<-sum(lambda.start)
rho<-c()
for (t in 2:N){
  sigma<-length(unique(lambda[[t-1]]))
  b<-sample(c(-1,1),1)
  if (b==1){
    i<-sample(sigma+1,1)
    rho<-add_y(lambda[[t-1]],i)$y
    if (add_y(lambda[[t-1]],i)$valid==TRUE){
      #accept with prob. min(alpha,1)  
      if (sum(lambda[[t-1]])>0) {
        alpha<-z*(sigma+1)/sigma
      } else {
        alpha <- z / 2
      }
      u<-runif(1,0,1)
      if (u<=alpha){
        lambda[[t]]<-rho
      }else{ #stay with prob. 1-alpha
        lambda[[t]]<-lambda[[t-1]]
      }
    }else{ #stay if after adding rho isn't valid
      lambda[[t]]<-lambda[[t-1]]
    }
    
  }else{     #when b==-1
    if(sum(lambda[[t-1]])>0){
      i<-sample(sigma,1)
      rho<-remove_y(lambda[[t-1]],i)$y
      if (remove_y(lambda[[t-1]],i)$valid==TRUE){
        #accept with prob. min(alpha,1)
        alpha<-z^(-1)*sigma/(sigma+1)
        u<-runif(1,0,1)
        if(u<=alpha){
          lambda[[t]]<-rho
        }else{
          lambda[[t]]<-lambda[[t-1]]
        } 
      }else{ #stay if after removing rho isn't valid
        lambda[[t]]<-lambda[[t-1]]
      }
    }else{ #if the partition goes 0 force it to be 1
      alpha <- z / 2
      u <- runif(1, 0, 1)
      if (u <= alpha) {
        lambda[[t]]=1
      } else {
        lambda[[t]]<-lambda[[t-1]]
      }
    }
  }
  lambda.abs[t]<-sum(lambda[[t]])
}

mean(lambda.abs)
plot.ts(lambda.abs)

hist(lambda.abs, probability = TRUE, 
     col="grey",
     breaks=seq(-0.5, max(lambda.abs)+1, by=1), 
     main = expression(paste ("Histogram of ", abs(lambda))),
     xlab = expression(abs(lambda)))

ll <- seq(0, max(lambda.abs)+1, by=1)
e <- numeric(length(ll))
for (i in 1:length(ll)) {
  l <- ll[i]
  e[i] <- z^l 
}
lines(ll, e/sum(e), col="red", lwd=1.5)
acf(lambda.abs)





pdf("MHYDTs.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
plot.ts(lambda.abs)
dev.off()

pdf("MHYDHistogram.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
hist(lambda.abs, probability = TRUE, 
     col="grey",
     breaks=seq(-0.5, max(lambda.abs)+1, by=1), 
     main = expression(paste ("Histogram of ", abs(lambda))),
     xlab = expression(abs(lambda)))
lines(ll, e/sum(e), col="red", lwd=1.5)
dev.off()


pdf("MHYDAcf.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
acf(lambda.abs)
dev.off()


(mean(lambda.abs))
(var(lambda.abs))












