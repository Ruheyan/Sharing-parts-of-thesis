#MCMC on Young Diagrams                          #integer of partitions
library(partitions)
set.seed(10) #6 gives good histogram

n <- 120 
z<- exp(-pi/sqrt(6*n)) #1- pi/sqrt(6*n)  #0.8796291   #or P(n-1)/P(n)  

N<-500000          #iterations of MCMC
lambda.start<-c(1)   #initial state (starting point)

#some functions 
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
  if (i > sigma){
    return(list(valid=FALSE))
  } else if (y[1]==0){
    return(list(valid=FALSE))}
  if (y[1]!=0){
    y[pos]<-y[pos]-1
    y=y
  }
  if (y[1]!=0){
    y <-y[y!=0]
  }else{ 
    y<-y
  } 
  return(list(valid=TRUE, y=y))
}

#MCMC method 
lambda<-vector("list",N)
lambda[[1]]<-lambda.start
lambda.abs<-matrix(0,nrow =N , ncol=1 )
lambda.abs[1,]<-sum(lambda.start)
lambda.parts<-matrix(0,nrow =N , ncol=1 )
lambda.parts[1,]<-length(lambda.start)
lambda.largest<-matrix(0,nrow =N , ncol=1 )
lambda.largest[1,]<-max(lambda.start)
rho<-c()
for (t in 2:N){
  s <- lambda[[t-1]]
  sigma<-length(unique(s))
  if (sum(s)>0){
    b<-sample(c(-1,1),1)
    if (b==1){
      i<-sample(sigma+1,1)
      rho<-add_y(s,i)$y
      alpha<-z*(sigma+1)/sigma
      u<-runif(1,0,1)
      if (u<=alpha){
        lambda[[t]]<-rho
      }else{ #stay with prob. 1-alpha
        lambda[[t]]<-s
      }
    }else{     #when b==-1
      i<-sample(sigma,1)
      rho<-remove_y(s,i)$y
      #accept with prob. min(alpha,1)
      alpha<-z^(-1)*sigma/(sigma+1)
      u<-runif(1,0,1)
      if(u<=alpha){
        lambda[[t]]<-rho
      }else{
        lambda[[t]]<-s
      }
    }
  } else {
    if (sum(s==0)){
      b<-sample(c(0,1),1)
      if (b==1){
        rho=c(1)
        alpha<-z
        u<-runif(1,0,1)
        if (u<=alpha){
          lambda[[t]]<-rho
        }else{ #stay with prob. 1-alpha
          lambda[[t]]<-s
        }
      } else{ # if b=0
        rho=c(0) 
        alpha<-1
        u<-runif(1,0,1)
        if (u<=alpha){
          lambda[[t]]<-rho
        }else{ #stay with prob. 1-alpha
          lambda[[t]]<-s
        }
      }
    }
  }
  lambda.abs[t]<-sum(lambda[[t]])
  lambda.parts[t]<-length(lambda[[t]])
  lambda.largest[t]<-max(lambda[[t]])
}

mean(lambda.abs)
plot.ts(lambda.abs,ylim=c(-20,300), ylab=expression(abs(lambda)))

hist(lambda.abs, probability = TRUE, 
     col="grey",
     breaks=seq(-0.5, max(lambda.abs)+1, by=1), 
     main = expression(paste ("Histogram of ", abs(lambda))),
     xlab = expression(abs(lambda)))

ll <- seq(0, max(lambda.abs)+1, by=1)
e <- numeric(length(ll))
for (i in 1:length(ll)) {
  l <- ll[i]
  e[i] <- z^l * P(l)
}
lines(ll, e/sum(e), col="red", lwd=1.5)


acf(lambda.abs, main="")
lag <- 1:N
lambda.abs.newlag<-lambda.abs[lag[which(lag%% 10==1)]]
acf(lambda.abs.newlag,main=" ", ylim=c(0,1))


pdf("MHYDTs.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
plot.ts(lambda.abs,ylim=c(-20,300), ylab=expression(abs(lambda)))
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
acf(lambda.abs.newlag,main=" ", ylim=c(0,1))
dev.off()


(mean(lambda.abs))
(var(lambda.abs))

#### for the length of partitions

plot.ts(lambda.parts, ylab="M", main=" ",ylim=c(-10,90))
hist(lambda.parts, probability = TRUE, col="grey",
     breaks=seq(-0.5, max(lambda.parts)+1, by=1), main="Histogram of M",xlab="M")
lambda.parts.newlag<-lambda.parts[lag[which(lag%% 10==1)]]
acf(lambda.parts)
acf(lambda.parts.newlag)


pdf("MHYDpartsTs.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
plot.ts(lambda.parts, ylab="M", main=" ",ylim=c(-10,90))
dev.off()

pdf("MHYDpartsHistogram.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
hist(lambda.parts, probability = TRUE, col="grey",
     breaks=seq(-0.5, max(lambda.parts)+1, by=1), main="Histogram of M",xlab="M")
dev.off()


pdf("MHYDpartsAcf.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
acf(lambda.parts.newlag)
dev.off()



#### for the largest part
plot.ts(lambda.largest,ylim=c(-5,75))
hist(lambda.largest,probability = TRUE, col="grey", 
     breaks = seq(-0.5,max(lambda.largest)+1,by=1))
acf(lambda.largest)




library(squash)
#hist2(lambda.parts, lambda.abs,zlab = "Counts",colFn =grayscale, xlab="M", ylab=expression(abs(lambda)),breaks = pretty,key=NULL) 

#colFn can be heat, jet, rainbow,grayscale
#breaks can be pretty or prettyLog

#lines(2.5651^(-1)*sqrt(lambda.abs)*log(lambda.abs)+0.2*sqrt(lambda.abs),lambda.abs,col="red")

hist2(lambda.abs,lambda.parts ,zlab = "Counts",colFn =grayscale, 
      ylab="M", xlab=expression(abs(lambda)),breaks = pretty,key=NULL) 
lines(lambda.abs,2.5651^(-1)*sqrt(lambda.abs)*log(lambda.abs)+0.1*sqrt(lambda.abs),col="red")


pdf("PartsandAbsHist.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
hist2(lambda.abs,lambda.parts ,zlab = "Counts",colFn =grayscale, 
      ylab="M", xlab=expression(abs(lambda)),breaks = pretty,key=NULL) 
lines(lambda.abs,2.5651^(-1)*sqrt(lambda.abs)*log(lambda.abs)+0.1*sqrt(lambda.abs),col="red")
dev.off()
