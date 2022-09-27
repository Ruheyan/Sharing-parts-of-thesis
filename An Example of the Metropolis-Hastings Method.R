### Basic MCMC example 
#s={0,1,2,3,4,...}
#pi(s)=z^s(1-z)   Stationary distribution (Target distribution)
#p(s,s+1)=1/2
#p(s,s-1)=1/2

set.seed(6)

z <- 0.85

N<-10000      #iterations of MCMC 10000
s0<-1          #initial state (starting point)

out<-numeric(N)
out[1]<-s0
for (t in 2:N){
  s<-out[t-1]
  if (s>0){
    b<-sample(c(1,-1),size=1)
  }else{
    b<-sample(c(0,1),size=1)
  }
  sp<-s+b
  alpha <- z^b
  u<-runif(1,0,1)
  if (u<=alpha){
    out[t] <- sp
  }else{
    out[t] <- s
  }
}
plot.ts(out,ylab="s",ylim=c(0,50))

hist(out,xlab = "s",main = " ",probability=TRUE,breaks=seq(-0.5,60.5,by=1),xlim=c(0,20))#target distribution
e <- numeric(101)
for (i in 0:100) {
  e[i+1] <- z^i*(1-z)
}
lines(0:100, e,col="red",lwd=1.5)

acf(out,main=" ")


pdf("MetropolisExampleTS.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
plot.ts(out,ylab="s",ylim=c(0,50))
dev.off()

pdf("MetropolisExampleHistogram.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
hist(out,xlab = "s",main = " ",probability=TRUE,breaks=seq(-0.5,60.5,by=1),xlim=c(0,20))#target distribution
lines(0:100, e,col="red",lwd=1.5)
dev.off()


pdf("MetropolisExampleACF.pdf",width = 5,height= 3.5)
par(mar = c(4, 4, 2, 2),xaxs= "r",yaxs  = "r",cex.axis = 1,cex.lab  = 1)
acf(out2,main=" ")
dev.off()


### to check the result
(mean.out<-mean(out2) )
(var.out<-var(out2) )
(mean.geo<-z/(1-z) )
(var.geo<-z/((1-z)^2) )





