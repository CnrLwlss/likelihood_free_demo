simDSLogistic=function(K,r,x0){
  if(x0>=K) return(data.frame(t=c(0,Inf),c=c(x0,x0)))
  unifs=runif(K-x0)
  clist=(x0+1):K
  dts=-log(1-unifs)/(r*clist*(1-clist/K))
  return(data.frame(t=c(0,cumsum(dts)),x=c(x0,clist)))
}

# Analytic solution of logistic model of population dynamics
logistic = function(K,r,x0,t) (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))

t = seq(0,10,0.5)
x = logistic(K=42.0, r=1.1, x0=1.0, t=t)
plot(t,x,type="b",main="Simulation from logistic model")

dat = simDSLogistic(K=42, r=1.1, x0=1)
plot(x~t,data=dat,type="s",main="Simulation from DSLM")


points(t,x,type="l",col="red",lwd=2)
points(x~t, data=simDSLogistic(K=42, r=1.1, x0=1),type="s")

t = seq(0,10,0.5)
x = logistic(K=42.0, r=1.1, x0=1.0, t=t)
plot(NULL,main="Simulations from logistic model & DSLM",xlim=range(t),ylim=c(0,42),xlab="t",ylab="x")
replicate(1000,points(x~t, data=simDSLogistic(K=42, r=1.1, x0=1),type="s",lwd=2,col=rgb(0,0,0,0.05)))
points(t,logistic(K=42.0, r=1.1, x0=1.0, t=t),col="red",lwd=3,type="l")



