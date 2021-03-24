library(MASS)
library(parallel)

# Analytic solution of logistic model of population dynamics
logistic = function(K,r,x0,t) (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))

# Alternative version of analytic solution of logistic model using Rcpp (~20% faster)
if(require(Rcpp)){
 cppFunction('NumericVector logistic(float K, float r, float x0, NumericVector t){
  NumericVector x = (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1.0));
  return x;
 }')
}

# Discrete stochastic simulation of logistic model adapted from http://lwlss.net/talks/discstoch/
simDSLogistic=function(K,r,N0){
  if(N0>=K) return(data.frame(t=c(0,Inf),c=c(N0,N0)))
  unifs=runif(K-N0)
  clist=(N0+1):K
  dts=-log(1-unifs)/(r*clist*(1-clist/K))
  return(data.frame(t=c(0,cumsum(dts)),c=c(N0,clist)))
}

stepLogistic = function(x0, t0, deltat, th) logistic(th[["K"]],th[["r"]],x0,deltat)

stepDSL = function(x0, t0, deltat, th){
  dsl = simDSLogistic(th[["K"]],th[["r"]],x0)
  if(deltat<=max(dsl$t)){
    dsl$c[which(dsl$t>deltat)[1]]
  }else{
    th[["K"]]
  }
}

# pfMLLik
# Particle filter for unbiased estimate of marginal likelihood (logged)
# Adapted from Stochastic Modelling for Systems Biology (3rd Ed.) by Darren Wilkinson

pfMLLik_gen = function (n, simx0, t0, stepFun, dataLik, data) {
  times = c(t0, as.numeric(rownames(data)))
  deltas = diff(times)
  return(function(th) {
    #xmat = as.matrix(simx0(n, t0))
    xmat = as.matrix(replicate(n,th[["x0"]]))
    ll = 0
    for (i in 1:length(deltas)) {
      xmat[] = apply(xmat[,drop=FALSE], 1, stepFun, t0 = times[i], deltat = deltas[i], th)
      lw = apply(xmat[,drop=FALSE], 1, dataLik, t = times[i + 1], y = data[i,], th)
      m = max(lw)
      rw = lw - m
      sw = exp(rw)
      ll = ll + m + log(mean(sw))
      rows = sample(1:n, n, replace = TRUE, prob = sw)
      xmat[] = xmat[rows,]
    }
    ll
  })
}

# Log likelihood calcuated directly from analytical solution (only possible because of simple model)
regularMLLik = function(dat){
  return(function(th)  sum(log(dnorm(dat$x,mean=logistic(th[["K"]],th[["r"]],th[["x0"]],as.numeric(rownames(dat))),sd=th[["stdev"]]))))
}

# Generate and visualise synthetic observations (and "true" dynamics) from logistic model
th = c(K=42,r=1.1,x0=1,stdev=2.5)

times = seq(0.2,7,length.out=5)
dat = data.frame(x = sapply(logistic(th[["K"]],th[["r"]],th[["x0"]],times),function(y) rnorm(1,y,th[["stdev"]])))
rownames(dat) = times

plot(dat$x~times,xlab="t",ylab="x(t)",main="Synthetic data from logistic model",ylim=c(0,max(c(th[["K"]],dat$x))),pch=16,col="red")
newf = function(x) logistic(th[["K"]],th[["r"]],th[["x0"]],x)
curve(newf,from=min(times),to=max(times),col="red",lwd=2,add=TRUE)

# Generate and visualise synthetic observations (and "true" dynamics) from DSLM

dsl = simDSLogistic(th[["K"]],th[["r"]],th[["x0"]])
dsl_curve = approxfun(dsl$t,dsl$c,method="constant")
dat_dsl = data.frame(x = sapply(dsl_curve(times),function(y) rnorm(1,y,th[["stdev"]])))
rownames(dat_dsl)=times

plot(dsl$t,dsl$c,type="s",col="blue",lwd=2,xlab="t",ylab="x(t)",main="Synthetic data from discrete stochastic logistic model",ylim=c(0,max(c(th[["K"]],dat_dsl$x))))
points(dat_dsl$x~times,pch=16,col="blue")
                                
png("compare.png",width=1000,height=1000,pointsize=24)
op = par(mar=c(5,5,3,3))
plot(NULL,xlim=c(0,10),ylim=c(0,40),xlab="Time (years)",ylab="Population size",cex.axis=1.5,cex.lab=2.55,main="Simulations from logistic model & DSLM",cex.main=1.55)
for(i in 1:500){
  dsl = simDSLogistic(th[["K"]],th[["r"]],th[["x0"]])
  points(dsl$t,dsl$c,col=rgb(0,0,0,0.15),lwd=2,type="s")
}
curve(newf,from=0,to=10,col="red",lwd=4,add=TRUE)
par(op)
dev.off()

simx0 = function(n, t0, th) rlnorm(n, meanlog=0,sdlog=2.5)
simx0 = function(n, t0, th) runif(n, 0, 10)
simx0 = function(n, t0, th) sample(1:4,n,replace=TRUE)

dataLik = function(x, t, y, th) sum(dnorm(y, x, th[["stdev"]],log=TRUE))

mLLik = pfMLLik_gen(1,simx0,0,stepLogistic,dataLik,dat)
#mLLik = regularMLLik(dat)
#mLLik = pfMLLik_gen(100,simx0,0,stepDSL,dataLik,dat_dsl)

priorlik = function(th){
  pK = dnorm(th[["K"]],mean=20,sd=20)
  pr = dnorm(th[["r"]],mean=2,sd=0.5)
  px0 = dunif(th[["x0"]],0,10)
  pstdev = dunif(th[["stdev"]],0,10)
  return(pK*pr*px0*pstdev)
}

failth = function(th) th[["K"]]<=0|th[["r"]]<=0|th[["x0"]]<=0|th[["stdev"]]<=0|th[["x0"]]>th[["K"]]

relprob = function(th){
  if(failth(th)){
   rprob = 0
  }else{
   rprob = exp(mLLik(th))*priorlik(th)
  }
  return(rprob)
}

nsamps = 210000
burnin = 10000
thin = 100
thc = c(30,2.0,2,4)
names(thc)=names(th)
ndim = length(thc)

trajectory = matrix(0, nrow=nsamps/thin, ncol = ndim)
trajectory[1,] = thc
colnames(trajectory)=c("K","r","x0","stdev")

nAccepted = 0
nRejected = 0

sd1 = 1.3; sd2 = 0.13; sd3 = 0.13; sd4 = 0.26
covarmat = matrix(c(sd1^2, 0,0,0,0,sd2^2,0,0,0,0,sd3^2,0,0,0,0,sd4^2), nrow=ndim, ncol=ndim)

rrow = 1
for(stepno in 1:(nsamps-1)){
  if((stepno-1)%%thin==0){
    trajectory[rrow,] = thc
    rrow = rrow + 1
  }
  proposedjump = mvrnorm(1,mu=rep(0,ndim),Sigma = covarmat)
  probaccept = min(1, relprob(thc+proposedjump)/relprob(thc))
  if(runif(1) < probaccept){
    thc = thc+proposedjump
    if(stepno>burnin) nAccepted = nAccepted + 1
  }else{
    thc = thc
    if(stepno>burnin) nRejected = nRejected + 1
  }
}

# Choose covarmat, aiming for acceptance ratio of about 0.2?
sprintf("Proposals accepted %i",nAccepted)
sprintf("Proposals rejected %i",nRejected)
nAccepted/(nAccepted+nRejected)

pretrajectory = trajectory[1:(burnin/thin),]
trajectory = trajectory[(burnin/thin+1):dim(trajectory)[1],]

traceplot = function(param,chain) {
  plot(chain[,param],type="l",ylab=param)
  abline(h=th[[param]],col="red",lwd=2)
}

pRngs = list(
K = c(0,60),
r = c(0,3),
x0 = c(0,10),
stdev = c(0,10)
)

densplot = function(param,chain,pRngs) {
  plot(density(chain[,param]),type="l",xlab=param,main="",xlim=pRngs[[param]])
  pfunc = pfuncs[[param]]
  curve(pfunc,from=pRngs[[param]][1],to=pRngs[[param]][2],col="green",lwd=2,add=TRUE)
  abline(v=th[[param]],col="red",lwd=2)
}

pfuncs = c(
K = function(x) dnorm(x,mean=20,sd=20),
r = function(x) dnorm(x,mean=2,sd=1.5),
x0 = function(x) dunif(x,0,10),
stdev = function(x) dunif(x,0,10)
)

op=par(mfrow=c(2,2))
traces = sapply(c("K","r","x0","stdev"),traceplot,trajectory)
par(op)

op=par(mfrow=c(2,2))
denses = sapply(c("K","r","x0","stdev"),densplot,trajectory,pRngs)
par(op)

pairs(trajectory,pch=16,col=rgb(0,0,0,0.3),cex=0.5)
