library(MASS)

set.seed(1234)

logist = function(K,r,x0,t) (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))

th = c(K=1,r=1.1,x0=0.01,stdev=0.05)

times = seq(0.2,10,1.5)
dat = data.frame(
  t = times,
  x = sapply(logist(th[["K"]],th[["r"]],th[["x0"]],times),function(y) rnorm(1,y,th[["stdev"]]))
)

plot(x~t,data=dat)
newf = function(x) logist(th[["K"]],th[["r"]],th[["x0"]],x)
curve(newf,from=min(dat$t),to=max(dat$t),col="red",lwd=2,add=TRUE)

lik = function(th,dat){
  sum(log(dnorm(dat$x,mean=logist(th[["K"]],th[["r"]],th[["x0"]],dat$t),sd=th[["stdev"]])))
}

priorlik = function(th){
  pK = dnorm(th[["K"]],mean=2,sd=0.5)
  pr = dnorm(th[["r"]],mean=2,sd=0.5)
  px0 = dlnorm(th[["x0"]],meanlog=0,sdlog=2.5)
  pstdev = dlnorm(th[["stdev"]],meanlog=0,sdlog=2.5)
  return(pK*pr*px0*pstdev)
}

pfuncs = c(
K = function(x) dnorm(x,mean=2,sd=2),
r = function(x) dnorm(x,mean=2,sd=2),
x0 = function(x) dlnorm(x,meanlog=0,sdlog=2.5),
stdev = function(x) dlnorm(x,meanlog=0,sdlog=2.5)
)

pRngs = list(
K = c(0,3),
r = c(0,3),
x0 = c(0,0.1),
stdev = c(0,0.1)
)

failth = function(th) th[["K"]]<0|th[["r"]]<0|th[["x0"]]<0|th[["stdev"]]<0|th[["x0"]]>th[["K"]]

relprob = function(th,dat){
  if(failth(th)){
   rprob = 0
  }else{
   rprob = exp(lik(th,dat))*priorlik(th)
  }
  return(rprob)
}

nsamps = 510000
burnin = 10000
ndim = 4

trajectory = matrix(0, nrow=nsamps, ncol = ndim)
trajectory[1,] = c(2,2,0.5,0.5)
colnames(trajectory)=c("K","r","x0","stdev")

nAccepted = 0
nRejected = 0

sd1 = 0.05; sd2 = 0.05; sd3 = 0.005; sd4 = 0.005
covarmat = matrix(c(sd1^2, 0,0,0,0,sd2^2,0,0,0,0,sd3^2,0,0,0,0,sd4^2), nrow=ndim, ncol=ndim)

for(stepno in 1:(nsamps-1)){
  thc = trajectory[stepno,]
  proposedjump = mvrnorm(1,mu=rep(0,ndim),Sigma = covarmat)
  probaccept = min(1, relprob(thc+proposedjump,dat)/relprob(thc,dat))
  if(runif(1) < probaccept){
    trajectory[stepno+1,] = thc+proposedjump
    if(stepno>burnin) nAccepted = nAccepted + 1
  }else{
    trajectory[stepno+1,] = thc
    if(stepno>burnin) nRejected = nRejected + 1
  }
}

accepted = trajectory[(burnin+1):dim(trajectory)[1],]

thinned = accepted[seq(1,dim(accepted)[1],length.out=1000),]

traceplot = function(param,chain) {
  plot(chain[,param],type="l",ylab=param)
  abline(h=th[[param]],col="red",lwd=2)
}

densplot = function(param,chain) {
  plot(density(chain[,param]),type="l",xlab=param,main="",xlim=pRngs[[param]])
  pfunc = pfuncs[[param]]
  curve(pfunc,from=pRngs[[param]][1],to=pRngs[[param]][2],col="green",lwd=2,add=TRUE)
  abline(v=th[[param]],col="red",lwd=2)
}

op=par(mfrow=c(2,2))
sapply(c("K","r","x0","stdev"),traceplot,thinned)
par(op)

op=par(mfrow=c(2,2))
sapply(c("K","r","x0","stdev"),densplot,thinned)
par(op)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(thinned,upper.panel=panel.cor,pch=16,cex=0.25,cex.cor=6,col=rgb(1,0,0,0.3))


