library(smfsb)
library(Rcpp)

# pfMLLik
# Particle filter for unbiased estimate of marginal likelihood (logged)
# Adapted from Stochastic Modelling for Systems Biology (3rd Ed.) by Darren Wilkinson

pfMLLik_gen = function (n, simx0, t0, stepFun, dataLik, data) {
  times = c(t0, as.numeric(rownames(data)))
  deltas = diff(times)
  return(function(...) {
    xmat = as.matrix(simx0(n, t0))
    ll = 0
    for (i in 1:length(deltas)) {
      xmat[] = apply(xmat[,drop=FALSE], 1, stepFun, t0 = times[i], deltat = deltas[i], ...)
      lw = apply(xmat[,drop=FALSE], 1, dataLik, t = times[i + 1], y = data[i,], ...)
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

# Logistic model of population dynamics
logistic = function(K,r,x0,t) (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1))

cppFunction('float logistic(float K, float r, float x0, float t){
  float x = (K*x0*exp(r*t))/(K+x0*(exp(r*t)-1.0));
  return x;
}')

# Generate and visualise synthetic observations (and "true" dynamics"
th = c(K=1,r=1.1,x0=0.01,stdev=0.05)

times = seq(0.2,10,1.5)
dat = data.frame(x = sapply(logistic(th[["K"]],th[["r"]],th[["x0"]],times),function(y) rnorm(1,y,th[["stdev"]])))
rownames(dat) = times

plot(dat$x~times,xlab="t",ylab="x(t)",main="Synthetic data from logistic model")
newf = function(x) logistic(th[["K"]],th[["r"]],th[["x0"]],x)
curve(newf,from=min(times),to=max(times),col="red",lwd=2,add=TRUE)

stepLogistic = function(x0, t0, deltat, th) logistic(th[["K"]],th[["r"]],th[["x0"]],deltat)
simx0 = function(n, t0, th) rlnorm(n, meanlog=0,sdlog=2.5)
dataLik = function(x, t, y, th) sum(dnorm(y, x, th[["stdev"]]))

mLLik=pfMLLik_gen(100,simx0,0,stepLogistic,dataLik,dat)

iters=1000
tune=0.01
thin=10

starting = Sys.time()

p=length(th)
ll=-1e99
thmat=matrix(0,nrow=iters,ncol=p)
colnames(thmat)=names(th)
# Main pMCMC loop
for (i in 1:iters) {
    message(paste(i,""),appendLF=FALSE)
    for (j in 1:thin) {
        thprop=th*exp(rnorm(p,0,tune))
        llprop=mLLik(thprop)
        if (log(runif(1)) < llprop - ll) {
            th=thprop
            ll=llprop
            }
        }
    thmat[i,]=th
    }
message("Done!")
# Compute and plot some basic summaries
mcmcSummary(thmat[,])

ending = Sys.time()

print(ending-starting)


