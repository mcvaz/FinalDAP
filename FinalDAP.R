#### FINAL DAP R SCRIPT

### PACKAGES
# Loading packages
library(R2jags)
library(lattice)
library(ggplot2)
library(xtable)
source("AddBurnin.R")


### DATA
dd = read.csv("HGEvsAS.csv",header=T,as.is=T) # data set
head(dd)
str(dd)

ss = length(unique(dd$sp)) # total number of species

is = table(dd$sp) # sampling size of each species

hist(log(dd$hge),main="",xlab="log(HGE)")
hist(dd$aSize,xlab="Adult size (m)",main="")
hist(aggregate(aSize~sp,dd,mean)$aSize,xlab="Adult size (m)",main="")

plot(hge~aSize,dd)
lines(predict(loess(hge~aSize,dd)),col=2,lwd=2)
summary(lm(hge~aSize,dd))

## Response variable: heigh gain efficiency (HGE; m/kg)
yy = matrix(,ncol=max(is),nrow=ss)
for(i in 1:nrow(yy)){
  yy[i,1:is[i]] = dd[dd$sp==names(is)[i],"hge"]
  }
head(yy)

## Predictor variable: tree adult height (m)
xx = matrix(,ncol=max(is),nrow=ss)
for(i in 1:nrow(xx)){
  xx[i,1:is[i]] = dd[dd$sp==names(is)[i],"aSize"]
  }
head(xx)


### JAGS MODEL
sink("model.txt")
  cat("
    model{

      # LIKELIHOOD
      for(i in 1:ss){
        for(j in 1:is[i]){
          y[i,j] =  beta0 + x[i,j]*beta1 + beta[i]+epsilon[i,j]
          epsilon[i,j] ~ dnorm(0,tau[i])
          }
        beta[i] ~ dnorm(0,taub)
        tau[i] ~ dunif(ai,bi)
        }

      # PRIORS
      beta0 ~ dnorm(m0,tau0)
      beta1 ~ dnorm(0,tau1)
      
    }    
    ", fill = TRUE)
sink()


### PRIORS
## Extreme cases: 
# inefficient tree: DSL=5cm, sH=150cm, cPos=1m, sDens=1.35g/cm^3
#maxHGE = 1/((pi*(5/2)^2*150/3)*1.35)
# efficient tree: DSL=2.5cm, sH=300cm, cPos=3m, sDens=0.1g/cm^3
#minHGE = 3/((pi*(2.5/2)^2*300/3)*0.1)
#dd[which(dd$hge==min(dd$hge)),]
#dd[which(dd$hge==max(dd$hge)),]

range(dd$hge) # range of HGE used to define the priors
range(dd$aSize) # range of adult sizes used to define priors

## beta0
aa0 = diff(range(dd$hge)) # range of possible values for the intercept is the same of the range for the data
tau0 = 4/aa0 # precision of beta0 according to the range method
m0 = aa0/2 # mean of beta0 assumed to be at the center of the range

## beta1
aa1 = diff(range(dd$hge))/diff(range(dd$aSize)) # maximum possible value for beta1
tau1 = 4/2*aa1 # precision of beta1 according to the range method

## beta[i]
taub = tau0 # precision of beta[i] is set to be the same of beta0

## tau[i] 
ai = 0.01 # minimum allowed value for tau[i]
bi = 100 # maximum allowed precision for tau[i]


### JAGS RUN
# Jags input info
data = list(n=nn, y=yy, N=nrow(nn), # raw data
            mm=0, mt=1/3.36, tm1=0, tm2=1/3.36, dm=0, dt=1/2.17, td1=0, td2=1/2.17) 
inits = rep(list(list(
  pie=matrix(rep(0.5,22), ncol=2, nrow=11, byrow=TRUE), 
  mu=rep(0,11), delta=rep(0,11),
  mu0=0, tau.mu=0.01, delta0=0, tau.delta=0.01 
)),3) # number of chains
params = c("pie[1:11,1:2]","mu[1:11]","delta[1:11]","mu0","delta0","tau.mu","tau.delta")

# Jags output
model.out = jags(data=data, inits=NULL, parameter=params, "model.txt", n.chains=3, n.iter=101000, n.burnin=0, n.thin=20, DIC=F)
outparams = dimnames(model.out$BUGSoutput$sims.array)[[3]] # all parameter names estimated by jags

# Take the first 1000 runs as burnin
Output = AddBurnin(model.out$BUGSoutput$sims.array, burnin=1000,n.thin=1)






