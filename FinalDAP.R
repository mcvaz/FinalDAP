#### FINAL DAP R SCRIPT

### PACKAGES
# Loading packages
library(R2jags)
library(lattice)
library(ggplot2)
library(xtable)

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
ai = 0.001 # minimum allowed value for tau[i]
bi = tau0 # maximum allowed precision for tau[i]





