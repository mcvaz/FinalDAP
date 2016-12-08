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
dd$hge = dd$hge/10 # scale correction
head(dd)
str(dd)

ss = length(unique(dd$sp)) # total number of species
is = table(dd$sp) # sampling size of each species

## Table 0. Data
tab0 = aggregate(cbind(hge,aSize)~sp,dd,mean)
tab0=cbind(tab0,as.numeric(is))
names(tab0) = c("Species","HGE","Size","n")
tab0a = tab0[1:ceiling(nrow(tab0)/2),]
tab0b = rbind(tab0[ceiling(nrow(tab0)/2+1):nrow(tab0),],rep(NA,4))
tab0f = cbind(tab0a,cbind(row.names(tab0b),tab0b))
head(tab0f)
tail(tab0f)
xtable(tab0f)


## Figure 0. Data
require(mvtnorm)
source("scatterBarNorm.R")
scatterBarNorm(dd[,c("aSize","hge")], xlab="Adult size (m)", ylab=expression(paste("Height gain efficency (m ",kg^-1,")")))
#scatterBarNorm(tab0[,c("Size","HGE")], xlab="Adult size (m)", ylab=expression(paste("Height gain efficency (m ",kg^-1,")")))
par(mfrow=c(1,1))


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
xx = xx[,1]


### JAGS MODEL
sink("model.txt")
  cat("
    model{

      # LIKELIHOOD
      for(i in 1:ss){
        for(j in 1:is[i]){
          y[i,j] ~ dnorm(mu[i],tau[i])
        }
        mu[i] = beta0 + x[i]*beta1 + beta[i]
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

## beta0  
aa0 = diff(c(0,10)) # range of possible values for the intercept
tau0 = (4/aa0)^2 # precision of beta0 according to the range method
m0 = aa0/2 # mean of beta0 assumed to be at the center of the range

## beta[i]
taub = tau0 # precision of beta[i] is set to be the same of beta0

## tau[i] 
ai = tau0 # minimum allowed value for tau[i] is the same of beta0
bi = 100 # maximum allowed precision for tau[i]

## beta1
range(dd$hge) # observed range of HGE
range(dd$aSize) # observed range of adult sizes
aa1 = diff(range(dd$hge))/diff(range(dd$aSize)) # maximum possible value for beta1
tau1 = (4/(2*aa1))^2 # precision of beta1 according to the range method

### JAGS RUN
# Jags input info
data = list(x=xx, y=yy, ss=ss, is=as.numeric(is), # raw data
            taub=taub,ai=ai,bi=bi,m0=m0,tau0=tau0,tau1=tau1) # constants 
inits = rep(list(list(
  beta=rep(0,ss), tau=rep(1,ss), # species-specific parameters
  beta0=0, beta1=0 # hyperparameters 
)),3) # number of chains
params = c("beta0","beta1","beta","mu","tau")

# Jags output
model.out = jags(data=data, inits=inits, parameter=params, "model.txt", 
                 n.chains=3, n.iter=501000, n.burnin=0, n.thin=100, DIC=F)
outparams = dimnames(model.out$BUGSoutput$sims.array)[[3]] # all parameter names estimated by jags

# Take the first 1000 runs as burnin
Output = AddBurnin(model.out$BUGSoutput$sims.array, burnin=1000,n.thin=1)

# Checking for MCMC convergence
pdf("acf.pdf",paper="letter")
par(mfrow=c(2,3))
for(i in outparams){
  acf(model.out$BUGSoutput$sims.array[1:5000, 1, i], lag.max= 160, main=i)
}
dev.off()
par(mfrow=c(1,1))

par(mfrow=c(2,3))
for(i in c("beta[1]","beta[10]","tau[1]","mu[1]","beta0","beta1")){
  acf(model.out$BUGSoutput$sims.array[1:5000, 1, i], lag.max= 160, main=i)
}
par(mfrow=c(1,1))

# Time series
cols = rainbow(3,alpha=0.7)
pdf("timeSeries.pdf",paper="letter")
par(mfrow=c(4,1))
for(i in outparams){
  plot(model.out$BUGSoutput$sims.array[1:1000, 1, i], type="l", col=cols[1], main=i, ylab="", xlab="Iteration")
  lines(model.out$BUGSoutput$sims.array[1:1000, 2, i], type="l", col=cols[2])
  lines(model.out$BUGSoutput$sims.array[1:1000, 3, i], type="l", col=cols[3])  
}
dev.off()
par(mfrow=c(1,1))

par(mfrow=c(3,2),mar=c(1, 4, 4, 2) + 0.1)
for(i in c("beta[1]","beta[10]","tau[1]","mu[1]","beta0","beta1")){
  plot(model.out$BUGSoutput$sims.array[1:1000, 1, i], type="l", col=cols[1], main=i, ylab="", xlab="Iteration")
  lines(model.out$BUGSoutput$sims.array[1:1000, 2, i], type="l", col=cols[2])
  lines(model.out$BUGSoutput$sims.array[1:1000, 3, i], type="l", col=cols[3])  
}
par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)


### RESULTS

## Main results
posts = Output$Burnin.sims.matrix
tab = round(Output$Burnin.Summary,3)[,1:4]
head(tab)

# Table1. Betas
tab1 = tab[c("beta0","beta1",paste("beta[",1:ss,"]",sep="")),]
tab1a = tab1[1:ceiling(nrow(tab1)/2),]
tab1b = rbind(tab1[ceiling(nrow(tab1)/2+1):nrow(tab1),],rep(NA,4))
tab1f = cbind(tab1a,cbind(row.names(tab1b),tab1b))
head(tab1f)
tail(tab1f)
xtable(tab1f)

# Table2. Mus
tab2 = tab[paste("mu[",1:ss,"]",sep=""),]
tab2a = tab2[1:ceiling(nrow(tab2)/2),]
tab2b = rbind(tab2[ceiling(nrow(tab2)/2+1):nrow(tab2),],rep(NA,4))
tab2f = cbind(tab2a,cbind(row.names(tab2b),tab2b))
head(tab2f)
tail(tab2f)
xtable(tab2f)

# Table3. Taus
tab3 = tab[paste("tau[",1:ss,"]",sep=""),]
tab3a = tab3[1:ceiling(nrow(tab3)/2),]
tab3b = rbind(tab3[ceiling(nrow(tab3)/2+1):nrow(tab3),],rep(NA,4))
tab3f = cbind(tab3a,cbind(row.names(tab3b),tab3b))
head(tab3f)
tail(tab3f)
xtable(tab3f)

## Regression parameters
par(mfrow=c(1,2))
plot(density(posts[,"beta0"]),xlab=expression(paste(beta[0])),main="",xlim=c(-1,10))
curve(dnorm(x,m0,sqrt(1/tau0)),add=T,lty=2)
legend(x=4,y=.5,lty=2:1,c("Prior","Posterior"),bty="n")
plot(density(posts[,"beta1"]),xlab=expression(paste(beta[1])),main="")
curve(dnorm(x,0,sqrt(1/tau1)),add=T,lty=2)
par(mfrow=c(1,1))

## Species-specific parameters
jpeg("spar.jpeg",height=800,quality=100,)
par(mfrow=c(3,1))
boxplot(posts[1:1000,paste("beta[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(beta[i])),xlab="Species")
abline(h=0,lty=3,col="red")

boxplot(posts[1:1000,paste("mu[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(mu[i])),xlab="Species")
#abline(h=mean(Output$Burnin.Summary[paste("mu[",1:ss,"]",sep=""),"mu.vect"]),lty=3,col="red")

boxplot(posts[1:1000,paste("tau[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(tau[i])),xlab="Species")
#abline(h=mean(Output$Burnin.Summary[paste("tau[",1:ss,"]",sep=""),"mu.vect"]),lty=3,col="red")
par(mfrow=c(1,1))
dev.off()

## Singletons
table(is)[1]/sum(table(is))

### SENSITIVITY
senan = function(ni=31000, m0=5, tau0=0.16, 
                 tau1=97.97005, ai=0.16, bi=100, taub=0.16){
  data = list(x=xx, y=yy, ss=ss, is=as.numeric(is),
              taub=taub,ai=ai,bi=bi,m0=m0,tau0=tau0,tau1=tau1) 
  model.out = jags(data=data, inits=NULL, parameter=params, "model.txt", n.chains=3, n.iter=ni, n.burnin=0, n.thin=20, DIC=F)
  Output = AddBurnin(model.out$BUGSoutput$sims.array, burnin=1000,n.thin=20)
  return(list(Output$Burnin.Summary,Output$Burnin.sims.matrix[,"beta1"]))
  }

## Set beta1 prior to be more or less variable
T1m = senan(tau1=tau1/2,ni=101000) # more variable
T1d = senan(tau1=tau1*2,ni=101000) # less variable

## Set beta0 prior to be centered around 0 or -5
M0 = senan(m0=0,ni=101000)
Mn5 = senan(m0=-5,ni=101000)

## Set beta[i] to be more or less variable
Tbm = senan(taub=taub/2,ni=101000)
Tbd = senan(taub=taub*2,ni=101000)

## Set lower limit of tau[i] to be 
A001 = senan(ai=0.01,ni=101000)
A50 = senan(ai=50,ni=101000)

## Plot sensitivities of priors parameters on beta1
par(mfrow=c(2,2))

# tau
rr = range(density(posts[,"beta1"])$y,
           density(T1m[[2]])$y,
           density(T1d[[2]])$y)
plot(density(posts[,"beta1"]),xlab="",main=expression(bold(paste(tau))),ylim=rr)
lines(density(T1m[[2]]),lty=2)
lines(density(T1d[[2]]),lty=3)
legend(x=.05,y=11,lty=1:3,c(expression(paste(tau," = 98*")),expression(paste(tau," = 49")),expression(paste(tau," = 196"))),bty="n",y.intersp=1.5,cex=.8)

# m0
rr = range(density(posts[,"beta1"])$y,
           density(M0[[2]])$y,
           density(Mn5[[2]])$y)
plot(density(posts[,"beta1"]),xlab="",main=expression(bold(paste(m[0]))),ylim=rr)
lines(density(M0[[2]]),lty=2)
lines(density(Mn5[[2]]),lty=3)
legend(x=.05,y=11,lty=1:3,c(expression(paste(m[0]," = 5*")),expression(paste(m[0]," = 0")),expression(paste(m[0]," = -5"))),bty="n",y.intersp=1.5,cex=.8)

# taub
rr = range(density(posts[,"beta1"])$y,
           density(Tbm[[2]])$y,
           density(Tbd[[2]])$y)
plot(density(posts[,"beta1"]),xlab="",main=expression(bold(paste(tau[b]))),ylim=rr)
lines(density(Tbm[[2]]),lty=2)
lines(density(Tbd[[2]]),lty=3)
legend(x=.03,y=14,lty=1:3,c(expression(paste(tau[b]," = 0.16*")),expression(paste(tau[b]," = 0.08")),expression(paste(tau[b]," = 0.32"))),bty="n",y.intersp=1.5,cex=.8)

# a
rr = range(density(posts[,"beta1"])$y,
           density(A001[[2]])$y,
           density(A50[[2]])$y)
plot(density(posts[,"beta1"]),xlab="",main=expression(bold(paste(a))),ylim=rr)
lines(density(A001[[2]]),lty=2)
lines(density(A50[[2]]),lty=3)
legend(x=.03,y=10,lty=1:3,c(expression(paste(a," = 0.16*")),expression(paste(a," = 0.01")),expression(paste(a," = 0.50"))),bty="n",y.intersp=1.5,cex=.8)

par(mfrow=c(1,1))








