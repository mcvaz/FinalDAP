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

## Table 0. Data



hist(dd$aSize,xlab="Adult size (m)",main="")
hist(aggregate(aSize~sp,dd,mean)$aSize,xlab="Adult size (m)",main="")

par(mfrow=c(2,2))
plot(hge~aSize,dd,xlab="Adult size (m)",ylab=expression(paste("Height gain efficency (m ",kg^-1,")")))
hist(dd$hge,main="",xlab="")
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
data = list(x=xx, y=yy, ss=ss, is=as.numeric(is), # raw data
            taub=taub,ai=ai,bi=bi,m0=m0,tau0=tau0,tau1=tau1) # constants 
inits = rep(list(list(
  beta=rep(0,ss), tau=rep(0.01,ss), # species-specific parameters
  beta0=0, beta1=0 # hyperparameters 
)),3) # number of chains
params = c("beta0","beta1","beta","mu","tau")

# Jags output
model.out = jags(data=data, inits=inits, parameter=params, "model.txt", 
                 n.chains=3, n.iter=101000, n.burnin=0, n.thin=20, DIC=F)
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
plot(density(posts[,"beta0"]),xlab=expression(paste(beta[0])),main="")
plot(density(posts[,"beta1"]),xlab=expression(paste(beta[1])),main="")
par(mfrow=c(1,1))

## Species-specific parameters
boxplot(posts[1:1000,paste("beta[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(beta[i])),xlab="Species")
abline(h=0,lty=3,col="red")

boxplot(posts[1:1000,paste("mu[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(mu[i])),xlab="Species")
abline(h=mean(Output$Burnin.Summary[paste("mu[",1:ss,"]",sep=""),"mu.vect"]),lty=3,col="red")

boxplot(posts[1:1000,paste("tau[",1:ss,"]",sep="")],range=0,xaxt="n",ylab=expression(paste(tau[i])),xlab="Species")
abline(h=mean(Output$Burnin.Summary[paste("tau[",1:ss,"]",sep=""),"mu.vect"]),lty=3,col="red")


### SENSITIVITY







