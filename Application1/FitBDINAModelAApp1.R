##-------------------------------------------------------------------------------
## This script implement the Chapter Bounded CDM model of
## the Minicourse
## Bazán, J. L. (2026). Regression and Latent Variable Models for Bounded Outcomes.
## XLVIII Jornadas Nacionales de Estadística 2026. Universidad de Valparaíso, 
## Octuber 20-23, 2026.
##--------------------------------------------------------------------------------


####Loading the packages 

library("rstudioapi")
library("coda")
library("readxl")
library("R2jags")
library("loo")
library("psych")

#### Setting the directory and loading the functions

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#### Reading Q-Matrix and data set

Q <- read.csv("QApp1.csv")
str(Q)
Q<-Q[,-1]
Q
dim(Q)
#19x3

datac <- read.csv("dataApp1.csv")
str(datac)
datac<-datac[,-1]
str(datac)
dim(datac)
#1668x19
describe(datac)

#### Preparing data set to JAGS

N=nrow(datac) 
J=ncol(datac) 
K=ncol(Q) 
C=2^K 
delta0 <- rep(1,2^K) #vector with size of possible profiles

#Matrix containing all possible profiles
As <- rep(0,K)
for(j in 1:K){
  temp <- combn(1:K,m=j)
  tempmat <- matrix(0,ncol(temp),K)
  for(j in 1:ncol(temp))
    tempmat[j,temp[,j]] <- 1
  As <- rbind(As,tempmat)
}
As <- as.matrix(As)


para <- c("u0", "phi0", "u1", "phi1", "pai","loglik","ppp")
data <- list(Y=datac, N=N, J=J, K=K, 
             C=C, q=Q, all.patterns=As, delta=delta0)

#### Fitting Model with JAGS

start_time <- proc.time()
jagsfit <- jags.parallel(data, inits=NULL, para, "BDINAmodelA.txt",  n.chains=2, 
                         n.iter = 5000, 
                         n.burnin = 1000, n.thin=1, jags.seed = 123)
end_time <- proc.time()
jags_time <- end_time-start_time

#### Saving results
results <- jagsfit$BUGSoutput$summary
results 
dim(results) 

write.csv(results[1:nrow(results),],"ResultsBDINAmodelAApp1.csv")

#### Post Analysis

jags_time
jagsfit
summary(jagsfit)
jagsfit$DIC
jagsfit$parameters.to.save


#Computing Item parameters
#deviance
results[1,]
#loglik 2 31692 (1668*19)
results[2:31693]

#pai
results[31694:31701]

#phio
results[31702:31720]

#phi1
results[31721:31739]

#ppp
results[31740]

#muo
results[31741:31759]

#mu1
results[31760:31778]

#model comparison criteria
library(loo)
#extract log-likelihood
loglik_m1<-jagsfit$BUGSoutput$sims.list$loglik
#waic
waic_m1<-waic(loglik_m1)
waic_m1
#looic
looic_m1<-loo(loglik_m1)
looic_m1

#### Saving Rdata

save.image("FitBDINAModelAApp1.RData")


