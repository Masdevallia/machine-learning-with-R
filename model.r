
# Clear the workspace
rm(list=ls())

# Packages loading
library(mclust) 
library(MASS)
library(ecodist) # distances
library(proxy) # similarities
library(gdata) # drop.levels
library(FD) # gowdis

# Working directory
# getwd()
setwd("~/../Ironhack/machine-learning-with-R/dataset")

# Input
data <- read.table("data.txt", header = TRUE, sep = "\t", dec=".")

# .........................................................................................................

# TRAITSPACE ON THE DATA:
# Traitspace Model (Chaitanya Joshi and Daniel C Laughlin, Univ. of Waikato, NZ).
# Predicting relative abundances of species given the environmental conditions.

# .........................................................................................................

# STEP 1a: Generalised Linear Model (GLM) to compute P(T/E)

# Dataset for prediction - values on which predictions are needed
env <- data[,6]
senv <- unique(env)
env2 <- data[,7]
senv2 <- unique(env2)
env_p <- data.frame(env=senv,env2=senv2)

# ...............................................................

# Fitting a linear model on trait 1: height

trait1 <- data[,2]
      
# lineal:
# lm1 = lm(log(trait1)~env+env2)
      
# cuadrática:
lm1 = lm(log(trait1)~env+I(env^2)+env2+I(env2^2)) # RAW POLYNOMIALS 
# lm1 = lm(log(trait1)~poly(env, 2,raw = TRUE)+ poly(env2, 2,raw = TRUE)) # RAW POLYNOMIALS 
# lm1 = lm(log(trait1)~poly(env, 2,raw = FALSE)+ poly(env2, 2,raw = FALSE)) # ORTHOGONAL POLYNOMIALS

# cúbica:
# lm1 = lm(log(trait1)~env+I(env^2)+I(env^3)+env2+I(env2^2)+I(env2^3)) # RAW POLYNOMIALS 
# lm1 = lm(log(trait1)~poly(env, 3,raw = TRUE)+ poly(env2, 3,raw = TRUE)) # RAW POLYNOMIALS 
# lm1 = lm(log(trait1)~poly(env, 3,raw = FALSE)+ poly(env2, 3,raw = FALSE)) # ORTHOGONAL POLYNOMIALS

summary(lm1)
std1=sd(lm1$residuals)
pred1=predict.lm(lm1,env_p,se.fit=TRUE,interval="prediction",level=0.95)
# pred1=predict(lm1,newdata=env_p,se.fit=TRUE,interval="prediction",level=0.95)

# ...............................................................

# Fitting a linear model on trait 2: area

trait2 <- data[,3]
lm2 = lm(log(trait2)~env+I(env^2)+env2+I(env2^2)) # RAW POLYNOMIALS
summary(lm2)
std2=sd(lm2$residuals)
pred2=predict.lm(lm2,env_p,se.fit=TRUE,interval="prediction",level=0.95)

# ...............................................................

# Fitting a linear model on trait 3: sla

trait3 <- data[,4]
lm3 = lm(log(trait3)~env+I(env^2)+env2+I(env2^2)) # RAW POLYNOMIALS
summary(lm3)
std3=sd(lm3$residuals)
pred3=predict.lm(lm3,env_p,se.fit=TRUE,interval="prediction",level=0.95)

# .........................................................................................................

# STEP 1b: Using Mclust to compute P(T/Sk)
# I'll do this step later together with the step 2b in a loop.

# .........................................................................................................

# STEP 2a: Drawing samples from P(T/E)

# Using the set of environmental conditions at which you are interested in predicting species abundances.
# For example, to test and validate the model, you would use environmental conditions measured at plot locations where
# species abundances are known.

# Using environmental conditions measured at plot locations where species abundances are known to test and validate the model:

N = 1000 #sample size per location on the environmental gradient

trait1_sample=rep(0,length(senv)*N)
trait2_sample=rep(0,length(senv)*N)
trait3_sample=rep(0,length(senv)*N)
temp_prd=rep(0,length(senv)*N) # temp array for plotting (env)
temp_prd_2=rep(0,length(senv)*N) # temp array for plotting (env2)

ct=1
for (j in 1:length(senv))
  { 
    ct=(j-1)*N+1
    trait1_sample[ct:(ct+(N-1))]= rnorm(N,pred1[[1]][j],std1)
    trait2_sample[ct:(ct+(N-1))]= rnorm(N,pred2[[1]][j],std2)
    trait3_sample[ct:(ct+(N-1))]= rnorm(N,pred3[[1]][j],std3)
    temp_prd[ct:(ct+(N-1))]= senv[j]
    temp_prd_2[ct:(ct+(N-1))]= senv2[j]
  }

trt_sample = cbind(trait1_sample,trait2_sample,trait3_sample)

# computing(P(T/E))

P_T1_E=rep(0,length(senv)*N)
P_T2_E=rep(0,length(senv)*N)
P_T3_E=rep(0,length(senv)*N)

ct=1
for (j in 1:length(senv))
  { 
    ct=(j-1)*N+1
     P_T1_E[ct:(ct+(N-1))]= dnorm(trait1_sample[ct:(ct+(N-1))],pred1[[1]][j],std1)
     P_T2_E[ct:(ct+(N-1))]= dnorm(trait2_sample[ct:(ct+(N-1))],pred2[[1]][j],std2)
     P_T3_E[ct:(ct+(N-1))]= dnorm(trait3_sample[ct:(ct+(N-1))],pred3[[1]][j],std3)
  }

P_T_E=exp(log(P_T1_E)+log(P_T2_E)+log(P_T3_E))

# .........................................................................................................

# Step 2b: Computing the likelihood P(T/Sk) using Mclust done earlier

# Applying Mclust together with the step 2b in a loop:

species <- unique(data$species)
par<-list()
P_T_S <- matrix(0,nrow(trt_sample),length(species))

for(i in 1:length(species)) # for each species
{
  datos.aux <- data[which(data$species==species[i]),2:4]
  pdf<- Mclust(datos.aux,warn=TRUE)
  par[[i]]<- pdf$parameters
  P_T_S[,i]<- dens(pdf$modelName,exp(trt_sample),parameters=par[[i]]) 
}

# .........................................................................................................

# Step 2c: Computing posterior P(Sk/T,E) using Bayes theorem

P_T_S_pr = P_T_S/length(species)  
# multiplying likelihood by flat prior - numerator in Bayes thm
# To put the same "prior" for all the species -> Divide P_T_S by the number of species (uninformative prior)                                    

P_T_S_pr_sum=apply(P_T_S_pr,1,sum) # denominator in Bayes thm.

P_S_T_E = matrix(0,dim(P_T_S)[1],length(species))

for (i in 1:dim(P_T_S)[1]) # for each location
  { P_S_T_E[i,] = exp(log(P_T_S_pr[i,]) - log(P_T_S_pr_sum[i]))} # using log

# .........................................................................................................

# Step 2d: Posterior P(Sk/T) by integrating out T's
      
P_S_E_all = matrix(0,length(senv)*N,length(species))
P_S_E_unnorm = matrix(0,length(senv),length(species)) # unnormalised P_S_E
P_S_E = matrix(0,length(senv),length(species))

# Computing the integrand (with log)

for (i in 1:dim(P_S_E_all)[1])
	{ P_S_E_all[i,]=exp(log(P_T_E[i])+log(P_S_T_E[i,]))}

# MC integration and normalisation

c=1
for (k in 1:length(senv))
    { 
      c=(k-1)*N+1
      P_S_E_unnorm[k,]=apply(P_S_E_all[c:(c+(N-1)),],2,mean) # MC
      P_S_E[k,]=P_S_E_unnorm[k,]/sum(P_S_E_unnorm[k,]) # normalisation
    }

  apply(P_S_E,1,sum) # should produce a vector of 1's indicating that valid posterior probability distributions have been obtained.

# Inspect predictions for each species (columns) at each location (rows)
P_S_E 

# Code to calculate probability distributions for each species:
    
# prob distribution for each species - normalising each column

P_S_E_distrn = matrix(0,nrow(env_p),length(species))

for (l in 1:length(species))
    {
      P_S_E_distrn[,l]= P_S_E_unnorm[,l]/sum(P_S_E_unnorm[,l])
    }
    
apply(P_S_E_distrn,2,sum) # abundances are relative to species (total abundance of each species along the gradient = 1)
apply(P_S_E,1,sum) # abundances are relative to plots (total abundance of each plot = 1)

# .........................................................................................................

# End of Program

# colnames(P_S_E) <- species
# write.table(P_S_E, "P_S_E.txt", sep = "\t", dec=".", row.names=T, col.names=T)
# colnames(P_S_E_distrn) <- colnames(P_S_E)  
# write.table(P_S_E_distrn, "P_S_E_distrn.txt", sep = "\t", dec=".", row.names=T, col.names=T)
