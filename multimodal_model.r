
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

# STEP 1a: Multimodal model: Uses Mlcust instead of regression: Sampling traits from mixture modelling.

# Dataset for prediction - values on which predictions are needed
env <- data[,6]
senv <- unique(env)
env2 <- data[,7]
senv2 <- unique(env2)
env_p <- data.frame(env=senv,env2=senv2)

# I'll do this step later together with the step 2a in a loop.

# .........................................................................................................

# STEP 1b: Using Mclust to compute P(T/Sk)
# I'll do this step later together with the step 2b in a loop.

# .........................................................................................................

# STEP 2a: Drawing samples using Mclust done earlier

# Drawing samples from mixture densities fitted in Step 1A
# Applying Mclust now together with the step 2a in a loop.

N = 1000 #sample size per location on the environmental gradient
trt_sample = matrix(0,length(senv)*N,3)
site <- list()
par_site <- list()
pdf_site <- list()

for (i in 1:length(senv))
  {
    site[[i]] <- log(data[which(data$env.var==senv[i]),2:4])  
    pdf_site[[i]] <- Mclust(site[[i]],warn=TRUE)
    par_site[[i]] = pdf_site[[i]]$parameters
    trt_sample[(((i-1)*N)+1):(i*N),] <- sim(pdf_site[[i]]$modelName, par_site[[i]],N)[,2:4]
  }

# Computing (P(T/E))

P_T_E = rep(0,length(senv)*N)

for (i in 1:length(senv))
  {
    P_T_E[(((i-1)*N)+1):(i*N)] <- dens(pdf_site[[i]]$modelName,trt_sample[(((i-1)*N)+1):(i*N),],parameters=par_site[[i]])
  }

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
{ P_S_T_E[i,] = exp(log(P_T_S_pr[i,]) - log(P_T_S_pr_sum[i]))} #using log

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

apply(P_S_E,1,sum)  # should produce a vector of 1's indicating that valid posterior probability distributions have been obtained.

# Inspect predictions for each species (columns) at each location (rows)

P_S_E

colnames(P_S_E) <- species 

# Code to calculate probability distributions for each species:
# prob distribution for each species - normalising each column

P_S_E_distrn = matrix(0,nrow(env_p),length(species))

for (l in 1:length(species))
  {
    P_S_E_distrn[,l]= P_S_E_unnorm[,l]/sum(P_S_E_unnorm[,l],na.rm=TRUE)
  }
  
apply(P_S_E_distrn,2,sum, na.rm=TRUE) # abundances are relative to species (total abundance of each species along the gradient = 1)
apply(P_S_E,1,sum) # abundances are relative to plots (total abundance of each plot = 1)

colnames(P_S_E_distrn) <- colnames(P_S_E)  

# .........................................................................................................

# End of Program

# write.table(P_S_E, "P_S_E.txt", sep = "\t", dec=".", row.names=T, col.names=T)
# write.csv(P_S_E,"P_S_E.csv")
# P_S_E <- read.table("P_S_E.txt", sep = "\t", dec=".") 

# write.table(P_S_E_distrn, "P_S_E_distrn.txt", sep = "\t", dec=".", row.names=T, col.names=T)
# write.csv(P_S_E_distrn,"P_S_E_distrn.csv")
# P_S_E_distrn <- read.table("P_S_E_distrn.txt", sep = "\t", dec=".")
