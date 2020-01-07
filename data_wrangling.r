
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
traits <- read.table("traits.txt", header = TRUE, sep = "\t", dec=".")
env <- read.table("enviro.txt", header = TRUE, sep = "\t", dec=".")

# .........................................................................................................

# PREPARING DATA FOR APPLYING TRAITSPACE MODEL:

# .........................................................................................................

# Environmental conditions: moisture, phosphorus, potassium, nitrogen -> PCA
pca <- princomp(env, cor = TRUE)
# summary(pca)
# loadings(pca)
# plot(pca)
# biplot(pca)
env.pca <- matrix(0,40,3)
env.pca[,1] <- unique(traits$plot)
env.pca[,2] <- pca$scores[,1] # Explains 55%
env.pca[,3] <- pca$scores[,2] # Explains 20% (Cumulative Proportion: 75%)
colnames(env.pca) <- c("plot","productivity_1","productivity_2")

# .........................................................................................................

# Traits -> Data wrangling:

plots <- traits[,1:3]

height.data <- traits[,4:13]
height.data <- cbind(plots, height.data) 
height <- height.data[,4:13]
height <- as.vector(t(height))

area.data <- traits[,34:43]
area.data <- cbind(plots, area.data)
area <- area.data[,4:13]
area <- as.vector(t(area))

sla.data <- traits[,44:53]
sla.data <- cbind(plots, sla.data)
sla <- sla.data[,4:13]
sla <- as.vector(t(sla))

# ...........................................................

data <- as.data.frame(matrix(0,length(height),7))
data[,1] <- rep(height.data$species, each=10)
data[,2] <- height
data[,3] <- area
data[,4] <- sla
data[,5] <- rep(height.data$plot, each=10)
data[,6] <- rep(height.data$plot, each=10)
data[,7] <- rep(height.data$plot, each=10)
colnames(data) <-c("species","height","area","sla","plot", "env.var","env.var.2")    

# ...........................................................

# Replacing plot by environmental variable 1 in column 6
for(i in 1:nrow(data))
  {
    for (j in 1:nrow(env.pca))
      {
      if(data$env.var[i]==env.pca[j,1]){data$env.var[i]=env.pca[j,2]}
      }
  }
      
# Replacing plot by environmental variable 2 in column 7
for(i in 1:nrow(data))
  {
    for (j in 1:nrow(env.pca))
      {
      if(data$env.var.2[i]==env.pca[j,1]){data$env.var.2[i]=env.pca[j,3]}
      }
  }

# ...........................................................

# Deleting rows that have NA on at least one of the traits:
data <- data[which(data[,2]!="NA"&data[,3]!="NA"&data[,4]!="NA"),]

# ...........................................................

write.table(data, "data.txt", sep = "\t", dec=".", row.names=T, col.names=T)
