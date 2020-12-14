source(url("http://bioinf.wehi.edu.au/voom/tspm.R"))

# Get distribution function of abundance proportions
# This distribution was generated from a real dataset
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))



######################################## initial set up ########################################

ngenes <- 2000 #How many genes (OTUs)

baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )

baselineprop <- baselineprop/sum(baselineprop) #For each gene, get their baseline proportion

nlibs <- 12 #How many libraries (samples)

ntrt <- 2 #How many treatment group

nbatch <- 2 #How many batch group

trt <- rep(1:ntrt,each=nlibs/ntrt) #Assign treatment, treatment group sizes are equal
#  trt <- sample(1:ntrt, nlibs, replace = TRUE) # Randomly assign treatment

#batch <- sample(1:nbatch, nlibs, replace = TRUE) # Randomly assign batch
batch <- rep(rep(1:nbatch, each = nlibs/ntrt/nbatch),2) # Blanced Block design

expected.lib.size <- 4e4 #Equal library sizes
# expected.lib.size <- 11e6* runif(nlibs,min = 0.1, max = 1) #Uneuqal library sizes



######################################## Differential expression ########################################

ndiff.trt <- round(runif(ntrt, ngenes/100, ngenes/50)) #How many genes are differentially expressed in each treatment group (1% to 2%)

dif.express <- list() #Empty list to store which genes are differentially expressed in each treatment group

fc.trt <- list() #Empty list to store Fold Change rate at which genes are differentially expressed in each treatment group

for(i in 1:ntrt){
  dif.express[[i]] <- sample(1:ngenes, ndiff.trt[i]) #Sample which genes are differentially expressed in each treatment group
  fc.trt[[i]] <- runif(ndiff.trt[i], 1, 3) #Sample FC of DEgenes in each treatment group
}

#Note that genes are differentially expressed relative to their baseline expression level.




######################################## Batch effect, biological ########################################

ndiff.bat <- round(runif(nbatch, ngenes/5, ngenes/3.3)) #How many genes are affected by batch effect in each batch group (20% to 30%)

dif.batch <- list() #Empty list to store which genes are affected by batch effect in each batch group

fc.bat <- list() #Empty list to store Fold Change rate at which genes are affected by batch effect in each batch group

for(i in 1:ntrt){
  dif.batch[[i]] <- sample(1:ngenes, ndiff.bat[i]) #Sample which genes are affected by batch effect in each batch group
  fc.bat[[i]] <- runif(ndiff.bat[i], 1, 3) #Sample FC of affected genes in each batch group
}

#Note that genes are differentially expressed relative to their baseline expression level.




######################################## counting data generation ########################################

prop.matrix <- matrix(rep(baselineprop,nlibs),ncol = nlibs) #set up the baseline proportion matrix

for(i in 1:ntrt){
  prop.matrix[dif.express[[i]],trt == i] <- prop.matrix[dif.express[[i]],trt == i] * fc.trt[[i]] #Differential expression by multiply the FC rates
}


for(i in 1:nbatch){
  prop.matrix[dif.batch[[i]],batch == i] <- prop.matrix[dif.batch[[i]],batch == i] * fc.bat[[i]] #Differential expression by multiply the FC rates
}


mu0 <- t(t(prop.matrix) * expected.lib.size) #set up the baseline counting matrix

BCV0 <- 0.2+1/sqrt(mu0) #Global trend of dispersion (See voom for details)

df.BCV <- 40

BCV <- BCV0*sqrt(df.BCV/rchisq(ngenes,df=df.BCV)) #Add biological noise

shape <- 1/BCV^2 #Shape of gamma distribution

scale <- mu0/shape #Scale of gamma distribution

mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs) # simulate negative binomial distribution with Gamma-Poisson mixture

counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)






plot(density(counts)) #Check for sparsity


